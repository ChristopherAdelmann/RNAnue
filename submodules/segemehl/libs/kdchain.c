/*
 *   segemehl - a read aligner
 *   Copyright (C) 2008-2017  Steve Hoffmann and Christian Otto
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */




/*
 *  kdchain.c
 *  implementation of kdchain
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 04/29/2008 07:01:30 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 85 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-11-18 15:34:44 +0100 (Tue, 18 Nov 2008) $
 *
 *  Id: $Id: kdchain.c 85 2008-11-18 14:34:44Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/kdchain.c $
 *  
 */

#include "manout.h"
#include "kdchain.h"
#include "mathematics.h"
#include "sufarray.h"
#include "container.h"
#include "kdseed.h"
#include "debug.h"
#include "karlin.h"
#include "bitvectoralg.h"
#include "iupac.h"
#include <assert.h>
#include <limits.h>
#include <unistd.h>
#include <float.h>
#include <math.h>


/*----------------------------- minDistfragment ------------------------------
 *    
 * @brief minimum distance of fragment hits
 * @author Steve Hoffmann 
 *   
 */
 
Uint
minDistFragmentHits (Suffixarray *arr, branchfragment_t *u, branchfragment_t *v)
{

  Uint i, j, d1idx, d2idx;
  Uint mindist = UINT_MAX;
  Uint d = UINT_MAX;

  for(i=u->branch->l; i <= u->branch->r; i++) {
    for(j=v->branch->l; j <= v->branch->r; j++) {
      d1idx = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[i]]);
      d2idx = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[j]]);
      if(d1idx == d2idx && (d=llabs((Lint)arr->suftab[i] - arr->suftab[j])) < mindist) {
        mindist = d;
      } 
    }
  }

  return d;
}


/*-------------------------------- wrapChains --------------------------------
 *    
 * @brief remove chains from heap
 * @author Steve Hoffmann 
 *   
 */
 
void
wrapChains(void *space, branchChain_t *chains, Uint noofchains) {
  Uint i,j;

  for(i=0; i < noofchains; i++) { 
    for(j=0; j < chains[i].nooffragments; j++) {
      FREEMEMORY(space, chains[i].f[j]);
    }
    FREEMEMORY(space, chains[i].f);
  }
  return;
}



/*---------------------------------- chain -----------------------------------
 *    
 * @brief add fragment to a chain
 * @author Steve Hoffmann 
 *   
 */
 
void
chain(branchChain_t *chain, branchfragment_t *f) {
  chain->score = chainscore(chain, f);
  chain->end = f->end;
  chain->nooffragments++;
  chain->f = ALLOCMEMORY(space, chain->f, branchfragment_t*, chain->nooffragments);
  chain->f[chain->nooffragments-1] = f;
}


/*------------------------------- fragmentovl --------------------------------
 *    
 * @brief fragment overlap
 * @author Steve Hoffmann 
 *                         f1->end
 *   f1 -----------------|
 *                  |--------------
 *                  f2->start
 */
 
int
fragmentovl (branchfragment_t *f1, branchfragment_t *f2)
{
  return (f1->end >= f2->start && f1->start <= f2->start) ? f1->end - f2->start : 0;
}

/*--------------------------------- chainovl ---------------------------------
 *    
 * @brief overlap of a fragment with chain on the query!
 * @author Steve Hoffmann 
 *   
 */
 

Lint
chainovl(branchChain_t *chain, branchfragment_t *f) {
  return ((Lint)chain->end - (Lint)f->start)+1;
}


/*-------------------------------- chainscore --------------------------------
 *    
 * @brief get score of a chain when a fragment is added to it
 * @author Steve Hoffmann 
 *   
 */
 

int
chainscore (branchChain_t *chain, branchfragment_t *f) {
  Lint ovl = 0;
  int score = 0;

  ovl = chainovl(chain, f); //v -- 0 or -ovl
  ovl = (ovl < 0) ? 0 : ovl;

  score = (chain->score + f->score) - ovl;
  return score;
}



/*------------------------------- chainscore2 --------------------------------
 *    
 * @brief get score of chain with two fragments chained
 * @author Steve Hoffmann 
 *   
 */
 
int
chainscore2 (branchChain_t *chain, branchfragment_t *f1, branchfragment_t *f2)
{
  Lint ovl =0;
  Lint ovl2=0;
  int score =0;

  ovl = chainovl(chain, f1);
  ovl = (ovl < 0) ? -ovl : ovl;

  score = (chain->score + f1->score) - ovl;
	
  ovl2 = (Lint)f1->end - (Lint)f2->start + 1;
  ovl2 = (ovl2 < 0) ? -ovl2 : ovl2;

  score += f2->score - ovl2; 
    
  return score;
}

/*----------------------------- cmp_chainscores ------------------------------
 *    
 * @brief compare the scores of a chain (qsort)
 * @author Steve Hoffmann 
 *   
 */
 
int
cmp_chainscores (const void *a, const void *b) {
  branchChain_t *first = (branchChain_t *) a;
  branchChain_t *second = (branchChain_t *) b;

  if(first->score < second->score) return 1;
  if(first->score == second->score) return 0;

  return -1;
}


/*---------------------------- cmp_chainlocality -----------------------------
 *    
 * @brief compare locality of chain
 * @author Steve Hoffmann 
 *   
 */
 
int
cmp_chainlocality (const void *a, const void *b)
{
  branchChain_t *first = (branchChain_t *) a;
  branchChain_t *second = (branchChain_t *) b;

  Uint i, swtch1=0, swtch2=0, strandswtch1=0, strandswtch2=0;

  for(i=1; i < first->nooffragments; i++) {
    if( first->f[i]->subidx != first->f[i-1]->subidx 
        || dist_uint(first->f[i]->substart, first->f[i-1]->substart) > 200000) 
      swtch1++;
    if(first->f[i]->strand != first->f[i-1]->strand) strandswtch1++;
  }

  for(i=1; i < second->nooffragments; i++) {
    if( second->f[i]->subidx != second->f[i-1]->subidx 
        || dist_uint(second->f[i]->substart, second->f[i-1]->substart) > 200000) 
      swtch2++;
    if(second->f[i]->strand != second->f[i-1]->strand) strandswtch2++; 
  }

  if(swtch1 > swtch2) return 1;
  if(swtch1 < swtch2) return -1;
  if(swtch1 == swtch2) {
    //fprintf(stderr, "strandswitch called (%d == %d); %d ? %d\n", swtch1, swtch2, strandswtch1, strandswtch2);
    //if(strandswtch1 > strandswtch2) return 1;
    //if(strandswtch1 < strandswtch2) return -1;
  }

  if (first->score < second->score) return 1;
  if (first->score > second->score) return -1;
  
  return 0;
}

/*--------------------------- cmp_branchfragments ----------------------------
 *    
 * @brief compare start positions of branch fragmentsi (qsort)
 * @author Steve Hoffmann 
 *   
 */

int
cmp_branchfragments (const void *a, const void *b) {
  branchfragment_t *first = (branchfragment_t*) a;
  branchfragment_t *second = (branchfragment_t*) b;

  if(first->start < second->start) return -1;
  if(first->start == second->start) return 0;

  return 1;
}

/*--------------------------- cmp_branchfragmentsptr ----------------------------
 *    
 * @brief compare start positions of branch fragmentsi (qsort)
 * @author Steve Hoffmann 
 *   
 */

int
cmp_branchfragmentsptr (const void *a, const void *b) {
  branchfragment_t **first = (branchfragment_t**) a;
  branchfragment_t **second = (branchfragment_t**) b;

  if(first[0]->start < second[0]->start) return -1;
  if(first[0]->start == second[0]->start) return 0;

  return 1;
}

/*--------------------------- cmp_branchfragmentssub ----------------------------
 *    
 * @brief compare start positions of branch fragmentsi (qsort)
 * @author Steve Hoffmann 
 *   
 */

int
cmp_branchfragmentssub (const void *a, const void *b) {
  branchfragment_t **first = (branchfragment_t**) a;
  branchfragment_t **second = (branchfragment_t**) b;

  if(first[0]->substart < second[0]->substart) return -1;
  if(first[0]->substart == second[0]->substart) return 0;

  return 1;
}

/*-------------------------------- initChains --------------------------------
 *    
 * @brief initalize chains
 * @author Steve Hoffmann 
 *   
 */
 
branchChain_t*
initChains (branchChain_t *chains, Uint k)
{
  Uint i;

  for(i=0; i < k; i++) {
    chains[i].start = 0;
    chains[i].end = 0;
    chains[i].f = NULL;
    chains[i].score= 0;
    chains[i].nooffragments=0;
  }
	
  return chains;
}




/*------------------------------ condenseChain -------------------------------
 *    
 * @brief merge fragments of the chain if they are too close on the reference 
 * practically [u.v][x,y] -> [u,y]; 
 * if(u<x &% u<y) 
 *  if(x<v || x-v < 50) 
 *      merge 
 *  else
 *      dont merge
 * @author Steve Hoffmann 
 *   
 */

  branchChain_t*
condenseChain (branchChain_t * chains, Uint noofchains, MultiCharSeq *seq, 
    Suffixarray *arr)
{
  Uint i, j, u, v, x, y, k, h, w, len1, len2, strand1, strand2, 
       chr1, chr2, nochain, d1idx=0, d2idx=0, d3idx=0, l, r, ll, rr, p=0, q=0;
  double  mindist = DBL_MAX, d1=0, d2=0, di, dj;
  branchChain_t *newchains = NULL;
  Uint sub_start, sub_end, subidx, ***bd, *beststart=NULL, *cnt=NULL;

  bd  = ALLOCMEMORY(space, NULL, Uint**, noofchains);
  cnt = ALLOCMEMORY(space, NULL, Uint, noofchains);
  beststart = ALLOCMEMORY(space, NULL, Uint, noofchains);



  for(k=0; k < noofchains; k++){ 

    l = chains[k].f[0]->branch->l;
    r = chains[k].f[0]->branch->r;
    mindist = DBL_MAX;

    bd[k] = ALLOCMEMORY(space, NULL, Uint*, r-l+1);
    cnt[k] = r-l+1;

    /******
     * first step: minimize distance of fragment hits within chain
     * for all possible start loci in [l,r] select
     * a chain of closest loci
     ******/
    for(i=0, p=l; p <= r; p++, i++) {
      bd[k][i] = calloc(chains[k].nooffragments, sizeof(Uint));
//      fprintf(stdout, "\n\n start chain: %d, iteration:%d, fragment:%d -> %u\n", k, i, 0, arr->suftab[p]);
      bd[k][i][0] = p;

      for(di=0, j=1; j < chains[k].nooffragments; j++) {
        ll = chains[k].f[j]->branch->l;
        rr = chains[k].f[j]->branch->r;

        for(dj=0, q=ll; q <= rr; q++) {
          d1idx = getMultiCharSeqIndex(seq, 
              &seq->sequences[arr->suftab[q]]);
          d2idx = getMultiCharSeqIndex(seq, 
              &seq->sequences[arr->suftab[bd[k][i][j-1]]]);

          if(d1idx != d2idx) {
            d1 = UINT_MAX;
          } else {
            d1 = llabs((LLint) arr->suftab[bd[k][i][j-1]] - arr->suftab[q]);
          }

          if(bd[k][i][j]) {
            d3idx = getMultiCharSeqIndex(seq, 
                &seq->sequences[arr->suftab[bd[k][i][j]]]);
            d2 = llabs((LLint) arr->suftab[bd[k][i][j-1]] - arr->suftab[bd[k][i][j]]);
          }

          if(d3idx != d2idx) {
            d2 = UINT_MAX;
          }

          if(!bd[k][i][j] || d1 < d2) {
 //           fprintf(stdout, "selected for chain: %d, iteration:%d, fragment:%d -> %u (dist:%f)\n", k, i, j, arr->suftab[q], d1);
            bd[k][i][j] = q;
            dj = d1;
          }       
        }
        di += dj;
      }

      if(di < mindist) {
 //       fprintf(stdout, "-> iteration %d selected for chain %d with new mindist %f\n", i, k, di);
        beststart[k] = i;
        mindist = di;
      }
    }
  }

  /*decode substarts to real coordinates*/
  for(i=0; i < noofchains; i++) {    
//    fprintf(stdout, "assigning for chain %d - beststart[%d]=%d\n", i, i, beststart[i]);
    for(j=0; j < chains[i].nooffragments; j++) {
//      fprintf(stdout, "assigning pos bd[%d][%d][%d] -> %u\n", i, beststart[i], j, arr->suftab[bd[i][beststart[i]][j]]);
      chains[i].f[j]->substart = arr->suftab[bd[i][beststart[i]][j]];
//      fprintf(stdout, "chains[%d].f[%d] = %u\n", i, j, chains[i].f[j]->substart);
    } 
//    fprintf(stdout, "\n");
  }

  /******
   *second step: merge fragments that are close
   ******/

  for(i=0; i < noofchains; i++) {    
    for(j=0; j < chains[i].nooffragments-1; j++) {
      for(k=j+1; k < chains[i].nooffragments; k++) { 
        //   for(k=j+1; k < j+2; k++) { 

        len1 = chains[i].f[j]->end - chains[i].f[j]->start;
        strand1 = chains[i].f[j]->strand;
        h = chains[i].f[j]->end;
        u = arr->suftab[bd[i][beststart[i]][j]];
        v = u + len1;
        chr1 = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[u]);

        len2 = chains[i].f[k]->end - chains[i].f[k]->start;
        strand2 = chains[i].f[k]->strand;
        w = chains[i].f[k]->start;
        x = arr->suftab[bd[i][beststart[i]][k]];
        y = x + len2;
        chr2 = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[x]);

        // read:           h)--(w
        // reference:   [u,v]  [x,y]
        // if the distance of h and w equals the distance of v and x we assume
        // that the fragments need to be merged
        if(u < x && v < y && chr1 == chr2 && strand1 == strand2 && strand1 == 0) {
          if(x < v || (x-v <= 20 && w-h <= 20) || (w >= h &&  x-v < w-h+20 && x-v+20 > w-h)) { 
            //merge
#ifdef DEBUGTRANSALIGN
            fprintf(stdout, "merging fragment %d with %d. [%d,%d;%d,%d] with [%d,%d;%d,%d]. u=%u : %u=y, u-y=%u \t h=%u : %u=w, w-h=%u\n, pass:%d,%d\n",
                j, k,
                chains[i].f[j]->start, chains[i].f[j]->end,
                u, v,
                chains[i].f[k]->start, chains[i].f[k]->end,
                x, y,
                u, y,
                u-y,
                h, w,
                w-h,
                chains[i].f[j]->pass,
                chains[i].f[k]->pass
                );
#endif
            chains[i].f[j]->start = MIN(chains[i].f[j]->start, chains[i].f[k]->start);
            chains[i].f[j]->end = MAX(chains[i].f[j]->end, chains[i].f[k]->end);
            chains[i].f[k]->start =  chains[i].f[j]->start;
            chains[i].f[k]->end = chains[i].f[j]->end;
            //chains[i].f[k]->pass = chains[i].f[j]->pass;
            bd[i][beststart[i]][j] = (u < x) ? bd[i][beststart[i]][j] : bd[i][beststart[i]][k];
            bd[i][beststart[i]][k] = bd[i][beststart[i]][j];
            chains[i].f[j]->substart = arr->suftab[bd[i][beststart[i]][j]];
            chains[i].f[k]->substart = arr->suftab[bd[i][beststart[i]][j]];

          } else {
            //dont merge
          }
        }
        //   h)--(w
        //[x,y]  [u,v]
        if(x < u && y < v && chr1 == chr2 && strand1 == strand2 && strand1 == 1) {
          //if(u < y || u-y < ((w >= h) ? w-h+20 : 20)) {
          if(u < y || (u-y <= 20 && w-h <= 20) || (w >= h &&  u-y < w-h+20 && u-y+20 > w-h)) { 
            //merge       
#ifdef DEBUGTRANSALIGN
            fprintf(stdout, "merging fragment %d with %d. [%d,%d;%d,%d] with [%d,%d;%d,%d]. u=%u : %u=y, u-y=%u \t h=%u : %u=w, w-h=%u, pass:%d,%d\n",
                j, k,
                chains[i].f[j]->start, chains[i].f[j]->end,
                u, v,
                chains[i].f[k]->start, chains[i].f[k]->end,
                x, y,
                u, y,
                u-y,
                h, w,
                w-h,
                chains[i].f[j]->pass,
                chains[i].f[k]->pass
                );
#endif
            chains[i].f[j]->start = MIN(chains[i].f[j]->start, chains[i].f[k]->start);
            chains[i].f[j]->end = MAX(chains[i].f[j]->end, chains[i].f[k]->end);
            chains[i].f[k]->start =  chains[i].f[j]->start;
            chains[i].f[k]->end = chains[i].f[j]->end;
            //chains[i].f[k]->pass = chains[i].f[j]->pass;
            bd[i][beststart[i]][j] = (u < x) ? bd[i][beststart[i]][j] : bd[i][beststart[i]][k];
            bd[i][beststart[i]][k] = bd[i][beststart[i]][j];
            chains[i].f[j]->substart = arr->suftab[bd[i][beststart[i]][j]];
            chains[i].f[k]->substart = arr->suftab[bd[i][beststart[i]][j]];

          } else {
            //dont merge
#ifdef DEBUGTRANSALIGN
            fprintf(stdout, "\t not merging fragment %d with %d. [%d,%d;%d,%d] with [%d,%d;%d,%d]. u=%u : %u=y, u-y=%u \t h=%u : %u=w, h-w=%u, pass:%d,%d\n",
                j, k,
                chains[i].f[j]->start, chains[i].f[j]->end,
                u, v,
                chains[i].f[k]->start, chains[i].f[k]->end,
                x, y,
                u, y,
                u-y,
                h, w,
                h-w,
                chains[i].f[j]->pass,
                chains[i].f[k]->pass
                );
#endif
          }
        }

        }
      }
      }

      newchains = ALLOCMEMORY(space, NULL, branchChain_t, noofchains);
      initChains (newchains, noofchains);

      for(i=0; i < noofchains; i++) {
        newchains[i].nooffragments = 0;
        newchains[i].f = NULL;
    
        for(j=0; j < chains[i].nooffragments; j++) {
          nochain = 0;
          len1 = chains[i].f[j]->end - chains[i].f[j]->start;
          for (k=j+1; k < chains[i].nooffragments; k++) {
            len2 = chains[i].f[k]->end - chains[i].f[k]->start;
            if(chains[i].f[j]->start  >= chains[i].f[k]->start  && 
                chains[i].f[j]->end    <= chains[i].f[k]->end    && 
                chains[i].f[j]->strand == chains[i].f[k]->strand &&
                chains[i].f[j]->substart >= chains[i].f[k]->substart &&
                chains[i].f[j]->substart+len1 <= chains[i].f[k]->substart+len2) {
              nochain = 1;
            } 
          }

          if(nochain) {
#ifdef DEBUGTRANSALIGN
            fprintf(stdout, "nochain at chain i:%d, f:%d\n", i, j);
#endif

          } else {
            subidx = getMultiCharSeqIndex(seq, &seq->sequences[chains[i].f[j]->substart]);
            getMultiCharSeqIdxBounds(seq, subidx, &sub_start, &sub_end);

#ifdef DEBUGTRANSALIGN
            fprintf(stdout, "for chain i:%d, f:%d: adding element [%d,%d]-%d to newchain -> %d (%u) (pass:%d)\n", i, j, chains[i].f[j]->start, chains[i].f[j]->end, chains[i].f[j]->score, chains[i].f[j]->substart-sub_start, chains[i].f[j]->substart, chains[i].f[j]->pass);
#endif
            chains[i].f[j]->subidx = subidx;

            if(newchains[i].f == NULL) {
              newchains[i].f = ALLOCMEMORY(NULL, NULL, branchfragment_t*, 1);

              branchfragment_t *copy = ALLOCMEMORY(NULL, NULL, branchfragment_t, 1);
              memmove(copy, chains[i].f[j], sizeof(branchfragment_t));
              newchains[i].f[0] = copy;
              newchains[i].nooffragments = 1;
              newchains[i].start = chains[i].f[j]->start; 
              newchains[i].end = chains[i].f[j]->end; 
              newchains[i].score = chains[i].f[j]->score;
            } else { 

              branchfragment_t *copy = ALLOCMEMORY(NULL, NULL, branchfragment_t, 1);
              memmove(copy, chains[i].f[j], sizeof(branchfragment_t));
              chain(&newchains[i], copy); 
            }

#ifdef DEBUGTRANSALIGN
            fprintf(stdout, "score:%d\n", newchains[i].score);
#endif
          }
        }
      }

      for(k=0; k < noofchains; k++) { 
        l = cnt[k];
        for (i=0; i < l; i++) {
          FREEMEMORY(space, bd[k][i]);
        }
        FREEMEMORY(space, bd[k]);
      }

      FREEMEMORY(space, bd);
      FREEMEMORY(space, cnt);
      FREEMEMORY(space, beststart);

      return newchains;
    }


/*----------------------------- appendFragments ------------------------------
 *    
 * @brief merge append list b to list a
 * @author Steve Hoffmann 
 *   
 */

  branchfragment_t *
appendFragments (branchfragment_t *a, Uint m, branchfragment_t *b, Uint n)
{
  Uint i;

  a = ALLOCMEMORY(space, a, branchfragment_t, m+n);

  for(i=0; i < n; i++) { 

    a[m+i].start = b[i].start;
    a[m+i].end = b[i].end;
    a[m+i].substart = b[i].substart;
    a[m+i].strand = b[i].strand;
    a[m+i].branchno = b[i].branchno;
    a[m+i].branch = b[i].branch;
    a[m+i].score = b[i].score;
    a[m+i].x = b[i].x;
    a[m+i].evalue = b[i].evalue;
    a[m+i].pass = b[i].pass;
    a[m+i].subidx = b[i].subidx;
  }
  
  return a;
}

/*----------------------------- filterFragments ------------------------------
 *    
 * @brief filter the fragments with respect to entropy, Evalue and maxocc
 * @author Steve Hoffmann 
 *   
 */

branchfragment_t *
filterFragments (void *space, Suffixarray *arr, matchstem_t **stems, char** seqs, 
    Uint len, karlin_t* stats, Uint maxocc, double maxevalue, double minentropy, Uint* nooffrags)
{

  Uint i, u, start, end, substart, k=0, l, r, s, x;
  double E, H;
  branch_t *branch;
  branchfragment_t *f = NULL;

  for (i = 0; i < len; i++) {
    for (u = 0; u < 2; u++) {
      x = (u == 0) ? i : len-1-i;

      for (s = 0; s < stems[u][x].noofbranches; s++) {     
        branch = &stems[u][x].branches[s];
        //        for (v=branch->l; v <= branch->r; v++) { 
        l = branch->l;
        r = branch->r;

        if (u == 0) {
          start = i;
          //CHANGED: end position one position too far
          //before: end = i + branch->mat;
          end = i + branch->mat - 1;
          substart = arr->suftab[l]; //l to v
        } else {
          start = i - branch->mat + 1;
          end = i;
          substart = arr->suftab[l] - branch->mat + 1; //l to v
        }

        E = kd_getBranchEvalue(stems[u], x, s, len, arr->numofsuffixes, stats);
        H = minshannonentropy(&seqs[0][start], end-start+1);

        Uint sub_idx, sub_start =0, sub_end=0;
        sub_idx = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[l]]);
        getMultiCharSeqIdxBounds(arr->seq, sub_idx, &sub_start, &sub_end);

#ifdef DEBUGTRANSALIGN       
        fprintf(stdout, "%d-[%d,%d] -> chr:%d-%d\tx:%d (strand:%d)\t", k, start, end, 
            getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[l]]),arr->suftab[l]-sub_start
            , x, u);
        fprintf(stdout, "Evalue: %f (max: %f), x:%u, s:%u, len:%u, H: %f (min %f), occ:%d, scr:%d\n", E, maxevalue, x, s, len, H, minentropy, r-l, kdscore(branch));
#endif
        if(E < maxevalue && H >= minentropy && l <= r && r - l < maxocc) {

#ifdef DEBUGTRANSALIGN          
          fprintf(stdout, "adding %d [%d,%d]\n", k, start, end);
#endif 
          k++;
          f = ALLOCMEMORY(space, f, branchfragment_t, k);
          f[k-1].start = start ;
          f[k-1].end = end;
          f[k-1].substart = substart;
          f[k-1].subidx = sub_idx;
          f[k-1].strand = (unsigned char) u;
          f[k-1].branchno = s;
          f[k-1].branch = branch;
          f[k-1].score = kdscore(branch);
          f[k-1].x = x;
          f[k-1].evalue = E;
          f[k-1].pass = 0;
        }
      }
      }
    }
    *nooffrags = k;
    return f;
  }


/*-------------------------- removeOverlapFragments --------------------------
 *    
 * @brief remove overlapping fragments
 * @author Steve Hoffmann 
 *   
 */

#define INCLUSIONMARGIN 10

branchfragment_t*
removeOverlapFragments (branchfragment_t *f, Uint nooffrags, Uint maxovl, Uint *newnooffrags)
{
      
  Uint i, j, k=0;
  branchfragment_t *g = NULL;

  for(i=0; i < nooffrags; i++) {
    for(j=0; j < k; j++){  

      Uint x1 = g[j].start;
      Uint y1 = g[j].end;
      Uint x2 = f[i].start;
      Uint y2 = f[i].end;
      Uint u1 = g[j].substart;
      Uint u2 = f[i].substart;
      //to remove: set s1 = s2;
      Uint s1 = g[j].strand;
      Uint s2 = f[i].strand;
      Uint c1 = g[j].subidx;
      Uint c2 = f[i].subidx;

#ifdef DEBUGTRANSALIGN            
      fprintf(stdout, "[%d,%d]->[%u] vs. [%d,%d]->[%u]; duint:%u, d1:%d, d2:%d\n", x1, y1, u1, x2, y2, u2,dist_uint(u1,u2), x2-x1, x1-x2);
#endif
      //inclusion
      if(x1 <= x2 && y2 <= y1 && s1 == s2 && c1 == c2 && dist_uint(u1,u2) < x2-x1+INCLUSIONMARGIN) //&& dist_uint(u1,u2) < 200000) 
      { 
#ifdef DEBUGTRANSALIGN
        fprintf(stdout, "inclusion consumption (1)\n");
#endif
        break;
      }
      if(x2 <= x1 && y1 <= y2 && s1 == s2 && c1 == c2 && dist_uint(u1,u2) < x1-x2+INCLUSIONMARGIN) //&& dist_uint(u1,u2) < 200000) 
      { 
#ifdef DEBUGTRANSALIGN
        fprintf(stdout, "inclusion consumption (2)\n");
#endif
        break;
      }
          //overlap
      if(y1 > x2 && y2 > y1 && y1-x2 >= maxovl){ 
#ifdef DEBUGTRANSALIGN
        fprintf(stdout, "overlap consumption (1)\n");
#endif
        //break;
      }
      if(y2 > x1 && y1 > y2 && y2-x1 >= maxovl) { 
#ifdef DEBUGTRANSALIGN
        fprintf(stdout, "overlap consumption (2)\n");
#endif
        //break;
      }
    }

    if(j < k) {
      if (kdscore(f[i].branch) > g[j].score) { 
#ifdef DEBUGTRANSALIGN        
        fprintf(stdout, "replacing %d [%d,%d] by [%d,%d]\n", k, g[i].start, g[i].end, f[i].start, f[i].end);
#endif
        g[j].start = f[i].start;  
        g[j].end = f[i].end;
        g[j].substart = f[i].substart;
        g[j].strand = (unsigned char) f[i].strand;
        g[j].branchno = f[i].branchno;
        g[j].branch = f[i].branch;
        g[j].score = kdscore(f[i].branch);
        g[j].x = f[i].x;
        g[j].evalue = f[i].evalue;
      }
    } else { 

#ifdef DEBUGTRANSALIGN          
      fprintf(stdout, "adding %d [%d,%d]\n", k, f[i].start, f[i].end);
#endif 
      k++;
      g = ALLOCMEMORY(space, g, branchfragment_t, k);
      g[k-1].start = f[i].start ;
      g[k-1].end = f[i].end;
      g[k-1].substart = f[i].substart;
      g[k-1].subidx = f[i].subidx;
      g[k-1].strand = (unsigned char)f[i].strand;
      g[k-1].branchno = f[i].branchno;
      g[k-1].branch = f[i].branch;
      g[k-1].score = kdscore(f[i].branch);
      g[k-1].x = f[i].x;
      g[k-1].evalue = f[i].evalue;
      g[k-1].pass = f[i].pass;
    }
  }

  *newnooffrags = k;
  return g;
}


/*------------------------------- branchChain --------------------------------
 *    
 * @brief find chain of branches
 * @author Steve Hoffmann 
 *   
 */
 
branchChain_t *
branchChain(void *space, Suffixarray *arr, matchstem_t **stems, char **seqs, 
    Uint len, karlin_t *stats, Uint *noofchains, branchfragment_t **fragments, 
    Uint maxocc, double maxevalue, double minentropy) {
  
  Uint i, j, k=0, l=0, c_prime=0, a=0; //q_prime = 0, c=0, q, v;
  int maxovl = 12, bestscr; //maxgap = 22;
  branchChain_t *chains = NULL, *extra = NULL;
  branchfragment_t *f = NULL, *g = NULL, *h = NULL;// *f_prime = NULL; 
  branchfragment_t *copy;
  
  g = filterFragments (space, arr, stems, seqs, len, stats, maxocc, maxevalue, minentropy, &l);
  f = removeOverlapFragments (g, l, maxovl, &k);
  FREEMEMORY(space, g);
  FREEMEMORY(space, h);

  qsort(f, k, sizeof(branchfragment_t), cmp_branchfragments);
  chains = ALLOCMEMORY(space, chains, branchChain_t, k);
  initChains(chains, k);

  for(i = 0; i < k; i++) {

    c_prime = i;
    bestscr = 0;

    // search for best precedessor chain
    for (j = 0; j < i; j++){

      // only allow short overlap
      if (chainovl(&chains[j], &f[i]) < maxovl || 
          (dist_uint(chains[j].f[0]->substart, f[i].substart) < 1000 && //allow local template switches explicitly
           chains[j].nooffragments == 1 && 
           chains[j].f[0]->strand != f[i].strand)) {

        //fprintf(stdout, "attempt to chain %d: [%d,%d] with [%d,%d]i - ovl:%lld, dist:%u\n", j, chains[j].start, chains[j].end, f[i].start, f[i].end, chainovl(&chains[j], &f[i]), dist_uint(chains[j].end, f[i].start));
        // update best precessor
        if (bestscr < chainscore(&chains[j], &f[i])){
          bestscr = chainscore(&chains[j], &f[i]);
          c_prime = j;
        }
        if (bestscr == chainscore(&chains[j], &f[i]) && (c_prime == i ||
            minDistFragmentHits(arr, chains[j].f[chains[j].nooffragments-1], &f[i]) <
            minDistFragmentHits(arr, chains[c_prime].f[chains[c_prime].nooffragments-1], &f[i]))) {          
          bestscr = chainscore(&chains[j], &f[i]);
          c_prime = j;
 
          Uint sub_start, sub_end;
          Uint subidx = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[f[i].branch->l]]);
          getMultiCharSeqIdxBounds(arr->seq, subidx, &sub_start, &sub_end);

        }
      }
    }

    if (c_prime != i){
      // TODO: add function to do the following
      chains[i].nooffragments = chains[c_prime].nooffragments + 1;
      chains[i].f = ALLOCMEMORY(space, NULL, branchfragment_t*, chains[i].nooffragments);

      for(Uint w=0; w < chains[c_prime].nooffragments; w++) {
        copy = ALLOCMEMORY(NULL, NULL, branchfragment_t, 1);
        memmove(copy, chains[c_prime].f[w], sizeof(branchfragment_t));
        chains[i].f[w] = copy;
      }
      
      copy = ALLOCMEMORY(NULL, NULL, branchfragment_t, 1);
      memmove(copy, &f[i], sizeof(branchfragment_t));

      chains[i].f[chains[i].nooffragments-1] = copy;
      chains[i].score = bestscr;
      chains[i].end = f[i].end;
      chains[i].start = chains[c_prime].start;
    }
    else {
      chains[i].nooffragments = 1;
      chains[i].f = 
        ALLOCMEMORY(space, NULL, branchfragment_t*, chains[i].nooffragments);

      copy = ALLOCMEMORY(NULL, NULL, branchfragment_t, 1);
      memmove(copy, &f[i], sizeof(branchfragment_t));

      chains[i].f[0] = copy;
      chains[i].score = f[i].score;
      chains[i].start = f[i].start;
      chains[i].end = f[i].end;
    }

    if(c_prime != i && chains[i].nooffragments > 1) { 

       Uint u = chains[i].nooffragments-1;
       Uint v = chains[i].nooffragments-2;
       Uint frag1 = chains[i].f[u]->branch->l;
       Uint frag2 = chains[i].f[v]->branch->l;
       Uint idx1 = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[frag1]]);
       Uint idx2 = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[frag2]]);
      
       if( 
          (chains[i].f[u]->strand == 0 && chains[i].f[u]->substart < chains[i].f[v]->substart) 
          || 
          (chains[i].f[u]->strand == 1 && chains[i].f[u]->substart > chains[i].f[v]->substart) 
          ||
          (chains[i].f[u]->strand != chains[i].f[v]->strand)
          ||
          (chains[i].f[u]->strand != chains[i].f[v]->strand)
          ||
          (idx1 != idx2)
        )
     
    {
       
      extra = ALLOCMEMORY(space, extra, branchChain_t, a+1);
      extra[a].start = f[i].start;
      extra[a].end = f[i].end;
      extra[a].f = ALLOCMEMORY(space, NULL, branchfragment_t*, 1);
      extra[a].score= f[i].score;
      extra[a].nooffragments=1;
      
      copy = ALLOCMEMORY(NULL, NULL, branchfragment_t, 1);
      memmove(copy, &f[i], sizeof(branchfragment_t));

      extra[a].f[0] = copy;
      a++;
    }
    }
  }
  

   if(a > 0) { 
   chains = ALLOCMEMORY(space, chains, branchChain_t, a+k);
   memmove(&chains[k], extra, sizeof(branchChain_t)*a);
   k += a;
   }

   FREEMEMORY(space, extra);

  (*noofchains) = k;
  (*fragments) = f;

  return chains;
}



/*-------------------------------- showChains --------------------------------
 *    
 * @brief dump the chains
 * @author Steve Hoffmann 
 *   
 */
 
void
showChains(branchChain_t *chains, Uint noofchains, Suffixarray *arr, 
    FILE *dev, char *seq, Uint len) {
  
  Uint i, j, q, subidx, sub_start, sub_end;
  double H;

  for(i=0; i < noofchains; i++) {
    fprintf(dev, "chain %d: %d-%d (%d)\n", i, chains[i].start, 
        chains[i].end, chains[i].score);
    
    for(j=0; j < chains[i].nooffragments; j++) {
     
      int score = chains[i].f[j]->score;
      int ovl = 0;
      int usedovl = 0;

      if(j>0) { 
        ovl = (chains[i].f[j-1]->end - chains[i].f[j]->start)+1;
        usedovl = (ovl < 0) ? 0 : ovl;
        score = (chains[i].f[j-1]->score + chains[i].f[j]->score) - usedovl;
      }


      fprintf(dev, "fragment %d: %d-%d (%d) (%d:%f); ovl: (%d,%d), cscore:%d; substart:", j, 
          chains[i].f[j]->start, chains[i].f[j]->end, 
          chains[i].f[j]->strand, chains[i].f[j]->score, 
          chains[i].f[j]->evalue, ovl, usedovl, score);
      
      for(q=chains[i].f[j]->branch->l; q <= chains[i].f[j]->branch->r; q++) {
        subidx = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[q]]);
        getMultiCharSeqIdxBounds(arr->seq, subidx, &sub_start, &sub_end);

        fprintf(dev,"%u (chr:%d) -> %u, ",arr->suftab[q], 
            getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[q]]),arr->suftab[q]-sub_start);
      }
      //CHANGED: from end-start to end-start+1 for length
      //before: H = shannonentropy(NULL, &seq[chains[i].f[j]->start], chains[i].f[j]->end - chains[i].f[j]->start, asize, tab);
      H = minshannonentropy(&seq[chains[i].f[j]->start], chains[i].f[j]->end - chains[i].f[j]->start + 1);

      fprintf(dev, "entropy: %f\n", H);
      fprintf(dev, "substart selected: %u\n", chains[i].f[j]->substart);
    }
    fprintf(dev, "\n");
  }
}
