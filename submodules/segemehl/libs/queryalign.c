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
 *  queryalign.c
 *  
 *  this is the replacement for kdmatch. It uses the mapfrag and matealign
 *  routines for alignment and than passes it to the replacment of the
 *  output manager which generically works on mapping_t and mappingset_t
 *  Usually the output manager uses sam.c to realise the output
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 12/19/2012 12:55:11 CET
 *  
 */

#include "debug.h"
#include "mapfrag.h"
#include "kdseed.h"
#include "locus.h"
#include "sufarray.h"
#include "string.h"
#include <limits.h>
#include <inttypes.h>
#include "mathematics.h"
#include "segemehl.h"
#include "manout.h"
#include "bitvectoralg.h"
#include "matealign.h"
#include "alignment.h"

/*----------------------------- bl_getGoodSeeds ------------------------------
 *    
 * @brief filter the branches and get the seeds
 * @author Steve Hoffmann 
 *   
 */

mapseedlist_t*
bl_getGoodSeeds (matchstem_t **stems, unsigned m, unsigned n, karlin_t *stats, 
    segemehl_t *nfo)
{

  unsigned int u, i, q;
  branch_t *b;
  double SPM;
  mapseedlist_t *l;

 
  l = ALLOCMEMORY(NULL, NULL, mapseedlist_t, 1);
  bl_initMapSeedList(l);

  SPM = spacemult(m, n, stats->H, stats->K);

  //get a list of suitable seeds first
  for(u = 0; u < 2; u++) {
    for(i = 0; i < m; i++) {
      for(q = 0; q < stems[u][i].noofbranches; q++) {
        b =  &stems[u][i].branches[q];
        l = bl_addMapSeedBranch(l, b, i, u, nfo->maxevalue, nfo->M, SPM, stats);
      }
    }
  }

  return l;
}


/*-------------------------- bl_seedAlign ---------------------------
 *    
 * @brief align the seeds (from getGoodSeeds) in the matchstem using bv 
 * algorithm
 * @author Steve Hoffmann 
 *   
 */

mappingset_t*
bl_seedAlign(Suffixarray *s, mappingset_t* set, MultiCharSeq *mseq, 
    char **seqs, char **qual, Uint len, char* qname, mapseedlist_t *l, 
    segemehl_t *nfo, Uint *enctab, bitvector* D, char ismate) {

  void *space = NULL;
  Uint pos, i, j, k, rc; 
  int maxedist, mrgn=0, ret;
  uint64_t seedstart;
  Uint *check=NULL;
  Uint checklen=0;
  Uint lclip = 0;
  Uint rclip = 0;
  bitvector *peq[2];
  MultiCharSeqAlignment mcsa, *copy;

  //if(l->n > nfo->maxmyers) return set;
  maxedist = getMaximumAlignmentEdist(len, nfo->accuracy);
 
  //get the extended length of the sequence (includes edist on both sides)
  Uint ext = getExtendedAlignmentLength(len, nfo->accuracy);

  mrgn = 40*((double)maxedist/100.);

  peq[0] = getpeq(NULL, seqs[0], len, mseq->map, mseq->mapsize, enctab);
  peq[1] = getpeq(NULL, seqs[1], len, mseq->map, mseq->mapsize, enctab);
 
  //now go and align the mapseedlist l
  for(i=0; i < l->n; i++) {
    if(l->l[i].good) { 
    //all interval on strand u
    rc = l->l[i].rc;
    seedstart = l->l[i].u;

    for(j = l->l[i].l; j <= l->l[i].r; j++) {
      pos = s->suftab[j];
      
      //get the ustart of the alignment on plus strand
      Uint beg = seedstart;
      //get the offset for each side and add the position of the seed in read
      Uint off = floor(((double)ext-len)/2.0) + beg;

      //init the multicharseq alignment
      if(l->l[i].mat != len || l->l[i].edist != 0) {
        initMultiCharSeqAlignment(space, &mcsa, mseq, pos, off, 
            ext, rc, qname, seqs[rc], qual[rc], len);
      } else {
        initMultiCharSeqAlignment(space, &mcsa, mseq, pos, 0, len, 
            rc, qname, seqs[rc], qual[rc], len);
      }
      //check if a similar alignment has been seen already
      for(k = 0; k < checklen; k++)  {
        if (check[k] >= mcsa.refstart-mrgn && check[k] <= mcsa.refstart+mrgn) {     
          break;	
        }
      }
      if (k < checklen) {
        wrapMultiCharSeqAlignment(space, &mcsa);
        continue;
      } else {
        check = ALLOCMEMORY(space, check, Uint, checklen+1);
        check[checklen++]= mcsa.refstart;
      }
      //fresh align de novo or save time if the seed is full length
      if(l->l[i].mat != len || l->l[i].edist != 0) {

#ifdef QUERYALIGNDEBUG
        DBG("call myers: seed=%"PRIu64", pos=%"PRIu64", loff=%d, maxedist=%d, ext=%d, len=%d, rc=%d\n",
            seedstart, pos, off, maxedist, ext, len, rc);
#endif

#ifndef LOCAL
        ret = bl_myersMCSA(&mcsa, mseq, enctab, peq[rc], D, maxedist);
#else
        int minscore = MAX(nfo->localparams[0], ((nfo->localparams[1]*len)/100)*nfo->scores[0]);
        ret = bl_localMCSA(&mcsa, mseq, nfo->scores, nfo->scores[2], minscore);
#endif

#ifdef QUERYALIGNDEBUG
        DBG("call myers: returned %d\n", ret);
#endif

      } else {

#ifdef QUERYALIGNDEBUG
        DBG("no call myers: full match of length %d\n", len);
#endif

        insertMeop(mcsa.al, Match, len);
             
        if(nfo->bisulfite) {
           //required to enforce C->T/G->A and exclude T->C/A->G
           reevalMultiCharSeqAlignment(&mcsa);
        }
       
        ret = 1;
      }

      //add if suitable 
      //change to orig: if mappingsetHasPos, comp edist, replace
      if(ret && !bl_mappingsetHasPos(set, mcsa.refstart+mcsa.al->voff)){
        //fprintf(stdout, "adding alignment: %d\n", ret);
        copy = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
        copy = bl_copyMCSA(copy, &mcsa);
        set->elem = ALLOCMEMORY(NULL, set->elem, mapping_t, set->n+1);
        bl_initMapping(&set->elem[set->n], seqs[rc], qual[rc], lclip, rclip); 
        bl_addMapFrag(&set->elem[set->n], copy, &l->l[i], ismate, 0);
        set->n += 1;
       }
      //clean up  
      wrapMultiCharSeqAlignment(NULL, &mcsa);
    }
    }
  }
  
  for(j=0; j < mseq->mapsize; j++) {
    FREEMEMORY(space, peq[0][j]);
    FREEMEMORY(space, peq[1][j]);
  }

  FREEMEMORY(space, peq[0]);
  FREEMEMORY(space, peq[1]);
  if(check) {
    FREEMEMORY(space, check);
  }

  return set;
}


