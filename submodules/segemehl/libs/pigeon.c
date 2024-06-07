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
 *  suffixmatch.c
 *  suffix matches
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 06.09.2016 13:50:21 CEST
 *  
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "segemehl.h"
#include "memory.h"
#include "fileio.h"
#include "stringutils.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "sufarray.h"
#include "mmchar.h"
#include "mathematics.h"
#include "biofiles.h"
#include "vtprogressbar.h"
#include "sort.h"
#include "bitArray.h"
#include "vqueue.h"
#include "vstack.h"
#include "container.h"
#include <pthread.h>
#include "kdseed.h"
#include "info.h"
#include "debug.h"
#include "mapfrag.h"
#include "alignment.h"
#include <assert.h>
#include "iupac.h"
#include "bitvectoralg.h"
#include "alignment.h"
#include "matealign.h"
#include "mapfrag.h"
#include "splitalign.h"
#include "segemehl_helper.h"

#define PIGEONHOLESIZE 25
#define MAXSEEDS 25


mappingset_t*
pigeon(void *space, mappingset_t *set, Suffixarray *s, MultiCharSeq *mseq, 
    seseq_t *query, Uint nsplits, char *qname,
    Uint *enctab, bitvector* D, char ismate, Uint nchars, segemehl_t *nfo) {
  int indel = -2;
  int transition = -10;
  int scores[]={1, -2};

  unsigned char trans;
  Uint totalcover=0;
  int totalscore=0;
  Uint m;
  m = query->len;

  MultiCharSeqAlignment *a;
#ifndef SEEDHASH
  PairUint *full[2];
  Uint fulllast[2];
#endif
  PairUint *hit[2];
  bitvector *peq[2];
  Uint rc, i, j, k; 
  int ret;
  Uint splitsize, rest;
  Uint pos, off, ext;
  Uint cursplitpos=0;
  MultiCharSeqAlignment mcsa, *copy;
  mappingset_t *tmp;
  //PairUint test;
  char fullmatch = 0;
  Uint *check=NULL;
  Uint checklen=0;
  char *seqs[2];
  char *quals[2];

  assert(nfo->hashsize <= PIGEONHOLESIZE);
  se_getData(space, query, seqs, quals, nfo->bisulfite, 1);


  //divide read into n splits
  splitsize = m/nsplits;
  //prepare array for storing the n intervals on 2 strands
  hit[0] = ALLOCMEMORY(space, NULL, PairUint, nsplits);
  hit[1] = ALLOCMEMORY(space, NULL, PairUint, nsplits);
  for(i=0; i < nsplits; i++) {
    hit[0][i].a = 1;
    hit[0][i].b = 0;
    hit[1][i].a = 1;
    hit[1][i].b = 0;
  }
  //prepare eq vector
  peq[0] = getpeq(NULL, seqs[0], m, mseq->map, mseq->mapsize, enctab);
  peq[1] = getpeq(NULL, seqs[1], m, mseq->map, mseq->mapsize, enctab);

  //first attempt full path match for full matches
#ifndef SEEDHASH
  full[0] = searchSuffixPath (space, s, seqs[0], m, &fulllast[0]);
  full[1] = searchSuffixPath (space, s, seqs[1], m, &fulllast[1]);
#else
  hit[0][0] = searchSuffixArrayHash(space, s, nfo->hash, nfo->hashsize, seqs[0], m);
  hit[1][0] = searchSuffixArrayHash(space, s, nfo->hash, nfo->hashsize, seqs[1], m);
#endif

  if(!nchars) { 
    for(rc=0; rc < 2; rc++) { 

#ifdef SEEDHASH
      if(hit[rc][0].a <= hit[rc][0].b && 
          hit[rc][0].b - hit[rc][0].a <= MAXSEEDS) {
#else
        if(full[rc][m-1].a <= full[rc][m-1].b && 
            full[rc][m-1].b - full[rc][m-1].a <= MAXSEEDS) {
#endif
      
        se_getData(space, query, seqs, quals, nfo->bisulfite, 0);

#ifdef SEEDHASH
          for(j=hit[rc][0].a; j <= hit[rc][0].b; j++){
#else
            for(j=full[rc][m-1].a; j <= full[rc][m-1].b; j++){
#endif
              pos = s->suftab[j];
              //add a full match to mappingset
              initMultiCharSeqAlignment(space, &mcsa, mseq, pos, 0, m, 
                  rc, qname, seqs[rc], quals[rc], m);

              insertMeop(mcsa.al, Match, m);

              if(nfo->bisulfite) {
                //required to enforce C->T/G->A and exclude T->C/A->G
                reevalMultiCharSeqAlignment(&mcsa);
              }

              copy = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
              copy = bl_copyMCSA(copy, &mcsa);
              set->elem = ALLOCMEMORY(NULL, set->elem, mapping_t, set->n+1);
              bl_initMapping(&set->elem[set->n], seqs[rc], quals[rc], 0, 0); 
              bl_addMapFrag(&set->elem[set->n], copy, NULL, ismate, 0);
              set->n += 1;

              wrapMultiCharSeqAlignment(NULL, &mcsa);
              fullmatch = 1;
            }
          }
        }
      }

  if(!fullmatch){ 

    
    for(i=0; i < nsplits; i++) {
      //substring to search has to have a minimum length
      if((rest = m - (i*splitsize)) >= PIGEONHOLESIZE) { 
#ifndef SEEDHASH
        if(i > 0) { 
#endif
          //note that first splitpos is zero
          cursplitpos = i*splitsize;
          //attempt match on forward strand
#ifndef SEEDHASH
          hit[0][i] = searchSuffix(space, s, &seqs[0][cursplitpos], PIGEONHOLESIZE);
#else
          hit[0][i] = searchSuffixArrayHash(space, s, nfo->hash, nfo->hashsize, &seqs[0][cursplitpos], PIGEONHOLESIZE);
#endif
          //fprintf(stdout, "hit[0][%d] = (%d,%d), test=(%d,%d)\n", i, hit[0][i].a, hit[0][i].b, test.a, test.b);
          //assert(hit[0][i].a == test.a && hit[0][i].b == test.b);
          //attempt match on backward strand
#ifndef SEEDHASH
          hit[1][i] = searchSuffix(space, s, &seqs[1][cursplitpos], PIGEONHOLESIZE);
#else
          hit[1][i] = searchSuffixArrayHash(space, s, nfo->hash, nfo->hashsize, &seqs[1][cursplitpos], PIGEONHOLESIZE);
#endif
          //fprintf(stdout, "hit[1][%d] = (%d,%d), test=(%d,%d)\n", i, hit[1][i].a, hit[1][i].b, test.a, test.b);
          //assert(hit[1][i].a == test.a && hit[1][i].b == test.b);
#ifndef SEEDHASH
        } else {
          hit[0][i] = full[0][PIGEONHOLESIZE-1];
          hit[1][i] = full[1][PIGEONHOLESIZE-1];
        }
#endif     
    
        se_getData(space, query, seqs, quals, nfo->bisulfite, 0);
        
        for(rc=0; rc <= 1; rc++) { 
          //seeds found, but not too many
          if(hit[rc][i].a <= hit[rc][i].b && hit[rc][i].b-hit[rc][i].a <= MAXSEEDS){
            //do alignments
            for(j=hit[rc][i].a; j <= hit[rc][i].b; j++ ) { 
              pos = s->suftab[j];
              off = cursplitpos + nsplits;
              ext = m + nsplits + 10;
              //inititalize alignment
              initMultiCharSeqAlignment(space, &mcsa, mseq, pos, off, 
                  ext, rc, qname, seqs[rc], quals[rc], m);

              //make sure that this alignment has not yet
              //been evaluated
              for(k = 0; k < checklen; k++)  {
                if (check[k] >= mcsa.refstart-10 && check[k] <= mcsa.refstart+10) {     
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

              //perform alignment
              ret = bl_myersMCSA(&mcsa, mseq, enctab, peq[rc], D, nsplits-1);

              if(ret) {
                copy = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
                copy = bl_copyMCSA(copy, &mcsa);
                set->elem = ALLOCMEMORY(NULL, set->elem, mapping_t, set->n+1);
                bl_initMapping(&set->elem[set->n], seqs[rc], quals[rc], 0, 0); 
                bl_addMapFrag(&set->elem[set->n], copy, NULL, ismate, 0);
                set->n += 1;
              }             

              wrapMultiCharSeqAlignment(NULL, &mcsa);
            }
          } 
        }
      } else {
        hit[0][i].a = 1;
        hit[0][i].b = 0;
        hit[1][i].a = 1;
        hit[1][i].b = 0;
      }   
    }
  } 


  if(!bl_hasQueryMapping(set) && nfo->split) { 
    for(rc=0; rc < 2; rc++) {
      //if the first and last seed are unique
      if(hit[rc][0].a == hit[rc][0].b 
          && hit[rc][nsplits-1].a == hit[rc][nsplits-1].b) {
        //if the distance is within a smaller range
        Uint l = s->suftab[hit[rc][0].a];
        Uint r = s->suftab[hit[rc][nsplits-1].b];
        if(DIST(l,r) > m && DIST(l,r) < 20000) {

          Uint *mystart = ALLOCMEMORY(NULL, NULL, Uint, 2);
          Uint *myends = ALLOCMEMORY(NULL, NULL, Uint, 2);
          uint64_t *mypos = ALLOCMEMORY(NULL, NULL, uint64_t, 2);
          char *myrc = ALLOCMEMORY(NULL, NULL, char, 2);

          if(rc) {
            mystart[0] = m-PIGEONHOLESIZE;
            myends[0] = m-1;
            mystart[1] = m - ((nsplits-1) *splitsize);
            myends[1] = mystart[1] + PIGEONHOLESIZE -1; 
            mypos[0] = r;
            mypos[1] = l;
          } else { 
            mystart[0] = 0;
            myends[0] = PIGEONHOLESIZE;
            mystart[1] = (nsplits-1) *splitsize;
            myends[1] = mystart[1] + PIGEONHOLESIZE -1; 
            mypos[0] = l;
            mypos[1] = r;
          }
          myrc[0] = rc;
          myrc[1] = rc;

          a= se_AlignSplitMap (mystart, myends, mypos, myrc, 2,  mseq, qname, seqs, quals, m, 
              scores, indel, transition);

          tmp = ALLOCMEMORY(NULL, NULL, mappingset_t, 1);
          bl_initMappingSet(tmp);

          se_kdAlignEvalSplitAlign (mseq, a,  tmp, &totalcover, &totalscore, 
              &trans, scores, indel, 2, ismate, seqs, quals, nfo);

          //fprintf(stdout, "totalcover %d, totalscore %d, m:%d, nsplits:%d, scores[0]:%d, scores[1]:%d\n", totalcover, totalscore, m, nsplits, scores[0], scores[1]);
          //bl_dumpMappingSet(stdout, tmp);

          if(totalcover >= m-nsplits && totalscore >= (m*scores[0])-abs(nsplits*scores[1])) {
            bl_removeMappingSet(set);
            FREEMEMORY(NULL, set);
            set = tmp;
            //fprintf(stdout, "accepted.\n");
          } else {
            bl_removeMappingSet(tmp);
            FREEMEMORY(NULL, tmp);
          }

          FREEMEMORY(NULL , mystart);
          FREEMEMORY(NULL , myends);
          FREEMEMORY(NULL , mypos);
          FREEMEMORY(NULL , myrc);

          for(i=0; i < 2; i++) { 
            wrapMultiCharSeqAlignment(NULL, &a[i]);
          }

          FREEMEMORY(NULL, a);


        }
      }
    }
  }

  for(j=0; j < mseq->mapsize; j++) {
    FREEMEMORY(space, peq[0][j]);
    FREEMEMORY(space, peq[1][j]);
  }

  FREEMEMORY(space, peq[0]);
  FREEMEMORY(space, peq[1]);
  FREEMEMORY(space, hit[0]);
  FREEMEMORY(space, hit[1]);
#ifndef SEEDHASH
  FREEMEMORY(space, full[0]);
  FREEMEMORY(space, full[1]);
#endif
  if(check) FREEMEMORY(space, check);

  return set;
}


