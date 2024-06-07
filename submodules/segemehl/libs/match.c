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
 *  match.c
 *  replacement for kdmatch
 *  it gets the seeds from kdseed.c,
 *  it uses queryalign, matealign and
 *  split align. the last three return a mapping_t
 *  data structure
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 12/21/2012 00:25:18 CET
 *  
 */

#include "match.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <pthread.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/times.h>
#include <time.h>

#include "basic-types.h"
#include "segemehl.h"
#include "mapfrag.h"
#include "bitVector.h"
#include "segemehl_helper.h"
#include "bitvectoralg.h"
#include "queryalign.h"
#include "kdseed.h"
#include "matealign.h"
#include "pigeon.h"
#include "splitalign.h"
#include "multicharseq.h"	//(skipped)
#include "biofiles.h"	//(skipped)
#include "iupac.h"	//(skipped)
#include "vtprogressbar.h"	//(skipped)
#include "manout.h"	//(skipped)
#include "memory.h"	//(skipped)
#include "samio.h"	//(skipped)

/*--------------------------- se_updateProgressBar ---------------------------
 *    
 * @brief keeping the user informed ... somehow
 * @author Steve Hoffmann 
 *   
 */

void
se_updateProgressBar(Uint k, segemehl_t *nfo) {

  if (!nfo->mute) {
    if (nfo->counter == NULL) {
      progressBarVT("reads matched.", nfo->reads->noofseqs, k, 25);
    } else {
      (*nfo->counter)++;
    }
  }
  return;
}

/*------------------------------ bl_fastaMaster ------------------------------
 *    
 * @brief a server of fasta chunks
 * @author Steve Hoffmann 
 *   
 */

  fasta_t*
se_fastaMaster (void *space, fasta_t *f, Uint size, segemehl_t *nfo)
{

  Uint k;
  fasta_t *piece = NULL;
  merge_t *m = nfo->merge;
  uint64_t did = (nfo->bisulfite + 1) % 2;

  //return full set if no threads are used
  if(nfo->threadno == 1) {
    if(nfo->nextfastaidx[0] ==0) {
      if(m->mapsz == 0) {
        m->map = ALLOCMEMORY(NULL, m->map, mergemap_t, 1);
        m->map[0].beg = 0;
        m->map[0].end = f->noofseqs-1;
        m->map[0].chk = 0;
        m->map[0].n = 0;
        m->map[0].binids = NULL; 
        m->map[0].queuepos = NULL; 
        m->map[0].processed = 0;
        m->mapsz=1;
      }
  
      //allocate space for a new bin id
      m->map[0].binids = ALLOCMEMORY(NULL, m->map[0].binids, Uint, m->map[0].n+1);
      m->map[0].queuepos = ALLOCMEMORY(NULL, m->map[0].queuepos, Uint, m->map[0].n+1);
      //the following calculation is the bin position in a consecutive
      //list of all bins accross all domains
      m->map[0].binids[m->map[0].n] = did;
      m->map[0].queuepos[m->map[0].n] = 0; 
      m->map[0].n += 1;

      nfo->nextfastaidx[0] = 1;
      return f;
    } else { 
      return NULL;
    }
  }

  //in case threads are used
  pthread_mutex_lock(nfo->mtx4);

  
  if(nfo->nextfastaidx[0] < f->chunkindex->size) {
    /*
     * get a copy of the fasta piece to be processed here
     *
     */
    piece = bl_fastxCopyIndex (space, f, nfo->nextfastaidx[0], size);
   /*
    * get map for merging
    *
    */
    uint64_t cumnoofseqs = f->chunkindex->ap[nfo->nextfastaidx[0]].cumnoofseqs;
    uint64_t noofseqs = f->chunkindex->ap[nfo->nextfastaidx[0]].noofseqs;
    uint64_t beg = cumnoofseqs-noofseqs;
    uint64_t end = cumnoofseqs-1;
    uint64_t chk = nfo->nextfastaidx[0];
    //domain id
    //thread id
    uint64_t tid = nfo->threadid;
    uint64_t nthreads = nfo->threadno;

    m = nfo->merge;

    for(k=0; k < m->mapsz; k++) {
      if(beg == m->map[k].beg && end == m->map[k].end)
        break;
    }
    /*
     * create a new map entry if the interval is seen
     * for the first time
     *
     */
    if(k == m->mapsz) {

      m->map = ALLOCMEMORY(NULL, m->map, mergemap_t, m->mapsz+1);
      m->map[k].beg = beg;
      m->map[k].end = end;
      m->map[k].chk = chk;
      m->map[k].n = 0;
      m->map[k].binids = NULL; 
      m->map[k].queuepos = NULL; 
      m->map[k].processed = 0;
      m->mapsz++;
    } 
     
    //allocate space for a new bin id
    m->map[k].binids = ALLOCMEMORY(NULL, m->map[k].binids, Uint, m->map[k].n+1);
    m->map[k].queuepos = ALLOCMEMORY(NULL, m->map[k].queuepos, Uint, m->map[k].n+1);
    //the following calculation is the bin position in a consecutive
    //list of all bins accross all domains
    m->map[k].binids[m->map[k].n] = (did*nthreads)+tid;
    m->map[k].queuepos[m->map[k].n] = 0; 
    m->map[k].n += 1;
    /*
     * update for next chunk reques update for next chunk requestt
     *
     */
    nfo->nextfastaidx[0] += piece->chunkindex->size;

  }   

  pthread_mutex_unlock(nfo->mtx4);

  return piece;
}


/*--------------------------------- se_clip ----------------------------------
 *    
 * @brief clipping sequences
 * @author Steve Hoffmann 
 *   
 */

  void
se_clip (void *space, fasta_t *reads, Uint elem, segemehl_t *nfo)
{

  if(nfo->hardclip3Prime || nfo->hardclip5Prime) {
    bl_fastaHardClip(space, reads, elem, nfo->hardclip5Prime, 
        nfo->hardclip3Prime);
    if(bl_fastaHasMate(reads)) {
      bl_fastaMateHardClip(space, reads, elem, nfo->hardclip5Prime, 
          nfo->hardclip3Prime);
    } 
  } 

  if(nfo->softclip3Prime || nfo->softclip5Prime) {
    bl_fastaSoftClip(space, reads, elem, 
        nfo->softclip5Prime, nfo->softclip5PrimeLen, nfo->minclipscr5,
        nfo->softclip3Prime, nfo->softclip3PrimeLen, nfo->clipacc, nfo->polyAlen);
    if(bl_fastaHasMate(reads)) {
      bl_fastaMateSoftClip(space, reads, elem, 
          nfo->softclip5Prime, nfo->softclip5PrimeLen, nfo->minclipscr5,
          nfo->softclip3Prime, nfo->softclip3PrimeLen, nfo->clipacc, nfo->polyAlen);
    }
  }

  return ;
}



/*--------------------------------- se_jump ----------------------------------
 *    
 * @brief a helper function for jump step calculation
 * @author Steve Hoffmann 
 *   
 */

  unsigned int
se_jump (unsigned int len, segemehl_t *nfo)
{
  unsigned int jump;

  //init jumps
  if (nfo->jump == 0) {
    jump = floor(len/75) * 2;
    jump = (jump > 0) ? jump : 1;
  } else {
    jump = nfo->jump;
  }

  return jump;
}


void
initStems(void *space, matchstem_t **stems, Uint len, Uint cond, segemehl_t *nfo) {
  Uint u,i;

  stems[0] = NULL; stems[1] = NULL; 

  /* restrict search to one strand */
  for (u = 0; u < 2; u++){
    /* nfo->strand == 1 : search only on plus strand
     * => init stems[1] as empty
     * nfo->strand == 2 : search only on minus strand
     * => init stems[0] as empty
     * Note: empty but initialized stems are ignored
     * in function kdbest
     */
    if (nfo->strand == ((cond==1)? 2-u : u+1)){	//2-u 
      stems[u] = ALLOCMEMORY(space, NULL, matchstem_t, len);
      for (i = 0; i < len; i++){
        stems[u][i].branches = NULL;
        stems[u][i].noofbranches = 0;
      }
    }
  }
}

void
getSeeds(void *space, Suffixarray *s, char **seqs, Uint len, Uint jump, 
    matchstem_t **stems, segemehl_t *nfo) {

  matchstem_t *b0[2]; 

  b0[0] = NULL; b0[1] = NULL;

  if (stems[0] == NULL){
    stems[0]=kdseeds(space, s, seqs[0], len, jump, nfo->s_ext, nfo->p_mis,
        nfo->Xoff, nfo->k_p, b0[0], nfo->nosuflinks);
  }

  if (stems[1] == NULL){
    stems[1]=kdseeds(space, s, seqs[1], len, jump, nfo->s_ext, nfo->p_mis,
        nfo->Xoff, nfo->k_p, b0[1], nfo->nosuflinks);
  }
}



/*------------------------------ se_matchlogic -------------------------------
 *    
 * @brief match logic
 * @author Steve Hoffmann 
 *   
 */

  void
se_matchlogic (Suffixarray *s,  
    bitvector *D, bitvector *Mv, bitvector *Mh, Uint *enctab, Uint maxL, 
    karlin_t *stats, fasta_t* reads, Uint k, segemehl_t* nfo)
{
  char *seqs[2], *mateseqs[2], *quals[2], *matequals[2];

  unsigned int jump; //, setmatches;
  Uint querysplitedist =0;
  Uint matesplitedist = 0;
  mappingset_t *tmp, *tmp2;
  mapseedlist_t *seeds=NULL, *mateseeds=NULL;
  matchstem_t *stems[2] = {NULL, NULL}, *matestems[2] = {NULL, NULL}; 
  void *space = NULL;
  char hasQuerySplit = 0;
  char hasMateSplit = 0;
  unsigned int len=0, matelen=0; //, setmatches;
  char *querydesc, *query, *queryqual=NULL; 
  char *mate=NULL, *matedesc=NULL, *matequal=NULL;
  Uint nchars;
  seseq_t myquery;
  seseq_t mymate;

  mappingset_t *maps, *matemaps, *querysplitmap, *matesplitmap;

  //inits
  seeds = NULL;
  mateseeds = NULL;

  stems[0] = NULL;
  stems[1] = NULL;
  matestems[0] = NULL;
  matestems[1] = NULL;
  querysplitmap = NULL;
  matesplitmap = NULL;
  hasQuerySplit = 0;
  hasMateSplit = 0;
  seqs[0] = NULL;
  seqs[1] = NULL;
  mateseqs[0] = NULL;
  mateseqs[1] = NULL;
  quals[0] = NULL;
  quals[1] = NULL;
  matequals[0] = NULL;
  matequals[1] = NULL;

  queryqual = NULL;
  mate = NULL;
  matedesc = NULL;
  matequal = NULL;

  querydesc = bl_fastaGetDescription(reads, k);
  query = bl_fastaGetSequence(reads, k);
  len = bl_fastaGetSequenceLength(reads, k);
  if(bl_fastaHasQuality(reads)) { 
    queryqual = bl_fastaGetQuality(reads, k);
  } 

  if(bl_fastaHasMate(reads)) {

    if(reads->checkid) {
      matedesc = bl_fastaGetMateDescription(reads, k);
    } else {
      matedesc = querydesc;
    }
    mate = bl_fastaGetMate(reads, k);
    matelen = bl_fastaGetMateLength(reads, k);
    if(bl_fastaHasQuality(reads)) {
      matequal = bl_fastaGetMateQuality(reads, k);
    }
    se_segemehlSeqInit(space, &mymate, mate, matequal, matelen); 
  } 

  se_segemehlSeqInit(space, &myquery, query, queryqual, len); 

  //fresh map
  maps = ALLOCMEMORY(space, NULL, mappingset_t, 1);
  bl_initMappingSet(maps);

  //clipping TODO
  //se_clip(space, reads, k, nfo);

  //attempt to split align poorly aligned reads
  querysplitedist = (((double)len)/100.0 * 2.0);
  //get data for the mate if avail
  matesplitedist = (((double)matelen)/100.0 * 2.0);

  /**************
   * STEP A: align first mate (aka query)
   *************/

  if(len >= nfo->minsize) {

    se_getData(space, &myquery, seqs, quals, nfo->bisulfite, 1);

    /*
     fprintf(stderr, "******* %s *******\n", querydesc); 
     fprintf(stderr, "checkpoint 1\n%s\n%s\n", seqs[0], seqs[1]);
    */
    
    /**************
     * STEP A1: if bestonly match read with pigeon approach with h holes 
     *************/

    if (!nfo->bisulfite && nfo->bestonly && (nchars = countNonMatchingChars(seqs[0], len)) <= 1){
       
      //get the pigeon map
      maps = pigeon(space, maps, s, nfo->seq, &myquery, 2, querydesc,
          enctab, D, 0, nchars, nfo);
           
      //fprintf(stderr, "A1.1: pigeon obtained map with %d matches.\n", maps->n);
    }

    
    /**************
     * STEP A2: if there was no (acceptable) hit -> get stems
     *************/

    if(!bl_hasQueryMapping(maps)) { 
      //init jumps, stems and get seeds
      jump = se_jump(len, nfo);
      initStems(space, stems, len, 1,  nfo);
      getSeeds(space, s, seqs, len, jump, stems, nfo);

      //filter the seeds
      seeds = bl_getGoodSeeds(stems, len, s->numofsuffixes, stats, nfo);

      //convert for alignment
      se_getData(space, &myquery, seqs, quals, nfo->bisulfite, 0);
      
      //fprintf(stderr, "checkpoint 2\n%s\n%s\n", seqs[0], seqs[1]);
      
      //align read to seeds
      bl_seedAlign(s, maps, nfo->seq, seqs, quals, len, querydesc, seeds, nfo, 
          enctab, D, 0);

      /*fprintf(stderr, "A2.1: seed align obtained map with %d matches.\n", maps->n);
        bl_dumpMappingSet(stderr, maps);
       */
    }

    /**************
     * STEP A3: if there was no (acceptable) or suboptimal hit -> split read
     *************/

    if(nfo->split && !bl_hasQueryMappingMaxEdist(maps, querysplitedist)) {
    //  fprintf(stdout, "A3.1: fixing query alignments.\n");

      bl_sortMappingSetByScore(maps, nfo->scores, nfo->scores[2]);

      bl_fixSplitAlignHoffmann (s, maps, nfo->seq, querydesc, 
          seqs, quals, len, 0, nfo);

      bl_fixSplitAlignBrendel(s, maps, nfo->seq, querydesc, 
          seqs, quals, len, 0, nfo);

      if(!bl_hasQueryMappingMaxEdist(maps, querysplitedist)) {

        if(stems[0] == NULL) {

          se_getData(space, &myquery, seqs, quals, nfo->bisulfite, 1);

          //init jumps, stems and get seeds
          jump = se_jump(len, nfo);
          initStems(space, stems, len, 1,  nfo);
          getSeeds(space, s, seqs, len, jump, stems, nfo);
        }

        //init a split map
        querysplitmap = ALLOCMEMORY(space, NULL, mappingset_t, 1);
        bl_initMappingSet(querysplitmap);

        //caluclate split map
        querysplitmap = bl_splitAlign(space, s,  nfo->seq, querysplitmap,
            querydesc, stems, seqs, quals, len, stats, enctab, D, 0, nfo);

        //fprintf(stdout, "A3.2: query split align obtained map with %d matches.\n", querysplitmap->n);
        //correct alignment
        if(bl_hasQueryMapping(querysplitmap)) { 

          bl_sortMappingSetByScore(querysplitmap, nfo->scores, nfo->scores[2]);
          
          bl_fixSplitAlignHoffmann (s, querysplitmap, nfo->seq, querydesc, 
              seqs, quals, len, 0, nfo);

          bl_fixSplitAlignBrendel(s, querysplitmap, nfo->seq, querydesc, 
              seqs, quals, len, 0, nfo);

          bl_concatMappingSet(maps, querysplitmap);
          hasQuerySplit =1;
        }

        //results have been concat to map, remove this map
        bl_removeMappingSet(querysplitmap);
        FREEMEMORY(space, querysplitmap);
      }
    }
  }

  /**************
   * STEP B: align the mate
   *************/
 

  if (matelen > 0 && matelen >= nfo->minsize) {
 
    //convert for alignment
    se_getData(space, &mymate, mateseqs, matequals, nfo->bisulfite, 0);
    //fprintf(stderr, "checkpoint 3 (mate)\n%s\n%s\n", mateseqs[0], mateseqs[1]);
    
    /**************
     * STEP B1: direct align to hits in map
     *************/

    if (!nfo->bisulfite && nfo->bestonly && (nchars = countNonMatchingChars(mateseqs[0], matelen)) <= 1){
      //get the pigeon map
      tmp = ALLOCMEMORY(space, NULL, mappingset_t, 1);
      bl_initMappingSet(tmp);
      
      tmp = pigeon(space, tmp, s, nfo->seq, &mymate, 2, matedesc, enctab, D, 1, nchars, nfo);

      if(bl_hasQueryMapping(maps) && bl_hasMateMapping(tmp)) {
        //tmp2 = bl_pairMateMapping(maps, tmp, 20000);
        tmp2 = bl_pairMateMapping(maps, tmp, nfo->maxpairinsertsize);
        if(bl_mappingsetHasPaired(tmp2)) { 
          bl_removeMappingSet(maps);
          FREEMEMORY(NULL, maps);
          maps = tmp2;
        } else {
          bl_removeMappingSet(tmp2);
          FREEMEMORY(NULL, tmp2);
        }
      }

      bl_removeMappingSet(tmp);
      FREEMEMORY(NULL, tmp);
      //fprintf(stdout, "A1.1: pigeon obtained map with %d matches.\n", maps->n);
    }
    
    if(!bl_mappingsetHasPaired(maps)) { 

      //find a direct mate for each of the hits
      tmp = bl_matealign(maps, nfo->seq, mateseqs, matequals, 
          matedesc, matelen, enctab, D, maxL, 1, s, nfo);

      bl_removeMappingSet(maps);
      FREEMEMORY(NULL, maps);
      maps = tmp;
    }

    /* 
     fprintf(stdout, "B1.1: direct align found paired matches? %d.\n", bl_mappingsetHasPaired(maps));
     bl_dumpMappingSet(stderr, maps);
    */
    
    /**************
     * STEP B2: if no or suboptimal pairs have been found
     *************/

    if(!bl_mappingsetHasPaired(maps)) {


      //convert for seeding
      se_getData(space, &mymate, mateseqs, matequals, nfo->bisulfite, 1);
      //fprintf(stderr, "checkpoint 4 (mate)\n%s\n%s\n", mateseqs[0], mateseqs[1]);
      
      //init mate jumps, stems and get seeds
      jump = se_jump(matelen, nfo);        
      initStems(space, matestems, matelen, 2,  nfo);
      getSeeds(space, s, mateseqs, matelen, jump, matestems, nfo);

      //convert for mate alignment
      se_getData(space, &mymate, mateseqs, matequals, nfo->bisulfite, 0);
      //fprintf(stderr, "checkpoint 5 (mate)\n%s\n%s\n", mateseqs[0], mateseqs[1]);
      
      //filter the seeds
      mateseeds = bl_getGoodSeeds(matestems, matelen, s->numofsuffixes, stats, nfo);

      //get a fresh mate map
      matemaps = ALLOCMEMORY(space, NULL, mappingset_t, 1);
      bl_initMappingSet(matemaps);

      //fprintf(stdout, "B2: seed align for mate triggered.\n");
      //align mate to the seeds
      bl_seedAlign(s, matemaps, nfo->seq, mateseqs, matequals, matelen, 
          matedesc,
          mateseeds, nfo, enctab, D, 1);

      //fprintf(stdout, "B2.1: mate seed align obtained map with %d matches.\n", matemaps->n);

      //now go check the query (bisulfite is already in alignment mode)
      if(len >= nfo->minsize) { 
        
        //convert for mate alignment
        se_getData(space, &myquery, seqs, quals, nfo->bisulfite, 0);

        tmp = bl_matealign(matemaps, nfo->seq, seqs, quals, 
            querydesc, len, enctab, D, maxL, 0, s, nfo);
        bl_removeMappingSet(matemaps);
        FREEMEMORY(NULL, matemaps);
        matemaps = tmp;
      }
      //fprintf(stdout, "B2.2: direct align (mate) found paired matches? %d.\n", bl_mappingsetHasPaired(matemaps));

      //pair the maps
      tmp = bl_pairMateMapping(maps, matemaps, nfo->maxpairinsertsize);
      bl_removeMappingSet(maps);
      bl_removeMappingSet(matemaps);
      FREEMEMORY(NULL, maps);
      FREEMEMORY(NULL, matemaps);
      maps = tmp;

      //fprintf(stdout, "B2.3: pairing (mate algin) found paired matches? %d.\n", bl_mappingsetHasPaired(maps));
    }

    /**************
     * STEP B3: attempt a cross correction       
     * *************/

    if(nfo->split && bl_mappingsetHasPaired(maps)) { 
      bl_sortMappingSetByScore(maps, nfo->scores, nfo->scores[2]);

      if(hasQuerySplit && len >= nfo->minsize && matelen >= nfo->minsize) {
        bl_crosscorrection(s, maps, nfo->seq, 
            matedesc, mateseqs, matequals, matelen, 
            len, matesplitedist, 1, nfo);

        //fprintf(stdout, "B3.1: cross correction of mate by query.\n");
      }
    }


    /**************
     * STEP B4: if no or suboptimal pairs have been found
     *************/

    if(nfo->split && (!bl_mappingsetHasPaired(maps) 
          || !bl_hasQueryMappingMaxEdist(maps, querysplitedist) 
          || !bl_hasMateMappingMaxEdist(maps, matesplitedist)
          || !bl_hasMappingPairedMaxEdist(maps, querysplitedist+matesplitedist))) {

     // fprintf(stdout, "B4.1: fixing query alignments.\n");

      bl_sortMappingSetByScore(maps, nfo->scores, nfo->scores[2]);


      bl_fixSplitAlignHoffmann (s, maps, nfo->seq, matedesc, 
          mateseqs, matequals, matelen, 1, nfo);

      bl_fixSplitAlignBrendel(s, maps, nfo->seq, matedesc, 
          mateseqs, matequals, matelen, 1, nfo);

      if (!bl_mappingsetHasPaired(maps) 
          || !bl_hasQueryMappingMaxEdist(maps, querysplitedist) 
          || !bl_hasMateMappingMaxEdist(maps, matesplitedist)
          || !bl_hasMappingPairedMaxEdist(maps, querysplitedist+matesplitedist)) { 

        if(matestems[0] == NULL) {
                  
          //convert for seeding
          se_getData(space, &mymate, mateseqs, matequals, nfo->bisulfite, 1);

          //init mate jumps, stems and get seeds
          jump = se_jump(matelen, nfo);        
          initStems(space, matestems, matelen, 2,  nfo);
          getSeeds(space, s, mateseqs, matelen, jump, matestems, nfo);
        }

        //init a split map
        matesplitmap = ALLOCMEMORY(space, NULL, mappingset_t, 1);
        bl_initMappingSet(matesplitmap);
        //caluclate split map
        matesplitmap = bl_splitAlign(space, s,  nfo->seq, matesplitmap,
            matedesc, 
            matestems, mateseqs, matequals, matelen, stats, enctab, D, 1, nfo);

        //fprintf(stdout, "B4.2: mate split align obtained map with %d matches.\n", matesplitmap->n);

        //correct alignment
        if(bl_hasMateMapping(matesplitmap)) { 

          bl_sortMappingSetByScore(matesplitmap, nfo->scores, nfo->scores[2]);
          
          bl_fixSplitAlignHoffmann (s, matesplitmap, nfo->seq, matedesc, 
              mateseqs, matequals, matelen, 1, nfo);

          //fprintf(stdout, "B4.3: fixing mate split alignments.\n");

          bl_fixSplitAlignBrendel(s, matesplitmap, nfo->seq, matedesc, 
              mateseqs, matequals, matelen, 1, nfo);

          hasMateSplit = 1;
        }

        //pair the maps
        tmp = bl_pairMateMapping(maps, matesplitmap, nfo->maxpairinsertsize);
        bl_removeMappingSet(maps);
        bl_removeMappingSet(matesplitmap);
        FREEMEMORY(NULL, maps);
        FREEMEMORY(NULL, matesplitmap);
        maps = tmp;

        //fprintf(stdout, "B3.3: pairing (matesplitmap) found paired matches? %d.\n", bl_mappingsetHasPaired(maps));
      }


      //there is a good query alignment, even a pair but it doesnt make a good pair.
      if(bl_mappingsetHasPaired(maps) && bl_hasQueryMappingMaxEdist(maps, querysplitedist) &&
          !bl_hasMappingPairedMaxEdist(maps, querysplitedist+matesplitedist)) {
        //fprintf(stdout, "CNT: split realigning query.\n");
      }
    }

    /**************
     * STEP B5: cross correct if a splitmap exists
     *************/

    bl_sortMappingSetByScore(maps, nfo->scores, nfo->scores[2]);

    if(hasQuerySplit && len >= nfo->minsize && matelen >= nfo->minsize) {
  //    fprintf(stdout, "B5.1: cross correction of mate by query.\n");
      bl_crosscorrection(s, maps, nfo->seq, 
          matedesc, mateseqs, matequals, matelen, 
          len, matesplitedist, 1, nfo);
    }

    if(hasMateSplit && matelen >= nfo->minsize && len >= nfo->minsize) {
//      fprintf(stdout, "B5.2: cross correction of query by mate.\n");
      bl_crosscorrection(s, maps, nfo->seq, 
          querydesc, seqs, quals, len, 
          matelen, querysplitedist, 0, nfo);
    }
  }

  bl_sortMappingSetByScore(maps, nfo->scores, nfo->scores[2]);

  //bl_dumpMappingSet(stdout, maps);
  //bl_removeBadMatesEdistQM(maps, querysplitedist, 0);
  //bl_removeBadMatesEdistQM(maps, matesplitedist, 1);

  mappingset_t* newmaps = bl_sammappingJoinFrags(maps, nfo);
  newmaps = bl_mappingsetRemoveDuplicates(newmaps, 10);

  /*GET THE BEST SEEDS (FOR UNMAPPED READS)*/ 
  mapseed_t *bestseed = NULL;
  if(len > 0 && seeds && seeds->n) { 
    bestseed = bl_getMapSeedListBest(seeds);
    bl_getMapSeedAdditionalInformation (bestseed, s, nfo->seq);
  }

  mapseed_t *bestmateseed = NULL;
  if(matelen > 0 && mateseeds && mateseeds->n) { 
    bestmateseed = bl_getMapSeedListBest(mateseeds);
    bl_getMapSeedAdditionalInformation (bestmateseed, s, nfo->seq);
  }

  /*
   fprintf(stderr, "writing mappingsets to output:\n");
   bl_dumpMappingSet(stderr, newmaps);
  */

  /*OUTPUT*/
  se_output(newmaps, reads, k, bestseed, bestmateseed, s, stats, nfo);


  /*LOOP CLEANUP BEGIN*/ 

  se_segemehlSeqDestruct(space, &myquery);
  if(matelen > 0) {
    se_segemehlSeqDestruct(space, &mymate);
  }


  if(seeds) { 
    bl_wrapSeedList(seeds);
    FREEMEMORY(space, seeds);
  }

  if(mateseeds) {
    bl_wrapSeedList(mateseeds);
    FREEMEMORY(space, mateseeds);
  }

  if(stems[0]) { 
    bl_kdMatchstemDestruct(space, stems[0], len);
    bl_kdMatchstemDestruct(space, stems[1], len);
    stems[0] = NULL;
    stems[1] = NULL;
  }

  bl_removeMappingSet(maps);
  FREEMEMORY(space, maps);
  bl_removeMappingSet(newmaps);
  FREEMEMORY(space, newmaps);

  if(matestems[0]) {
    bl_kdMatchstemDestruct(space, matestems[0], matelen);
    bl_kdMatchstemDestruct(space, matestems[1], matelen);
    matestems[0] = NULL;
    matestems[1] = NULL;
  }
 
  return ;  
}



void
match(void *space, Suffixarray *s, fasta_t *reads, segemehl_t *nfo) {

  Uint i;
  fasta_t* readchunk;  
  unsigned int *enctab, dim, maxL, wordno; 
  karlin_t stats;
  bitvector *D, *Mv, *Mh;

  //initialize data for stats and bv
  enctab = encodetab(nfo->seq->map, nfo->seq->mapsize);
  dim = getExtendedAlignmentLength(reads->maxlen+1000, nfo->accuracy);
  maxL = nfo->maxaligninsertsize;

  //add the maximum insert size if we use the myers for looking up mate
  if(bl_fastaHasMate(reads)) {
    dim += nfo->maxaligninsertsize;
  }

  wordno = (reads->maxlen/BITVECTOR_WORDSIZE)+1; 

  D = ALLOCMEMORY(space, NULL, bitvector, 3*(dim+1));
  Mv = &D[dim+1];
  Mh = &D[2*(dim+1)];

  for(i=0; i <= dim; i++) {
    D[i] = initbitvector(space, wordno*BITVECTOR_WORDSIZE);
    Mv[i] = initbitvector(space, wordno*BITVECTOR_WORDSIZE);
    Mh[i] = initbitvector(space, wordno*BITVECTOR_WORDSIZE);
  }  

  stats.lambda = 1.098612;
  stats.H = 0.549306;
  stats.K = 0.333333;

  while((readchunk = se_fastaMaster(space, reads, 1, nfo))) { 

    for(i=0; i < readchunk->noofseqs; i++) { 

      //progress bar   
      if(!nfo->mute) se_updateProgressBar(i, nfo);
 
      se_matchlogic (s, D, Mv, Mh, enctab, maxL, &stats, readchunk, i, nfo);
    }

    if(nfo->threadno > 1) { 
     bl_fastxDestructSequence(space, readchunk);
     bl_fastxDestructChunkIndex(space, readchunk);
     FREEMEMORY(space, readchunk);
    }
  }

  wrapBitmatrix(space, D, 3*(dim+1));
  FREEMEMORY(space, D);
  FREEMEMORY(space, enctab);
}

