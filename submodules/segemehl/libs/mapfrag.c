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
 *  mapfrag.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 09.12.2012 16:21:44 CET
 *  
 */

#include "mapfrag.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <inttypes.h>
#include <pthread.h>

#include "basic-types.h"

#include "mathematics.h"
#include "segemehl.h"
#include "karlin.h"
#include "locus.h"
#include "memory.h"	//(skipped)
#include "filebuffer.h"	//(skipped)
#include "stringutils.h"	//(skipped)
#include "multicharseq.h"	//(skipped)
#include "alignment.h"	//(skipped)


/*
 *  functions to manage alignment/map
 *  fragments on reference genomes.
 *
 *
 */

/*------------------------------- bl_copyMCSA --------------------------------
 *    
 * @brief copy a multi char seq alignment
 * @author Steve Hoffmann 
 *   
 */

MultiCharSeqAlignment *
bl_copyMCSA(MultiCharSeqAlignment *dest, MultiCharSeqAlignment *src) {

  dest->subidx = src->subidx;
  dest->substart = src->substart;
  dest->subend = src->subend;

  dest->refdesc = src->refdesc;    
  dest->refseq = src->refseq;
  dest->refstart = src->refstart;
  dest->reflen = src->reflen;   
  dest->floff = src->floff;

  dest->qrydesc = src->qrydesc;
  dest->query = src->query;
  dest->qual = src->qual;
  dest->qrystart = src->qrystart;
  dest->strand = src->strand;
  dest->qrylen = src->qrylen;
  dest->maxalignlen = src->maxalignlen;
  dest->pass = src->pass;

  dest->al = ALLOCMEMORY(space, NULL, Alignment, 1);
  copyAlignment(dest->al, src->al);

  return dest;
}

/*----------------------------- bl_wrapSeedList ------------------------------
 *    
 * @brief destruct list of seeds
 * @author Steve Hoffmann 
 *   
 */

  void
bl_wrapSeedList (mapseedlist_t *l)
{
  FREEMEMORY(NULL, l->l);
  return ;
}


/*----------------------------- bl_initMapFrag ------------------------------
 *    
 * @brief add a mapping fragment from multicharseqalignment 
 * @author Steve Hoffmann 
 *   
 */

  void
bl_initMapFrag (mapfrag_t *f, char* seq, char *qual, 
    MultiCharSeqAlignment *mcsa, mapseed_t *seed,
    unsigned char mate, unsigned char issplit)
{
  //Uint p, q, u, v, e;
  //Uint mat, mis, ins, del;

  f->seed = seed;
  f->mcsa = mcsa;
  f->mate = mate;
  f->seq = seq;
  f->qual = qual;
  f->issplit = issplit;
  f->leftgap = 0;
  f->rightgap = 0;
  f->lclip = 0;
  f->rclip = 0;
  f->mapq = 0;
  f->nextnoncollinear =0;
  f->prevnoncollinear =0;
  f->mapq_dbl = -500.0;
  countEops(f->mcsa->al, &f->mat, &f->mis, &f->ins, &f->del, &f->lmat);

  return ;
}

/*--------------------------- bl_getMapFragQryDesc ---------------------------
 *    
 * @brief get the query name
 * @author Steve Hoffmann 
 *   
 */

  char*
bl_getMapFragQryDesc (mapfrag_t *f)
{
  return f->mcsa->qrydesc;
}

/*----------------------------- bl_getMapFragQry -----------------------------
 *    
 * @brief get the map frag query
 * @author Steve Hoffmann 
 *   
 */

  char*
bl_getMapFragQry (mapfrag_t *f)
{
  char *seq;
  Uint uallen=0;

  uallen = getUalignlen(f->mcsa->al);

  seq = ALLOCMEMORY(NULL, NULL, char, uallen+1);
  memmove(seq, &f->mcsa->al->u[f->mcsa->al->uoff], uallen);

  /*  if(bl_getMapFragStrand(f)) { ASSUMPTION THAT STRANDEDNESS IS TAKEN CARE OF
      tmp = charIUPACcomplement(NULL, seq, uallen);
      FREEMEMORY(NULL, seq);
      seq = tmp;
      }
      */
  seq[uallen] = 0;

  return seq;
}

/*---------------------------- bl_getMapFragQual -----------------------------
 *    
 * @brief get the quality string. as the quality is not part of the MCSA,
 * it has to taken from the parent mapping_t
 * @author Steve Hoffmann 
 *   
 */

  char*
bl_getMapFragQual (mapfrag_t *f)
{
  char *qual = NULL;
  Uint uallen = 0;

  if(f->mcsa->qual) { 
    uallen = getUalignlen(f->mcsa->al);

    qual = ALLOCMEMORY(NULL, NULL, char, uallen+1);
    memmove(qual, &f->mcsa->qual[f->mcsa->al->uoff], uallen);
    qual[uallen] = 0;

    if(bl_getMapFragStrand(f)) {
      qual = strrev(qual, uallen);
    }
  } 

  return qual;
}

/*--------------------------- bl_setMapFragMapQual ----------------------------
 *    
 * @brief set MapFragMapQual
 * @author Steve Hoffmann 
 *   
 */

  void
bl_setMapFragMapQual (mapfrag_t *f, double mapqual)
{

  /*  if(mapqual > f->mapq_dbl) {
      fprintf(stdout, "changeing maxqual from %f to %f\n", f->mapq_dbl, mapqual);
      } else {
      fprintf(stdout, "keeping maxqual %f instead of %f\n", f->mapq_dbl, mapqual);
      }
      */
  f->mapq_dbl = MAX(mapqual, f->mapq_dbl);


  return;
}

/*--------------------------- bl_getMapFragMapQual ----------------------------
 *    
 * @brief set MapFragMapQual
 * @author Steve Hoffmann 
 *   
 */

  double
bl_getMapFragMapQual (mapfrag_t *f)
{
  return f->mapq_dbl;
}



/*----------------------------- getMapFragIsMate -----------------------------
 *    
 * @brief return the mate status
 * @author Steve Hoffmann 
 *   
 */

  unsigned char
bl_getMapFragIsMate (mapfrag_t *f)
{
  return f->mate;
}

/*--------------------------- bl_getMapFragStrand ----------------------------
 *    
 * @brief get strand
 * @author Steve Hoffmann 
 *   
 */

  unsigned char
bl_getMapFragStrand (mapfrag_t *f)
{
  return f->mcsa->strand;
}

/*------------------------------ bl_getMapFragRight ------------------------------
 *    
 * @brief get the fragment end wrt query
 * @author Steve Hoffmann 
 *   
 */

unsigned int
bl_getMapFragRight(mapfrag_t *f) {
  Uint uallen = 0;

  uallen = getUalignlen(f->mcsa->al);
  return f->mcsa->al->uoff + uallen - 1;

}

/*------------------------------ bl_getMapFragLeft ------------------------------
 *    
 * @brief get the fragment start wrt to query
 * @author Steve Hoffmann 
 *   
 */

unsigned int
bl_getMapFragLeft(mapfrag_t *f) {

  return f->mcsa->al->uoff;

}
/*------------------------------ bl_getMapFragV ------------------------------
 *    
 * @brief get the fragment end wrt query
 * @author Steve Hoffmann 
 *   
 */

unsigned int
bl_getMapFragV(mapfrag_t *f) {
  Uint uallen = 0;

  if(bl_getMapFragStrand(f)) {
    return f->mcsa->al->ulen - f->mcsa->al->uoff - 1;
  } else {
    uallen = getUalignlen(f->mcsa->al);
    return f->mcsa->al->uoff + uallen - 1;
  }

  return 0;
}

/*------------------------------ bl_getMapFragU ------------------------------
 *    
 * @brief get the fragment start wrt to query
 * @author Steve Hoffmann 
 *   
 */

unsigned int
bl_getMapFragU(mapfrag_t *f) {
  unsigned int uallen = 0;

  if(bl_getMapFragStrand(f)) {
    uallen = getUalignlen(f->mcsa->al);

    if(f->mcsa->al->ulen < f->mcsa->al->uoff + uallen) {
      fprintf(stderr, "uoff: %d, uallen: %d, qryllen: %d, voff:%d\n'%s'\n", 
          f->mcsa->al->uoff, uallen, f->mcsa->al->ulen, f->mcsa->al->voff, f->seq);
      showAlign(f->mcsa->al, stderr);
      exit(0);
    }

    return f->mcsa->al->ulen - (f->mcsa->al->uoff + uallen);
  } else {
    return f->mcsa->al->uoff;
  }
  return 0;
}

/*-------------------------- bl_getMapFragSubstart ---------------------------
 *    
 * @brief get the start position of the mapfrags chromosome
 * @author Steve Hoffmann 
 *   
 */

  uint64_t
bl_getMapFragSubstart (mapfrag_t *f)
{
  return f->mcsa->substart;
}
/*-------------------------- bl_getMapFragAlignment --------------------------
 *    
 * @brief get the alignment
 * @author Steve Hoffmann 
 *   
 */

  Alignment*
bl_getMapFragAlignment (mapfrag_t *f)
{
  return f->mcsa->al;
}

/*------------------------------ bl_getMapFragP ------------------------------
 *    
 * @brief get the fragment start wrt reference
 * @author Steve Hoffmann 
 *   
 */

unsigned int
bl_getMapFragP(mapfrag_t *f) {
  return f->mcsa->refstart + f->mcsa->al->voff;
}

/*------------------------------ bl_getMapFragQ ------------------------------
 *    
 * @brief get the fragment end wrt reference
 * @author Steve Hoffmann 
 *   
 */


unsigned int
bl_getMapFragQ(mapfrag_t *f) {

  //vlen was mb.a directly from bitmatrix. Is that the same? NO!
  return bl_getMapFragP(f) + getValignlenAndSkipped(f->mcsa->al) -1; //getValignlen(f->mcsa->al) -1;
}

/*----------------------------- bl_getMapFragULen ----------------------------
 *    
 * @brief get the frag alignment length of the query
 * @author Steve Hoffmann 
 *   
 */

  unsigned int
bl_getMapFragGetUlen (mapfrag_t *f)
{
  return getUalignlen(f->mcsa->al);
}
/*--------------------------- bl_getMapFragULenNoClip -------------------------
 *    
 * @brief get the frag alignment length of the query w/o clipping
 * @author Steve Hoffmann 
 *   
 */

  unsigned int
bl_getMapFragGetUlenNoClip (mapfrag_t *f)
{
  return getUalignlenNoClip(f->mcsa->al);
}

/*----------------------------- bl_getMapFragVLen ----------------------------
 *    
 * @brief get the frag alignment length of the reference
 * @author Steve Hoffmann 
 *   
 */

  unsigned int
bl_getMapFragGetVlen (mapfrag_t *f)
{
  return getValignlen(f->mcsa->al);
}

/*---------------------------- bl_getMapFragEdist ----------------------------
 *    
 * @brief get the edit distance of a map frag
 * @author Steve Hoffmann 
 *   
 */

unsigned int
bl_getMapFragEdist(mapfrag_t *f) {
  return f->mis + f->ins + f->del;
}
/*------------------------- bl_getMapFragLongestMatch ------------------------
 *    
 * @brief get the longest match
 * @author Steve Hoffmann 
 *   
 */

unsigned int
bl_getMapFragLongestMatch(mapfrag_t *f) {
  return f->lmat;
}


/*---------------------------- bl_getMapFragScore ----------------------------
 *    
 * @brief get the score of a mapping fragment
 * @author Steve Hoffmann 
 *   
 */

  unsigned int
bl_getMapFragScore (mapfrag_t *f, int *scores, int indel)
{
  return getAlignScore(f->mcsa->al, scores, indel);
}

/*---------------------------- bl_getMapFragMCSA -----------------------------
 *    
 * @brief get the multicharseq alignment of a map frag
 * @author Steve Hoffmann 
 *   
 */ 

MultiCharSeqAlignment*
bl_getMapFragMCSA(mapfrag_t *f) {
  return f->mcsa;
}

/*---------------------------- bl_getmapFragSeed -----------------------------
 *    
 * @brief get the seed of a mapping frag. Could be NULL
 * @author Steve Hoffmann 
 *   
 */

mapseed_t *
bl_getMapFragSeed(mapfrag_t *f) {
  return f->seed;
}

/*--------------------------- bl_getMapFragChrIdx ----------------------------
 *    
 * @brief get the index number of the reference chr the mapfrag is aligned to
 * @author Steve Hoffmann 
 *   
 */

unsigned int
bl_getMapFragChrIdx(mapfrag_t *f) {
  return f->mcsa->subidx;
}

/*--------------------------- bl_getMapFragChrDesc ---------------------------
 *    
 * @brief get the description of the reference chr of the mapfrag
 * @author Steve Hoffmann 
 *   
 */

char*
bl_getMapFragRefDesc(mapfrag_t *f) {
  return f->mcsa->refdesc;
}

/*---------------------------- bl_getMapFragSplit ----------------------------
 *    
 * @brief get the issplit field
 * @author Steve Hoffmann 
 *   
 */

  unsigned int
bl_getMapFragSplit (mapfrag_t *f)
{

  return f->issplit;
}
/*---------------------------- bl_getMappingUlen -----------------------------
 *    
 * @brief get the length of the alignment on the query (mate and query)
 * @author Steve Hoffmann
 *   
 */

  unsigned int
bl_getMappingUlen (mapping_t *l)
{
  unsigned int i, ulen = 0;

  for(i=0; i < l->n; i++) {
    ulen += bl_getMapFragGetUlen(&l->f[i]);
  }
  return ulen;
}

/*------------------------- bl_getMappingMinFragDist -------------------------
 *    
 * @brief 
 * @author Steve Hoffmann 
 *   
 */

  Uint 
bl_getMappingMinFragDist (mapping_t *l)
{

  unsigned int i, prevstart = 0, prevend =0, curstart =0, curend =0, dist=-1;

  if(l->n == 0) return 0;

  prevstart = bl_getMapFragP(&l->f[0]);
  prevend = bl_getMapFragQ(&l->f[0]);

  for(i=1; i < l->n; i++) {

    curstart = bl_getMapFragP(&l->f[i]);
    curend = bl_getMapFragQ(&l->f[i]);

    dist = MIN(dist, DIST(prevstart, curend));
    dist = MIN(dist, DIST(prevend, curstart));

    prevstart = curstart;
    prevend = curend;
  }

  return dist;
}

/*-------------------------- bl_getMappingFastaSet ---------------------------
 *    
 * @brief get the fasta set 
 * @author Steve Hoffmann 
 *   
 */

  char*
bl_getMappingSeq (mapping_t *f)
{
  return f->seq;
}



/*---------------------------- bl_getMappingQual -----------------------------
 *    
 * @brief get the quality string of the mapped sequence
 * @author Steve Hoffmann 
 *   
 */

  char*
bl_getMappingQual (mapping_t *f)
{
  return f->qual;
}



/*------------------------------ bl_initMapping ------------------------------
 *    
 * @brief initialize a mapping structure
 * @author Steve Hoffmann 
 *   
 */

void
bl_initMapping(mapping_t *l, char* seq, char* qual, Uint lclip, Uint rclip) {

  l->seq = seq;
  l->qual = qual;
  //Replace the following to avoid multiple calculations
  l->seqlen = strlen(seq);
  l->lclip = lclip;
  l->rclip = rclip;
  l->n = 0;
  l->f = NULL;
  l->P = -1;
  l->Q = 0;
  l->matestatus= 0;
  l->consecutive = 1;
  l->maxqual = -100.0;
  l->scr = INT_MIN;
  return;
}

/*------------------------ bl_getMappingIsConsecutive ------------------------
 *    
 * @brief is mapping consecutive
 * @author Steve Hoffmann 
 *   
 */

  char
bl_getMappingIsConsecutive(mapping_t *m)
{
  return m->consecutive;
}

/*---------------------------- bl_getMappingScore ----------------------------
 *    
 * @brief get the score of the mapping (query and mate)
 * @author Steve Hoffmann 
 *   
 */

  int
bl_getMappingScore (mapping_t *l, int * scores, int indel)
{
  unsigned int i;
  int scr = 0;

  for(i=0; i < l->n; i++) {
    scr += bl_getMapFragScore(&l->f[i], scores, indel);
  }
  return scr;
}

/*---------------------------- bl_getMappingScore ----------------------------
 *    
 * @brief get the score of the mapping (query and mate)
 * @author Steve Hoffmann 
 *   
 */

  int
bl_getMappingScoreQM (mapping_t *l, int * scores, int indel, int *qscore, int *mscore)
{
  unsigned int i;
  int queryscore=0, matescore=0, score;

  for(i=0; i < l->n; i++) {
    if(bl_getMapFragIsMate(&l->f[i])) { 
      matescore += bl_getMapFragScore(&l->f[i], scores, indel);
    } else {
      queryscore += bl_getMapFragScore(&l->f[i], scores, indel);
    }
  }

  score = queryscore + matescore;

  *qscore = queryscore;
  *mscore = matescore;

  return score;
}



/*---------------------------- bl_getMappingEdist ----------------------------
 *    
 * @brief get the edist of the mapping (query and mate)
 * @author Steve Hoffmann 
 *   
 */

  unsigned int
bl_getMappingEdist (mapping_t *l)
{
  unsigned int i, edist = 0;

  for(i=0; i < l->n; i++) {
    edist += bl_getMapFragEdist(&l->f[i]);
  }
  return edist;
}


/*--------------------------- bl_getMappingEdistQM ---------------------------
 *    
 * @brief get mapping edists for query and mate
 * @author Steve Hoffmann 
 *   
 */

  unsigned int
bl_getMappingEdistQM (mapping_t *l, Uint *qedist, Uint *medist)
{
  unsigned int i, queryedist = 0, mateedist = 0, edist;

  for(i=0; i < l->n; i++) {
    if(bl_getMapFragIsMate(&l->f[i])) { 
      mateedist += bl_getMapFragEdist(&l->f[i]);
    } else {
      queryedist += bl_getMapFragEdist(&l->f[i]);
    }
  }

  edist = queryedist + mateedist;

  *qedist = queryedist;
  *medist = mateedist;

  return edist;
}


/*----------------------------- bl_getMappingCovQM ---------------------------
 *    
 * @brief get mapping edists for query and mate
 * @author Steve Hoffmann 
 *   
 */

  unsigned int
bl_getMappingCovQM (mapping_t *l, Uint *qcover, Uint *mcover)
{
  unsigned int i, qcov = 0, mcov = 0, cov;

  for(i=0; i < l->n; i++) {
    if(bl_getMapFragIsMate(&l->f[i])) { 
      mcov += bl_getMapFragGetUlenNoClip(&l->f[i]);
    } else {
      qcov+= bl_getMapFragGetUlenNoClip(&l->f[i]);
    }
  }

  cov = qcov + mcov;

  *qcover = qcov;
  *mcover = mcov;

  return cov;
}


  int
bl_getMappingFragNoQM (mapping_t *l, Uint *qno, Uint *mno)
{
  unsigned int i;
  int queryno = 0, mateno = 0, no;

  for(i=0; i < l->n; i++) {
    if(bl_getMapFragIsMate(&l->f[i])) { 
      mateno++;
    } else {
      queryno++; 
    }
  }

  no = queryno + mateno;

  *qno = queryno;
  *mno = mateno;

  return no;
}


/*----------------------------- bl_removeMapping -----------------------------
 *    
 * @brief removes a mapping from the heap. ATTENTION: seeds are NOT removed!
 * @author Steve Hoffmann 
 *   
 */

  void
bl_removeMapping (mapping_t *l)
{
  unsigned int i;

  if (l == NULL) return ;

  for(i=0; i < l->n; i++) {
    wrapMultiCharSeqAlignment(NULL, l->f[i].mcsa);
    FREEMEMORY(NULL, l->f[i].mcsa);
  }
  FREEMEMORY(NULL, l->f);
  l->n = 0;
  l->matestatus = 0;

  return ;
}



/*----------------------------- bl_removeMapping -----------------------------
 *    
 * @brief removes a mapping from the heap. ATTENTION: seeds are NOT removed!
 * @author Steve Hoffmann 
 *   
 */

  void
bl_removeMappingQM (mapping_t *l, char ismate)
{
  unsigned int i, newsize=0, P=-1, Q=0;
  char consecutive=1;
  MultiCharSeqAlignment *mcsa;

  if (l == NULL) return ;

  mapfrag_t *myfrags = NULL;

  l->P = -1;
  l->Q = 0;

  for(i=0; i < l->n; i++) {
    if(bl_getMapFragIsMate(&l->f[i]) == ismate) {
      wrapMultiCharSeqAlignment(NULL, l->f[i].mcsa);
      FREEMEMORY(NULL, l->f[i].mcsa);
    } else {

      myfrags = ALLOCMEMORY(space, myfrags, mapfrag_t, newsize+1);
      memmove(&myfrags[newsize],&l->f[i],sizeof(mapfrag_t));

      if(newsize && bl_getMapFragChrIdx(&myfrags[0]) != bl_getMapFragChrIdx(&l->f[i])) {
        consecutive =0;
      }

      mcsa = l->f[i].mcsa;
      //keep the list sorted

      P = mcsa->refstart+mcsa->al->voff;
      Q = P + getValignlen(mcsa->al) - 1;

      if(consecutive) {
        //update first and last position
        l->P = (l->P < P) ? l->P : P;
        l->Q = (l->Q > Q) ? l->Q : Q;
      } else { 
        l->P = 0;
        l->Q = 0;
      }

      newsize++;
    }
  }

  int mask = (1 << ismate);
  l->matestatus &= ~mask;

  FREEMEMORY(NULL, l->f);
  l->consecutive = consecutive;
  l->n = newsize;
  l->f = myfrags;


  return ;
}

/*---------------------------- bl_getMappingRange ----------------------------
 *    
 * @brief get the range of a mapping; 0 if not consecutive
 * @author Steve Hoffmann 
 *   
 */

  Uint
bl_getMappingRange (mapping_t *m)
{
  if(!m->consecutive || m->P == -1) return 0;

  return m->Q - m->P + 1;
}

/*--------------------------- bl_removeMappingSet ----------------------------
 *    
 * @brief remove a mapping set
 * @author Steve Hoffmann 
 *   
 */

  void
bl_removeMappingSet (mappingset_t *s)
{
  unsigned int i;

  if (s == NULL) return ;
  for(i=0; i < s->n; i++) {
    bl_removeMapping(&s->elem[i]);
  }

  FREEMEMORY(NULL, s->elem);
  s->n = 0;

  return ;
}


/*----------------------------- bl_copyMapping -------------------------------
 *    
 * @brief 
 * @author Steve Hoffmann 
 *   
 */

  mapping_t*
bl_copyMapping (mapping_t* dest, mapping_t* source)
{

  Uint i;
  MultiCharSeqAlignment *copy;

  bl_initMapping(dest, source->seq, source->qual, source->lclip, source->rclip);

  for(i=0; i < source->n; i++) {
    copy = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
    copy = bl_copyMCSA(copy, source->f[i].mcsa);
    bl_addMapFrag(dest, copy, NULL, source->f[i].mate, source->f[i].issplit);
  }

  return dest;
}


/*------------------------------ bl_getMappings ------------------------------
 *    
 * @brief get the mappings from mapping set
 * @author Steve Hoffmann 
 *   
 */

mapping_t*
bl_getMappings(mappingset_t *s, Uint *size) {
  Uint i, k=0;
  mapping_t *mymappings = NULL;

  for(i=0; i < s->n; i++) {
    mymappings = ALLOCMEMORY(NULL, mymappings, mapping_t, k+1);
    memmove(&mymappings[k], &s->elem[i], sizeof(mapping_t));
    k++;
  }

  *size =k;
  return mymappings;
}

/*------------------------ bl_getMapFragIsChimeric ------------------------
 *    
 * @brief map frag is chimeric
 * @author Steve Hoffmann 
 *   
 */

  char
bl_getMapFragIsChimeric(mapfrag_t *f)
{
  return f->nextnoncollinear || f->prevnoncollinear;
}


/*------------------------ bl_getMappingIsChimericQM -------------------------
 *    
 * @brief 
 * @author Steve Hoffmann 
 *   
 */

  char
bl_getMappingIsChimericQM (mapping_t *m, char ismate)
{
  Uint i;

  for(i=0; i < m->n; i++) {
    if(bl_getMapFragIsMate(&m->f[i]) == ismate) { 
      if(bl_getMapFragIsChimeric(&m->f[i])) 
        return 1;
    } 
  }

  return 0;
}

char
bl_mappingHasPairedEndOrientation(mapping_t *mapping) {
  Uint i, rc=-1;
  char mate;

  for(i=0; i < mapping->n; i++) {
    if(rc == -1) { 
      rc = bl_getMapFragStrand(&mapping->f[i]); 
      mate = bl_getMapFragIsMate(&mapping->f[i]);
    } 
    if(mate == bl_getMapFragIsMate(&mapping->f[i]) && 
        rc != bl_getMapFragStrand(&mapping->f[i]))
      return 0;
    if(mate != bl_getMapFragIsMate(&mapping->f[i]) &&
        rc == bl_getMapFragStrand(&mapping->f[i]))
      return 0;
  }

  return 1;
}


/*------------------------------ bl_getMapFrags ------------------------------
 *    
 * @brief get the mapping frags from mapping
 * @author Steve Hoffmann 
 *   
 */

  mapfrag_t*
bl_getMapFrags (mapping_t *m, Uint *size)
{
  Uint i, k=0;
  mapfrag_t* myfrags = NULL;

  for(i=0; i < m->n; i++) {
    myfrags = ALLOCMEMORY(NULL, myfrags, mapfrag_t, k+1);
    memmove(&myfrags[k], &m->f[i], sizeof(mapfrag_t));
    k++;
  }

  *size = k;
  return myfrags ;
}

/*----------------------------- bl_getMapFragsQM -----------------------------
 *    
 * @brief get the mapping frags from mapping
 * @author Steve Hoffmann 
 *   
 */

  mapfrag_t*
bl_getMapFragsQM (mapping_t *m, Uint *size, char ismate)
{
  Uint i, k=0;
  mapfrag_t* myfrags = NULL;

  for(i=0; i < m->n; i++) {
    if( bl_getMapFragIsMate(&m->f[i]) == ismate) { 
      myfrags = ALLOCMEMORY(NULL, myfrags, mapfrag_t, k+1);
      memmove(&myfrags[k], &m->f[i], sizeof(mapfrag_t));
      k++;
    }
  }

  *size = k;
  return myfrags ;
}

/*------------------------------ bl_addMapFrag -------------------------------
 *    
 * @brief add a mapping fragment to a mapping maintains the order of the 
 * fragments on the query. First query, second mate!
 * @author Steve Hoffmann 
 *   
 */

Uint
bl_addMapFrag(mapping_t *l, MultiCharSeqAlignment *mcsa, 
    mapseed_t *seed, unsigned char mate, unsigned char issplit) {

  unsigned int i=0, U=0, P=0, Q, curidx = mcsa->subidx;
  int gap=0;

  l->n++;
  l->f = ALLOCMEMORY(NULL, l->f, mapfrag_t, l->n); 

  //keep the list sorted
  P = mcsa->refstart+mcsa->al->voff;
  Q = P + getValignlen(mcsa->al) - 1;

  Uint uallen = getUalignlen(mcsa->al);

  if(mcsa->strand) {
    U = mcsa->al->ulen - (mcsa->al->uoff + uallen);
  } else {
    U = mcsa->al->uoff;
  }

  //fast forward to mate positions
  while(mate && i+1 < l->n && !bl_getMapFragIsMate(&l->f[i])) i++;
  //keep order wrt U
  while (i+1 < l->n && bl_getMapFragIsMate(&l->f[i]) == mate && 
      bl_getMapFragU(&l->f[i]) < U) i++;


  //fragment is inserted
  if(i+1 < l->n) { 
    curidx = bl_getMapFragChrIdx(&l->f[i]);
    memmove(&l->f[i+1], &l->f[i], sizeof(mapfrag_t)*((l->n-1)-i));
  } else if (i > 0){
    curidx = bl_getMapFragChrIdx(&l->f[i-1]);
  }

  bl_initMapFrag(&l->f[i], l->seq, l->qual, mcsa, seed, mate, issplit);


  //update the left gap
  if(i > 0 && bl_getMapFragIsMate(&l->f[i-1]) == mate) {

    Uint thisstart = bl_mcsaGet5PrimeU(mcsa);
    MultiCharSeqAlignment *prevmcsa = bl_getMapFragMCSA(&l->f[i-1]);
    Uint prevend = bl_mcsaGet3PrimeU(prevmcsa);
    gap = (thisstart-prevend)-1;
    l->f[i-1].rightgap = gap;
    l->f[i].leftgap = gap;

    if(bl_getMapFragStrand(&l->f[i-1]) != bl_getMapFragStrand(&l->f[i])) { 
      l->f[i-1].nextnoncollinear |= 1 << 0;
      l->f[i].prevnoncollinear |= 1 << 0;
    }
    if(bl_getMapFragChrIdx(&l->f[i-1]) != bl_getMapFragChrIdx(&l->f[i])) {
      l->f[i-1].nextnoncollinear |= 1 << 1;
      l->f[i].prevnoncollinear |= 1 << 1;
    }

    if(bl_getMapFragQ(&l->f[i-1]) >= bl_getMapFragP(&l->f[i]) && !bl_getMapFragStrand(&l->f[i])) {
      l->f[i-1].nextnoncollinear |= 1 << 2;
      l->f[i].prevnoncollinear |= 1 << 2;
    }

    if(bl_getMapFragP(&l->f[i-1]) <=  bl_getMapFragQ(&l->f[i]) && bl_getMapFragStrand(&l->f[i])) {
      l->f[i-1].nextnoncollinear |= 1 << 3;
      l->f[i].prevnoncollinear |= 1 << 3;
    }

  } else {
    l->f[i].leftgap = bl_mcsaGet5PrimeU(mcsa);
    if(i > 0) { 
      MultiCharSeqAlignment *prevmcsa = bl_getMapFragMCSA(&l->f[i-1]);
      Uint prevend = bl_mcsaGet3PrimeU(prevmcsa);
      Uint len = prevmcsa->qrylen;
      l->f[i-1].rightgap = (len - prevend) - 1;
    }
  }

  //update the right gap
  if(i+1 < l->n && bl_getMapFragIsMate(&l->f[i+1]) == mate) {
    MultiCharSeqAlignment *nextmcsa = bl_getMapFragMCSA(&l->f[i+1]);
    Uint nextstart = bl_mcsaGet5PrimeU(nextmcsa);
    Uint thisend = bl_mcsaGet3PrimeU(mcsa);
    gap = (nextstart-thisend)-1;
    l->f[i].rightgap = gap;
    l->f[i+1].leftgap = gap;

    if(bl_getMapFragStrand(&l->f[i+1]) != bl_getMapFragStrand(&l->f[i])) { 
      l->f[i+1].nextnoncollinear |= 1 << 0;
      l->f[i].prevnoncollinear |= 1 << 0;
    }
    if(bl_getMapFragChrIdx(&l->f[i+1]) != bl_getMapFragChrIdx(&l->f[i])) {
      l->f[i+1].nextnoncollinear |= 1 << 1;
      l->f[i].prevnoncollinear |= 1 << 1;
    }

    if(bl_getMapFragQ(&l->f[i]) >= bl_getMapFragP(&l->f[i+1]) && !bl_getMapFragStrand(&l->f[i])) {
      l->f[i+1].nextnoncollinear |= 1 << 2;
      l->f[i].prevnoncollinear |= 1 << 2;
    }

    if(bl_getMapFragP(&l->f[i]) <=  bl_getMapFragQ(&l->f[i+1]) && bl_getMapFragStrand(&l->f[i])) {
      l->f[i+1].nextnoncollinear |= 1 << 3;
      l->f[i].prevnoncollinear |= 1 << 3;
    }


  } else {
    Uint len = mcsa->qrylen;
    Uint thisend = bl_mcsaGet3PrimeU(mcsa);
    l->f[i].rightgap = (len - thisend) - 1;
    if(i+1 < l->n) {
      MultiCharSeqAlignment *nextmcsa = bl_getMapFragMCSA(&l->f[i+1]);
      Uint nextstart = bl_mcsaGet5PrimeU(nextmcsa);
      l->f[i+1].leftgap = nextstart;
    }
  }



  if((l->n > 1 && bl_getMapFragChrIdx(&l->f[i]) != curidx)) {
    l->consecutive = 0;
    l->P = 0;
    l->Q = 0;
  }

  l->matestatus |= 1 << mate;

  if(l->consecutive) {
    //update first and last position
    l->P = (l->P < P) ? l->P : P;
    l->Q = (l->Q > Q) ? l->Q : Q;
  } else { 
    l->P = 0;
    l->Q = 0;
  }

  if(l->n == 1) {
    l->P = P;
    l->Q = Q;
  }

  return i;
}


/*----------------------------- bl_concatMapping -----------------------------
 *    
 * @brief concat two mappings to one. Need to take care of the mem yourself!
 * @author Steve Hoffmann 
 *   
 */

mapping_t*
bl_concatMapping(mapping_t* l, mapping_t *r) {
  unsigned int i;
  mapping_t *join;

  join = ALLOCMEMORY(NULL, NULL, mapping_t, 1);
  bl_initMapping(join, l->seq, l->qual, l->lclip, l->rclip); 

  for(i=0; i < l->n; i++) {
    bl_addMapFrag(join, l->f[i].mcsa, l->f[i].seed, l->f[i].mate, l->f[i].issplit);
  }

  for(i=0; i < r->n; i++) {
    bl_addMapFrag(join, r->f[i].mcsa, r->f[i].seed, r->f[i].mate, r->f[i].issplit);
  }

  return join;
}


/*------------------------------ bl_distMapping ------------------------------
 *    
 * @brief get the distance of two mappings on the reference, check
 * consecutiveness of the alignments and mapping to the same chrom
 * @author Steve Hoffmann 
 *   
 */

long long int
bl_distMapping(mapping_t *l, mapping_t *r) {

  long long int dist = LLONG_MAX;

  if(l->n == 0 || r->n == 0) {
    return LLONG_MAX;
  } 

  if(bl_getMappingIsConsecutive(l) && bl_getMappingIsConsecutive(r)) { 
    if(bl_getMapFragChrIdx(&l->f[0]) == bl_getMapFragChrIdx(&r->f[0])) {
      if(l->P > r->Q) { 
        return l->P - r->Q; 
      } else {
        return r->Q - l->P;
      }
    }
  }

  if(bl_getMapFragChrIdx(&l->f[l->n-1]) == bl_getMapFragChrIdx(&r->f[0])) {
    dist = bl_getMapFragP(&r->f[0]) - bl_getMapFragQ(&l->f[l->n-1]);
  } 

  if(bl_getMapFragChrIdx(&r->f[r->n-1]) == bl_getMapFragChrIdx(&l->f[0])) {
    dist = MIN(dist, bl_getMapFragP(&l->f[0]) - bl_getMapFragQ(&r->f[r->n-1]));
  }

  return dist;
}


/*---------------------------- bl_isPairedMapping ----------------------------
 *    
 * @brief find out if the mapping was paired with a mate
 * @author Steve Hoffmann 
 *   
 */

char
bl_isPairedMapping (mapping_t *l) {
  return (l->matestatus == 3);
}

char
bl_isQueryMapping (mapping_t *l) {
  return (l->matestatus & 1);
}

char
bl_isMateMapping (mapping_t *l) {
  return (l->matestatus & 2);
}

/*----------------------- bl_hasMultipleQueryMappings ------------------------
 *    
 * @brief has multiple query mappings
 * @author Steve Hoffmann 
 *   
 */

  char
bl_hasMultipleQueryMappings (mappingset_t *s)
{
  Uint i;
  char ret = 0;

  for(i=0; i < s->n; i++) {
    ret += bl_isQueryMapping(&s->elem[i]);
  }

  return (ret > 1) ;
}

/*----------------------- bl_hasMultipleMateMappings ------------------------
 *    
 * @brief has multiple mate mappings
 * @author Steve Hoffmann 
 *   
 */

  char
bl_hasMultipleMateMappings (mappingset_t *s)
{
  Uint i;
  char ret = 0;

  for(i=0; i < s->n; i++) {
    ret += bl_isPairedMapping(&s->elem[i]);
  }

  return (ret > 1) ;
}

/*----------------------- bl_hasMultipleParedMappings ------------------------
 *    
 * @brief has multiple paired mappings
 * @author Steve Hoffmann 
 *   
 */

  char
bl_hasMultiplePairedMappings (mappingset_t *s)
{
  Uint i;
  char ret = 0;

  for(i=0; i < s->n; i++) {
    ret += bl_isPairedMapping(&s->elem[i]);
  }

  return (ret > 1) ;
}


/*---------------------------- bl_getMateStartIdx ----------------------------
 *    
 * @brief returns the pointer to the first mate fragment in mapping
 * @author Steve Hoffmann 
 *   
 */

  Uint
bl_getMateStartIdx (mapping_t *l)
{
  Uint i;

  for(i=0; i < l->n; i++) {
    if(bl_getMapFragIsMate(&l->f[i])) break;
  }

  if (i < l->n) return i;

  return -1;
}


/*---------------------------- bl_getMateStartPos ----------------------------
 *    
 * @brief get the start position of the mate alginment on reference
 * @author Steve Hoffmann 
 *   
 */

  Uint
bl_getMateStartPos (mapping_t *l)
{
  Uint idx;

  idx = bl_getMateStartIdx(l);

  if(idx != -1) {
    return bl_getMapFragP(&l->f[idx]);
  }

  return 0;
}

/*---------------------------- bl_getQueryStartIdx ----------------------------
 *    
 * @brief returns the pointer to the first mate fragment in mapping
 * @author Steve Hoffmann 
 *   
 */

  Uint
bl_getQueryStartIdx (mapping_t *l)
{
  Uint i;

  for(i=0; i < l->n; i++) {
    if(!bl_getMapFragIsMate(&l->f[i])) break;
  }

  if (i < l->n) return i;

  return -1;
}

/*---------------------------- bl_getQueryStartPos ----------------------------
 *    
 * @brief get the start position of the query alginment on reference
 * @author Steve Hoffmann 
 *   
 */

  Uint
bl_getQueryStartPos (mapping_t *l)
{
  Uint idx;

  idx = bl_getQueryStartIdx(l);

  if(idx != -1) {
    return bl_getMapFragP(&l->f[idx]);
  }

  return 0;
}


/*---------------------------- bl_hasQueryMapping ----------------------------
 *    
 * @brief check if mapping set contains query map
 * @author Steve Hoffmann 
 *   
 */

  char
bl_hasQueryMapping (mappingset_t *s)
{
  Uint i;

  if (s == NULL) return 0;

  for(i=0; i < s->n; i++) {
    if(bl_isQueryMapping(&s->elem[i])){  
      return 1;
    }
  }

  return 0;
}

/*---------------------------- bl_hasQueryMappingMaxEdist ----------------------------
 *    
 * @brief check if mapping set contains query map with max edist
 * @author Steve Hoffmann 
 *   
 */

  char
bl_hasQueryMappingMaxEdist (mappingset_t *s, Uint maxedist)
{
  Uint i, qedist=0, medist=0;

  if (s == NULL) return 0;


  for(i=0; i < s->n; i++) {
    if(bl_isQueryMapping(&s->elem[i])){  
      bl_getMappingEdistQM (&s->elem[i], &qedist, &medist);
      if(qedist <= maxedist)
        return 1;
    }
  }

  return 0;
}

/*---------------------------- bl_hasMateMapping -----------------------------
 *    
 * @brief check if mapping set has mate mapping
 * @author Steve Hoffmann 
 *   
 */

  char
bl_hasMateMapping (mappingset_t *s)
{
  Uint i;

  if(s == NULL) return 0;

  for(i=0; i < s->n; i++) {
    if(bl_isMateMapping(&s->elem[i])){ 
      return 1;
    } 
  }

  return 0;
}

/*---------------------------- bl_hasMateMappingMaxEdist ----------------------------
 *    
 * @brief check if mapping set contains mate map with max edist
 * @author Steve Hoffmann 
 *   
 */

  char
bl_hasMateMappingMaxEdist (mappingset_t *s, Uint maxedist)
{
  Uint i, qedist=0, medist=0;

  if (s == NULL) return 0;


  for(i=0; i < s->n; i++) {
    if(bl_isMateMapping(&s->elem[i])){  
      bl_getMappingEdistQM (&s->elem[i], &qedist, &medist);
      if(medist <= maxedist)
        return 1;
    }
  }

  return 0;
}

/*---------------------------- bl_hasMappingMaxEdist ----------------------------
 *    
 * @brief check if mapping set contains mate map with max edist
 * @author Steve Hoffmann 
 *   
 */

  char
bl_hasMappingPairedMaxEdist (mappingset_t *s, Uint maxedist)
{
  Uint i, qedist=0, medist=0;

  if (s == NULL) return 0;


  for(i=0; i < s->n; i++) {
    if(bl_isMateMapping(&s->elem[i]) && bl_isQueryMapping(&s->elem[i])){  
      bl_getMappingEdistQM (&s->elem[i], &qedist, &medist);
      if(medist+qedist <= maxedist)
        return 1;
    }
  }

  return 0;
}


/*--------------------------- bl_isSplitMappingQM ----------------------------
 *    
 * @brief 
 * @author Steve Hoffmann 
 *   
 */

  char
bl_isSplitMappingQM (mapping_t *mapping, char ismate)
{
  Uint i, ret=0;

  for(i=0; i < mapping->n; i++) {
    if(bl_getMapFragIsMate(&mapping->f[i]) == ismate || ismate ==2) {
      if(bl_getMapFragSplit(&mapping->f[i])) return 1;
    }
  }

  return ret;
}

/*--------------------------- bl_isSplitMappingQM ----------------------------
 *    
 * @brief 
 * @author Steve Hoffmann 
 *   
 */

  char
bl_getMappingStrandQM (mapping_t *mapping, char ismate)
{
  Uint i, ret=0;

  for(i=0; i < mapping->n; i++) {
    if(bl_getMapFragIsMate(&mapping->f[i]) == ismate || ismate ==2) {
      return bl_getMapFragStrand(&mapping->f[i]);
    }
  }

  return ret;
}

/*---------------------------- bl_initMappingSet -----------------------------
 *    
 * @brief init a mapping set
 * @author Steve Hoffmann 
 *   
 */

  void
bl_initMappingSet(mappingset_t *set)
{
  set->n =0;
  set->elem = NULL;

  return ;
}


/*------------------------------ bl_addMapping -------------------------------
 *    
 * @brief add a mapping to a mappingset
 * @author Steve Hoffmann 
 *   
 */

  void
bl_addMapping (mappingset_t *s, mapping_t *l)
{

  unsigned int i;
  MultiCharSeqAlignment *copy;
  s->elem = ALLOCMEMORY(NULL, s->elem, mapping_t, s->n+1);

  s->elem[s->n].f = NULL;
  s->elem[s->n].matestatus = l->matestatus;
  s->elem[s->n].consecutive = l->consecutive;
  s->elem[s->n].P = l->P;
  s->elem[s->n].Q = l->Q;
  s->elem[s->n].seq = l->seq;
  s->elem[s->n].qual = l->qual;
  s->elem[s->n].lclip = l->lclip;
  s->elem[s->n].rclip = l->rclip;
  s->elem[s->n].n = 0;
  s->elem[s->n].mapqual = l->mapqual;
  s->elem[s->n].sigma = l->sigma;
  s->elem[s->n].maxqual = l->maxqual;

  for(i=0; i < l->n; i++) {
    copy = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
    copy = bl_copyMCSA(copy, l->f[i].mcsa);
    bl_addMapFrag(&s->elem[s->n], copy, l->f[i].seed, l->f[i].mate, 
        l->f[i].issplit);
  }

  s->n++;

  return ;
}


/*--------------------------- bl_concatMappingSet ----------------------------
 *    
 * @brief concat two mapping sets into one
 * @author Steve Hoffmann 
 *   
 */

  void
bl_concatMappingSet (mappingset_t *dest, mappingset_t*source)
{
  Uint i;

  assert(dest);

  for(i=0; i < source->n; i++) {
    bl_addMapping(dest, &source->elem[i]);
  }

  return ;
}

/*------------------------------ bl_copyMappingSet ------------------------------
 *    
 * @brief get the mappings from mapping set
 * @author Steve Hoffmann 
 *   
 */

mappingset_t*
bl_copyMappingSet(mappingset_t *s) {
  Uint i;
  mappingset_t *myset;

  myset = ALLOCMEMORY(space, NULL, mappingset_t, 1);
  bl_initMappingSet(myset);

  for(i=0; i < s->n; i++) {
    bl_addMapping(myset, &s->elem[i]);
  }

  return myset;
}

/*-------------------------- bl_getMappingMinEdist ---------------------------
 *    
 * @brief get the minimum edit distance present in the set of mappings
 * @author Steve Hoffmann 
 *   
 */

unsigned int
bl_getMappingMinEdist(mappingset_t *s) {
  unsigned int i, e=0;

  for(i=0; i < s->n; i++) {
    if(i==0) {
      e = bl_getMappingEdist(&s->elem[0]);
    } else {
      e = MIN(e, bl_getMappingEdist(&s->elem[i]));
    }
  }

  return e;
}
/*-------------------------- bl_getMappingMaxScore ---------------------------
 *    
 * @brief get the minimum edit distance present in the set of mappings
 * @author Steve Hoffmann 
 *   
 */

int
bl_getMappingMaxScore(mappingset_t *s, int *scores, int indel) {
  unsigned int i;
  int e=0;

  for(i=0; i < s->n; i++) {
    if(i==0) {
      e = bl_getMappingScore(&s->elem[0], scores, indel);
    } else {
      e = MAX(e, bl_getMappingScore(&s->elem[i], scores, indel));
    }
  }

  return e;
}

/*-------------------------- bl_getMappingMaxScoreQM ---------------------------
 *    
 * @brief get the minimum edit distance present in the set of mappings
 * @author Steve Hoffmann 
 *   
 */

int
bl_getMappingMaxScoreQM(mappingset_t *s, int *scores, int indel, char ismate) {
  unsigned int i,k=0;
  int e=0, qscore=0, mscore=0;

  for(i=0; i < s->n; i++) {
    if(ismate == 2 || (!ismate && bl_isQueryMapping(&s->elem[i])) ||
        (ismate && bl_isMateMapping(&s->elem[i]))) { 
      bl_getMappingScoreQM(&s->elem[i], scores, indel, &qscore, &mscore);
      if(k==0) {
        if(ismate == 2) { 
          e = mscore + qscore;
        } else { 
          e = (ismate) ? mscore : qscore;
        }
        k++;
      } else {
        if(ismate == 2) {
          e = MAX(e, (qscore+mscore));
        } else { 
          e = MAX(e, ((ismate)? mscore : qscore));
        }
      }
    }
  }

  return e;
}

/*------------------------- bl_removeUnpairedMapping -------------------------
 *    
 * @brief remove all mappings that do not have a query/mate pair
 * @author Steve Hoffmann 
 *   
 */

void
bl_removeUnpairedMapping(mappingset_t *s) {

  unsigned int i;
  for(i=0; i < s->n; i++) {
    if(!bl_isPairedMapping(&s->elem[i])) {
      bl_removeMapping(&s->elem[i]);
      if(i+1 < s->n) { 
        memmove(&s->elem[i], &s->elem[i+1], sizeof(mapping_t)*(s->n-(i+1)));
      }
      if(s->n > 1) { 
        s->elem = ALLOCMEMORY(NULL, s->elem, mapping_t, s->n-1);
      } else {
        FREEMEMORY(NULL, s->elem);
        s->n =0;
        break;
      }
      s->n--;
      i--;
    }
  }

  return;
}

/*------------------------ bl_removeSuboptimalMapping ------------------------
 *    
 * @brief remove all mappings that are suboptimal
 * @author Steve Hoffmann 
 *   
 */

void
bl_removeSuboptimalMapping (mappingset_t *s, int *scores, int indel) {
  unsigned int i;
  int r,q;

  //e = bl_getMappingMinEdist(s);
  r = bl_getMappingMaxScore(s, scores, indel);

  for(i=0; i < s->n; i++) {
    if((q=bl_getMappingScore(&s->elem[i], scores, indel)) < r) {
      //fprintf(stdout, "removeSuboptimalMapping: removing mapping %d w score:%d (max:%d)\n",i, q, r);
      bl_removeMapping(&s->elem[i]);
      if(i+1 < s->n) { 
        memmove(&s->elem[i], &s->elem[i+1], sizeof(mapping_t)*(s->n-(i+1)));
      }
      if(s->n > 1) { 
        s->elem = ALLOCMEMORY(NULL, s->elem, mapping_t, s->n-1);
      } else {
        FREEMEMORY(NULL, s->elem);
        s->n =0;
        break;
      }
      s->n--;
      i--;
    }
  }

  return ;
}


/*------------------------ bl_removeSuboptimalMappingQM ------------------------
 *    
 * @brief remove all mappings of a certain type that are suboptimal
 * @author Steve Hoffmann 
 *   
 */

void
bl_removeSuboptimalMappingQM (mappingset_t *s, int *scores, int indel, char ismate) {
  unsigned int i;
  int r,q,qscore=0,mscore=0;

  //e = bl_getMappingMinEdist(s);
  r = bl_getMappingMaxScoreQM(s, scores, indel, ismate);

  for(i=0; i < s->n; i++) {

    bl_getMappingScoreQM(&s->elem[i], scores, indel, &qscore, &mscore);
    //fprintf(stdout, "mscore:%d, qscore:%d\n", mscore, qscore);
    q = (ismate) ? mscore : qscore;

    if(q < r) {
      //fprintf(stdout, "removeSuboptimalMappingQM: ismate:%d removing mapping %d w score:%d (max:%d)\n", ismate, i, q, r);
      bl_removeMappingQM(&s->elem[i], ismate);    
      if(s->elem[i].n == 0) { 
        if(i+1 < s->n) { 
          memmove(&s->elem[i], &s->elem[i+1], sizeof(mapping_t)*(s->n-(i+1)));
        }
        if(s->n > 1) { 
          s->elem = ALLOCMEMORY(NULL, s->elem, mapping_t, s->n-1);
        } else {
          FREEMEMORY(NULL, s->elem);
          s->n =0;
          break;
        }
        s->n--;
        i--;
      }
    }
  }

  return ;
}




/*----------------------- bl_removeErrorneousMapping ----------------------
 *    
 * @brief if one mate alignment is erroneous both are removed
 * @author Steve Hoffmann 
 *   
 */

  void
bl_removeBadMappings (mappingset_t *s, Uint querylen, Uint matelen, 
    Uint acc)
{ 

  Uint i, queryedist =0, mateedist = 0;
  Uint maxedist = querylen-ceil((acc*querylen)/100);
  Uint maxmateedist = matelen-ceil((acc*matelen)/100);

  for(i=0; i < s->n; i++) {
    bl_getMappingEdistQM(&s->elem[i], &queryedist, &mateedist);
    if(queryedist > maxedist || 
        mateedist > maxmateedist) {

      bl_removeMapping(&s->elem[i]);
      if(i+1 < s->n) { 
        memmove(&s->elem[i], &s->elem[i+1], sizeof(mapping_t)*(s->n-(i+1)));
      }
      if(s->n > 1) { 
        s->elem = ALLOCMEMORY(NULL, s->elem, mapping_t, s->n-1);
      } else {
        FREEMEMORY(NULL, s->elem);
        s->n = 0;
        break;
      }
      s->n--;
      i--;
    }
  }

  return ;
}

/*----------------------- bl_removeBadMatesAcc ----------------------
 *    
 * @brief removes only erroneous mate or both
 * @author Steve Hoffmann 
 *   
 */

  void
bl_removeBadMatesAcc (mappingset_t *s, Uint querylen, Uint matelen, Uint acc)
{ 

  Uint i, queryedist =0, mateedist = 0;
  Uint maxedist = querylen-ceil((acc*querylen)/100);
  Uint maxmateedist = matelen-ceil((acc*matelen)/100);

  for(i=0; i < s->n; i++) {
    bl_getMappingEdistQM(&s->elem[i], &queryedist, &mateedist);
    //fprintf(stdout, "queryedist:%d (max:%d), mateedist:%d (max:%d)\n", queryedist, maxedist, mateedist, maxmateedist);

    if( (queryedist > maxedist && mateedist > maxmateedist) 
        || (mateedist > maxmateedist && !bl_isQueryMapping(&s->elem[i])) 
        || (queryedist > maxedist && !bl_isMateMapping(&s->elem[i]))) {

      //fprintf(stdout, "BadmatesAcc: removing mapping %d: qedist:%d, mateedist:%d, maxedist:%d\n", i,queryedist, mateedist, maxedist);
      bl_removeMapping(&s->elem[i]);
      if(i+1 < s->n) { 
        memmove(&s->elem[i], &s->elem[i+1], sizeof(mapping_t)*(s->n-(i+1)));
      }
      if(s->n > 1) { 
        s->elem = ALLOCMEMORY(NULL, s->elem, mapping_t, s->n-1);
      } else {
        FREEMEMORY(NULL, s->elem);
        s->n = 0;
        break;
      }
      s->n--;
      i--;
    } else

      if(queryedist > maxedist) {
        bl_removeMappingQM(&s->elem[i], 0);    
      } else

        if(mateedist > maxmateedist) {
          bl_removeMappingQM(&s->elem[i], 1);    
        }
  }

  return;
}


/*----------------------- bl_removeBadMatesScr ----------------------
 *    
 * @brief removes only low scoring mate or both
 * @author Steve Hoffmann 
 *   
 */

  void
bl_removeBadMatesScr (mappingset_t *s, Uint querylen, Uint matelen, int *scores, int indel, int minscr, int minmatescr)
{ 

  Uint i;
  int queryscr =0, matescr = 0;

  for(i=0; i < s->n; i++) {
    bl_getMappingScoreQM(&s->elem[i], scores, indel, &queryscr, &matescr);

    if( (queryscr < minscr && matescr <  minmatescr) 
        || (matescr <  minmatescr && !bl_isQueryMapping(&s->elem[i])) 
        || (queryscr <  minscr && !bl_isMateMapping(&s->elem[i]))) {

      bl_removeMapping(&s->elem[i]);
      if(i+1 < s->n) { 
        memmove(&s->elem[i], &s->elem[i+1], sizeof(mapping_t)*(s->n-(i+1)));
      }
      if(s->n > 1) { 
        s->elem = ALLOCMEMORY(NULL, s->elem, mapping_t, s->n-1);
      } else {
        FREEMEMORY(NULL, s->elem);
        s->n = 0;
        break;
      }
      s->n--;
      i--;
    } else

      if(queryscr < minscr) {
        bl_removeMappingQM(&s->elem[i], 0);    
      } else

        if(matescr <  minmatescr) {
          bl_removeMappingQM(&s->elem[i], 1);    
        }
  }

  return;
}



/*----------------------- bl_removeBadMatesCov ----------------------
 *    
 * @brief removes only erroneous mate or both
 * @author Steve Hoffmann 
 *   
 */

  void
bl_removeBadMatesCov (mappingset_t *s, Uint querylen, Uint matelen, Uint mc)
{ 

  Uint i, querycover =0, matecover = 0;
  Uint minmc = ceil((mc*querylen)/100);
  Uint minmatemc = ceil((mc*matelen)/100);

  for(i=0; i < s->n; i++) {

    bl_getMappingCovQM(&s->elem[i], &querycover, &matecover);

    if ( (querycover < minmc && matecover < minmatemc) 
        || (matecover < minmatemc && !bl_isQueryMapping(&s->elem[i])) 
        || (querycover < minmc && !bl_isMateMapping(&s->elem[i]))) {

      //fprintf(stdout, "BadmatesCov: removing mapping %d: querycover:%d, matecover:%d, mincov:%d minmatecov:%d\n", 
      //     i, querycover, matecover, minmc, minmatemc);
      bl_removeMapping(&s->elem[i]);
      if(i+1 < s->n) { 
        memmove(&s->elem[i], &s->elem[i+1], sizeof(mapping_t)*(s->n-(i+1)));
      }
      if(s->n > 1) { 
        s->elem = ALLOCMEMORY(NULL, s->elem, mapping_t, s->n-1);
      } else {
        FREEMEMORY(NULL, s->elem);
        s->n = 0;
        break;
      }
      s->n--;
      i--;
    } else

      if(querycover < minmc) {
        bl_removeMappingQM(&s->elem[i], 0);    
      } else

        if(matecover < minmatemc) {
          bl_removeMappingQM(&s->elem[i], 1);    
        }
  }

  return;
}


/*--------------------- bl_removeBadMatesQMEdist ---------------------
 *    
 * @brief removes only erroneous mate or both
 * @author Steve Hoffmann 
 *   
 */

  void
bl_removeBadMatesEdistQM (mappingset_t *s, Uint maxedist, Uint ismate)
{ 

  Uint i, queryedist =0, mateedist = 0;

  for(i=0; i < s->n; i++) {
    bl_getMappingEdistQM(&s->elem[i], &queryedist, &mateedist);
    if((!ismate && mateedist > maxedist && !bl_isQueryMapping(&s->elem[i])) 
        ||(ismate && queryedist > maxedist && !bl_isMateMapping(&s->elem[i]))) {
      //fprintf(stderr, "removing mapping %d: qedist:%d, mateedist:%d, maxedist:%d\n", i,queryedist, mateedist, maxedist);
      bl_removeMapping(&s->elem[i]);
      if(i+1 < s->n) { 
        memmove(&s->elem[i], &s->elem[i+1], sizeof(mapping_t)*(s->n-(i+1)));
      }
      if(s->n > 1) { 
        s->elem = ALLOCMEMORY(NULL, s->elem, mapping_t, s->n-1);
      } else {
        FREEMEMORY(NULL, s->elem);
        s->n = 0;
        break;
      }
      s->n--;
      i--;
    } else

      if(!ismate && queryedist > maxedist) {
        bl_removeMappingQM(&s->elem[i], 0);    
      } else

        if(ismate && mateedist > maxedist) {
          bl_removeMappingQM(&s->elem[i], 1);    
        }
  }

  return;
}





/*------------------------- bl_countMultipleMappings -------------------------
 *    
 * @brief count the number of querys and mates in the alignment
 * @author Steve Hoffmann 
 *   
 */

  void
bl_countMultipleMappings (mappingset_t *s, Uint *nqueries, Uint *nmates)
{
  Uint i, noofmates =0, noofqueries = 0;

  for(i=0; i < s->n; i++) {
    if(bl_isQueryMapping(&s->elem[i])) noofqueries++; 
    if(bl_isMateMapping(&s->elem[i])) noofmates++;
  }

  *nqueries = noofqueries;
  *nmates = noofmates;


  return ;
}

/*---------------------------- bl_initMapSeedList ----------------------------
 *    
 * @brief initialize the list of mapseed_t
 * @author Steve Hoffmann 
 *   
 */

  mapseedlist_t*
bl_initMapSeedList (mapseedlist_t *l)
{
  l->n = 0;
  l->l = NULL;

  return l;
}


/*------------------------------ bl_getBestSeed ------------------------------
 *    
 * @brief get the best seed from the seed list
 * @author Steve Hoffmann 
 *   
 */

  mapseed_t*
bl_getMapSeedListBest (mapseedlist_t *l)
{
  Uint i;
  mapseed_t *best = NULL;
  Uint max = 0;

  for(i=0; i < l->n; i++) { 
    if(i == 0 || l->l[i].mat > max) { 
      best = &l->l[i];
      max = l->l[i].mat;
    }
  }

  return best;
}


/*--------------------------- bl_getMapSeedChrIdx ----------------------------
 *    
 * @brief get the chromosome/contig index for mapseed i and locus j
 * @author Steve Hoffmann 
 *   
 */

  locus_t *
bl_getMapSeedLocus (mapseed_t *seed, Uint j, MultiCharSeq *mseq, Suffixarray *s)
{

  locus_t *locus;
  char rc;
  uint64_t pos, len;


  assert(seed->l + j <= seed->r);

  pos = s->suftab[seed->l+j];
  rc = seed->u;
  len = seed->mat + seed->mis + seed->ins;

  locus = ALLOCMEMORY(NULL, NULL, locus_t, 1);
  locus = bl_initLocus(locus, mseq, pos, rc, len, 0);

  return locus;
}

/*-------------------- bl_getMapSeedAdditionalInformation --------------------
 *    
 * @brief stores additional map frag information for seeds (for output)
 * @author Steve Hoffmann 
 *   
 */

  void
bl_getMapSeedAdditionalInformation (mapseed_t *seed, Suffixarray *s, MultiCharSeq *mseq)
{

  /*best query seed*/

  locus_t *locus = bl_getMapSeedLocus(seed, 0, mseq, s);
  seed->refpos = bl_getLocusChromPos(locus);
  seed->refname = bl_getLocusChromDesc(mseq, locus);
  seed->refidx = bl_getLocusChromIdx(locus); 
  seed->seedlen = seed->mat+seed->mis+seed->ins;


  FREEMEMORY(NULL, locus);
  return ;
}
/*--------------------------- bl_branch2mapseed ---------------------------
 *    
 * @brief conversion of types. matchstem branch to seed. Should be replaced
 * @author Steve Hoffmann 
 *   
 */

  mapseedlist_t*
bl_addMapSeedBranch (mapseedlist_t *s, branch_t *b, unsigned int u, unsigned int strand, 
    double maxE, unsigned int maxM, double SPM, karlin_t *stats)
{
  int scr; double E;

  if(b->l > b->r) return s;

  scr = b->mat - (b->mis + b->ins + b->del);
  E = evalue(stats->lambda, stats->K, SPM, scr);

  s->l = ALLOCMEMORY(NULL, s->l, mapseed_t, s->n+1);
  s->l[s->n].good = 1;
  s->l[s->n].maxinterval = 0;
  s->l[s->n].maxevalue = 0;


  if((b->r - b->l) > maxM) {
    s->l[s->n].good = 0;
    s->l[s->n].maxinterval = 1;
  } 

  if(E > maxE) {
    s->l[s->n].good = 0; 
    s->l[s->n].maxevalue = 1;
  }


  s->l[s->n].l = b->l;
  s->l[s->n].r = b->r;
  s->l[s->n].mat = b->mat;
  s->l[s->n].mis = b->mis;
  s->l[s->n].ins = b->ins;
  s->l[s->n].del = b->del;
  s->l[s->n].edist = b->mis + b->ins + b->del;
  s->l[s->n].evalue = E;
  s->l[s->n].score = scr;
  s->l[s->n].u = u;
  s->l[s->n].rc = strand; 

  s->n++;

  return s;
}


/*---------------------------- bl_mappingFindPos -----------------------------
 *    
 * @brief find a mapping that starts at a given position
 * @author Steve Hoffmann 
 *   
 */

  unsigned int
bl_mappingsetHasPos (mappingset_t *s, unsigned int pos)
{
  unsigned int i, j;

  for(i=0; i < s->n; i++) {
    for(j=0; j < s->elem[i].n; j++) {
      if(bl_getMapFragP(&s->elem[i].f[j]) == pos) {
        return 1;
      }
    }
  }

  return 0;
}


/*-------------------------- bl_mappingsetHasPaired --------------------------
 *    
 * @brief check if there are paired mappnigs in a set
 * @author Steve Hoffmann 
 *   
 */

  char
bl_mappingsetHasPaired (mappingset_t *s)
{

  unsigned int i;
  for(i=0; i < s->n; i++) {
    if(bl_isPairedMapping (&s->elem[i]))
      return 1;
  }

  return 0;
}


/*----------------------- bl_mappingsetHasQueryMatches -----------------------
 *    
 * @brief check if there are query matches in the mapping set
 * @author Steve Hoffmann 
 *   
 */

  char
bl_mappingsetHasQueryMatches (mappingset_t *s)
{

  unsigned int i;
  for(i=0; i < s->n; i++) {
    if(bl_isQueryMapping (&s->elem[i])){
      return 1;
    } 
  }

  return 0;
}

/*----------------------- bl_mappingsetHasMateMatches -----------------------
 *    
 * @brief check if there are mate matches in the mapping set
 * @author Steve Hoffmann 
 *   
 */

  char
bl_mappingsetHasMateMatches (mappingset_t *s)
{

  unsigned int i;
  for(i=0; i < s->n; i++) {
    if(bl_isMateMapping (&s->elem[i]))
      return 1;
  }

  return 0;
}

/*--------------------------- bl_findNextMapFragU ---------------------------
 *    
 * @brief given fragment k, find the next fragment wrt to query
 * @author Steve Hoffmann 
 *              |------------------|
 *     |--------------------|
 */

  Uint
bl_getNextMapFragU (mapfrag_t *f, mapping_t *m, char ismate, Uint *lgap, Uint *rgap)
{
  unsigned int i, fstart, fend, len, gstart, gend, tmp=-1;
  unsigned int curstart = -1, curend = 0;
  char endset = 0;

  MultiCharSeqAlignment *mcsa;

  mcsa = bl_getMapFragMCSA(f);
  fstart = bl_mcsaGet5PrimeU(mcsa);
  fend = bl_mcsaGet3PrimeU(mcsa);
  len = mcsa->qrylen;
  *rgap = 0;
  *lgap = fstart;


  for(i=0; i < m->n; i++) {

    mcsa = bl_getMapFragMCSA(&m->f[i]);
    gstart = bl_mcsaGet5PrimeU(mcsa);
    gend = bl_mcsaGet3PrimeU(mcsa);

    if (gstart > fstart &&  curstart > gstart && bl_getMapFragIsMate(&m->f[i]) == ismate){ 
      tmp = i;
      curstart = gstart;
      endset = 1;
      if(gstart > fend) {
        *rgap = (gstart-fend)-1;
      }
    }

    if(fstart > gstart && gend >= curend && bl_getMapFragIsMate(&m->f[i]) == ismate) {
      curend = gend;
      if(fstart > gend) {
        *lgap = (fstart-gend)-1;
      }
    }
  }

  //sorted order
  if(!endset && len > fend) { 
    *rgap = (len-fend)-1;
  }


  return tmp;
}

locuslist_t*
bl_getMappingLocusList(mapping_t *mapping, MultiCharSeq* mseq, char ismate) {
  Uint i, k=0;
  locuslist_t *mylist = NULL;
  MultiCharSeqAlignment *mcsa=NULL;

  for(i=0; i < mapping->n; i++) {
    if(bl_getMapFragIsMate(&mapping->f[i]) == ismate || ismate == 2) { 
      mcsa = ALLOCMEMORY(NULL, mcsa, MultiCharSeqAlignment, k+1);
      memmove(&mcsa[k], mapping->f[i].mcsa, sizeof(MultiCharSeqAlignment));
      k++;
    }
  }

  //fprintf(stderr, "initialized %d mcsas in getMappingLocusList\n", k);
  mylist = bl_getLocusList (mcsa, mseq, k);

  FREEMEMORY(space, mcsa);

  return mylist;
}

char
bl_isCollinearMapping(mapping_t *m, char ismate) {
  Uint i=0;
  Uint result =1;

  if(!m || m->f == NULL || m->n == 0) 
    return 1;

  //fast forward to mate positions
  while(ismate && i < m->n-1 && !bl_getMapFragIsMate(&m->f[i])) i++;

  Uint rc = bl_getMapFragStrand(&m->f[i]);
  Uint chr =  bl_getMapFragChrIdx(&m->f[i]);
  Uint pos = bl_getMapFragP(&m->f[i]);
  i++;

  while (i < m->n-1 && bl_getMapFragIsMate(&m->f[i]) == ismate) { 
    //same mate different strand
    if(rc != bl_getMapFragStrand(&m->f[i])){ 
      return 0;
    }
    //different chromosome
    if(chr != bl_getMapFragChrIdx(&m->f[i])){ 
      return 0;
    }
    //cirucular
    if((rc && pos <  bl_getMapFragP(&m->f[i])) || 
        (!rc && pos > bl_getMapFragP(&m->f[i]))) { 
      return 0;
    }

    rc = bl_getMapFragStrand(&m->f[i]);
    chr =  bl_getMapFragChrIdx(&m->f[i]);
    pos = bl_getMapFragP(&m->f[i]);
    i++;
  }

  return result;
}


/*----------------------------- bl_dumpMapFrags ------------------------------
 *    
 * @brief show the mapfrags
 * @author Steve Hoffmann 
 *   
 */

  void
bl_dumpMapFrags (FILE *dev, mapfrag_t *f)
{

  MultiCharSeqAlignment *mcsa;
  Uint left, right, chridx, lgap, rgap;
  uint64_t pos, end;
  char rc,mate;

  mcsa  = bl_getMapFragMCSA(f);
  left  = bl_getMapFragLeft(f);
  right = bl_getMapFragRight(f);
  rc = bl_getMapFragStrand(f);
  chridx = bl_getMapFragChrIdx(f);
  pos = bl_getMapFragP(f);
  end = bl_getMapFragQ(f);
  mate = bl_getMapFragIsMate(f);
  lgap = f->leftgap;
  rgap = f->rightgap;


  fprintf(dev, "]-%d-[%d,%d]-%d-[ (mate:%d) -> %d-[%"PRIu64",%"PRIu64"](rc:%d)\n", lgap, left, right, rgap, (int)mate, chridx, pos, end, (int)rc);
  showAlign(mcsa->al, dev);

  return ;
}


/*------------------------------ bl_dumpMapping ------------------------------
 *    
 * @brief dump a mapping to a device
 * @author Steve Hoffmann 
 *   
 */

  void
bl_dumpMapping (FILE *dev, mapping_t *m)
{
  Uint i;

  fprintf(dev, "fragments:%d, consecutive:%d, matestatus:%d, [%d,%d] (range:%d)\n",
      m->n, m->consecutive, (int) m->matestatus, m->P, m->Q, bl_getMappingRange(m));

  for(i=0; i < m->n; i++) {
    fprintf(dev, "\t\t%d\t", i);
    bl_dumpMapFrags(dev, &m->f[i]);
  }
  return ;
}


/*---------------------------- bl_dumpMappingSet -----------------------------
 *    
 * @brief dump the mapping set to a device
 * @author Steve Hoffmann 
 *   
 */

  void
bl_dumpMappingSet (FILE *dev, mappingset_t *set)
{
  Uint i;

  fprintf(dev, "dumping mapping set with %d mappings.\n", set->n);

  for(i=0; i < set->n; i++) {
    fprintf(dev, "mapping %d\t", i);
    bl_dumpMapping(dev, &set->elem[i]);
  }
  return ;
}


/*--------------------------- bl_cmpMappingScores ----------------------------
 *    
 * @brief compare mapping scores
 * @author Steve Hoffmann 
 *   
 */

  int
bl_cmpMappingScores (const void *a, const void *b)
{
  mapping_t *first = (mapping_t*) a;
  mapping_t *secnd = (mapping_t*) b;

  if(first->scr > secnd->scr) return -1;
  if(first->scr <  secnd->scr) return 1;

  return 0;
}

/*------------------------- bl_sortMappingSetByScore -------------------------
 *    
 * @brief sort the mapping set by score
 * @author Steve Hoffmann 
 *   
 */

  mappingset_t*
bl_sortMappingSetByScore (mappingset_t *set, int* scores, int indel)
{
  Uint i;
  for(i=0; i < set->n; i++) {
    set->elem[i].scr = bl_getMappingScore (&set->elem[i], scores, indel);
  }

  qsort(set->elem, set->n, sizeof(mapping_t), bl_cmpMappingScores);

  return set;
}


/*----------------------------- bl_mergeMappings -----------------------------
 *    
 * @brief merge mappings if there is exactly one mate and one query
 * @author Steve Hoffmann 
 *   
 */

  void
bl_mergeMappings (mappingset_t *set)
{
  int i;
  
  MultiCharSeqAlignment *copy;


  //make sure that above condition is fulfilled
  if(set->n != 2) return;
  if(set->elem[0].matestatus == set->elem[1].matestatus) return;
  if(set->elem[0].matestatus == 3 || set->elem[1].matestatus == 3) return;


  /*assert the proper pairing in +/-*/
 /* 
  unsigned char strand=255;
  for(i=0; i < set->elem[0].n; i++) {
    if(i == 0) {
      strand = set->elem[0].f[i].mcsa->strand;
    } else {
      if(strand != set->elem[0].f[i].mcsa->strand) return;
    }
  }

  if(strand != 255) {
 
    for(i=0; i < set->elem[1].n; i++) {
        if(strand == set->elem[1].f[i].mcsa->strand) return;
    }
  }
 */
  
  /*copy elems*/
  for(i=0; i < set->elem[1].n; i++) {
    copy = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
    copy = bl_copyMCSA(copy, set->elem[1].f[i].mcsa);
    bl_addMapFrag(&set->elem[0], copy, NULL, set->elem[1].f[i].mate, 
        set->elem[1].f[i].issplit);
  }

  bl_removeMapping (&set->elem[1]);
  set->elem = ALLOCMEMORY(NULL, set->elem, mapping_t, 1);
  set->n = 1;

  return ;
}



mappingset_t* 
bl_mappingsetRemoveDuplicates(mappingset_t *set, Uint maxdist) {
  Uint i,k,n=0;
  char remove = 0, hasquery, hasmate;
  Uint querystart, matestart;
  mappingset_t *new;
  PairUint *poslist;

  new = ALLOCMEMORY(NULL, NULL, mappingset_t, 1);
  bl_initMappingSet(new);

  poslist = ALLOCMEMORY(NULL, NULL, PairUint, set->n);
  memset(poslist, 0, sizeof(PairUint)*set->n);

  for(i=0; i < set->n; i++) {
    hasquery = bl_isQueryMapping(&set->elem[i]);
    hasmate = bl_isMateMapping(&set->elem[i]);
    querystart = bl_getQueryStartPos(&set->elem[i]);
    matestart = bl_getMateStartPos(&set->elem[i]);
    remove = 0;

    for(k=0; k < n; k++) { 
      if((DIST(poslist[k].a, querystart) < maxdist || !hasquery)
          && (DIST(poslist[k].b, matestart) < maxdist || !hasmate)) {
        remove =1;
      } 
    }

    if(!remove) {
      bl_addMapping (new, &set->elem[i]); 
    }

    poslist[n].a = querystart;
    poslist[n].b = matestart;
    n++;
  }

  bl_removeMappingSet(set);
  FREEMEMORY(space, set);
  FREEMEMORY(NULL, poslist);

  return new;
}



/*---------------------- bl_dumpSpliceJunctions -------------------------
 *    
 * @brief return the regular splice junctions
 * @author Steve Hoffmann 
 *   
 */

  void
bl_dumpSpliceJunctions( mapping_t *m, char ismate, MultiCharSeq *mseq, 
    char *basename, char *queryname, segemehl_t *nfo)
{

  mapfrag_t *frags;
  MultiCharSeqAlignment *mcsa=NULL;
  locuslist_t *list = NULL;
  mapping_t* newmap;
  multilocus_t *mult, *nvrt;
  Uint i, j, nooffrags=0;
  Uint noofparts=0;
  Uint noofmultiloci = 0;
  char *item;
  char *queryid;
  Uint queryidlen=strlen(queryname);
  uint8_t score;
  double phred = 0;
  char sep = ';';

  //needs to take care of joint (N separted and single mapfrags)
  queryid = strclip(NULL, queryname, &queryidlen);
  item = ALLOCMEMORY(NULL, NULL, char, strlen(basename)+4+queryidlen+3);
  memmove(item, basename, strlen(basename));
  item[strlen(basename)] = sep;

  //  list = ALLOCMEMORY(NULL, NULL, locuslist_t, 1);
  //  bl_initLocuslist(list);
  frags = bl_getMapFragsQM(m, &nooffrags, ismate);
  newmap= ALLOCMEMORY(NULL, NULL, mapping_t, 1);
  bl_initMapping(newmap, m->seq, m->qual, m->lclip, m->rclip);


  for(i=0; i < nooffrags; i++) {
    MultiCharSeqAlignment *tmpfrag = bl_getMapFragMCSA(&frags[i]);
    mcsa = bl_getPartialMultiCharSeqAlignments(tmpfrag, mseq, &noofparts);

    if(i==0) { 
      phred = MAX(m->mapqual, frags[i].mapq_dbl) * -4.34294481903252;
    } else { 
      phred = MIN(phred, (MAX(m->mapqual, frags[i].mapq_dbl) * -4.34294481903252)); 
    }

    for(j=0; j < noofparts; j++) {
      MultiCharSeqAlignment *copy = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
      copy = bl_copyMCSA(copy, &mcsa[j]);
      bl_addMapFrag(newmap, copy, frags[i].seed, frags[i].mate, 
          frags[i].issplit);
      wrapMultiCharSeqAlignment(NULL, &mcsa[j]);
    }
    FREEMEMORY(NULL, mcsa);
  }

  FREEMEMORY(NULL, frags);

  if(phred < 1) { 
    score = 1;
  } else if (phred >= 60) {
    score = 60;
  } else {
    score = (uint8_t) phred;
  }

  list = bl_getMappingLocusList(newmap, mseq, ismate);
  mult = bl_createMultiLocusFromLocusList(list, &noofmultiloci);

  for(i=0; i < noofmultiloci; i++) {

    bl_setMultiLocusScore(&mult[i], score);
    bl_setMultiLocusName(&mult[i], basename);

    if(mult[i].noofloci > 1) { 

      nvrt = bl_invertMultiLocus (&mult[i]);


      item[strlen(basename)+1] = 'R';
      item[strlen(basename)+2] = sep;
      memmove(&item[strlen(basename)+3], queryid, queryidlen);
      item[strlen(basename)+3+queryidlen] = sep;
      item[strlen(basename)+3+queryidlen+1] = (ismate) ? '2' : '1' ;
      item[strlen(basename)+3+queryidlen+2] = '\0';




      if (nfo->threadno > 1) {
        pthread_mutex_lock(nfo->mtx6);
      }

      bl_dumpMultiLocusSingle(nfo->singlesplitdev, nvrt, mseq, item);

      if (nfo->threadno > 1) {
        pthread_mutex_unlock(nfo->mtx6);
      }


      if (nfo->threadno > 1) {
        pthread_mutex_lock(nfo->mtx7);
      }

      bl_dumpMultiLocusJoint(nfo->multisplitdev, nvrt, mseq, item);

      if (nfo->threadno > 1) {
        pthread_mutex_unlock(nfo->mtx7);
      }



      bl_wrapMultiLocus(nvrt);
      FREEMEMORY(NULL, nvrt);
    } else {

      if(i > 0) {


        if(mult[i].strand == mult[i-1].strand && mult[i].idx == mult[i-1].idx 
            && DIST(mult[i].pos,mult[i-1].pos) < 20000) {

          if(!mult[i-1].strand && mult[i-1].pos + mult[i-1].len - 1 > mult[i].pos) {

            if(mult[i].pos+mult[i].len-1 < mult[i-1].pos) { 
              item[strlen(basename)+1] = 'C';
            } else { 
              item[strlen(basename)+1] = 'B';
            }

            item[strlen(basename)+2] = sep;
            memmove(&item[strlen(basename)+3], queryid, queryidlen);
            item[strlen(basename)+3+queryidlen] = sep;
            item[strlen(basename)+3+queryidlen+1] = (ismate) ? '2' : '1' ;
            item[strlen(basename)+3+queryidlen+2] = '\0';

            char *tmp;

            bl_asprintf(&tmp, "%s\t%"PRIu64"\t%"PRIu64"\t%s\t%d\t%c\n", 
                ((CharSequence*)mseq->ref[mult[i].idx].ref)->description,
                mult[i].pos-mult[i].chrstart,                              //start and end reversed here
                mult[i-1].pos+mult[i-1].len-1 - mult[i-1].chrstart+1, //+2
                item, MIN(mult[i-1].score, mult[i].score),
                ((mult[i-1].strand) ? '-' : '+')
                );



            if (nfo->threadno > 1) {
              pthread_mutex_lock(nfo->mtx6);
            }

            fprintf(nfo->singlesplitdev, "%s", tmp);

            if (nfo->threadno > 1) {
              pthread_mutex_unlock(nfo->mtx6);
            }




            FREEMEMORY(NULL, tmp);

          }

          if(mult[i-1].strand && mult[i-1].pos < mult[i].pos + mult[i].len - 1) {


            if(mult[i-1].pos+mult[i-1].len-1 < mult[i].pos) { 
              item[strlen(basename)+1] = 'C';
            } else { 
              item[strlen(basename)+1] = 'B';
            }

            item[strlen(basename)+2] = sep;
            memmove(&item[strlen(basename)+3], queryid, queryidlen);
            item[strlen(basename)+3+queryidlen] = sep;
            item[strlen(basename)+3+queryidlen+1] = (ismate) ? '2' : '1' ;
            item[strlen(basename)+3+queryidlen+2] = '\0';

            char *tmp;

            bl_asprintf(&tmp, "%s\t%"PRIu64"\t%"PRIu64"\t%s\t%d\t%c\n", 
                ((CharSequence*)mseq->ref[mult[i].idx].ref)->description,
                mult[i-1].pos-mult[i-1].chrstart,                    //+2
                mult[i].pos+mult[i].len-1 - mult[i-1].chrstart+1,        //+0
                item, MIN(mult[i-1].score, mult[i].score),
                ((mult[i-1].strand) ? '-' : '+') 
                );


            if (nfo->threadno > 1) {
              pthread_mutex_lock(nfo->mtx6);
            }

            fprintf(nfo->singlesplitdev, "%s", tmp);

            if (nfo->threadno > 1) {
              pthread_mutex_unlock(nfo->mtx6);
            }





            FREEMEMORY(NULL, tmp);

          }
        } else {
          char *tmp = NULL;

          bl_bsprintf(&tmp, "%s,%"PRIu64",%c,%u,%"PRIu64",%d,%d\t", 
              ((CharSequence*)mseq->ref[mult[i-1].idx].ref)->description, 
              mult[i-1].pos - mult[i-1].chrstart+1,
              ((mult[i-1].strand) ? '-' : '+'), 
              mult[i-1].readstart+1,
              mult[i-1].len, 
              mult[i-1].edist,
              mult[i-1].score);

          bl_bsprintf(&tmp, "%s,%"PRIu64",%c,%u,%"PRIu64",%d,%d\t", 
              ((CharSequence*)mseq->ref[mult[i].idx].ref)->description, 
              mult[i].pos - mult[i].chrstart+1,
              ((mult[i].strand) ? '-' : '+'), 
              mult[i].readstart+1,
              mult[i].len, 
              mult[i].edist,
              mult[i].score);

          bl_bsprintf(&tmp, "%s\n", queryid);  




          if (nfo->threadno > 1) {
            pthread_mutex_lock(nfo->mtx8);
          }

          fprintf(nfo->transsplitdev, "%s", tmp);

          if (nfo->threadno > 1) {
            pthread_mutex_unlock(nfo->mtx8);
          }





          FREEMEMORY(NULL, tmp);
        }
      }
    }
  }

  for(i=0; i < noofmultiloci; i++) {
    bl_wrapMultiLocus(&mult[i]);
  }

  FREEMEMORY(NULL, mult);
  bl_removeMapping(newmap);
  FREEMEMORY(NULL, newmap);
  wrapLocuslist(list);
  FREEMEMORY(NULL, list);
  FREEMEMORY(NULL, item);
  FREEMEMORY(NULL, queryid);

  return;
}


/*--------------------------- bl_isCircularMapping ---------------------------
 *    
 * @brief check if mapping contains a circle
 * @author Steve Hoffmann 
 *   
 */

  char
bl_mappingGetType (mapping_t *m, char ismate)
{

  Uint i, idx0, idx1;
  char frc0, frc1;
  uint64_t pos0, pos1;
  uint64_t len0, len1;
  mapfrag_t* frags;
  Uint nooffrags = 0;
  char type = 0;

  frags = bl_getMapFragsQM (m, &nooffrags, ismate);

  for(i=1; i < nooffrags; i++) {

    idx0 = bl_getMapFragChrIdx(&frags[i-1]);
    idx1 = bl_getMapFragChrIdx(&frags[i]);
    frc0 = bl_getMapFragStrand(&frags[i-1]);
    frc1 = bl_getMapFragStrand(&frags[i]);
    pos0 = bl_getMapFragP(&frags[i-1]);
    pos1 = bl_getMapFragP(&frags[i]);
    len0 = bl_getMapFragGetVlen(&frags[i-1]);
    len1 = bl_getMapFragGetVlen(&frags[i]);

    if(idx0 == idx1 && frc0 == frc1) { 

      if(!frc0 && pos0 + len0 -1 > pos1) {
        if (pos1+len1-1 < pos0) {
          type |= 1 << 0; //circular        
        } else {
          type |= 1 << 1;
        }
      }

      if(frc0 && pos0 < pos1 + len1 -1) {
        if(pos0+len0-1 < pos1) { 
          type |= 1 << 0;
        } else {
          type |= 1 << 1;
        }

        type |= 1 << 0; //circular
      }

    } else {
      if(idx0 == idx1) { 
        type |= 1 << 2;
      } else {
        type |= 1 << 3;
      }
    }
  }

  FREEMEMORY(NULL, frags);
  return type;
}
