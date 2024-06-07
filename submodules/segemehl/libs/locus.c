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
 *  locus.c
 *  a genomic locus given a reference
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 09.12.2012 14:15:53 CET
 *  
 */

#include "locus.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>

#include "basic-types.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "memory.h"
#include "stringutils.h"
#include "alignment.h"
#include "nw.h"
#include "sw.h"


/*------------------------------- bl_initLocus ------------------------------
 *    
 * @brief inits a locus structure on a multicharseq with global(!) position
 * @author Steve Hoffmann 
 *   
 */ 

locus_t*
bl_initLocus (locus_t* loc, MultiCharSeq *mseq, uint64_t pos,
    unsigned char strand, uint64_t length, uint64_t loff) 
{
  unsigned int chrstart, chrend;// lcut, rcut;

  loc->idx = getMultiCharSeqIndex(mseq, &mseq->sequences[pos]);
  getMultiCharSeqIdxBounds(mseq, loc->idx, &chrstart, &chrend);
  assert(length > 0);

  //left boundary check
  if(pos < loff+chrstart) {
    loc->pos = chrstart;
//    lcut = loff+chrstart-pos;
  } else {
    loc->pos = pos-loff;
//    lcut = 0;
  }

  //right boundary check
  if(loc->pos + length >= chrend) {
    loc->len = chrend-loc->pos+1;
//    rcut = loc->pos+length-chrend+1;
  } else {
    loc->len = length;
  }

  loc->chrstart = chrstart;
  loc->chrend = chrend;
  loc->strand = strand;
  loc->score = 0;
  loc->name = NULL;
  loc->readstart = 0;
  loc->readend = 0;
  loc->edist = 0;

  return loc;
}


/*----------------------------- bl_locusSetName ------------------------------
 *    
 * @brief set the name of a locus
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_locusSetName (locus_t *loc, char *name)
{
    loc->name = bl_strdup(name);
	return ;
}


/*----------------------------- bl_locusSetScore -----------------------------
 *    
 * @brief set the score
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_locusSetScore (locus_t *loc, uint8_t score)
{
   loc->score = score;
   return ;
}

/*------------------------------ bl_getLocusPos ------------------------------
 *    
 * @brief important: this function returns the global(!) position
 * @author Steve Hoffmann 
 *   
 */
 
uint64_t
bl_getLocusPos (locus_t *loc)
{
	return loc->pos;
}

/*------------------------------ bl_getMultiLocusPos ------------------------------
 *    
 * @brief important: this function returns the global(!) position
 * @author Steve Hoffmann 
 *   
 */
 
uint64_t
bl_getMultiLocusPos (multilocus_t *m)
{
	return m->pos;
}



/*--------------------------- bl_getLocusChromPos ----------------------------
 *    
 * @brief important: this function returns the local(!) position, 
 * ie. the pos on the chrom
 * @author Steve Hoffmann 
 *   
 */
 
uint64_t
bl_getLocusChromPos (locus_t *loc)
{
	return loc->pos-loc->chrstart;
}

/*--------------------------- bl_getMultiLocusChromPos ----------------------------
 *    
 * @brief important: this function returns the local(!) position, 
 * ie. the pos on the chrom
 * @author Steve Hoffmann 
 *   
 */
 
uint64_t
bl_getMultiLocusChromPos (multilocus_t *m)
{
	return m->pos-m->chrstart;
}

/*------------------------------ bl_getLocusPos ------------------------------
 *    
 * @brief important: this function returns the global(!) position
 * @author Steve Hoffmann 
 *   
 */
 
uint64_t
bl_getLocusEndPos (locus_t *loc)
{
	return loc->pos+loc->len-1;
}

/*-------------------------- bl_getLocusChromEndPos --------------------------
 *    
 * @brief return the end of the locus (locally!)
 * @author Steve Hoffmann 
 *   
 */
 
uint64_t
bl_getLocusChromEndPos (locus_t *loc)
{
	return (loc->pos-loc->chrstart)+loc->len-1;
}

/*-------------------------- bl_getMultiLocusChromEndPos --------------------------
 *    
 * @brief return the end of the locus (locally!)
 * @author Steve Hoffmann 
 *   
 */
 
uint64_t
bl_getMultiLocusChromEndPos (multilocus_t *m)
{
	return (m->pos-m->chrstart)+m->len-1;
}


/*--------------------------- bl_getLocusChromDesc ---------------------------
 *    
 * @brief get the description of chrom of a locus
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_getLocusChromDesc (MultiCharSeq *mseq, locus_t *loc)
{
	return ((CharSequence*)mseq->ref[loc->idx].ref)->description;
}

/*------------------------------ bl_getLocusSeq ------------------------------
 *    
 * @brief get the sequence of a locus
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_getLocusSeq (MultiCharSeq *mseq, locus_t *loc)
{
  return &mseq->sequences[loc->pos];
}

/*------------------------- bl_getLocusSeqOffset ----------------------------
 *    
 * @brief get the sequence of a locus
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_getLocusSeqOffset (MultiCharSeq *mseq, locus_t *loc, Uint off)
{
  if(loc->pos > loc->chrstart+off) return &mseq->sequences[loc->pos - off];
  if(loc->chrstart > 0) return &mseq->sequences[loc->chrstart];

  return &mseq->sequences[0];
}



/*------------------------ bl_getLocusChromPosOffset -------------------------
 *    
 * @brief get postion with save offset
 * @author Steve Hoffmann 
 *   
 */
 
uint64_t
bl_getLocusChromPosOffset (locus_t *loc, Uint loff)
{
	
  if(loc->pos > loc->chrstart+loff) return loc->pos - loff;
  if(loc->chrstart > 0) return loc->chrstart;

  return 0;
}

/*--------------------------- bl_getLocusChromIdx ----------------------------
 *    
 * @brief get idx
 * @author Steve Hoffmann 
 *   
 */
 
uint64_t
bl_getLocusChromIdx (locus_t *loc)
{
	return loc->idx;
}

/*--------------------------- bl_getLocusLenOffset ---------------------------
 *    
 * @brief just access the length of a locus
 * @author Steve Hoffmann 
 *   
 */
 
unsigned int
bl_getLocusLenOffset(MultiCharSeq *mseq, locus_t *loc, Uint loff, Uint roff)
{
  Uint right, left;

  assert(loc->pos+loc->len-1 < loc->chrend);

  if(loc->pos + loc->len + roff >= loc->chrend) 
    right = loc->chrend - loc->pos;
  else
    right = loc->len+roff;

  if(loc->pos <= loc->chrstart + loff) 
    left = loc->pos - loc->chrstart;
  else
    left = loff;

  return right+left;
}

/*------------------------------ bl_getLocusLen ------------------------------
 *    
 * @brief just access the length of a locus
 * @author Steve Hoffmann 
 *   
 */
 
unsigned int
bl_getLocusLen(locus_t *loc)
{

  return loc->len;
}


/*----------------------------- bl_initLocuslist -----------------------------
 *    
 * @brief initalize a list of loci
 * @author Steve Hoffmann 
 *   
 */

locuslist_t*
bl_initLocuslist(locuslist_t *list) {
  list->noofloci = 0;
  list->loci = NULL;

  return list;
}


/*-------------------------- bl_addLocustoLocusList --------------------------
 *    
 * @brief at a locus to a list
 * @author Steve Hoffmann 
 *   
 */
 

void
bl_addLocusToLocuslist(locuslist_t *l, locus_t *loc) {

  l->loci = ALLOCMEMORY(NULL, l->loci, locus_t, l->noofloci+1);
  memmove(&l->loci[l->noofloci], loc, sizeof(locus_t));

  l->noofloci++;
}



/*------------------- bl_getLocusFromMultiCharSeqAlignment -------------------
 *    
 * @brief extract locus information from a MSCA
 * @author Steve Hoffmann 
 *   
 */
 
locus_t*
bl_getLocusFromMultiCharSeqAlignment(MultiCharSeqAlignment *mcsa, MultiCharSeq *mseq, locus_t *loc) {
  
  Uint chrstart, chrend;

  loc->idx = mcsa->subidx;
  getMultiCharSeqIdxBounds(mseq, loc->idx, &chrstart, &chrend);
 
  loc->strand = mcsa->strand;
  loc->chrstart = mcsa->substart;
  loc->pos = mcsa->refstart+mcsa->al->voff;
  loc->len = getValignlen(mcsa->al);
  loc->chrend = chrend;
  loc->readstart = bl_alignGet5PrimeU(mcsa->al, mcsa->strand); 
  loc->readend = bl_alignGet3PrimeU(mcsa->al, mcsa->strand);
  loc->edist = getEdist(mcsa->al);
  loc->score = 0; //TODO

  return loc;
}


/*------------------------------- bl_cmpLocus --------------------------------
 *    
 * @brief compare loci wrt to position for sorting with qsort 
 * (numerical/ascending)
 * @author Steve Hoffmann 
 *   
 */

int
bl_cmpLocusPos(const void *a, const void *b) {
  const locus_t *da = (const locus_t *) a;
  const locus_t *db = (const locus_t *) b;

  if(da->idx > db->idx) return 1;
  if(da->idx < db->idx) return -1;

  if(da->strand > db->strand) return 1;
  if(da->strand < db->strand) return -1;

  if(da->pos > db->pos) return 1;
  if(da->pos < db->pos) return -1;

  return 0;
}
/*------------------------------- bl_cmpLocus --------------------------------
 *    
 * @brief compare loci wrt to position for sorting with qsort 
 * (numerical/ascending)
 * @author Steve Hoffmann 
 *   
 */

int
bl_cmpLocusPosNoStrand(const void *a, const void *b) {
  const locus_t *da = (const locus_t *) a;
  const locus_t *db = (const locus_t *) b;

  if(da->idx > db->idx) return 1;
  if(da->idx < db->idx) return -1;

  if(da->pos > db->pos) return 1;
  if(da->pos < db->pos) return -1;

  return 0;
}

/*------------------------------- bl_cmpLocus --------------------------------
 *    
 * @brief compare loci wrt to length for sorting with qsort 
 * (descending)
 * @author Steve Hoffmann 
 *   
 */

int
bl_cmpLocusLen(const void *a, const void *b) {
  const locus_t *da = (const locus_t *) a;
  const locus_t *db = (const locus_t *) b;

  if(da->len > db->len) return -1;
  if(da->len < db->len) return 1;

  if(da->idx > db->idx) return 1;
  if(da->idx < db->idx) return -1;

  if(da->strand > db->strand) return 1;
  if(da->strand < db->strand) return -1;

  if(da->pos > db->pos) return 1;
  if(da->pos < db->pos) return -1;

  return 0;
}


/*------------------------- bl_sortLocuslistByLength -------------------------
 *    
 * @brief sort a list of loci by length
 * @author Steve Hoffmann 
 *   
 */
 
locuslist_t *
bl_sortLocuslistByLength (locuslist_t *list)
{

  qsort(list->loci, list->noofloci, sizeof(locus_t), bl_cmpLocusLen);
	
  return list;
}

/*-------------------------- bl_mergeLocibyDistanceNoStrand --------------------------
 *    
 * @brief merge loci in a locuslist by distance d
 * @author Steve Hoffmann 
 *   
 */
 
locuslist_t*
bl_mergeLocibyDistanceNoStrand (locuslist_t *list, uint64_t d)
{

  Uint i, idx;
  unsigned char strand;
  uint64_t pos, chrstart, chrend, len;
  locuslist_t *mylist;
  locus_t mylocus;



/*  
  fprintf(stderr, "merging loci by distances:\n");

  for(i=0; i < list->noofloci; i++) {
    fprintf(stderr, "%d - [%"PRIu64",%"PRIu64"]\n", i, list->loci[i].pos, list->loci[i].pos+list->loci[i].len-1);
  }
*/

  idx = list->loci[0].idx;
  pos = list->loci[0].pos;
  strand = list->loci[0].strand;
  len = list->loci[0].len;
  chrstart = list->loci[0].chrstart;
  chrend = list->loci[0].chrend;
  
  mylist = ALLOCMEMORY(NULL, NULL, locuslist_t, 1);
  bl_initLocuslist(mylist);

  for(i=1; i < list->noofloci; i++) {
    if(idx == list->loci[i].idx && 
       //strand == list->loci[i].strand && 
       pos + len + d > list->loci[i].pos) {
      //merge
      assert(idx == list->loci[i].idx);
      len = (list->loci[i].pos + list->loci[i].len -1) - pos + 1; 

    } else {
      //add locus
      mylocus.idx = idx;
      mylocus.pos = pos;
      mylocus.strand = strand;
      mylocus.chrstart = chrstart;
      mylocus.chrend = chrend;
      mylocus.len = len;

      bl_addLocusToLocuslist(mylist, &mylocus); 

      pos = list->loci[i].pos;
      chrstart = list->loci[i].chrstart;
      chrend = list->loci[i].chrend;
      idx = list->loci[i].idx;
      len = list->loci[i].len;
      strand = list->loci[i].strand;
    }
  }

  mylocus.idx = idx;
  mylocus.pos = pos;
  mylocus.strand = strand;
  mylocus.chrstart = chrstart;
  mylocus.chrend = chrend;
  mylocus.len = len;

  bl_addLocusToLocuslist(mylist, &mylocus); 

  
  return mylist;
}


/*-------------------------- bl_mergeLocibyDistance --------------------------
 *    
 * @brief merge loci in a locuslist by distance d
 * @author Steve Hoffmann 
 *   
 */
 
locuslist_t*
bl_mergeLocibyDistance (locuslist_t *list, uint64_t d)
{

  Uint i, idx;
  unsigned char strand;
  uint64_t pos, chrstart, chrend, len;
  locuslist_t *mylist;
  locus_t mylocus;


  qsort(list->loci, list->noofloci, sizeof(locus_t), bl_cmpLocusPos);
 

  idx = list->loci[0].idx;
  pos = list->loci[0].pos;
  strand = list->loci[0].strand;
  len = list->loci[0].len;
  chrstart = list->loci[0].chrstart;
  chrend = list->loci[0].chrend;
  
  mylist = ALLOCMEMORY(NULL, NULL, locuslist_t, 1);
  bl_initLocuslist(mylist);

  for(i=1; i < list->noofloci; i++) {
    if(idx == list->loci[i].idx && 
       strand == list->loci[i].strand && 
       pos + len + d > list->loci[i].pos) {
      //merge
      assert(idx == list->loci[i].idx);
      len = (list->loci[i].pos + list->loci[i].len -1) - pos + 1; 

    } else {
      //add locus
      mylocus.idx = idx;
      mylocus.pos = pos;
      mylocus.strand = strand;
      mylocus.chrstart = chrstart;
      mylocus.chrend = chrend;
      mylocus.len = len;

      bl_addLocusToLocuslist(mylist, &mylocus); 

      pos = list->loci[i].pos;
      chrstart = list->loci[i].chrstart;
      chrend = list->loci[i].chrend;
      idx = list->loci[i].idx;
      len = list->loci[i].len;
      strand = list->loci[i].strand;
    }
  }

  mylocus.idx = idx;
  mylocus.pos = pos;
  mylocus.strand = strand;
  mylocus.chrstart = chrstart;
  mylocus.chrend = chrend;
  mylocus.len = len;

  bl_addLocusToLocuslist(mylist, &mylocus); 

  
  return mylist;
}

void
wrapLocuslist(locuslist_t *list) {
  FREEMEMORY(NULL, list->loci);
}


/*---------------------------- bl_getLocusStrand -----------------------------
 *    
 * @brief get the strand of the locus
 * @author Steve Hoffmann 
 *   
 */
 
unsigned char
bl_getLocusStrand (locus_t *loc)
{
	return loc->strand;
}



/*--------------------------- bl_getLociList ---------------------------
 *    
 * @brief get the locus list from mcsa
 * @author Steve Hoffmann 
 *   
 */
 
locuslist_t*
bl_getLocusList (MultiCharSeqAlignment *a, MultiCharSeq *mseq, unsigned int noofaligns)
{
	
  Uint i;
  locuslist_t *mylist;
  locus_t myloc;

  mylist = ALLOCMEMORY(NULL, NULL, locuslist_t, 1);

  bl_initLocuslist(mylist);

  for(i=0; i < noofaligns; i++) {
    bl_getLocusFromMultiCharSeqAlignment(&a[i], mseq, &myloc);
    bl_addLocusToLocuslist(mylist, &myloc);
  }
 

  return mylist;
}
/*--------------------------- bl_mergeLocuList ---------------------------
 *    
 * @brief add locus list b to locus list a 
 * @author Steve Hoffmann 
 *   
 */
 
locuslist_t*
bl_mergeLocusList (locuslist_t *a, locuslist_t *b)
{
	
  Uint i;
  
  for(i=0; i < b->noofloci; i++) {
    bl_addLocusToLocuslist(a, &b->loci[i]);
  }
 

  return a;
}



/*-------------------------- bl_locusListAddOffset ---------------------------
 *    
 * @brief add offset to leftmost and/or rightmost locus of list
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_locusListAddOffset ( MultiCharSeq *mseq, locuslist_t *list, Uint loff, Uint roff)
{
  uint64_t left, right, d, n;

  n = list->noofloci;

  if(n == 0) return;

  left   = bl_getLocusChromPosOffset(&list->loci[0], loff);
  right  = bl_getLocusLenOffset(mseq, &list->loci[n-1], 0, roff);
  

  d = list->loci[0].pos - left;

  list->loci[0].pos = left;
  list->loci[0].len += d;
  list->loci[n-1].len = right;


  return ;
}


/*--------------------------- bl_getLocusSequence ----------------------------
 *    
 * @brief return a concat sequence of loci (i.e. concat exons)
 * @input sorted! list of loci
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_getLocusListSequence (MultiCharSeq *mseq, locuslist_t *list, Uint loff, Uint roff, Uint *seqlen)
    //Uint *start, Uint **splitpos)
{
  Uint i, len=0, tmplen, curloff=0, curroff=0;
  char *myseq = NULL, *tmp;

  for(i=0; i < list->noofloci; i++) {

    curloff = 0;
    curroff = 0;
    
    if(i==0) curloff = loff;
    if(i+1 == list->noofloci) curroff = roff;

    tmp = bl_getLocusSeqOffset (mseq, &list->loci[i], curloff); 
    tmplen = bl_getLocusLenOffset(mseq, &list->loci[i], curloff, curroff);
    
    myseq = ALLOCMEMORY(NULL, myseq, char, len+tmplen+1);
    memmove(&myseq[len], tmp, sizeof(char)*tmplen);
    myseq[len+tmplen] = '\0';

    len += tmplen;
  }

  *seqlen = len;

  return myseq;
}

Uint
bl_getLocusDist(locus_t *a, locus_t *b) {
  
  return b->pos - (a->pos + a->len);
}



char
bl_getLocusListCheck(locuslist_t *list, Uint maxdist) {
  Uint i;

  if(list->noofloci == 0 || list == NULL) return 0;

  for(i=0; i < list->noofloci-1; i++) {
    if(list->loci[i].idx != list->loci[i+1].idx) return 0;
    if(list->loci[i].pos+list->loci[i].len+maxdist < list->loci[i+1].pos) return 0;
  }

  return 1;
}


Alignment**
expandAlignmentDisjoint(Alignment *al, locuslist_t* list, MultiCharSeq *mseq) { 
  Uint i, j, p=0, q=0, k=0, vlen, voff = 0, x = 0, uoff =0,  transcnt=0, o =0;
  char *v;
  Alignment **new = NULL;

  new = ALLOCMEMORY(NULL, NULL, Alignment*, list->noofloci);
  memset(new, 0, sizeof(Alignment*)*list->noofloci);

  voff = al->voff;
  x = list->loci[0].len -1; 
  
  //advancing end of new alignment my expandpos if it begins afterwards
  if(al->voff > x) { 

    while(al->voff > x) {
      o = x;
      x += list->loci[k+1].len;
      k++;
    } 
    voff = al->voff - o - 1;
  } 

  v = bl_getLocusSeq (mseq, &list->loci[k]);
  vlen = bl_getLocusLen(&list->loci[k]);
  uoff = al->uoff;
  
  //fprintf(stdout, "original alignment: voff:%d, uoff:%d\nnew alignment starts with locus %d, voff:%d, uoff:%d, vbreak:%d\n", 
  //    al->voff, al->uoff, k, voff, uoff, x);

  new[k]  = ALLOCMEMORY(NULL, NULL, Alignment, 1);
  initAlignment(new[k], al->u, al->ulen, uoff, v, vlen, voff); 


  //fprintf(stdout, "first locuslen %d-[%"PRIu64"]\n", voff, list->loci[k].len);
  
  for(i=0; i < al->numofmeops; i++) {
    for(j=0; j < al->meops[i].steps; j++) { 

      if(k < list->noofloci-1 && q + al->voff > x) {
 
    //    fprintf(stdout, "transition to next locus\n" );
        transcnt++;
        x += list->loci[k+1].len;
        v = bl_getLocusSeq (mseq, &list->loci[k+1]);
        vlen = bl_getLocusLen(&list->loci[k+1]);
        uoff = al->uoff + p;

        //make sure that previous alignment has ullen >0 in case
        //the alignment already included the query (q>0)
        //otherwise discard previous alignment and exit the loops.

        if(q >0 && getUalignlen(new[k]) ==0) {
      //    fprintf(stdout, "removing alignment as the alignment was unsuccessful");
          wrapAlignment(new[k]);
          FREEMEMORY(NULL, new[k]);
          return new;
        }

      //  fprintf(stdout, "init alignment at u offset %d\n", uoff);
        new[k+1] = ALLOCMEMORY(NULL, NULL, Alignment, 1);
        initAlignment(new[k+1], al->u, al->ulen, uoff, v, vlen, 0);
             
        k++;
      }
         
      insertEop(new[k], al->meops[i].eop);     
      
      if(al->meops[i].eop == Replacement 
          || al->meops[i].eop == Mismatch 
          || al->meops[i].eop == Match) {
        p++;
        q++;
      }

      //here is the problem [uoff,uoff+alignmentlength-1] -> [98,97] for alignment length 0 
      if (al->meops[i].eop == Deletion) {
        q++;
      } 

      if (al->meops[i].eop == Insertion){
        p++;
      }

     if (al->meops[i].eop == Softclip){
        p++;
      }

    }
  }


  if(getUalignlen(new[k]) ==0) {
    wrapAlignment(new[k]);
    FREEMEMORY(NULL, new[k]);
  }

  return new;
}

Alignment*
bl_locusListAlign(MultiCharSeq *mseq, locuslist_t *list, char *qry, Uint qrylen, 
    int* scores, int indel) {

  char *myseq;
  Uint seqlen = 0;
  Alignment *al;
  int *M;

  al = ALLOCMEMORY(NULL, NULL, Alignment, 1);

  //get the sequence
  myseq = bl_getLocusListSequence(mseq, list, 0, 0, &seqlen);
  //apply the semiglobal aligment of mate 1 to new sequence
  M = sgmatrix (NULL, qry, qrylen, myseq, seqlen, indel, constscr, scores);
  //init the alignment
  initAlignment(al, qry, qrylen, 0, myseq, seqlen, 0);
  //trace it back
  sgtraceback(NULL, M, qry, qrylen, myseq, seqlen, indel, constscr, scores, al);

  FREEMEMORY(NULL, M);
  
  return al;

}

void
bl_showLocusList(FILE *dev, locuslist_t *list) {
  Uint i;
  
  fprintf(dev, "showing list with %d loci\n", list->noofloci);
  for(i=0; i < list->noofloci; i++) {
    fprintf(dev, "locus %d: %"PRIu64" (len:%"PRIu64")(rc:%d); readstart:%d\n", i, bl_getLocusPos(&list->loci[i]), list->loci[i].len, list->loci[i].strand, list->loci[i].readstart);
  }

  return;
}


/*-------------------------- bl_locusListGetStartPos -------------------------
 *    
 * @brief get an array of reference (global) start poistions
 * @author Steve Hoffmann 
 *   
 */
 
uint64_t*
bl_locusListGetStartPos (locuslist_t *list)
{
  Uint i;
  uint64_t *arr;

  arr = ALLOCMEMORY(NULL, NULL, uint64_t, list->noofloci);

  for(i=0; i < list->noofloci; i++) {
    arr[i] = bl_getLocusPos(&list->loci[i]); 
  }
	
  return arr;
}

/*-------------------------- bl_locusListGetEndPos -------------------------
 *    
 * @brief get an array of reference (global) end positions
 * @author Steve Hoffmann 
 *   
 */
 
uint64_t*
bl_locusListGetEndPos (locuslist_t *list)
{
  Uint i;
  uint64_t *arr;

  arr = ALLOCMEMORY(NULL, NULL, uint64_t, list->noofloci);

  for(i=0; i < list->noofloci; i++) {
    arr[i] = bl_getLocusEndPos(&list->loci[i]); 
  }
	
  return arr;
}

/*-------------------------- bl_locusListStrand -------------------------
 *    
 * @brief get an array of reference (global) strand
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_locusListStrand (locuslist_t *list)
{
  Uint i;
  char *arr;

  arr = ALLOCMEMORY(NULL, NULL, uint64_t, list->noofloci);

  for(i=0; i < list->noofloci; i++) {
    arr[i] = bl_getLocusStrand(&list->loci[i]); 
  }
	
  return arr;
}

/*------------------------ bl_locusListReadStart ----------------------
 *    
 * @brief get an array of read start positions
 * @author Steve Hoffmann 
 *   
 */
 
Uint*
bl_locusListGetReadStart (locuslist_t *list)
{
  Uint i;
  Uint *arr;

  arr = ALLOCMEMORY(NULL, NULL, uint64_t, list->noofloci);

  for(i=0; i < list->noofloci; i++) {
    arr[i] = list->loci[i].readstart; 
  }
	
  return arr;
}

/*------------------------ bl_locusListReadEnd ----------------------
 *    
 * @brief get an array of read start positions
 * @author Steve Hoffmann 
 *   
 */
 
Uint*
bl_locusListGetReadEnd (locuslist_t *list)
{
  Uint i;
  Uint *arr;

  arr = ALLOCMEMORY(NULL, NULL, uint64_t, list->noofloci);

  for(i=0; i < list->noofloci; i++) {
    arr[i] = list->loci[i].readend; 
  }
	
  return arr;
}


/*------------------------ bl_locusListIsConsecutive -------------------------
 *    
 * @brief check if all loci in list are on the same strand 
 * @author Steve Hoffmann 
 *   
 */
 
char
bl_locusListIsConsecutive (locuslist_t *list)
{
  Uint i;
  char ret=1;

  for(i=1; i < list->noofloci; i++) {
    if(bl_getLocusStrand(&list->loci[i]) != bl_getLocusStrand(&list->loci[i-1])) ret=0;
  }
	
  return ret;
}



/*---------------------------- bl_initMultiLocus -----------------------------
 *    
 * @brief 
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_initMultiLocus (multilocus_t *m)
{

  m->idx = 0;
  m->pos = 0;
  m->len = 0;
  m->strand = 0;
  m->noofloci=0;
  m->loci = NULL;
  m->name = NULL;
  m->score = 0;
  m->readstart = 0;
  m->edist = 0;

  return ;
}

/*--------------------------- bl_setMultiLocusName ---------------------------
 *    
 * @brief set the name of a multilocus
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_setMultiLocusName (multilocus_t *m, char *name)
{
  m->name = bl_strdup(name); 
  return ;
}

/*--------------------------- bl_setMultiLocusScore --------------------------
 *    
 * @brief set the score of a multilocus
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_setMultiLocusScore (multilocus_t *m, uint8_t score)
{
  m->score = score;
  return ;
}


/*------------------------- bl_addLocusToMultiLocus --------------------------
 *    
 * @brief add a locus to a multilocus
 * @author Steve Hoffmann 
 *   
 */
 
char
bl_addLocusToMultiLocus (multilocus_t *m, locus_t *l)
{

 // uint64_t start = 0, len = 0;;

  if(m->noofloci > 0) {
    if(m->idx != l->idx || m->strand != l->strand || 
        (!m->strand && m->loci[m->noofloci-1].pos + m->loci[m->noofloci-1].len -1 > l->pos) ||
        ( m->strand && m->loci[m->noofloci-1].pos < l->pos + l->len -1))
 //     (!m->strand && m->pos > l->pos && m->loci[0].readstart < l->readstart) ||
 //     ( m->strand && m->pos < l->pos && m->loci[0].readstart < l->readstart)) 
    {
      return 0;
    } 

    if(m->pos > l->pos) { 
      m->len = (m->pos + m->len -1) - l->pos +1;
      m->pos = l->pos;
    }
    
    if(m->pos+m->len-1 < l->pos+l->len-1) { 
      m->len = (l->pos+l->len-1) - m->pos + 1;
    }

    //len = end - pos + 1 -> end = pos + len -1;
    //if(m->pos+m->len-1 < l->pos+l->len-1) len = (l->pos+l->len-1) - m->pos + 1; 
    //if(m->pos+m->len-1 < l->pos+l->len-1) len = (l->pos+l->len-1) - m->pos + 1; 
    //if(m->pos > l->pos) m->pos = start;
    //if(m->pos+m->len-1 < l->pos+l->len-1) m->len = len;


    if(m->score > l->score) m->score = l->score;
  } else {
    m->idx = l->idx;
    m->pos = l->pos;
    m->len = l->len;
    m->strand = l->strand;
    m->chrstart = l->chrstart;
    m->chrend = l->chrend;
    m->score = l->score;
    m->readstart = l->readstart;
  } 

  m->loci = ALLOCMEMORY(NULL, m->loci, locus_t, m->noofloci+1);
  memmove(&m->loci[m->noofloci], l, sizeof(locus_t));
  m->noofloci++;
  
  return 1;
}


/*---------------------------- bl_sortMultiLocus -----------------------------
 *    
 * @brief sort a multi locus list
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_sortMultiLocus (multilocus_t *m)
{
	  
  qsort(m->loci, m->noofloci, sizeof(locus_t), bl_cmpLocusPos);

  return ;
}


/*----------------------- bl_addLocusListToMultiLocus ------------------------
 *    
 * @brief convert a list to a multi locus
 * @author Steve Hoffmann 
 *   
 */
 
multilocus_t*
bl_createMultiLocusFromLocusList (locuslist_t *list, Uint *noofmultiloci)
{
  Uint i, k=0;
  multilocus_t *m;

  m = ALLOCMEMORY(NULL, NULL, multilocus_t, 1);
  bl_initMultiLocus(&m[0]);
  
  for(i=0; i < list->noofloci; i++) {
    
    if(!bl_addLocusToMultiLocus(&m[k], &list->loci[i])) {

      //copy last element of existing multilocus with at least two
      //elems to preserve circles
      if(m[k].noofloci > 1 && i > 1) {
        k++; 
        m = ALLOCMEMORY(NULL, m, multilocus_t, k+1);
        bl_initMultiLocus(&m[k]);
        bl_addLocusToMultiLocus(&m[k], &list->loci[i-1]);
      }

      k++;
      m = ALLOCMEMORY(NULL, m, multilocus_t, k+1);
      bl_initMultiLocus(&m[k]);
      bl_addLocusToMultiLocus(&m[k], &list->loci[i]);
    }
  }
	
  *noofmultiloci = k+1;
  return m;
}



/*---------------------------- bl_invertMultiLocusList ----------------------------
 *    
 * @brief this function "inverts" a list [a,b],[c,d],[e,f] -> ]b,c[ ]d,e[ 
 * @author Steve Hoffmann 
 *   
 */
 
multilocus_t*
bl_invertMultiLocus (multilocus_t *m)
{
  Uint j;
  multilocus_t* invert;
  locus_t loc;

  invert = ALLOCMEMORY(NULL, NULL, multilocus_t ,1);
  bl_initMultiLocus (invert);
  

  for(j=1; j < m->noofloci; j++) {
      loc.idx = m->loci[j].idx;
      loc.strand = m->loci[j].strand;
      loc.chrstart = m->loci[j].chrstart;
      loc.chrend = m->loci[j].chrend;
      loc.readstart = m->loci[j].readstart;
      loc.score = m->loci[j].score;

      if(!loc.strand) { 
       loc.pos = m->loci[j-1].pos + m->loci[j-1].len; //-1
       loc.len = m->loci[j].pos - loc.pos;            //+1  
      } else {
       loc.pos = m->loci[j].pos + m->loci[j].len;    //-1
       loc.len = m->loci[j-1].pos - loc.pos;         //+1
      }

      bl_addLocusToMultiLocus (invert, &loc);
  }

  bl_setMultiLocusScore(invert, m->score);
  bl_setMultiLocusName(invert, m->name);

  qsort(invert->loci, invert->noofloci, sizeof(locus_t), bl_cmpLocusPos);

  return invert;
}

/*---------------------------- bl_printMultiLocusSingle -----------------------------
 *    
 * @brief print the multi locus to a string
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_printMultiLocusSingle (multilocus_t *m, MultiCharSeq *mseq, char *name)
{
  Uint i;
  char *string = NULL;

  for(i=0; i < m->noofloci; i++) {

    bl_bsprintf(&string, "%s\t%"PRIu64"\t%"PRIu64"\t%s\t%d\t%c\n", 
       ((CharSequence*)mseq->ref[m->loci[i].idx].ref)->description, 
       bl_getLocusChromPos(&m->loci[i]), bl_getLocusChromEndPos(&m->loci[i])+1, //bed -offset ist 0
       name, m->score, ((m->strand) ? '-' : '+')); 
  }
  
  return string;
}




/*---------------------------- bl_dumpMultiLocusSingle -----------------------------
 *    
 * @brief dump the multi locus to a device
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_dumpMultiLocusSingle (FILE *dev, multilocus_t *m, MultiCharSeq *mseq, char *name)
{
  Uint i;


  for(i=0; i < m->noofloci; i++) {

    fprintf(dev, "%s\t%"PRIu64"\t%"PRIu64"\t%s\t%d\t%c\n", 
       ((CharSequence*)mseq->ref[m->loci[i].idx].ref)->description, 
       bl_getLocusChromPos(&m->loci[i]), bl_getLocusChromEndPos(&m->loci[i])+1, //bed -offset ist 0
       name, m->score, ((m->strand) ? '-' : '+')); 
  }
  
  return ;
}

/*---------------------------- bl_printMultiLocusJoint -----------------------------
 *    
 * @brief print the multi locus to a device
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_printMultiLocusJoint (multilocus_t *m, MultiCharSeq *mseq, char *name)
{
  Uint i;
  char *colorstring;
  char *string = NULL;

  if(m->noofloci == 0) return NULL;

  if(m->noofloci > 1) {
    colorstring = "255,108,0";
  } else { 
    colorstring = "0,229,255";
  }

  bl_bsprintf(&string, "%s\t%"PRIu64"\t%"PRIu64"\t%s\t%d\t%c\t%"PRIu64"\t%"PRIu64"\t%s", 
       ((CharSequence*)mseq->ref[m->idx].ref)->description, 
       bl_getMultiLocusChromPos(m), bl_getMultiLocusChromEndPos(m)+1, //0-offset, last not part of display 
       name, m->score,
       ((m->strand) ? '-' : '+'), 
       bl_getMultiLocusChromPos(m), bl_getMultiLocusChromEndPos(m)+1, colorstring); 

  bl_bsprintf(&string,"\t%d\t", m->noofloci);


  for(i=0; i < m->noofloci; i++) {
    if(i==0) 
      bl_bsprintf(&string,  "%"PRIu64"", m->loci[i].len);
    else  
      bl_bsprintf(&string, ",%"PRIu64"", m->loci[i].len);
  }

  bl_bsprintf(&string, "\t");
  
  for(i=0; i < m->noofloci; i++) {
    if(i==0) 
      bl_bsprintf(&string,  "%"PRIu64"", bl_getLocusChromPos(&m->loci[i])-bl_getMultiLocusChromPos(m));
      //fprintf(dev,  "%"PRIu64"", bl_getLocusChromPos(&m->loci[i]));
    else
      bl_bsprintf(&string, ",%"PRIu64"", bl_getLocusChromPos(&m->loci[i])-bl_getMultiLocusChromPos(m));
      //fprintf(dev, ",%"PRIu64"", bl_getLocusChromPos(&m->loci[i]));
  }

  bl_bsprintf(&string, "\n");

  
  return string;
}


/*---------------------------- bl_dumpMultiLocusJoint -----------------------------
 *    
 * @brief print the multi locus to a device
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_dumpMultiLocusJoint (FILE *dev, multilocus_t *m, MultiCharSeq *mseq, char *name)
{
  Uint i;
  char *colorstring;

  if(m->noofloci == 0) return;

  if(m->noofloci > 1) {
    colorstring = "255,108,0";
  } else { 
    colorstring = "0,229,255";
  }

  fprintf(dev, "%s\t%"PRIu64"\t%"PRIu64"\t%s\t%d\t%c\t%"PRIu64"\t%"PRIu64"\t%s", 
       ((CharSequence*)mseq->ref[m->idx].ref)->description, 
       bl_getMultiLocusChromPos(m), bl_getMultiLocusChromEndPos(m)+1, //0-offset, last not part of display 
       name, m->score,
       ((m->strand) ? '-' : '+'), 
       bl_getMultiLocusChromPos(m), bl_getMultiLocusChromEndPos(m)+1, colorstring); 

  fprintf(dev,"\t%d\t", m->noofloci);


  for(i=0; i < m->noofloci; i++) {
    if(i==0) 
      fprintf(dev,  "%"PRIu64"", m->loci[i].len);
    else  
      fprintf(dev, ",%"PRIu64"", m->loci[i].len);
  }

  fprintf(dev, "\t");
  
  for(i=0; i < m->noofloci; i++) {
    if(i==0) 
      fprintf(dev,  "%"PRIu64"", bl_getLocusChromPos(&m->loci[i])-bl_getMultiLocusChromPos(m));
      //fprintf(dev,  "%"PRIu64"", bl_getLocusChromPos(&m->loci[i]));
    else
      fprintf(dev, ",%"PRIu64"", bl_getLocusChromPos(&m->loci[i])-bl_getMultiLocusChromPos(m));
      //fprintf(dev, ",%"PRIu64"", bl_getLocusChromPos(&m->loci[i]));
  }

  fprintf(dev, "\n");

  
  return ;
}


/*---------------------------- bl_wrapMultiLocus -----------------------------
 *    
 * @brief destroy multi locus
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_wrapMultiLocus (multilocus_t *m)
{
  m->noofloci = 0;
  FREEMEMORY(NULL, m->loci);
  FREEMEMORY(NULL, m->name);
 return ;
}



/*-------------------------- bl_initMultiLocusList ---------------------------
 *    
 * @brief init the list of multiloci
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_initMultiLocusList (multilocuslist_t *l)
{

  l->noofmultiloci = 0;
  l->loci = NULL;
  return ;
}


/*-------------------------- bl_addMultiLocusToList --------------------------
 *    
 * @brief add a multilocus_t to a multi locus list
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_addMultiLocusToList (multilocuslist_t *l, multilocus_t *m)
{
  Uint j;

  l->loci = ALLOCMEMORY(NULL, l->loci, multilocus_t, l->noofmultiloci+1);
  bl_initMultiLocus (&l->loci[l->noofmultiloci]);
  memmove(&l->loci[l->noofmultiloci], m, sizeof(multilocus_t));

  for(j=1; j < m->noofloci; j++) {
      bl_addLocusToMultiLocus (&l->loci[l->noofmultiloci], &m->loci[j]);
  }

  l->noofmultiloci += 1;
  return ;
}


/*----------------------------- bl_cmpMultiLoci ------------------------------
 *    
 * @brief 
 * @author Steve Hoffmann 
 *   
 */
 
int
bl_cmpMultiLocusPos ( const void *a, const void *b) {
  const multilocus_t *da = (const multilocus_t *) a;
  const multilocus_t *db = (const multilocus_t *) b;

  if(da->idx > db->idx) return 1;
  if(da->idx < db->idx) return -1;

  if(da->strand > db->strand) return 1;
  if(da->strand < db->strand) return -1;

  if(da->pos > db->pos) return 1;
  if(da->pos < db->pos) return -1;

  return 0;

}

/*-------------------------- bl_sortMultiLocusList ---------------------------
 *    
 * @brief sort the multi locus list
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_sortMultiLocusList (multilocuslist_t *l)
{

  qsort(l->loci, l->noofmultiloci, sizeof(multilocus_t), bl_cmpMultiLocusPos);
  return ;
}


/*------------------------- bl_isIdenticalMultiLocus -------------------------
 *    
 * @brief return 1 if two loci are identical
 * @author Steve Hoffmann 
 *   
 */
 
int
bl_cmpMultiLocusParts (void *a, void *b)
{
  int val;
  Uint j;

  const multilocus_t *da = (const multilocus_t *) a;
  const multilocus_t *db = (const multilocus_t *) b;

  val = bl_cmpMultiLocusPos(a, b);


  if(val == 0) { 

    if(da->noofloci > db->noofloci) return -1;
    if(da->noofloci < db->noofloci) return 1;

    for(j=1; j < da->noofloci; j++) {

    }
  }

  return val;
}
