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
 *  multiseq.c
 *  some functions to handle multiseqs (type char)
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/15/06 11:42:53 CET
 *  
 *  SVN
 *  Revision of last commit: $Rev: 66 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-10-02 13:38:05 +0200 (Thu, 02 Oct 2008) $
 *
 *  Id: $Id: multicharseq.c 66 2008-10-02 11:38:05Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sufarray/multicharseq.c $
 */
#include "multicharseq.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <pthread.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/times.h>
#include <time.h>

#include "basic-types.h"
#include "mathematics.h"
#include "charsequence.h"
#include "sort.h"
#include "info.h"	//(skipped)
#include "memory.h"	//(skipped)
#include "alignment.h"	//(skipped)


/*----------------------------- bl_mcsaGet5primeU -----------------------------
 *    
 * @brief get 5'-end of msca on query
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_mcsaGet5PrimeU (MultiCharSeqAlignment *al)
{

	return bl_alignGet5PrimeU (al->al, al->strand);
}



/*----------------------------- bl_mcsaGet5primeV -----------------------------
 *    
 * @brief get 5'-end of msca on reference
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_mcsaGet5primeV (MultiCharSeqAlignment *al)
{

	return al->refstart + bl_alignGet5PrimeV(al->al, al->strand);
}



/*----------------------------- bl_mcsaGet3primeU -----------------------------
 *    
 * @brief get 3'-end of msca on query
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_mcsaGet3PrimeU (MultiCharSeqAlignment *al)
{

	return bl_alignGet3PrimeU(al->al, al->strand);
}




/*---------------------------- concatCharSequences ----------------------------
 *    
 * concatenates CharSequences using a given Uint delimiter
 * and stores them in a MultiCharSeq container.
 * 
 */
 
MultiCharSeq *
concatCharSequences (void *space, CharSequence **s, Uint len, 
					char delim, char sentinel)
{
    char *buf=NULL;
    char *map = NULL;
    Uint i, j, k=0, 
		 totallength=0, 
		 *markpos;
	MultiCharSeq *mseq;

	mseq = ALLOCMEMORY(space, NULL, MultiCharSeq, 1);
	markpos = ALLOCMEMORY(space, NULL, Uint, len);
	mseq->ref = ALLOCMEMORY(space, NULL, SeqReference, len);
    map = ALLOCMEMORY(space, NULL, char, 257);
    memset(map, 0, 256);
    mseq->delim = delim;

	for(i=0; i < len; i++) {

        mseq->ref[i].ref = s[i];
	  
        totallength += (s[i]->length+1);
		buf = ALLOCMEMORY(space, buf, char, totallength+1);
		if (buf==NULL) {
          NFO("allocation of %d bytes failed: exiting\n", totallength);
          exit(-1);
        }

		for(j=0; j < s[i]->length; j++) {
			buf[k] = s[i]->sequence[j];
                        if ((Uint)buf[k] == 0){
                          NFO("invalid character (NUL) in database sequences. Exit forced\n", NULL);
                          exit(-1);
                        }
            map[(Uint)buf[k]]=buf[k];
            k++;
		}
		/*separate sequences or finalize*/
		if (i == (len-1)) {
		  buf[k] = sentinel;
          map[(Uint)buf[k]]=buf[k];
		  markpos[i] = k;
		  k++;
          buf[k]='\0';

		} else {
		  buf[k] = delim;
		  map[(Uint)buf[k]]=buf[k];
          markpos[i] = k;
		  k++;
		}
        
        /*FREEMEMORY(space, s[i]->sequence);*/
	}
	mseq->totallength = totallength;
	mseq->numofsequences = len;
	mseq->sequences = buf;
	mseq->markpos = markpos;

    for(i=0; i < 256; i++) {
        if(map[i]==0) {
           j=i+1;
           while(j<256 && map[j]==0) j++;
           if (j < 256) {
             map[i]=map[j];
             map[j]=0;
           } else {
             break;
           }
        }
    }

    map = ALLOCMEMORY(space, map, char, i+1);
    mseq->map = map;
    mseq->mapsize = i;


	return mseq;
}



/*----------------------------- destructMultiSeq -----------------------------
 *    
 * destructs a MultiSeq structure
 * 
 */

void
destructMultiCharSeq (void *space, MultiCharSeq *mseq)
{
    
	FREEMEMORY(space, mseq->sequences);
	if (mseq->markpos != NULL)
      FREEMEMORY(space, mseq->markpos);
	if (mseq->map != NULL)
      FREEMEMORY(space, mseq->map);
    if (mseq->ref != NULL)
      FREEMEMORY(space, mseq->ref);
    FREEMEMORY(space, mseq);
	return ;
}


/*------------------------------- cmp_markpos --------------------------------
 *    
 * compare function for getMultiSeqIndex
 * 
 */
 
Uint
cmp_markpos (Uint a, void *data, void *key, void *info)
{
    Uint *d = (Uint*) data;
	Uint *k = (Uint*) key;
	
	if (d[a] > *k) {
		if (a > 0) {
			if (d[a-1] < *k) {
				return 0;
			} else {
				return 1;
			}
		} else {
			return 0;
		}
	}
	
    if (d[a] < *k) return 2;
	return 0;
}

/*-------------------------- getMultiSeqIndex --------------------------
 *    
 * returns index of a sequence in multiseq addressed by a pointer
 * 
 */
 
Uint
getMultiCharSeqIndex (MultiCharSeq *mseq, char *ptr)
{	
	Uint pos, i;
	
	if (mseq->numofsequences == 1){
	  return 0;
	}
	pos = (ptr - mseq->sequences); 
    if (mseq->numofsequences < MSEQ_BSEARCH_THRESHOLD) {
      i=binarySearch(mseq->markpos, mseq->numofsequences, &pos, 
          cmp_markpos, NULL);
    } else {
      for (i=0; i < mseq->numofsequences; i++) {
        if (mseq->markpos[i] > pos) break;
      }
    }

	return i;
}

/*------------------------- getMultiCharSeqIdxBounds -------------------------
 *    
 * @brief return start and end of idx 
 * @author Steve Hoffmann 
 *   
 */
 
void
getMultiCharSeqIdxBounds(MultiCharSeq *mseq, Uint idx, Uint *start, Uint *end)
{

  *start = (idx > 0) ? mseq->markpos[idx-1]+1 : 0;
  *end = mseq->markpos[idx]; 
  return ;
}


/*---------------------------- nextMultiSeqDelim -----------------------------
 *    
 * returns positions of next delimiter in multiseq (ie. end of current seq)
 * 
 */

Uint
nextMultiSeqDelim (MultiCharSeq *mseq, char *ptr)
{	
  return mseq->markpos[getMultiCharSeqIndex(mseq,ptr)];
}

 

/*---------------------------- getMultiSeqRelPos -----------------------------
 *    
 * returns the relative position of a pointer to multiseq
 * with respect to the addressed sequence.
 * 
 */
 
Uint
getMultiCharSeqRelPos (MultiCharSeq *mseq, char *ptr)
{
  Uint idx;
  CharSequence *seq;
  idx = getMultiCharSeqIndex(mseq, ptr);
  seq = getCharSequence(mseq, idx);
  return (ptr - seq->sequence);
}


/*------------------------------- dumpMultiSeq -------------------------------
 *    
 * dumps a multiseq to the screen
 * 
 */

void
dumpMultiCharSeq (MultiCharSeq *mseq)
{
  	Uint i;

	for(i=0; i < mseq->totallength; i++) {
		printf("%c-", mseq->sequences[i]);	
	}

	printf("\n");
	return ;
}


/*----------------------------- getCharSequence ------------------------------
 *    
 * @brief return the CharSequence at index idx
 * @author Steve Hoffmann 
 *   
 */
 
CharSequence*
getCharSequence (MultiCharSeq *mseq, Uint idx)
{
	return (CharSequence*) mseq->ref[idx].ref;
}


/*------------------------- getMaximumAlignmentEdist --------------------------
 *    
 * @brief get the maximum edit distances for a sequence length and an accuracy
 * @author Steve Hoffmann 
 *   
 */
 
Uint
getMaximumAlignmentEdist(Uint seqlength, float accuracy )
{
  Uint edist;

  edist = seqlength-floor((accuracy*seqlength)/100);
  return edist;
}


/*------------------------ getExtendedAlignmentLength ------------------------
 *    
 * @brief returns the length of the region needed for an errneous alignment
 * @author Steve Hoffmann 
 *   
 */
 
Uint
getExtendedAlignmentLength (Uint seqlength, float accuracy)
{
  Uint len;
  Uint maxedist;

  //get the maximum edist
  maxedist = getMaximumAlignmentEdist(seqlength, accuracy);
  //add the seqlength and the maximum edit distance for each side of length
  len  = seqlength + 2*(maxedist+1) ;
  //allow an additional margin of 10%
  len += ceil(((float)len/100.0)*10.0);

  return len;
}

/*------------------------ initMultiCharSeqAlignment -------------------------
 *    
 * @brief initalize an alignment in the multichar seq
 *
 * @author Steve Hoffmann 
 *   
 */
 
int
initMultiCharSeqAlignment(
    void *space, MultiCharSeqAlignment* a, MultiCharSeq *seq, Uint pos, 
    Uint loff, Uint len, unsigned char strand, 
    char *qrydesc, char *query, char *qual, Uint qrylen) 
{
  Uint sub_start, 
       sub_end;
  //Uint     i;

  a->subidx = getMultiCharSeqIndex(seq, &seq->sequences[pos]);
  getMultiCharSeqIdxBounds(seq, a->subidx, &sub_start, &sub_end);
  a->substart = sub_start;
  a->subend = sub_end;
   
  a->refstart = MAX(sub_start, (Lint)pos-loff);
  
  if(a->refstart > sub_end) {
    fprintf(stderr, "refstart > substart: skipping MultiCharSeqAlignment\n");
    return 0;
  }

  a->reflen = (sub_end > (Lint)a->refstart + len - 1) ? len : (sub_end - a->refstart); //changed from refstart+len ? ... (sub_end - a->refstart) +1
  a->refseq = &seq->sequences[a->refstart];
  a->refdesc = ((CharSequence*)seq->ref[a->subidx].ref)->description;
  a->qrydesc = qrydesc;
  a->query = query;
  a->strand = strand;
  a->qrylen = qrylen;
  a->qual = qual;
/*
  for(i=0; i < a->reflen; i++) {
    fprintf(stdout, "%c", a->refseq[i]);
  }
  fprintf(stdout, "\n");
*/  
  a->al = ALLOCMEMORY(space, NULL, Alignment, 1);
  initAlignment(a->al, query, qrylen, 0, a->refseq, a->reflen, 0);
	
  return 1;
}

void
wrapMultiCharSeqAlignment(void *space, MultiCharSeqAlignment *a) {
  wrapAlignment(a->al);
  FREEMEMORY(space, a->al);
}
/*----------------------- initMultiCharSeqAlignmentOpt -----------------------
 *    
 * @brief init a mcsa with query bounds
 * @author Steve Hoffmann 
 *   
 */
 
int
initMultiCharSeqAlignmentOpt(
    void *space, MultiCharSeqAlignment* a, MultiCharSeq *seq, Uint pos, 
    char *qrydesc, char *query, char *qual, Uint start, Uint end, 
    Uint qrylen, Uint floff, Uint flen, Uint uloff, Uint uroff, Uint maxoff, unsigned char strand) 
{

   Uint sub_start, 
       sub_end,
       rlen, qstart, qend, qlen;
   
  //get bounds and length of reference sequence
  a->subidx = getMultiCharSeqIndex(seq, &seq->sequences[pos]);
  getMultiCharSeqIdxBounds(seq, a->subidx, &sub_start, &sub_end);
  a->substart = sub_start;
  a->subend = sub_end;
  a->refstart = MAX(sub_start, (Lint)pos-floff); //maxoff
  a->floff = pos - a->refstart;

  //this should not happen  
  if(a->refstart > sub_end) {
    fprintf(stderr, "refstart > substart: skipping MultiCharSeqAlignment\n");
    return 0;
  }

  rlen = flen;
  a->refseq = &seq->sequences[a->refstart];
//  a->reflen = (sub_end > (Lint)a->refstart + rlen) ? rlen : (sub_end - a->refstart)+1;
  a->reflen = (sub_end > (Lint)a->refstart + rlen - 1) ? rlen : (sub_end - a->refstart); //changed from refstart+len ? ... (sub_end - a->refstart) +1
 
  //get bounds and length of query sequence
  qstart = (start > maxoff+uloff) ? start-(maxoff+uloff) : 0;
  qend = (end + uroff + maxoff < qrylen) ? end + uroff + maxoff : qrylen;
  qlen = qend - qstart;
  a->query = query;
  a->qrystart = qstart;
//  fprintf(stderr, "storing mcsa qrylen with qlen:%d, original qrylen was: %d\n", qlen, qrylen);
  a->maxalignlen = qlen;
  //a->qrylen = qlen;
  a->qrylen = qrylen;
  a->strand = strand;
  a->qual = qual;

  //descriptions
  a->refdesc = ((CharSequence*)seq->ref[a->subidx].ref)->description;
  a->qrydesc = qrydesc;
  //init alignment and return
  a->al = ALLOCMEMORY(space, NULL, Alignment, 1);
  initAlignment(a->al, query, qrylen, 0, a->refseq, a->reflen, 0);
	
  return 1;
}

//#define JOINDEBUG

MultiCharSeqAlignment*
joinalignments(MultiCharSeqAlignment **al, Uint noofaligns, unsigned char rev, Uint head, Uint tail, Uint lsize, Uint rsize, Eoptype cliptype) {

  Uint k,i,j,l, u, vlen, cur = 0;
  Alignment *new;
  MultiCharSeqAlignment *mcsa;
  char *vseq=NULL, *vptr;

  mcsa = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
  new = ALLOCMEMORY(NULL, NULL, Alignment, 1); 
  //caution: this is a non standard Alignment init!
  new->u = al[0]->al->u;
  new->ulen = al[0]->al->ulen;
  new->vlen = al[0]->al->vlen; 
  new->voff = 0;
  new->numofmeops = 0;
  new->rmvseq = 1;
  new->rmuseq = 0;

  u = (rev) ?  noofaligns-1 : 0;
  new->uoff = al[u]->al->uoff;
  new->meops = ALLOCMEMORY(NULL, NULL, Multieop, (al[u]->al->numofmeops+2)*sizeof(Multieop)); 

#ifdef JOINDEBUG
  Uint minlen = MIN(new->ulen, new->vlen);
  fprintf(stdout, "-> segment ends are: %d-%d (rev:%d)\n", al[u]->refstart+al[u]->al->voff, 
      al[u]->refstart+al[u]->al->voff+getValignlen(al[u]->al), rev);
  fprintf(stdout, "2: minimum length: %d, this fragment is head=%d and tail=%d, lsize:%d, rsize:%d, uoff:%d\n", 
      minlen, head, tail, lsize, rsize, new->uoff);
  showAlign(al[u]->al, stdout);
#endif

  if(lsize) {
    new->uoff -= lsize;
#ifdef JOINDEBUG
    char *dseq = ALLOCMEMORY(NULL, NULL, char, lsize+1);
    for(i=0; i < lsize; i++) {
      dseq[i] = al[u]->al->u[new->uoff+i];
    }
    dseq[i]=0;
    fprintf(stdout, "(rev) sequence deleted from reference: %s\n", dseq);
    FREEMEMORY(NULL, dseq);
#endif

    for(i=0; i < lsize; i++) { 
      insertEop(new, cliptype);
    }
  }
  
  memmove(&new->meops[lsize>0], al[u]->al->meops, al[u]->al->numofmeops*sizeof(Multieop));
  new->numofmeops += al[u]->al->numofmeops;
  vptr = al[u]->al->v;
  vlen = getValignlen(al[u]->al);
  vseq = ALLOCMEMORY(NULL, vseq, char, cur+vlen+1);
  memmove(&vseq[cur], &vptr[al[u]->al->voff], vlen);
  cur+= vlen;
  vseq[cur] = 0;

  
  for(k=1; k < noofaligns; k++) { 
    if(k > 0) assert(al[k]->al->u == al[k-1]->al->u);
    u = (rev) ? noofaligns - k -1 : k;


    vptr = al[u]->al->v;
    vlen = getValignlen(al[u]->al);

    vseq = ALLOCMEMORY(NULL, vseq, char, cur+vlen+1);
    memmove(&vseq[cur], &vptr[al[u]->al->voff], vlen);
    cur+= vlen;
    vseq[cur] = 0;
    Uint dsize = 0;
    new->meops = ALLOCMEMORY(NULL, new->meops, Multieop, (al[u]->al->numofmeops + new->numofmeops +2)*sizeof(Multieop)); 


#ifdef JOINDEBUG    
        Uint ulen = getUalignlen(al[u]->al);
    fprintf(stdout, "-> segment ends are: %d-%d (rev:%d)\n", al[u]->refstart+al[u]->al->voff, 
      al[u]->refstart+al[u]->al->voff+getValignlen(al[u]->al), rev);
    fprintf(stdout, "processing seq:%p with offset:%u with len:%d (total:%d)\n", al[u]->al->v, al[u]->al->voff, vlen, cur);
    fprintf(stdout, "query uoff:%d, ulen:%d\n", al[u]->al->uoff, ulen);
    fprintf(stdout, "partial alignment\n");
    showAlign(al[u]->al, stdout);
#endif

    if(!rev && k > 0 && al[u-1]->al->uoff + getUalignlen(al[u-1]->al) < al[u]->al->uoff) {
      Uint dstart = (al[u-1]->al->uoff + getUalignlen(al[u-1]->al));
      dsize =  al[u]->al->uoff - dstart;
#ifdef JOINDEBUG      
      char *dseq = ALLOCMEMORY(NULL, NULL, char, dsize+1);
      for(i=0; i < dsize; i++) {
        dseq[i] = al[u]->al->u[dstart+i];
      }
      dseq[i]=0;
      fprintf(stdout, "(fwd) sequence deleted from reference: %s\n", dseq);
      FREEMEMORY(NULL, dseq);
#endif
      for(j=0; j < dsize; j++) {
        insertEop(new, Insertion);
      }
    }

    if(!rev) { 
      vlen = getValignlen(al[u-1]->al);
      Uint voff = al[u-1]->al->voff;
      Uint vstart = voff + al[u-1]->refstart;
      Uint vend = vstart + vlen -1;
      Uint vnext = al[u]->refstart + al[u]->al->voff;
      Uint gap = (vnext - vend - 1); //- dsize?
      for(j=0; j < gap; j++) { 
        insertEop(new, Skipped);
      }
    }

    //if alignment in reverse,
    if(rev && u < noofaligns-1 && al[u]->al->uoff > al[u+1]->al->uoff+getUalignlen(al[u+1]->al)) { 
      Uint dstart = (al[u+1]->al->uoff+getUalignlen(al[u+1]->al));
      Uint dsize =  al[u]->al->uoff - dstart;
#ifdef JOINDEBUG
      char *dseq = ALLOCMEMORY(NULL, NULL, char, dsize+1);
      for(i=0; i < dsize; i++) {
        dseq[i] = al[u]->al->u[dstart+i];
      }
      dseq[i]=0;
      fprintf(stdout, "(rev) sequence deleted from reference: %s\n", dseq);
      FREEMEMORY(NULL, dseq);
#endif
      for(j=0; j < dsize; j++) {
        insertEop(new, Insertion);
      }
    }

    if(rev) { 
      vlen = getValignlen(al[u+1]->al);
      Uint voff = al[u+1]->al->voff;
      Uint vstart = voff + al[u+1]->refstart;
      Uint vend = vstart + vlen -1;
      Uint vnext = al[u]->refstart + al[u]->al->voff;
      Uint gap = (vnext - vend - 1); //- dsize?
      for(j=0; j < gap; j++) {
        insertEop(new, Skipped);
      }
    }

    //insert silent insertions
      for(j=0; j < al[u]->al->numofmeops; j++) {
        for(l=0; l < al[u]->al->meops[j].steps; l++) {
          insertEop(new, al[u]->al->meops[j].eop);
        }
      }
  }  
 
  if(rsize) {
    for(i=0; i < rsize; i++) { 
      insertEop(new, cliptype);
    }
  }
  
  /*register new in mcsa*/
  
  u = (rev) ?  noofaligns-1 : 0;
  new->v = vseq;
  new->vlen = cur;
  new->voff = 0;

#ifdef JOINDEBUG
  fprintf(stdout, "voff: %d, vlen:%d, v:%s, u:%s, ulen:%d\n", new->voff, new->vlen, vseq, new->u, new->ulen);
  char *mymeop =  multieopstring(new, 0, 0, 0);
  char *mymdstr = mdstring(new, 0);
  char *mycigar = cigarstring(new, 0, 0, 0, 0);

  fprintf(stdout, "multieop: %s\n", mymeop);
  fprintf(stdout, "mdstring: %s\n", mymdstr);
  fprintf(stdout, "cigarstr: %s\n", mycigar);
  showAlign(new, stdout);

  FREEMEMORY(NULL, mymdstr);
  FREEMEMORY(NULL, mymeop);
  FREEMEMORY(NULL, mycigar);

  fprintf(stdout, "*** final ends are: %d-%d (rev:%d)\n", al[u]->refstart+al[u]->al->voff, 
      al[u]->refstart+al[u]->al->voff+getValignlen(new), rev);
#endif

  mcsa->subidx = al[u]->subidx;
  mcsa->substart = al[u]->substart;
  mcsa->subend = al[u]->subend;
  //resetting the offsets in the concat alignment requires to correct the refstart by the 
  //first partial voff
  mcsa->refstart = al[u]->refstart+al[u]->al->voff;
  mcsa->floff = al[u]->floff;
  mcsa->refseq = al[u]->refseq;
  mcsa->reflen = al[u]->reflen;
  mcsa->query = new->u;
  mcsa->qrystart = new->uoff;
  mcsa->maxalignlen = MAX(getUalignlen(new), getValignlen(new));
  mcsa->qrylen = al[u]->qrylen;
  mcsa->strand = rev;
  mcsa->qual = al[u]->qual;
  mcsa->refdesc = al[u]->refdesc;
  mcsa->qrydesc = al[u]->qrydesc;
  mcsa->al = new;


  return mcsa;
}


/*------------------- bl_getPartialMultiCharSeqAlignments --------------------
 *    
 * @brief get partial multi char seq alignments from a single alignment with skips.
 * @author Steve Hoffmann 
 *   
 */

MultiCharSeqAlignment*
bl_getPartialMultiCharSeqAlignments (MultiCharSeqAlignment *mcsa, MultiCharSeq *mseq, Uint *noofaligns)
{

  Uint i, noofparts=0;
  Alignment *al, *copy;
  MultiCharSeqAlignment *new=NULL;

  al = bl_getPartialAlignments(mcsa->al, &mseq->sequences[mcsa->refstart], &noofparts);
  new = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, noofparts);

  for(i=0; i < noofparts; i++) {

    new[i].subidx = mcsa->subidx;
    new[i].substart = mcsa->substart;
    new[i].subend = mcsa->subend;
    // all alignments have the same offset = 0
    new[i].refstart = mcsa->refstart;
    new[i].floff = mcsa->floff;
    new[i].refseq = mcsa->refseq;
    new[i].reflen = mcsa->reflen;
    new[i].query = mcsa->query;
    new[i].qrylen = mcsa->qrylen;
    new[i].qrystart = al[i].uoff;
    new[i].maxalignlen = MAX(getUalignlen(&al[i]), getValignlen(&al[i]));
    new[i].strand = mcsa->strand;
    new[i].qual = mcsa->qual;
    new[i].refdesc = mcsa->refdesc;
    new[i].qrydesc = mcsa->qrydesc;
    copy = ALLOCMEMORY(NULL, NULL, Alignment, 1);
    copyAlignment(copy, &al[i]);
    new[i].al = copy;
    wrapAlignment(&al[i]);
  }

  FREEMEMORY(NULL, al);
  *noofaligns = noofparts;
  return new;
}

void
reevalMultiCharSeqAlignment(MultiCharSeqAlignment *mcsa) {
  Alignment* new;
  
  new=reevalAlignment(mcsa->al);
  
  wrapAlignment(mcsa->al);
  FREEMEMORY(NULL, mcsa->al);
                
  mcsa->al=new; 
  return;
}
