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


#ifndef MULTI_SEQ_H
#define MULTI_SEQ_H

/*
 *	multiseq.h
 *  declarations for a datastructure containing
 *  multiple integer sequences
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 12/11/06 15:09:15 CET  
 *
 *  SVN
 *  Revision of last commit: $Rev: 65 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-09-22 00:14:54 +0200 (Mon, 22 Sep 2008) $
 *
 *  Id: $Id: multicharseq.h 65 2008-09-21 22:14:54Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sufarray/multicharseq.h $
 */


#include "basic-types.h" 
#include "charsequence.h"
#include "alignment.h"

#define MSEQ_BSEARCH_THRESHOLD 10

typedef struct {
  void *ref; 	
} SeqReference;


typedef struct {
  Uint numofsequences;
  Uint totallength;
  Uint *markpos;		/*markpos[i] is the position of a*/
  /*separator character between S_i and S_i+1*/
  char *sequences; 	/*array of concatenated sequences*/
  SeqReference *ref;  /*ref[i] points to the original sequence*/
  /*that starts at position markpos[i]*/
  char *map;
  Uint mapsize;
  char delim;

} MultiCharSeq;


typedef struct {
  Uint subidx;
  Uint substart;
  Uint subend;

  char *refdesc;
  char *refseq;
  Uint refstart;
  Uint reflen;
  Uint floff;

  char *qrydesc;
  char *query;
  char *qual;
  Uint qrystart;
  Uint qrylen;
  Uint maxalignlen;
  Alignment *al;
  unsigned char strand;
  unsigned char pass;

} MultiCharSeqAlignment;

void dumpMultiCharSeq (MultiCharSeq *);
MultiCharSeq* concatCharSequences(void *, CharSequence **, Uint, char, char);
void destructMultiCharSeq(void*, MultiCharSeq *);
Uint getMultiCharSeqIndex(MultiCharSeq *, char *);
Uint getMultiCharSeqRelPos(MultiCharSeq *, char *);
CharSequence* getCharSequence(MultiCharSeq *, Uint idx);
void getMultiCharSeqIdxBounds(MultiCharSeq *mseq, Uint idx, Uint *start, Uint *end);
int initMultiCharSeqAlignment(
    void *space, MultiCharSeqAlignment* a, MultiCharSeq *seq, Uint pos, Uint loff, Uint len, 
    unsigned char strand, char *querydesc, char *query, char* qual, Uint qlen);
void wrapMultiCharSeqAlignment(void *space, MultiCharSeqAlignment *a);

int initMultiCharSeqAlignmentOpt(
    void *space, MultiCharSeqAlignment* a, MultiCharSeq *seq, Uint pos, 
    char *qrydesc, char *query, char*qual, Uint start, Uint end, 
    Uint qrylen, Uint floff, Uint flen, Uint uloff, Uint uroff, Uint maxoff, unsigned char strand) ;

MultiCharSeqAlignment*
joinalignments(MultiCharSeqAlignment **al, Uint noofaligns, unsigned char rev, Uint head, Uint tail, Uint lsize, Uint rsize, Eoptype cliptype);
Uint getExtendedAlignmentLength (Uint readlength, float accuracy);
Uint getMaximumAlignmentEdist(Uint seqlength, float accuracy );
Uint bl_mcsaGet5PrimeU (MultiCharSeqAlignment *al);
Uint bl_mcsaGet3PrimeU (MultiCharSeqAlignment *al);
Uint bl_mcsaGet5primeV (MultiCharSeqAlignment *al);
MultiCharSeqAlignment* bl_getPartialMultiCharSeqAlignments (MultiCharSeqAlignment *mcsa, MultiCharSeq *mseq, Uint *noofaligns);
void reevalMultiCharSeqAlignment(MultiCharSeqAlignment *mcsa);

#endif
