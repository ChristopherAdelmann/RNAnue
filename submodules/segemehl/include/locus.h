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


#ifndef LOCUS_H
#define LOCUS_H

/*
 *
 *	locus.h
 *  a genomic locus given a reference
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 09.12.2012 14:16:18 CET  
 *
 */

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



typedef struct locus_s {
  unsigned int idx;
  uint64_t pos;
  uint64_t len;
  uint64_t chrstart;
  uint64_t chrend;
  char strand; 
  char *name;

  uint8_t score;
  Uint readstart;
  Uint readend;
  Uint edist;

} locus_t;


typedef struct locuslist_s {
  unsigned int noofloci;
  locus_t *loci;
} locuslist_t;

typedef struct multilocus_s {
  unsigned int idx;
  uint64_t pos;
  uint64_t len;
  char strand;
  uint64_t chrstart;
  uint64_t chrend;
  char *name;

  uint8_t score;
  Uint readstart;
  Uint edist;
  
  Uint noofloci;

  locus_t *loci;
} multilocus_t;

typedef struct multilocuslist_s {
  Uint noofmultiloci;
  multilocus_t *loci;
} multilocuslist_t;

  uint64_t bl_getLocusChromPos ( locus_t *loc);
  uint64_t bl_getLocusChromEndPos ( locus_t *loc);
  char* bl_getLocusChromDesc (MultiCharSeq *mseq, locus_t *loc);
  locus_t* bl_initLocus (locus_t* loc, MultiCharSeq *mseq, uint64_t pos,
    unsigned char strand, uint64_t length, uint64_t loff) ;
  uint64_t bl_getLocusChromIdx (locus_t *loc);
  locus_t* bl_getLocusFromMultiCharSeqAlignment(MultiCharSeqAlignment *mcsa, MultiCharSeq *mseq, locus_t *loc);
  void bl_addLocusToLocuslist(locuslist_t *l, locus_t *loc);
  locuslist_t* bl_initLocuslist(locuslist_t *list);
  int bl_cmpLocusLen(const void *a, const void *b);
  locuslist_t * bl_sortLocuslistByLength (locuslist_t *list);
  locuslist_t* bl_mergeLocibyDistance (locuslist_t *list, uint64_t d);
  void wrapLocuslist(locuslist_t *list);
  char* bl_getLocusSeq (MultiCharSeq *mseq, locus_t *loc);
  unsigned int bl_getLocusLen(locus_t *loc);
  unsigned int bl_getLocusLenOffset(MultiCharSeq *mseq, locus_t *loc, Uint loff, Uint roff);
  char* bl_getLocusSeqOffset (MultiCharSeq *mseq, locus_t *loc, Uint off);
  unsigned char bl_getLocusStrand (locus_t *loc);
  uint64_t bl_getLocusChromPosOffset (locus_t *loc, Uint off);
  char* bl_getLocusListSequence (MultiCharSeq *mseq, locuslist_t *list, Uint loff, Uint roff, Uint *seqlen);//Uint *start, Uint **splitpos);
  locuslist_t* bl_getLocusList (MultiCharSeqAlignment *a, MultiCharSeq *mseq, unsigned int noofaligns);
  int bl_cmpLocusPos(const void *a, const void *b);
  void bl_locusListAddOffset ( MultiCharSeq *mseq, locuslist_t *list, Uint loff, Uint roff);
  Alignment* bl_getLocusListExpandedAlignment(Alignment *al, locuslist_t *list);
  char bl_getLocusListCheck(locuslist_t *list, Uint maxdist);
  Alignment** expandAlignmentDisjoint(Alignment *al, locuslist_t* list, MultiCharSeq *mseq);
  Alignment* bl_locusListAlign(MultiCharSeq *mseq, locuslist_t *list, char *qry, Uint qrylen, int* scores, int indel);
 locuslist_t* bl_mergeLocibyDistanceNoStrand (locuslist_t *list, uint64_t d);
 int bl_cmpLocusPosNoStrand(const void *a, const void *b);
 void bl_showLocusList(FILE *, locuslist_t *list);
 uint64_t* bl_locusListGetStartPos (locuslist_t *list);
 uint64_t* bl_locusListGetEndPos (locuslist_t *list);
 char* bl_locusListStrand (locuslist_t *list);
 Uint* bl_locusListGetReadStart (locuslist_t *list);
 Uint* bl_locusListGetReadEnd (locuslist_t *list);
 char bl_locusListIsConsecutive (locuslist_t *list);
 uint64_t bl_getLocusEndPos (locus_t *loc);
 multilocus_t* bl_createMultiLocusFromLocusList (locuslist_t *list, Uint *noofmultiloci);
 void bl_wrapMultiLocus (multilocus_t *m);
 void bl_dumpMultiLocusSingle (FILE *dev, multilocus_t *m, MultiCharSeq *mseq, char *name);
 void bl_dumpMultiLocusJoint (FILE *dev, multilocus_t *m, MultiCharSeq *mseq, char *name);
 multilocus_t* bl_invertMultiLocus (multilocus_t *m);
 void bl_setMultiLocusScore (multilocus_t *m, uint8_t score);
 void bl_setMultiLocusName (multilocus_t *m, char *name);
 locuslist_t* bl_mergeLocusList (locuslist_t *a, locuslist_t *b);
 char* bl_printMultiLocusJoint (multilocus_t *m, MultiCharSeq *mseq, char *name);
 char* bl_printMultiLocusSingle (multilocus_t *m, MultiCharSeq *mseq, char *name);


#endif
