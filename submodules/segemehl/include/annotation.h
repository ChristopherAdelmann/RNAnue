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


#ifndef _ANNOTATION_
#define _ANNOTATION_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stringutils.h"
#include "basic-types.h"
#include "charsequence.h"
#include "gzidx.h"

typedef struct {
  int64_t dir5prime;
  int64_t dir3prime; 
  int64_t left;
  int64_t right;
  int8_t select;
} annotationoffs_t;


typedef struct {

  unsigned char type;
  Uint trackid;
  /*
   *
   * chromname is gff seqname
   *
   */
  char *chromname; 
  Uint chromnamelen;

  /*
   * BED und GFF: 1-offset
   * for personalSNP: start with 0-offset
   * end base is not part of the feature ie. if 
   * end = 100 last feature base is 99. see below.
   *
   */

  uint64_t start;
  uint64_t end;

  /* 
   * GFF: name is the feature key
   * BED: name is the name of the acutal feature
   * for personalSNP track name is alleles A,C,T,G separated by '/'. 
   * Leading '-' is indel: insertion if start-end=0 
   *
   * */
  
  char *name;     
  Uint namelen;
  double score;
  unsigned char strand;

  /*GFF fields*/
  unsigned char frame; 
  char *source; 
  Uint sourcelen;
  Uint noofattributes;
  char **attributes;
  Uint *attributelen;

  /*BED fields*/
  uint64_t thickStart;
  uint64_t thickEnd;
  Uint *itemRgb;
  Uint blockCount;
  uint64_t* blockSizes;
  uint64_t* blockStarts;
  Uint noofovl;
  Uint firstovl;
  Uint level;
  /*extension*/
  char **blockRefseqs;
  char *blockStrands;

  /*SNPitem*/
  Uint alleleCount;//number of alleles in name
  Uint *alleleFreq;//from comma separated list of number of observed alleles - if unkowns 0
  Uint *alleleScores;//from a comma separated list - if unkown 0

} annotationitem_t;

typedef struct {
    uint64_t noofchr;
    char **chrname;
    Uint *chrnamelen;
    Uint sz;
    uint64_t *noofbins;
    uint64_t **last; //stores the last interval strictly left of indexed interval [i*sz,i+1*sz];

} annotationindex_t;

typedef struct {
 
  uint64_t init;
  char sorted;
  char *trackname;
  Uint tracknamelen;
  char *description;
  Uint descriptionlen;
  Uint noofitems;
  annotationitem_t *items;
  char *filename;
  Uint filenamelen;
} annotationtrack_t;


typedef struct {
    uint64_t init;
    char sorted;
    Uint nooftracks;
    char **trackname;
    Uint *tracknamelen;
    char **description;
    Uint *descriptionlen;
    char **filename;
    Uint *filenamelen;
    Uint noofitems;
    annotationitem_t *items;
} annotationmultitrack_t;

void bl_annotationtrackDestruct (void *space, annotationtrack_t *track);
Uint bl_annotationtrackGetStats (void *space, annotationtrack_t *track);
void bl_annotationitemInit (annotationitem_t *item, unsigned char type);
void bl_annotationtrackInit (annotationtrack_t *track);
void bl_annotationitemDestruct (void *space, annotationitem_t *item);
Uint bl_annotationitem_cmp_track (Uint item, void *track, void *elem, void *nfo);
int bl_annotationitem_nostrand_cmp (void const *a, void const *b);
int bl_annotationitem_cmp (void const *a, void const *b);
void bl_annotationtrackDestruct (void *space, annotationtrack_t *track);
void bl_annotationtrackAssignTrackLevel(annotationtrack_t *track);
annotationmultitrack_t* bl_annotationtrackJoin(void *space, annotationmultitrack_t *dest, annotationtrack_t *src);
annotationitem_t* bl_annotationitemCopy(annotationitem_t *dest, annotationitem_t *src);
void bl_annotationmultitrackInit (annotationmultitrack_t *track);
Uint bl_annotationitem_cmp_multitrack (Uint item, void *track, void *elem, void *nfo);
void bl_annotationtrackSetItems(annotationtrack_t* track, annotationitem_t* items, Uint n);
annotationindex_t* bl_annotationIndex(annotationtrack_t *t);
void bl_annotationtrackDumpIndex(annotationindex_t *idx,annotationtrack_t *track); 
void bl_annotationmultitrackDestruct (void *space, annotationmultitrack_t *track);
void bl_annotationitemApplyOffset(annotationitem_t *item, int64_t off3, int64_t off5, int64_t left, int64_t right, int8_t select);
void  bl_annotationitemDump(FILE *dev, annotationitem_t *item);
void bl_annotationmultitrackApplyOffset(annotationmultitrack_t *mtrack, annotationoffs_t *off);
#endif

