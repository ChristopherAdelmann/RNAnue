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

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <errno.h>
#include <inttypes.h>
#include "debug.h"
#include "zlib.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "annotation.h"
#include "fileio.h"
#include "seqclip.h"
#include "charsequence.h"
#include "assert.h"
#include "gzidx.h"
#include "info.h"
#include "bitVector.h"


/*-------------------------- bl_annotationtrackInit --------------------------
 *
 * @brief init
 * @author Steve Hoffmann
 *
 */

void
bl_annotationtrackInit (annotationtrack_t *track)
{

  track->init = MAGIC_INIT;
  track->sorted = 0;
  track->trackname = NULL;
  track->tracknamelen = 0;
  track->description = NULL;
  track->descriptionlen = 0;
  track->noofitems = 0;
  track->items = NULL;
  track->filename = NULL;
  track->filenamelen = 0;
  return ;
}

/*----------------------- bl_annotationmultitrackInit -----------------------
 *
 * @brief init
 * @author Steve Hoffmann
 *
 */

void
bl_annotationmultitrackInit (annotationmultitrack_t *track)
{

  track->init = MAGIC_INIT;
  track->sorted = 0;
  track->nooftracks = 0;
  track->trackname = NULL;
  track->tracknamelen = NULL;
  track->description = NULL;
  track->descriptionlen = NULL;
  track->noofitems = 0;
  track->items = NULL;
  track->filename =NULL;
  track->filenamelen=NULL;

  return ;
}


/*------------------------ bl_annotationtrackitemInit ------------------------
 *
 * @brief init item
 * @author Steve Hoffmann
 *
 */

void
bl_annotationitemInit (annotationitem_t *item, unsigned char type)
{

  item->trackid = 0;
  item->type = type;
  item->chromname = NULL;
  item->chromnamelen = 0;
  item->source=NULL;
  item->sourcelen=0;
  item->start = 0;
  item->end = 0;
  item->name = NULL;
  item->namelen = 0;
  item->score = .0;
  item->strand = '0';
  item->thickStart = 0;
  item->thickEnd = 0;
  item->itemRgb = NULL;
  item->blockCount = 0;
  item->blockSizes = NULL;
  item->blockStarts = NULL;
  item->blockStrands = NULL;
  item->blockRefseqs = NULL;
  item->noofovl = 0;
  item->firstovl = -1;
  item->level = 0;
  item->source = NULL;
  item->sourcelen = 0;
  item->noofattributes=0;
  item->attributes = NULL;
  item->attributelen = NULL;

  return ;
}

/*------------------------- bl_annotationitemDestruct -----------------------
 *
 * @brief destruct annotation item
 * @author Steve Hoffmann
 *
 */

void
bl_annotationitemDestruct (void *space, annotationitem_t *item)
{

  Uint i;

  if(item->chromname) FREEMEMORY(space, item->chromname);
  if(item->name) FREEMEMORY(space, item->name);
  if(item->itemRgb) FREEMEMORY(space, item->itemRgb);
  if(item->blockSizes) FREEMEMORY(space, item->blockSizes);
  if(item->blockStarts) FREEMEMORY(space, item->blockStarts);
  if(item->blockRefseqs) FREEMEMORY(space, item->blockRefseqs);
  if(item->blockStrands) FREEMEMORY(space, item->blockStrands);
  if(item->source) FREEMEMORY(space, item->source);

  if(item->noofattributes) {
    for(i=0; i < item->noofattributes; i++) {
      FREEMEMORY(space, item->attributes[i]);
    }
    FREEMEMORY(space, item->attributes);
    FREEMEMORY(space, item->attributelen);
  }

  return ;
}

/*--------------------------- bl_annotationitemCopy --------------------------
 *
 * @brief cpy annotation item
 * @author Steve Hoffmann
 *
 */

annotationitem_t*
bl_annotationitemCopy(annotationitem_t *dest, annotationitem_t *src) {

    Uint i;


    memmove(dest, src, sizeof(annotationitem_t));

    dest->chromname=bl_strdup(src->chromname);
    dest->name=bl_strdup(dest->name);
    dest->source=bl_strdup(dest->source);

    if(src->noofattributes) {
        dest->attributes=ALLOCMEMORY(NULL, NULL, char*, src->noofattributes);
        dest->attributelen=ALLOCMEMORY(NULL, NULL, Uint, src->noofattributes);
        for(i=0; i < src->noofattributes; i++) {
            dest->attributes[i] = ALLOCMEMORY(NULL, NULL, char, src->attributelen[i]);
            memmove(dest->attributes[i], src->attributes[i], sizeof(char)*(src->attributelen[i]+1));
            dest->attributelen[i] = src->attributelen[i];
        }
    }

    if(src->blockCount) {
        dest->blockSizes=ALLOCMEMORY(NULL, NULL, uint64_t, src->blockCount);
        dest->blockStarts=ALLOCMEMORY(NULL, NULL, uint64_t, src->blockCount);

        if(src->blockRefseqs) {
            dest->blockRefseqs = ALLOCMEMORY(NULL, NULL, char*, src->blockCount);
        }

        if(src->blockStrands) {
            dest->blockRefseqs = ALLOCMEMORY(NULL, NULL, char, src->blockCount);
        }

        for(i=0; i < src->blockCount; i++) {
            dest->blockSizes[i] = src->blockSizes[i];
            dest->blockStarts[i]= src->blockStarts[i];

            if(src->blockRefseqs) {
                dest->blockRefseqs[i] = bl_strdup(src->blockRefseqs[i]);
            }

            if(src->blockStrands) {
                dest->blockStrands[i]= src->blockStrands[i];
            }
        }
    }

   if(src->alleleCount) {
        dest->alleleFreq=ALLOCMEMORY(NULL, NULL, Uint, src->alleleCount);
        dest->alleleScores=ALLOCMEMORY(NULL, NULL, Uint, src->alleleCount);
        memmove(dest->alleleFreq, src->alleleFreq, sizeof(Uint)*src->alleleCount);
        memmove(dest->alleleScores, src->alleleScores, sizeof(Uint)*src->alleleCount);
    }

   return dest;
}


/*------------------------ bl_annotationitem_cmp_track ------------------------
 *
 * @brief find annotation item in track
 * @author Steve Hoffmann
 *
 */

Uint
bl_annotationitem_cmp_track (Uint item, void *track, void *elem, void *nfo)
{
  annotationitem_t *l, *r;
  annotationtrack_t *t;
  int chr;

  t = (annotationtrack_t*) track;

  l = (annotationitem_t*) &t->items[item];
  r = (annotationitem_t*) elem;


  if ((chr = strcmp(l->chromname, r->chromname))) {
    if(chr < 0) return 2;
    if(chr > 0) return 1;
  }

  if(l->end < r->start) {
    return 2;
  }

  if(l->end > r->start) {
    return 1;
  }

  return 0;
}

/*--------------------- bl_annotationitem_cmp_multitrack ---------------------
 *
 * @brief find annotation item in track
 * @author Steve Hoffmann
 *
 */

Uint
bl_annotationitem_cmp_multitrack (Uint item, void *track, void *elem, void *nfo)
{
  annotationitem_t *l, *r;
  annotationmultitrack_t *t;
  int chr;

  t = (annotationmultitrack_t*) track;

  l = (annotationitem_t*) &t->items[item];
  r = (annotationitem_t*) elem;


  if ((chr = strcmp(l->chromname, r->chromname))) {
    if(chr < 0) return 2;
    if(chr > 0) return 1;
  }

  if(l->end < r->start) {
    return 2;
  }

  if(l->end > r->start) {
    return 1;
  }

  return 0;
}

/*--------------------------- bl_annotationitem_nostrand_cmp --------------------------
 *
 * @brief compare annotation lines w/o strandedness
 * @author Steve Hoffmann
 *
 */

int
bl_annotationitem_nostrand_cmp (void const *a, void const *b)
{
  annotationitem_t *l, *r;
  int chr;

  l = (annotationitem_t*) a;
  r = (annotationitem_t*) b;

  if ((chr = strcmp(l->chromname, r->chromname))) {
    return chr;
  }

  if(l->start < r->start) {
    return -1;
  }

  if(l->start > r->start) {
    return 1;
  }

  if(l->end < r->end) {
    return -1;
  }

  if(l->end > r->end) {
    return 1;
  }


  return 0;
}



/*--------------------------- bl_annotationitem_cmp --------------------------
 *
 * @brief compare annotation lines
 * @author Steve Hoffmann
 *
 */

int
bl_annotationitem_cmp (void const *a, void const *b)
{
  annotationitem_t *l, *r;
  int chr;

  l = (annotationitem_t*) a;
  r = (annotationitem_t*) b;

  if ((chr = strcmp(l->chromname, r->chromname))) {
    return chr;
  }

  if(l->start < r->start) {
    return -1;
  }

  if(l->start > r->start) {
    return 1;
  }

  if(l->end < r->end) {
    return -1;
  }

  if(l->end > r->end) {
    return 1;
  }

  if(l->strand < r->strand) {
    return -1;
  }

  if(l->strand > r->strand) {
    return 1;
  }

  return 0;
}

/*------------------------ bl_annotationtrackDestruct -------------------------
 *
 * @brief init
 * @author Steve Hoffmann
 *
 */

void
bl_annotationtrackDestruct (void *space, annotationtrack_t *track)
{

  Uint i;

  if(track->trackname) FREEMEMORY(space, track->trackname);
  if(track->filename) FREEMEMORY(space, track->filename);
  if(track->description) FREEMEMORY(space, track->description);

  for(i=0; i < track->noofitems; i++) {
    bl_annotationitemDestruct(space, &track->items[i]);
  }

  track->noofitems = 0;
  if(track->items) FREEMEMORY(space, track->items) ;

  return ;
}

/*---------------------- bl_annotationmultitrackDestuct ---------------------
 *
 * @brief destruct
 * @author Steve Hoffmann
 *
 */

void
bl_annotationmultitrackDestruct (void *space, annotationmultitrack_t *track)
{
  uint32_t i;

  for(i=0; i < track->nooftracks; i++) {
    FREEMEMORY(space, track->trackname[i]);
    FREEMEMORY(space, track->description[i]);
    FREEMEMORY(space, track->filename[i]);
  }
    FREEMEMORY(space, track->filenamelen);
    FREEMEMORY(space, track->tracknamelen);
    FREEMEMORY(space, track->descriptionlen);

    FREEMEMORY(space, track->filename);
    FREEMEMORY(space, track->trackname);
    FREEMEMORY(space, track->description);


  //ITEMS STILL BELONG TO INDIVIDUAL TRACKS
/*  for(i=0; i < track->noofitems; i++) {
    bl_annotationitemDestruct(space, &track->items[i]);
  }
*/
  track->noofitems = 0;
  if(track->items) FREEMEMORY(space, track->items) ;

  return ;
}

/*-------------------------- bl_annotationtrackJoin ---------------------------
 *
 * @brief init
 * @author Steve Hoffmann
 *
 */

annotationmultitrack_t*
bl_annotationtrackJoin(void *space, annotationmultitrack_t *dest,
        annotationtrack_t *src) {

  uint64_t i, j, k, n;

  assert(dest->init == MAGIC_INIT && src->init == MAGIC_INIT);

  i = dest->noofitems;
  j = src->noofitems;
  k = dest->nooftracks;


  n = i + j;

  dest->items = ALLOCMEMORY(NULL, dest->items, annotationitem_t, n);
  memmove(&dest->items[i], src->items, j*sizeof(annotationitem_t));
  dest->noofitems = n;

  dest->trackname = ALLOCMEMORY(NULL, dest->trackname, char*, k+1);
  dest->tracknamelen = ALLOCMEMORY(NULL, dest->tracknamelen, Uint, k+1);

  dest->description = ALLOCMEMORY(NULL, dest->description, char*, k+1);
  dest->descriptionlen = ALLOCMEMORY(NULL, dest->descriptionlen, Uint, k+1);

  dest->filename = ALLOCMEMORY(NULL, dest->filename, char*, k+1);
  dest->filenamelen = ALLOCMEMORY(NULL, dest->filenamelen, Uint, k+1);

  dest->trackname[k] = bl_strdup(src->trackname);
  dest->description[k] = bl_strdup(src->description);
  dest->filename[k] = bl_strdup(src->filename);

  dest->tracknamelen[k]=src->tracknamelen;
  dest->descriptionlen[k]=src->descriptionlen;
  dest->filenamelen[k]=src->filenamelen;

  for(; i < n; i++) {
      dest->items[i].trackid = k;
  }

  dest->nooftracks++;
  return dest;
}


/*------------------- bl_annotationtrackAssignTrackLevel -------------------
 *
 * @brief assign the track numbers to annota
 * @author Steve Hoffmann
 *
 */

void
bl_annotationtrackAssignTrackLevel(annotationtrack_t *track)
{

  Uint i, k, p;
  bitvector a;

  for(i=0; i< track->noofitems; i++ ) {
    for(k=i+1; k < track->noofitems; k++) {
      if(track->items[i].end+1 <= track->items[k].start ||
         strcmp(track->items[i].chromname, track->items[k].chromname))
        break;
      if(track->items[k].firstovl == -1)
        track->items[k].firstovl = i;
      track->items[k].noofovl++;
    }
  }

  for(i=0; i < track->noofitems; i++) {
    if(track->items[i].noofovl < 2) {
      track->items[i].level = track->items[i].noofovl;
    } else {
      a=initbitvector(NULL, 255);
      for(k=track->items[i].firstovl; k < i; k++) {
        if(track->items[k].end+1 >= track->items[i].start &&
           !strcmp(track->items[i].chromname, track->items[k].chromname)) {
          bitvector_setbit(a, track->items[k].level, 1);
        }
      }
      for(p=0; p < 255; p++) {
        if (bitvector_getbit(a,p) == 0) break;
      }
      track->items[i].level = p;
      FREEMEMORY(space, a);
    }
  }

  return ;
}


/*------------------------ bl_annotationtrackGetStats ------------------------
 *
 * @brief get number of different loci in annotation track
 * @author Steve Hoffmann
 *
 */

Uint
bl_annotationtrackGetStats (void *space, annotationtrack_t *track)
{

  Uint i=0,j,noofdups, len;
  char *attr;
  annotationitem_t *a, *b;

  while(i < track->noofitems) {
    a = &track->items[i];
    for(j=i+1; j < track->noofitems; j++) {
      b = &track->items[j];
      if(a->start != b->start || a->end != b->end || a->strand != b->strand) {
        break;
      }
    }

    noofdups = j-i;

    for(j=i; j < i+noofdups; j++) {
      b = &track->items[j];

      len = snprintf(NULL, 0, "loci_cnt %d %d", j-i+1, noofdups);
      attr = ALLOCMEMORY(space, NULL, char, len+1);
      snprintf(attr, len+1,"loci_cnt %d %d", j-i+1, noofdups);
      bl_GFFAddAttribute(space, b, attr, len);
    }

    i+=noofdups;

  }

  return 0;
}

/*-------------------------- bl_annotationitemDump ---------------------------
 *
 * @brief dump item
 * @author Steve Hoffmann
 *
 */

void
bl_annotationitemDump(FILE *dev, annotationitem_t *item) {
    fprintf(dev, "%s\t%"PRIu64"\t%"PRIu64"\n", item->chromname, item->start,
            item->end);
}

/*----
 *
 *
 *
 */

void
bl_annotationtrackSetItems(annotationtrack_t* track, annotationitem_t* items, Uint n) {
  assert(track->init == MAGIC_INIT);
  qsort(items, n, sizeof(annotationitem_t), bl_annotationitem_cmp);
  track->sorted = 1;
  track->noofitems = n;
  track->items = items;
}

/*----
 *
 *
 *
 */

annotationindex_t*
bl_annotationIndex(annotationtrack_t *t) {

    //assert that track is sorted
    assert(t->sorted);

    Uint sz = 100000;
    char *curchrom = NULL;
    annotationindex_t *idx;

    idx = ALLOCMEMORY(NULL, NULL, annotationindex_t, 1);
    idx->noofchr = 0;
    idx->chrname = NULL;
    idx->chrnamelen = NULL;
    idx->sz = sz;
    idx->noofbins = NULL;
    idx->last = NULL;

    for(uint64_t i=0; i < t->noofitems; i++) {

        if(!curchrom || strcmp(t->items[i].chromname,curchrom)) {
            idx->chrname = ALLOCMEMORY(NULL, idx->chrname, char*, idx->noofchr+1);
            idx->chrname[idx->noofchr] = bl_strdup(t->items[i].chromname);
            idx->chrnamelen = ALLOCMEMORY(NULL, idx->chrnamelen, char*, idx->noofchr+1);
            idx->chrnamelen[idx->noofchr] = t->items[i].chromnamelen;
            idx->noofbins = ALLOCMEMORY(NULL, idx->noofbins, uint64_t, idx->noofchr+1);
            idx->noofbins[idx->noofchr] = 0;
            idx->last = ALLOCMEMORY(NULL, idx->last, uint64_t, idx->noofchr+1);
            idx->last[idx->noofchr] = NULL;
            idx->noofchr++;
            curchrom = t->items[i].chromname;
        }

        uint64_t bin = t->items[i].end/sz;

        if(bin >= idx->noofbins[idx->noofchr-1]){
            idx->last[idx->noofchr-1] = ALLOCMEMORY(NULL, idx->last[idx->noofchr-1], uint64_t, bin+1);

            for(uint64_t k=idx->noofbins[idx->noofchr-1]; k < bin; k++) {
                idx->last[idx->noofchr-1][k] = i;
            }

            idx->noofbins[idx->noofchr-1] = bin+1;
            idx->last[idx->noofchr-1][bin] = i;
        }

        if(t->items[idx->last[idx->noofchr-1][bin]].end < t->items[i].end) {
            idx->last[idx->noofchr-1][bin] =i;
        }
    }

    return idx;
}

/*----
 *
 *
 *
 */


void
bl_annotationtrackDumpIndex(annotationindex_t *idx, annotationtrack_t *track) {
    uint64_t i, k, sz;

    sz = idx->sz;

    for(i=0; i < idx->noofchr; i++) {
        fprintf(stderr, "index chromosome %"PRIu64" ('%s') of '%"PRIu64"'\n",
                i, idx->chrname[i], idx->noofchr);
        for(k=0; k < idx->noofbins[i]-1; k++) {
            fprintf(stderr, "\t%"PRIu64"[%"PRIu64",%"PRIu64"]=%"PRIu64"\n",
                    k, k*sz, (k+1)*sz, idx->last[i][k]);
            fprintf(stderr, "\t%"PRIu64" | %"PRIu64"\n", track->items[idx->last[i][k]].end,
                    track->items[idx->last[i][k+1]].end);
        }
    }

    return;
}

/*----
 *
 *
 *
 */

void
bl_annotationitemApplyOffset(annotationitem_t *item, int64_t off3,
    int64_t off5, int64_t left, int64_t right, int8_t select) {

  Uint i;
  unsigned char err =0;
  int64_t l, r, d, s, t, e;

  if(left != LONG_MIN) {
    s = item->start;
    l = left;
  } else {
    s = item->end;
    l = 0;
  }

  if(right != LONG_MIN) {
    e = item->end;
    r = right;
  } else {
    e = item->start;
    r = 0;
  }

  if(select == 1) {
    if(item->strand == '+') {
      s = item->start;
      e = item->start;
    } else if(item->strand == '-') {
      s = item->end;
      e = item->end;
    }
  } else if (select == 2) {
     if(item->strand == '+') {
      s = item->end;
      e = item->end;
    } else if(item->strand == '-') {
      s = item->start;
      e = item->start;
    }
  }

  if(item->strand == '+') {

    if(off5 != LONG_MIN) {
      l += off5;
    } else {
      s = item->end;
    }
    if(off3 != LONG_MIN) {
      r += off3;
    } else {
      e = item->start;
    }

  } else if (item->strand == '-') {

    if(off5 != LONG_MIN) {
      r -= off5;
    } else {
      e = item->start;
    }
    if(off3 != LONG_MIN) {
      l -= off3;
    } else {
      s = item->end;
    }
  }

  if(l >= 0 || s > labs(l)) {
    s += l;
  } else {
    s = 0;
  }

  if(r >= 0 || e > labs(r)) {
    e += r;
  } else {
    e = 0;
  }

  if(e < s) {
    t = e;
    e = s;
    s = t;
  }

  d = s - item->start;

  if(d && item->blockCount) {
    for(i=0; i < item->blockCount; i++){
      if(item->blockStarts[i] + d < s ||
          item->blockStarts[i] + item->blockSizes[i] + d > e) {
        /*violation of BED convention*/
        err = 1;
      }
      item->blockStarts[i] += d;
    }
  }

  if(err) {
    NFO("violation of BED chromStart <= [block] <= chromEnd convention", NULL);
  }

  //fprintf(stderr, "[%lu,%lu]->[%ld,%ld]\n", item->start, item->end, s, e);

  item->start = s;
  item->end = e;

  return;
}

void
bl_annotationtrackApplyOffset(annotationtrack_t *track,  annotationoffs_t *off) {
  Uint i;

  for(i=0; i < track->noofitems; i++) {
    bl_annotationitemApplyOffset(&track->items[i], off->dir3prime,
        off->dir5prime, off->left, off->right, off->select);
  }
}

void
bl_annotationmultitrackApplyOffset(annotationmultitrack_t *mtrack, annotationoffs_t *off) {
  Uint i;

  for(i=0; i < mtrack->noofitems; i++) {
    bl_annotationitemApplyOffset(&mtrack->items[i], off->dir3prime,
        off->dir5prime, off->left, off->right, off->select);
  }
}

void
bl_annotationitemInitOffset(annotationoffs_t *off) {
  off->dir3prime = 0;
  off->dir5prime = 0;
  off->left = 0;
  off->right = 0;
  off->select = 0;
}
