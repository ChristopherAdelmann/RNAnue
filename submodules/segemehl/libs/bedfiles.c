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
 *  bedfiles.c
 *  routines for reading and writing biofiles
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 10/30/16 17:38:47 CET
 *  
 */
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
#include "biofiles.h"
#include "fileio.h"
#include "seqclip.h"
#include "charsequence.h"
#include "assert.h"
#include "gzidx.h"
#include "info.h"
#include "bitVector.h"
#include "biofiles.h"
#include "annotation.h"

/*-------------------------------- bl_BEDread --------------------------------
 *    
 * @brief read a bed file
 * @author Steve Hoffmann 
 *   
 */
 
annotationtrack_t*
bl_BEDread (void *space, char *filename)
{
  stringset_t **set;
  annotationtrack_t *track;
  annotationitem_t *item;
  char *str, *pch, *tmp;
  Uint linecount, i, j, k, u, v, noofstr, len, ulen;
 
  track = ALLOCMEMORY(space, NULL, annotationtrack_t, 1);
  bl_annotationtrackInit(track);
  track->filename = bl_strdup(filename);
  track->filenamelen = strlen(filename);
  

  set = readcsv(space, filename, "\t", &linecount);
  track->items = ALLOCMEMORY(space, NULL, annotationitem_t, linecount);

  for(i=0; i < linecount; i++) {
    noofstr = set[i]->noofstrings;

    if(noofstr) { 
      str = set[i]->strings[0].str;
      len = strlen(str);
  
      //comment line
      if(strncmp(str, "#", 1) == 0) { 
       destructStringset(space, set[i]);
        continue;
      }

      //track description
      if(len >= 5 && !strncmp(str, "track", 5)) {
        for(j=1; j < noofstr; j++) {
          str = set[i]->strings[j].str;
          len = strlen(str);

          if(len > 5 && !strncmp(str, "name=", 5)) {
            track->tracknamelen = len-5;
            track->trackname = ALLOCMEMORY(space, NULL, char, len-4);
            memmove(track->trackname, &str[5], len-5);
            track->trackname[len-5] = '\0';
            fprintf(stderr, "storing trackname %s\n",track->trackname);
          }
   
          if(len > 12 && !strncmp(str, "description=", 12)) {
            track->descriptionlen = len-12;
            track->description = ALLOCMEMORY(space, NULL, char, len-11);
            memmove(track->description, &str[5], len-12);
            track->description[len-12] = '\0';
          }
        }
        destructStringset(space, set[i]);
        continue;
      }

      //real data
      if(noofstr >= 3) { 
        item = &track->items[track->noofitems];
        bl_annotationitemInit(item, BEDITEM);

        for(j=0; j < noofstr; j++) {
          str = set[i]->strings[j].str;
          len = strlen(str);

          switch(j) {
            case 0:             
              item->chromnamelen = len;
              item->chromname = ALLOCMEMORY(space, NULL, char, len+1);
              memmove(item->chromname, str, len);
              item->chromname[len] = '\0';
              break;
            case 1:
              item->start = atoi(str);
              if(!item->start && str[0] != '0') {
                DBG("BED '%s' %d:%d: atoi failed", filename, i, j);
                exit(-1);
              }
              break;
            case 2:
              item->end = atoi(str);
              if(!item->end && str[0] != '0') {
                DBG("BED '%s' %d:%d: atoi failed", filename, i, j);
                exit(-1);
              }
              break;
            case 3:
              item->namelen = len;
              item->name = ALLOCMEMORY(space, NULL, char, len+1);
              memmove(item->name, str, len);
              item->name[len] = '\0';
              break;
            case 4:
              item->score = atof(str);
              if(item->score == 0.0 && str[0] != '0' && str[0] != '.') {
                DBG("BED '%s' %d:%d: %f(%s) :atof failed", filename, i, j, item->score, str);
                exit(-1);
              }
              break;
            case 5:
              if(str[0] != '-' && str[0] != '+' && str[0] != '.') { 
                DBG("BED '%s' %d:%d: atof failed", filename, i, j);
                exit(-1);
              }
              item->strand = str[0];
              break;
            case 6:
              item->thickStart = atoi(str);
              if(!item->thickStart && str[0] != '0') {
                DBG("BED '%s' %d:%d: %s:atoi failed", filename, i, j, str);
                exit(-1);
              }
              break;
            case 7:
              item->thickEnd = atoi(str);
              if(!item->thickEnd && str[0] != '0') {
                DBG("BED '%s' %d:%d: atoi failed", filename, i, j);
                exit(-1);
              }
              break;
            case 8:
              k = 0;
              pch = strtok(str, ",");
              while(pch != NULL) { 
                item->itemRgb = ALLOCMEMORY(space, NULL, Uint, k+1);
                item->itemRgb[k] = atoi(pch); 
                if(!item->itemRgb[k] && pch[0] != '0' && k > 2) {
                  DBG("BED '%s' %d:%d: atoi failed", filename, i, j);
                  exit(-1);
                }
                k++;
                pch = strtok(NULL, ",");
              }
              if(k == 1) {
                FREEMEMORY(space, item->itemRgb);
                item->itemRgb = NULL;
              }
              if(k != 1 &&  k != 3) {
                DBG("BED '%s' %d:%d: wrong igb code", filename, i, j);
                exit(-1);
              }
              break;
            case 9:
              item->blockCount = atoi(str);
              if(!item->blockCount && str[0] != '0') {
                DBG("BED '%s' %d:%d: %s: atoi failed", filename, i, j, str);
                exit(-1);
              }
              break;
            case 10:
              k = 0;
              pch = strtok(str, ",");
              while(pch != NULL) { 
                item->blockSizes = ALLOCMEMORY(space, item->blockSizes, Uint, k+1);
                item->blockSizes[k] = atoi(pch);
                if(!item->blockSizes[k] && pch[0] != '0') {
                  DBG("BED '%s' %d:%d: %s: atoi failed", filename, i, j, pch);
                  exit(-1);
                } 
                 k++;
                pch = strtok(NULL, ",");
              }
              if(k != item->blockCount) {
                DBG("BED '%s' %d:%d: %d!=%d: wrong block count", filename, i, j, k, item->blockCount);
                exit(-1);
              }
              break;
            case 11:
              k = 0;
              pch = strtok(str, ",");
              while(pch != NULL) { 
                item->blockStarts = ALLOCMEMORY(space, item->blockStarts, Uint, k+1);
                item->blockRefseqs = ALLOCMEMORY(space, item->blockRefseqs, char*, k+1);
                item->blockStrands = ALLOCMEMORY(space, item->blockStrands, char, k+1);
                ulen = strlen(pch);
                for(u=0; u < ulen; u++) {
                  if(pch[u] == ':') break;
                }
                if(u < ulen) {
                  assert(u>0);

                  item->blockRefseqs[k] = ALLOCMEMORY(space, NULL, char, u+1);
                  memmove(item->blockRefseqs[k], pch, u);
                  item->blockRefseqs[k][u] = 0;
                  v = u+1;
                  
                  for(u=v; u < ulen; u++) {
                    if(pch[u]==':') break;
                  }
                  assert(u>v);
                  
                  tmp = ALLOCMEMORY(space, NULL, char, u-v+1);
                  memmove(tmp, &pch[v], u-v);
                  tmp[u-v] = 0;
                  item->blockStarts[k] = atoi(tmp);
                  
                  if(!item->blockStarts[k] && tmp[0] != '0') {
                    DBG("BED '%s' %d:%d: atoi failed while reading extension", filename, i, j);
                    exit(-1);
                  }
                  assert(pch[u+1]=='-' || pch[u+1] == '+');
                  item->blockStrands[k] = pch[u+1];
                  
                } else { 
                  item->blockStarts[k] = atoi(pch);
                  item->blockRefseqs[k] = NULL;
                  item->blockStrands[k] = item->strand;
                  if(!item->blockStarts[k] && pch[0] != '0') {
                    DBG("BED '%s' %d:%d: atoi failed", filename, i, j);
                    exit(-1);
                  }
                }
                k++;
                pch = strtok(NULL, ",");
              }
              if(k != item->blockCount) {
                DBG("BED '%s' %d:%d: wrong block count", filename, i, j);
                exit(-1);
              }
              break;
            default:
              DBG("'%s' not in BED format\n", filename);
              exit(-1);
              break;
          }
        }

        track->noofitems++;
      }

    }
    destructStringset(space, set[i]);
  }

  /* qsort(track->items, track->noofitems, sizeof(annotationitem_t), */ 
  /*     bl_annotationitem_cmp); */

  /* bl_annotationtrackAssignTrackLevel(track); */

  FREEMEMORY(space,set);
  return track;
}

/*------------------------------- bl_BEDwrite --------------------------------
 *    
 * @brief write a bed to a file
 * @author Steve Hoffmann 
 *   
 */

  void
bl_BEDwrite (annotationitem_t *items, Uint n,  FILE *fp)
{
  Uint i,j;
  annotationitem_t *b;


  for(i=0; i < n; i++) {
    b = &items[i];  
    fprintf(fp,"%s\t%"PRIu64"\t%"PRIu64"\t%u\t", b->chromname, b->start, b->end, b->trackid);
    if(b->name) { 
      fprintf(fp,"%s\t", b->name);
      if(b->score >=0) { 
        fprintf(fp, "%f\t", b->score);
        if(b->strand) {
          fprintf(fp, "%c\t", b->strand);
          if(b->thickStart) {
            fprintf(fp, "%"PRIu64"\t", b->thickStart);
            if(b->thickEnd) {
              fprintf(fp, "%"PRIu64"\t", b->thickEnd);
              if(b->itemRgb) { 
                fprintf(fp, "%d,%d,%d\t", b->itemRgb[0], b->itemRgb[1], b->itemRgb[2]);
              } else { 
                fprintf(fp, "0\t");
              }
              if(b->blockCount) {
                fprintf(fp, "%d\t", b->blockCount);
                if(b->blockSizes) {
                  for(j=0; j < b->blockCount; j++) { 
                    fprintf(fp, "%"PRIu64"", b->blockSizes[j]);
                    if (j < b->blockCount-1) fprintf(fp, ",");
                    else fprintf(fp, "\t");
                  }
                  if(b->blockStarts) {
                    for(j=0; j < b->blockCount; j++) { 
                      if(b->blockRefseqs && b->blockRefseqs[j]) {
                        fprintf(fp, "%s:%"PRIu64":%c", b->blockRefseqs[j], 
                            b->blockStarts[j], b->blockStrands[j]);
                      } else { 
                        fprintf(fp, "%"PRIu64"", b->blockStarts[j]);
                      }
                      
                      if (j < b->blockCount-1) fprintf(fp, ",");
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    fprintf(fp, "\n");
  } 
  return ;
}


