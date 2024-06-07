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


#ifndef MERGE_H
#define MERGE_H

/*
 * merge.h
 * functions to merge matches
 *
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fileBins.h"
#include "biofiles.h"
#include "samio.h"
#include "segemehl.h"

#define NOT_PAIRED 0
#define PAIRED 1

typedef struct {
  char* qname;
  int matchid;
  samrec_t *read;
  samrec_t *mate;
} bl_mergefilematch_t;

typedef struct {
  /* file pointer */
  FILE *fp;
  /* EOF read */
  unsigned char eof;
  /* current read alignment */
  bl_mergefilematch_t *entry;
  /* current entry complete? */
  unsigned char complete;
  char *buffer;
  circbuffer_t cb;
  pthread_mutex_t *mtx;
} bl_mergefile_t;

typedef struct {
  Uint nooffiles;
  bl_mergefile_t *f;
} bl_mergefiles_t;


typedef struct {

  bl_fileBinDomains_t *dms; 
  fasta_t *reads;
  samheader_t *head;
  FILE *dev; 
  bl_fileBinDomains_t *odms; 
  Uchar remove;
  Uint bestonly; 
  mappingstats_t *stats; 
  segemehl_t *nfo;
  mergemap_t *map;
  Uint mapsz;
  bl_mergefiles_t *files; 
  Uint cur;
  pthread_mutex_t *mastermtx;
  Uint *queuebeg;
  Uint processed;
} mergeinfo_t;


void bl_mergefilesInit(bl_mergefiles_t *files, Uint nooffiles);
void bl_mergefilesDestruct(bl_mergefiles_t *files);
void bl_mergefileInit(bl_mergefile_t *file, FILE *fp);
void bl_mergefileDestruct(bl_mergefile_t *file);
void bl_mergefilematchInit(bl_mergefilematch_t *match);
int bl_mergeCompareUnmapped(samrec_t *i, samrec_t *j);
int bl_mergefilematchComparePairingState(bl_mergefilematch_t *i, bl_mergefilematch_t *j, Uint *pairingstate);
void bl_mergefilematchRemoveSuboptimalPairs(bl_mergefilematch_t **list, Uint n);
void bl_mergefilematchRemoveSuboptimalSingletons(bl_mergefilematch_t **list, Uint n);
unsigned char bl_mergefileFastaIDCompare(char *desc1, Uint desc1len, char *desc2, Uint desc2len);
void bl_mergefilematchDestruct(bl_mergefilematch_t *match);
unsigned char bl_mergeParseLine(samheader_t* head, bl_mergefilematch_t *match, char *line, Uint *len);
void bl_mergeReadNext(samheader_t *head, bl_mergefile_t *file);
void bl_mergeUpdateTag(samrec_t *rec, Uint matchid, Uint noofmatches);
void bl_mergeCountAlignedMappings(bl_mergefilematch_t **list, int n, Uint *nqueries, Uint *nmates);
void bl_mergeCountMappings(bl_mergefilematch_t **list, int n, Uint *nqueries, Uint *nmates);
void se_mergeBisulfiteBins (bl_fileBinDomains_t *bsdomains, fasta_t *reads, samheader_t *head,
			    FILE *dev, bl_fileBinDomains_t *chrdomains, unsigned char remove,
                            Uint bestonly, mappingstats_t *stats, segemehl_t *nfo);
void se_mergeBisulfiteBinsThreaded(bl_fileBinDomains_t *dms, fasta_t *reads, 
    samheader_t *head, FILE *dev, bl_fileBinDomains_t *odms, Uchar remove, 
    Uint bestonly, mappingstats_t *stats, segemehl_t *nfo); 
void bl_mergefilematchRemoveRedundantUnmapped(bl_mergefilematch_t **list, Uint n);
void bl_mergefilematchRemoveNextInformation(bl_mergefilematch_t **list, Uint n);
Uint bl_mergeGetReadDevices(bl_fileBinDomains_t *dms, Uint k, Uint *list, segemehl_t *nfo);
void se_mergeStatsUpdate( bl_mergefilematch_t **list, Uint n, Uint pairingstate, 
    char hasMate, mappingstats_t *stats, pthread_mutex_t *stats_mtx);
void se_mergeWriteOutput( bl_mergefilematch_t **best, Uint noofbest, 
    mappingstats_t *stats, bl_fileBinDomains_t *odms, segemehl_t* nfo);
int cmp_mergemap_qsort(const void *a, const void *b); 
void se_mergeBisulfiteBinsNew(bl_fileBinDomains_t *dms, fasta_t *reads, 
    samheader_t *head, FILE *dev, bl_fileBinDomains_t *odms, Uchar remove, 
    Uint bestonly, mappingstats_t *stats, segemehl_t *nfo); 
void se_simulate(bl_fileBinDomains_t *dms, fasta_t *reads, 
    samheader_t *head, FILE *dev, bl_fileBinDomains_t *odms, Uchar remove, 
    Uint bestonly, mappingstats_t *stats, segemehl_t *nfo);
Uint se_mergeComplexMaster(mergeinfo_t *mi); 
mergeinfo_t * se_initMergeInfo(bl_fileBinDomains_t *dms, fasta_t *reads, 
    samheader_t *head, FILE *dev, bl_fileBinDomains_t *odms, Uchar remove, 
    Uint bestonly, mappingstats_t *stats, segemehl_t *nfo, mergemap_t *map, 
    Uint mapsz, bl_mergefiles_t *files, pthread_mutex_t *mtx); 
#endif /* MERGE_H */
