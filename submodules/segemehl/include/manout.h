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


#ifndef MANOUT_H
#define MANOUT_H

/*
 * manout.h
 * attempt for flexible output of genome mapping w/ SEGEMEHL
 *
 * @author Christian Otto
 * @email christan@bioinf.uni-leipzig.de
 * @date Wed Sep 24 10:56:23 CEST 2008
 *
 */

#include "basic-types.h"
#include "biofiles.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "alignment.h"
#include "kdseed.h"
#include "mapfrag.h"
#include "segemehl.h"
#include "karlin.h"
#include "samio.h"

#define MINUSSTRAND 0
#define PLUSSTRAND 1
#define SPLIT_NEXT_PLUS    ((unsigned char) (1 << 5))
#define SPLIT_PREV_PLUS    ((unsigned char) (1 << 6))



typedef enum matchstatus_e {
  QUERY, 
  MATE, 
  PAIR, 
  PAIR_REV, 
  PAIR_INS, 
  QUERY_SPL_NO_MATE, 
  QUERY_SPL_FULL_MATE,
  MATE_SPL_NO_QUERY,
  MATE_SPL_FULL_QUERY,
  PAIR_SPL
} matchstatus_t;

typedef struct gmate_s {

  unsigned char isset;
  Uint p;
  Uint q;
  int scr;
  int mat;
  int mis;
  int ins;
  int del;
  int edist;
  char *materefdesc;
  Uint materefdesclen;
  char *materefseq;

  Alignment *al;
  Uint subject;
 
} gmate_t;

typedef struct gmatch_s{
  Uint subject;
  unsigned char rc;
  Uint i;
  Uint j;
  Uint p;
  Uint q;
  int scr; 
  int mat;
  int mis;
  int ins;
  int del;
  int edist;
  Alignment *al;
  double evalue;

  Uint noofmatematches;
  Uint mateminedist;
  gmate_t mates[4];

  Uint fragno;
  Uint previdx;
  Uint prevpos;
  char prevflags;
  Uint nextidx;
  Uint nextpos;
  char nextflags;

  Uint prevseqstart;
  char *prevseqrefdesc;
  Uint nextseqstart;
  char *nextseqrefdesc;
  char *refdesc;
  Uint refdesclen;
  char *refseq;


  unsigned char skip;
} gmatch_t;


typedef struct gmatchlist_s{

  Uint minedist;
  Uint mateminedist;
  Uint pairminedist;

  Uint *n;
  gmatch_t **matches;

} gmatchlist_t;

typedef struct gread_s{
  Uint id;
  Uint noofmatepairs;
  Uint noofmatches;

  Uint n[2];
  gmatch_t* matches[2];

} gread_t;


typedef struct Gmap_s{

  MultiCharSeq *mseq;
  Uint mapoffset;
  Uint noofreads;
  gread_t *reads;

} Gmap;

void
se_destructMatches(void *space, gread_t *read);

unsigned char
se_hasMatches(gmatchlist_t *list);

unsigned char
se_hasMateMatches(gmatchlist_t *list);

Uint
se_kdSetMate(void *space, gmatch_t *match, 
    Uint chr_idx, Uint chr_start, Uint chr_end, Uint edist,
    Alignment *al, unsigned char downstream, unsigned char rc);

gmatchlist_t*
se_kdMatchListAdd(gmatchlist_t *list, 
    Uint chr_idx, 
    Uint chr_start, 
    Uint chr_end, 
    Uint edist,
    int scr,
    Uint start,
    Uint end, 
    double evalue, Alignment *al, Uint u, 
    Uint previdx, Uint prevpos, char prevstrand, 
    Uint nextidx, Uint nextpos, char nextstrand, Uint fragno);

gmatchlist_t*
se_kdMatchListSet(void *space,
    gmatchlist_t *list, 
    Uint chr_idx, 
    Uint chr_start, 
    Uint chr_end, 
    Uint edist,
    int scr,
    Uint start,
    Uint end, 
    double evalue, Alignment *al, Uint u, Uint n);
  
void reportSplicedMatch(void *space, char *qrydesc, 
    MultiCharSeqAlignment *mcsa, Uint noofaligns,
    Uint coverage, Uint edist,  int score, segemehl_t *nfo);
Uint se_kdMatchListLength(gmatchlist_t *list, unsigned char strand);
unsigned char se_kdMatchListhasMatches(gmatchlist_t *list);
unsigned char se_kdMatchListhasMates(gmatchlist_t *list);
Uint se_kdMatchListLength(gmatchlist_t *list, unsigned char strand);
Uint se_kdMatchListScore(gmatchlist_t *list);
gmatch_t* se_kdMatchListGet(gmatchlist_t *list, unsigned char strand, 
    Uint elem);
gmate_t* se_kdMatchGetMates(gmatch_t *match);
Uint se_kdMatchGetSubject(gmatch_t *match);
Uint se_kdMatchGetRefStart(gmatch_t *match);
extern void reportMap(FILE*, Gmap *map, Uint level);
extern void initMatch(gmatch_t *);
void initRead(gread_t *, Uint);
void initGmap(Gmap *, MultiCharSeq *, Uint);
extern void setMatches(gread_t*, gmatch_t *, Uint, unsigned char, 
    unsigned char);
extern void setReads(Gmap *, gread_t *, Uint);
extern Uint reportMatch (void *, Gmap *, fasta_t *, segemehl_t *, 
    matchstatus_t pairStatus, unsigned char mate);
Uint se_setMatches(void *space, gread_t *read, gmatchlist_t *list, Uint maxedist, segemehl_t *nfo, char rep);
void matchHeader(FILE* dev, Uint level);
void genericOutput (FILE *dev, char **list, Uint rep_type, char);
void bl_gmatchlistDestruct(void *space, gmatchlist_t *list);
gmatchlist_t* bl_gmatchlistInit(void *space, int maxedist, int matemaxedist);
void se_openOutputDevices(void *space, segemehl_t *info);
void se_closeOutputDevices(void *space, segemehl_t *info);
bl_fileBins_t* se_createChromBins (void *space, fasta_t *f, int maxbins, char 
    *template, Uint tmplen);
bl_fileBinDomains_t* se_createChromDomains (void *space, fasta_t *f, 
    Uint minbins, Uint maxbins, char *filetemplate, Uint tmplen);
bl_fileBinDomains_t*
se_createBisulfiteBins (void *space, Uint noofdomains, Uint threadno, char *filetemplate, Uint tmplen);

char * se_defaultHeader (void *space, segemehl_t *info, char, char);
void se_storeHeader(void *space, char *filename, char **header, Uint *headerlen);
void se_output(mappingset_t *s, fasta_t * reads, unsigned int k, mapseed_t *, mapseed_t *, Suffixarray *arr, karlin_t *stats, segemehl_t *nfo);
void se_printMappingStats(FILE *device, segemehl_t *nfo);
void se_printStats(FILE *dev, segemehl_t *nfo); 
samheader_t* se_getSamHeaderStruct(segemehl_t *info);
samheader_t* se_header(segemehl_t *nfo, Uint binno) ;
#endif /* MANOUT_H */

