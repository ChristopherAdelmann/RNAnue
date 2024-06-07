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


#ifndef MAPFRAG_H
#define MAPFRAG_H

/*
 *
 *	mapfrag.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 10.12.2012 10:25:35 CET  
 *
 */


#include <string.h>
#include <limits.h>

#include "basic-types.h"
#include "locus.h"
#include "multicharseq.h"
#include "sufarray.h"
#include "kdseed.h"
#include "biofiles.h"
#include "segemehl.h"

typedef struct mapseed_s {
  uint64_t u;
  char rc; //strand TODO
  double evalue;
  double score;
  unsigned int len;
  unsigned int mat;
  unsigned int mis;
  unsigned int ins;
  unsigned int del;
  unsigned int edist;
  unsigned int l;
  unsigned int r;
  char good;
  char maxinterval;
  char maxevalue;
  locus_t *locus;
  /*these variables are not always set*/
  uint64_t seedlen;
  uint64_t refidx;
  uint64_t refpos;
  char *refname;
} mapseed_t;

typedef struct mapseedlist_s {
  unsigned int n;
  mapseed_t *l;
} mapseedlist_t;

typedef struct mapfrag_s {
  mapseed_t *seed;
  char *seq;
  char *qual;
  MultiCharSeqAlignment *mcsa;
  unsigned int lclip;
  unsigned int rclip;
  unsigned int mat;
  unsigned int mis;
  unsigned int ins;
  unsigned int del;
  unsigned int lmat;
  unsigned char mate;
  unsigned char mapq;
  int leftgap;
  int rightgap;
  unsigned char nextnoncollinear;
  unsigned char prevnoncollinear;
  unsigned char issplit;
  double mapq_dbl;
} mapfrag_t;

typedef struct mapping_s {
  char *seq;
  char *qual;
  Uint seqlen;
  Uint lclip;
  Uint rclip;
  int scr;
  unsigned int n;
  unsigned int P;
  unsigned int Q;
  char consecutive;
  char matestatus;
  mapfrag_t *f;
  double mapqual;
  double sigma;
  double maxqual;
} mapping_t;

typedef struct mappingset_s {
  unsigned int n;
  mapping_t* elem;
} mappingset_t;

void bl_removeMapping (mapping_t *l);
void bl_removeMappingSet (mappingset_t *s);
void bl_removeUnpairedMapping(mappingset_t *s);
char bl_mappingsetHasPaired (mappingset_t *s);
long long int bl_distMapping(mapping_t *l, mapping_t *r);
mapping_t* bl_concatMapping(mapping_t* l, mapping_t *r); 
void bl_initMappingSet(mappingset_t *set);
void bl_addMapping (mappingset_t *s, mapping_t *l);
Uint bl_addMapFrag(mapping_t *l, MultiCharSeqAlignment *mcsa, 
    mapseed_t *seed, unsigned char mate, unsigned char issplit) ;
mapseedlist_t*
bl_addMapSeedBranch (mapseedlist_t *s, branch_t *b, unsigned int u, unsigned int strand, 
    double maxE, unsigned int maxM, double SPM, karlin_t *stats);
unsigned int
bl_mappingsetHasPos (mappingset_t *s, unsigned int pos);
void bl_removeSuboptimalMapping (mappingset_t *s, int *scores, int indel);
char* bl_getMapFragQryDesc (mapfrag_t *f);
char* bl_getMapFragQry (mapfrag_t *f);
char* bl_getMapFragQual (mapfrag_t *f);
unsigned char bl_getMapFragStrand (mapfrag_t *f);
unsigned int bl_getMapFragV(mapfrag_t *f);
unsigned int bl_getMapFragU(mapfrag_t *f);
unsigned int bl_getMapFragP(mapfrag_t *f);
unsigned int bl_getMapFragQ(mapfrag_t *f);
char* bl_getMapFragRefDesc(mapfrag_t *f);
Alignment* bl_getMapFragAlignment (mapfrag_t *f);
unsigned int bl_getMappingEdist (mapping_t *l);
unsigned int bl_getMapFragEdist(mapfrag_t *f);
mapseedlist_t* bl_initMapSeedList (mapseedlist_t *l);
void bl_removeMappingSet (mappingset_t *s);
void bl_wrapSeedList (mapseedlist_t *l);
void bl_initMapping(mapping_t *l, char*, char*, Uint, Uint);
MultiCharSeqAlignment * bl_copyMCSA(MultiCharSeqAlignment *dest, MultiCharSeqAlignment *src);
char bl_isPairedMapping (mapping_t *l);
char bl_isQueryMapping (mapping_t *l);
char bl_isMateMapping (mapping_t *l);
char bl_hasMateMapping (mappingset_t *s);
char bl_hasQueryMapping (mappingset_t *s);
unsigned int bl_getMapFragChrIdx(mapfrag_t *f);
unsigned int bl_getMapFragSplit (mapfrag_t *f);
unsigned char bl_getMapFragIsMate (mapfrag_t *f);
MultiCharSeqAlignment* bl_getMapFragMCSA(mapfrag_t *f);
Uint bl_getMateStartIdx (mapping_t *l);
Uint bl_getQueryStartIdx (mapping_t *l);
mapseed_t *bl_getMapSeedListBest (mapseedlist_t *l);
locus_t * bl_getMapSeedLocus (mapseed_t *seed, Uint j, MultiCharSeq *mseq, Suffixarray *s);
void bl_getMapSeedAdditionalInformation(mapseed_t *seed, Suffixarray *s, MultiCharSeq *mseq);
void bl_addMapFragPartialEdist (mapfrag_t *frag, Uint edist);
unsigned int bl_getMapFragGetUlen (mapfrag_t *f);
uint64_t bl_getMapFragSubstart (mapfrag_t *f);
void bl_removeBadMappings (mappingset_t *s, Uint querylen, Uint matelen, Uint acc);
void bl_removeBadMates (mappingset_t *s, Uint querylen, Uint matelen, Uint acc);
void bl_countMultipleMappings (mappingset_t *s, Uint *nqueries, Uint *nmates);
unsigned int bl_getMappingEdistQM (mapping_t *l, Uint *qedist, Uint *medist);
void bl_setMapFragMapQual (mapfrag_t *f, double mapqual);
double bl_getMapFragMapQual (mapfrag_t *f);
unsigned int bl_getMapFragLongestMatch(mapfrag_t *f);
mapfrag_t* bl_getMapFrags (mapping_t *m, Uint *size);
mapping_t* bl_getMappings(mappingset_t *s, Uint *size);
unsigned int bl_getMapFragScore (mapfrag_t *f, int *scores, int indel);
void bl_removeMappingQM (mapping_t *l, char ismate);
Uint bl_getNextMapFragU (mapfrag_t *f, mapping_t *m, char ismate, Uint *lgap, Uint *rgap);
int bl_getMappingScore (mapping_t *l, int * scores, int indel);
mappingset_t* bl_copyMappingSet(mappingset_t *s);
locuslist_t* bl_getMappingLocusList(mapping_t *mapping, MultiCharSeq* mseq, char ismate);
char bl_isSplitMappingQM (mapping_t *mapping, char ismate);
char bl_getMappingStrandQM (mapping_t *mapping, char ismate);
char bl_isCollinearMapping(mapping_t *m, char ismate);
mapfrag_t* bl_getMapFragsQM (mapping_t *m, Uint *size, char ismate);
unsigned int bl_getMapFragLeft(mapfrag_t *f);
unsigned int bl_getMapFragRight(mapfrag_t *f);
int bl_getMappingScoreQM (mapping_t *l, int *scores, int indel, int *qscore, int *mscore);
int bl_getMappingFragNoQM (mapping_t *l, Uint *qno, Uint *mno);
char bl_hasQueryMappingMaxEdist (mappingset_t *s, Uint maxedist);
char bl_hasMateMappingMaxEdist (mappingset_t *s, Uint maxedist);
Uint bl_getMappingRange(mapping_t *m);
char bl_getMappingIsConsecutive(mapping_t *m);
void bl_dumpMappingSet (FILE *dev, mappingset_t *set);
int bl_getMappingMaxScore(mappingset_t *s, int *scores, int indel);
void bl_removeBadMatesEdistQM (mappingset_t *s, Uint maxedist, Uint ismate);
mappingset_t* bl_sortMappingSetByScore (mappingset_t *set, int* scores, int indel);
void bl_removeBadMatesAcc (mappingset_t *s, Uint querylen, Uint matelen, Uint acc); 
int bl_getMappingMaxScoreQM(mappingset_t *s, int *scores, int indel, char ismate);
void bl_removeSuboptimalMappingQM (mappingset_t *s, int *scores, int indel, char ismate);
void bl_mergeMappings (mappingset_t *set);
void bl_removeBadMatesCov (mappingset_t *s, Uint querylen, Uint matelen, Uint mc);
void bl_concatMappingSet (mappingset_t *dest, mappingset_t*source);
char bl_hasMappingPairedMaxEdist (mappingset_t *s, Uint maxedist);
mappingset_t* bl_mappingsetRemoveDuplicates(mappingset_t *set, Uint maxdist);
void bl_dumpMapping (FILE *dev, mapping_t *m);
mapping_t* bl_copyMapping (mapping_t* dest, mapping_t* source);
char bl_hasMultipleQueryMappings (mappingset_t *s);
char bl_hasMultipleMateMappings (mappingset_t *s);
char bl_hasMultiplePairedMappings (mappingset_t *s);
char bl_isCircularMappingQM (mapping_t *m, char ismate);
void bl_dumpSpliceJunctions(mapping_t *m, char ismate, MultiCharSeq *mseq, char *basename, char *name, segemehl_t *nfo);
char bl_mappingGetType (mapping_t *m, char ismate);
Uint  bl_getMappingMinFragDist (mapping_t *l);
char bl_getMapFragIsChimeric(mapfrag_t *f);
char bl_getMappingIsChimericQM (mapping_t *m, char ismate);
void bl_removeBadMatesScr (mappingset_t *s, Uint querylen, Uint matelen, int *scores, int indel, int minscr, int minmatescr);
char bl_mappingHasPairedEndOrientation(mapping_t *mapping);
#endif

