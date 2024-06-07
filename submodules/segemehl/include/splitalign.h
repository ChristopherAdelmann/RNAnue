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


#ifndef SPLITALIGN_H
#define SPLITALIGN_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

#include "basic-types.h"
#include "bitvectoralg.h"
#include "karlin.h"
#include "mapfrag.h"

/*
 *
 *	splitalign.h
 *  splitaligns
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/23/14 14:36:19 CEST  
 *
 */


/*
typedef struct split_s {
  Uint subidx;
  char strand;
  Uint start;
  Uint end;
  uint16_t i;
  uint16_t j;
} split_t;
*/

typedef struct{
  //the coords of the read
  uint64_t uoff;
  uint64_t ulen; 
  //the coords of the reference
  uint64_t voff;
  uint64_t vlen; 
  //the alignment
  char *cigar;

} split_t;


typedef struct{
  char *u;
  uint64_t ulen;
  uint64_t upos;
  uint64_t vpos;
  Uint vidx;
  char rc;
 
  //trans alignments pointing to here
  uint64_t vdonpos;
  Uint vdonidx;
  char vdonrc;
  //the linear splits of this alignment
  Uint noofsplits;
  split_t *splits;
  //trans alignments pointing out
  uint64_t vaccpos;
  Uint vaccidx;
  char vaccrc;
  
} splitalignment_t;



typedef struct spliceevent_s {
  uint64_t vdonpos;
  Uint vdonidx;
  char vdonrc;
  char donerr; //in percent NOT YET USED!
  char donins; //trailing insertions max:254 NOT YET USED!
  uint64_t vaccpos;
  Uint vaccidx;
  char vaccrc;
  char accerr; //in percent NOT YET USED!
  char accins; //trailing insertions max:254 NOT YET USED!
} spliceevent_t;

typedef struct spliceevents_s {
  Uint noofevents;
  spliceevent_t *event;
} spliceevents_t;


typedef struct spliceventmapelem_s{
  unsigned char type; 
  uint8_t site;
  spliceevent_t *ptr;
} spliceeventmapelem_t;

typedef struct spliceeventmap_s{
  Uint size;
  spliceeventmapelem_t *map;
}spliceeventmap_t;

mappingset_t* bl_splitAlign(void *space, Suffixarray *arr, MultiCharSeq *seq, mappingset_t *set, 
    char *querydesc, matchstem_t **stems, char **seqs, char **quals, Uint len, 
    karlin_t *stats, Uint *enctab, bitvector *D, unsigned char ismate, segemehl_t *nfo);
char* bl_cigarSplitAlignment(char *cigar, uint64_t **, uint64_t **, Uint *nsplits);
splitalignment_t * bl_cigarGetSplitAlignment(
    char* u, uint64_t ulen, uint64_t upos,
    uint64_t vpos, Uint idx, char rc, uint64_t vdonpos, Uint vdonidx, char donrc, 
    uint64_t vaccpos, Uint vaccidx, char accrc, char *cigar) ;
spliceevents_t* bl_splitAlignmentGetSpliceEvents (splitalignment_t *aln);
void bl_dumpSpliceEvents (spliceevents_t *events);
void bl_destructSplitAlignment (splitalignment_t *aln);
void bl_destructSpliceEvents (spliceevents_t *events);
void se_kdHMMsplit (char **seqs, char **quals, Uint qrylen, char *qrydesc, MultiCharSeq *mseq, 
    mappingset_t *set, Uint minlen, Uint minscore);
void bl_fixSplitAlignBrendel (Suffixarray *arr, mappingset_t *set, MultiCharSeq *mseq, char *qrydesc, 
    char **seqs, char **quals, Uint qrylen, char ismate, segemehl_t *nfo);
char bl_crosscorrection(Suffixarray *arr, mappingset_t *set, MultiCharSeq *mseq, 
    char *qrydesc, char **seqs, char **quals, Uint qrylen, Uint matelen, Uint checkedist, char ismate, segemehl_t *nfo);


void
se_kdAlignEvalSplitAlign (MultiCharSeq *seq, MultiCharSeqAlignment *a,  
    mappingset_t* set, Uint *totalcover, int *totalscore, 
    unsigned char *trans, int *scores, int indel, unsigned int noofaligns, 
    char ismate, char **seqs, char **quals, segemehl_t *nfo);

MultiCharSeqAlignment*
se_AlignSplitMap (Uint *mystarts, Uint *myends, uint64_t *mypos, char *myrc, Uint nooffrags, MultiCharSeq *seq, char *querydesc, 
    char **seqs, char **quals, Uint qrylen, int *scores, int indel, int transition);
  
mappingset_t*
bl_fixSplitAlignHoffmann (Suffixarray *arr, mappingset_t *set, MultiCharSeq *mseq, char *qrydesc, 
    char **seqs, char **quals, Uint qrylen, char ismate, segemehl_t *nfo);

#endif
