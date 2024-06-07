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


#ifndef SAMOUT_H
#define SAMOUT_H

#include "multicharseq.h"
#include "mapfrag.h"
#include "segemehl.h"
#include "biofiles.h"
#include "samheader.h"

/*
 *
 *	samout.h
 *  samout
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 21.03.2012 08:41:55 CET  
 *
 */

typedef struct samtag_s {
  char *tag;
  char *key;
  char *type;
  char *val;
} samtag_t;

typedef struct samrec_s {
  char *qname;
  unsigned flag;
  char *rname;
  uint64_t pos;
  uint8_t mapq;
  char *cigar;
  char *rnext;
  uint64_t pnext;
  int64_t tlen;
  char *seq;
  char *qual;
  Uint nooftags;
  samtag_t *tags;
} samrec_t;

typedef struct samlist_s {
  unsigned int noofrecs;
  samrec_t *recs;
} samlist_t;

/*
typedef struct mapfragment_s {

  char *id;
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
  char *seq;
  char *qual;
  char *rname;

  char issplit;
  char usenextfieldforsplits;

  Uint fragno;
  Uint previdx;
  Uint prevpos;
  char prevflags;
  Uint nextidx;
  Uint nextpos;
  char nextflags;

  unsigned char skip;

} mapfrag_t;

typedef struct {
  Uint nooffrags;
  mapfrag_t *frags;

} mapfraglist_t;*/


samlist_t *bl_samgetSamList (fasta_t *reads, Uint id, mapping_t *l, MultiCharSeq *mseq, char ismultiple, Uint noofqueries, Uint noofmates, Uint queryid, Uint mateid, char unmappedtornext, segemehl_t *nfo); 
void bl_samprintSamrec (FILE *dev, samrec_t* r, char lf);
char* bl_samprintSamrec2Buffer (samrec_t* r, char lf);
void bl_samprintSamlist (samlist_t *l, mapping_t *f, segemehl_t*);
void bl_samdestructSamList (samlist_t *list);
mappingset_t* bl_sammappingJoinFrags (mappingset_t *s, segemehl_t *nfo);
void bl_samprintEmptyAlign (char *desc, char* seq, char* qual, char hasPaired, char isQuery, char nomatemapped, char *nextchr, int64_t nextrpos, char nextrc, char ismultiple, char ischimeric, mapseed_t *seed, segemehl_t *nfo);
samrec_t * bl_samline2rec(char *line, Uint len, samheader_t *head);
void bl_samDestruct(samrec_t *samrec);

samtag_t* bl_samgetTag (samrec_t *rec, char* key);
#endif
