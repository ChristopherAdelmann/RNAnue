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


#ifndef SEGEMEHL_HELPER_H
#define SEGEMEHL_HELPER_H

/*
 *
 *	segemehl_helper.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 12/01/16 23:55:07 CET  
 *
 */
#include "segemehl.h"


typedef struct seseq_s {
  
  char *sequence;
  char *sequence_rc;
  char *quality;
  char *quality_re;
  char *convert;
  char *convert_rc;
  Uint len;
  char converttype;
  char converttype_rc;

} seseq_t;

void
se_getData(void *space, seseq_t *seq, char **seqs, char **quals, char bisulfite, char type);

void
se_segemehlSeqInit(void *space, seseq_t *seq, char *orig, char *qual, Uint len); 

void
getqualandseq(void *space, char *seq, char *qual, char **seqs, char **quals, Uint len);

void
convert(void *space, char **seqs, char *orig, Uint len, Uint phase, Uint type, segemehl_t *nfo);

void
se_segemehlSeqDestruct(void *space, seseq_t *seq);

#endif
