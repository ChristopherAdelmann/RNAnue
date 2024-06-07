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


#ifndef SAMHEADER_H
#define SAMHEADER_H

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

typedef struct samheader_s{
  char *version;
  char **rnames;
  uint64_t *rlens;
  uint64_t nrnames;
  char **rgroups;
  char **rgroupsinfo;
  Uint nrgroups;
  char *cmd;
} samheader_t;

samheader_t* bl_samreadHeader(char* filename);
    //unsigned char gzip, struct access *index,)

void bl_samaddReadGroup (samheader_t *head, char *id, char *info);
void bl_saminitHeader (samheader_t *head);
void bl_samdumpHeader (samheader_t  *head);
void bl_samdestructHeader (samheader_t *head);
  char* bl_samHeader (void *space, segemehl_t *info, Uint binno);
char* bl_samwriteHeader(samheader_t *hdr, Uint binno, char order, char sep, char lf) ;
samheader_t* bl_samparseHeaderLine (samheader_t *head, char *line);
void bl_samaddReferenceSequence (samheader_t *head, char *name, uint64_t len);

#endif
