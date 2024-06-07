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


#ifndef ALIGNMENT_H
#define ALIGNMENT_H

/*
 *
 *	alignment.h
 *  alignment representation
 *  
 *  idea: 
 *  Stephan Kurtz, Gordon Gremme. Foundations
 *  of sequence analysis, University Hamburg, 2005
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 02/03/2009 11:56:27 AM CET  
 *
 */

#include "basic-types.h"
#include "biofiles.h"

typedef enum 
{
  Replacement, Deletion, Insertion, Skipped, Softclip, Match, Mismatch
} Eoptype;

typedef struct
{
  Eoptype eop;
  Uint steps;
} Multieop;

typedef struct {
  char *u;
  char *v;
  Uint ulen;
  Uint vlen;

  /*start of aligment (use in approx string matching, local align)*/
  Uint uoff;
  Uint voff;
  Multieop *meops;
  Uint numofmeops;
  unsigned char rmvseq;
  unsigned char rmuseq;
} Alignment;



void copyAlignment(Alignment *to, Alignment *from);
void showmultieoplist(FILE *dev, Alignment *al);
void showDynmultieoplist(Alignment* al, int size);
void showAlign(Alignment* al, FILE *dev);
void showAlignLF(Alignment* al, FILE *dev, char);
void initAlignment(Alignment *al, char *u, Uint ulen, Uint uoff, char *v, Uint vlen, Uint voff);
void insertEop(Alignment *al, Eoptype eop);
void revMeops(Alignment *al);
void wrapAlignment(Alignment *al);
Uint getEdist(Alignment *al);
Uint getBisulfiteMismatches(Alignment *al, Uint bisulfite);
Uint getWrongStrandBisulfiteMismatches(Alignment *al, Uint bisulfite);
void countEops(Alignment *al, Uint *mat, Uint *mis, Uint *ins, Uint *del, Uint *lmat);
char * multieopstring(Alignment *al, Uint leftpad, Uint rightpad, unsigned char rev);
Uint getUalignlen(Alignment *al);
Uint getValignlen(Alignment *al);
Uint getValignlenAndSkipped(Alignment *al);
int getSubstringEdist(Alignment *al, Uint u, Uint v);
int getAlignScore(Alignment *al, int *scores, int indel);
char* cigarstring(Alignment *al, Uint leftpad, Uint rightpad, char clipch, unsigned char rev, char brief);
char* mdstring(Alignment *al, unsigned char rev);
Uint bl_cigarGetAlignLen(char *cigar);
char* bl_mdGetDiffString(char *MD);
char* getNTcodekey(void *space);
void getSoftClipScores(Alignment *al, int polyAlen, int *scores, int indel, int *pAscr, int *adscr, int *adlen) ;
char* mdstrings(Alignment **al, Uint noofaligns, unsigned char rev);
char* bl_cigarGetAlignString(char *cigar, uint64_t **, uint64_t **, Uint *nsplits);
Uint* getSplitEdist(Alignment *al, Uint *noofsplits);
Uint* getUPartialAlignlen(Alignment *al, Uint *noofparts);
char* getEopString (Alignment *al);
char isMatch (Alignment *al, Uint i, Uint j);
void  showAlignModel(Alignment* al, FILE *dev, gene_t *model);
Alignment* getSubAlignment (Alignment *al, Uint l, Uint r);
Uint bl_alignGet5PrimeV (Alignment *al, char rc);
Uint bl_alignGet3PrimeU (Alignment *al, char rc);
Uint bl_alignGet5PrimeU (Alignment *al, char rc);
void insertMeop(Alignment *al, Eoptype eop, Uint steps);
Alignment* expandAlignment(Alignment *al, Uint* expandpos, Uint *expandlen);
Uint getPartialAlignNumber(Alignment *al);
int getUalignlenNoClip(Alignment *al);
Alignment* bl_getPartialAlignments (Alignment *al, char *vseq, Uint *noofparts);
char* getAlignString(Alignment* al, char lf);
Alignment* reevalAlignment(Alignment *al);
char predictStrand(Alignment *al, char *vseq); 
#endif
