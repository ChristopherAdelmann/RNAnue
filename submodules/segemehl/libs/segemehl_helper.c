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
 *  segemehl_helper.c
 *  helper functions
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 12/01/16 23:54:00 CET
 *  
 */

#include "segemehl.h"
#include "karlin.h"
#include "mapfrag.h"
#include "bitvectoralg.h"
#include "sufarray.h"
#include "kdseed.h"
#include "mathematics.h"
#include "iupac.h"
#include "queryalign.h"
#include "matealign.h"
#include "splitalign.h"
#include "manout.h"
#include "vtprogressbar.h"
#include "samio.h"
#include "segemehl_helper.h"
void
se_segemehlSeqInit(void *space, seseq_t *seq, char *orig, char *qual, Uint len) {

  seq->convert = NULL;
  seq->convert_rc = NULL;
  seq->converttype = 0;
  seq->converttype_rc =0;
  seq->quality = NULL;
  seq->quality_re = NULL;

  //init sequence
  seq->sequence = ALLOCMEMORY(space, NULL, char, len+1);
  memmove(seq->sequence, orig, len);
  seq->sequence[len] = 0;
  seq->sequence_rc = charIUPACcomplement(space, seq->sequence, len);

  //init quality
  if(qual) { 
  seq->quality = ALLOCMEMORY(space, NULL, char, len+1);
  memmove(seq->quality, qual, len);
  seq->quality[len] = 0; 
  
  seq->quality_re = ALLOCMEMORY(space, NULL, char, len+1);
  memmove(seq->quality_re, qual, len);
  seq->quality_re[len] = 0;
  seq->quality_re = strrev(seq->quality_re, len);
  }

  seq->len = len;
}


void
se_segemehlSeqDestruct(void *space, seseq_t *seq) {
  if(seq->sequence) {
    FREEMEMORY(space, seq->sequence);
  }
  if(seq->sequence_rc) {
    FREEMEMORY(space, seq->sequence_rc);
  }
  if(seq->quality) {
    FREEMEMORY(space, seq->quality);
  }
  if(seq->quality_re) {
    FREEMEMORY(space, seq->quality_re);
  }
  if(seq->convert) {
    FREEMEMORY(space, seq->convert);
  }
  if(seq->convert_rc) {
    FREEMEMORY(space, seq->convert_rc);
  }

  seq->len = 0;
  seq->converttype = 0;
  seq->converttype_rc = 0;
}

char*
se_segemehlGetSequence(void *space, seseq_t *seq, char rc, char bisulfite, char type) {

  if(!bisulfite) {
    if(!rc) { 
      return seq->sequence;
    } else {
      return seq->sequence_rc;
    }
  } else {
    if(!rc) {
      if(!seq->convert || seq->converttype != type) {
        if(!seq->convert) {
          seq->convert = ALLOCMEMORY(space, NULL, char, seq->len+1);
        }
        memmove(seq->convert, seq->sequence, seq->len+1);
        bl_convertBisulfite(seq->convert, seq->len, bisulfite, type);
        seq->converttype = type;
      } 
      return seq->convert;
    } else {
       if(!seq->convert_rc || seq->converttype_rc != type) {
        if(!seq->convert_rc) {
          seq->convert_rc = ALLOCMEMORY(space, NULL, char, seq->len+1);
        }
        memmove(seq->convert_rc, seq->sequence_rc, seq->len+1);
        bl_convertBisulfite(seq->convert_rc, seq->len, bisulfite, type);
        seq->converttype_rc = type;
      } 
      return seq->convert_rc;
    }
  }

  return NULL;
}

char*
se_segemehlGetQuals(void *space, seseq_t *seq, char rc) {

  if(!rc) {
    return seq->quality;
  } else {
    return seq->quality_re;
  }

  return NULL;
}


void
se_getData(void *space, seseq_t *seq, char **seqs, char **quals, char bisulfite, char type) {
  seqs[0] = se_segemehlGetSequence(space, seq, 0, bisulfite, type);
  seqs[1] = se_segemehlGetSequence(space, seq, 1, bisulfite, type);
  quals[0] = se_segemehlGetQuals(space, seq, 0);
  quals[1] = se_segemehlGetQuals(space, seq, 1);

  return;
}


void
getqualandseq(void *space, char *seq, char *qual, char **seqs, char **quals, Uint len) {

  //get the query (and mate) and convert for bisulfite if necessary
  seqs[0] = seq;
  seqs[1] = charIUPACcomplement(space, seqs[0], len);
  //get the quals and reverse them
  if(qual) { 
    quals[0] = qual;
    quals[1] = ALLOCMEMORY(NULL, NULL, char, len+1);
    memmove(quals[1], quals[0], len+1);
    quals[1] = strrev(quals[1], len);
  } else {
    quals[0] = NULL;
    quals[1] = NULL;
  }
}


void
convert(void *space, char **seqs, char *orig, Uint len, Uint phase, Uint type, segemehl_t *nfo) {

  /* convert for alignment */
  //se_convert(seqs, bl_fastaGetSequence(reads,k), len, nfo, 0);
  /* convert for alignment */

  if(phase == 1 && nfo->bisulfite){    
    seqs[0] = ALLOCMEMORY(space, NULL, char, len+1);
    memmove(seqs[0], orig, len+1);
    bl_convertBisulfite(seqs[0], len, nfo->bisulfite, type);
    bl_convertBisulfite(seqs[1], len, nfo->bisulfite, type);
  }

  if (phase == 2 && nfo->bisulfite){
    FREEMEMORY(space, seqs[1]);
    memmove(seqs[0], orig, len+1);
    seqs[1] = charIUPACcomplement(space, seqs[0], len);
    bl_convertBisulfite(seqs[0], len, nfo->bisulfite, type); 
    bl_convertBisulfite(seqs[1], len, nfo->bisulfite, type); 
  }


}
