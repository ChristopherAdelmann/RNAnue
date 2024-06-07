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
 *  junctions.c
 *  routines to report and store split read junctions
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 10/30/16 19:57:19 CET
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "biofiles.h"
#include "mapfrag.h"
#include "mathematics.h"
#include "segemehl.h"


void
bl_getSpliceJunctionsFromMappingSet(mappingset_t *set, MultiCharSeq *mseq, char *readname, segemehl_t *nfo) {

  Uint i,k;

  for(i=0; i < set->n; i++) {
    for(k=0; k < 2; k++) {

     
    //  bl_dumpMappingSet(nfo->multisplitdev, set);
      bl_dumpSpliceJunctions(&set->elem[i], k, mseq, nfo->readgroupid, readname, nfo);
    
       }
  }

  return;
}


