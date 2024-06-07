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
 *
 *	pigeon.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 10.09.2016 19:16:49 CEST  
 *
 */
#ifndef PIGEON_H
#define PIGEON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "segemehl.h"
#include "memory.h"
#include "fileio.h"
#include "stringutils.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "sufarray.h"
#include "mmchar.h"
#include "mathematics.h"
#include "biofiles.h"
#include "vtprogressbar.h"
#include "sort.h"
#include "bitArray.h"
#include "vqueue.h"
#include "vstack.h"
#include "container.h"
#include <pthread.h>
#include "kdseed.h"
#include "info.h"
#include "debug.h"
#include "mapfrag.h"
#include "alignment.h"
#include <assert.h>
#include "iupac.h"
#include "bitvectoralg.h"
#include "alignment.h"
#include "matealign.h"
#include "segemehl_helper.h"

mappingset_t*
pigeon(void *space, mappingset_t *set, Suffixarray *s, MultiCharSeq *mseq, 
    seseq_t *data, Uint nsplits, char *qname,
    Uint *enctab, bitvector* D, char ismate, Uint nchars, segemehl_t *nfo);


#endif
