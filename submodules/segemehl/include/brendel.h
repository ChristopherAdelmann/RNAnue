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
 *	brendel.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/01/2016 03:25:21 PM CEST  
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include "basic-types.h"
#include "alignment.h"
#include "memory.h"
#include "mathematics.h"
#include "iupac.h"
#include "biofiles.h"

Alignment* splicedaligndp (char *read, unsigned int m, char *genome, unsigned int n, gene_t **model);
mapping_t* bl_dpsplicealign2map(Alignment *al, gene_t *model, MultiCharSeq *mseq, Uint vpos, Uint vlen, char strand, char *querydesc, char *query, char *qual, Uint ulen, char ismate) ;
char bl_checkSpliceAlign(mapping_t *m);
Alignment* splicedaligndpopt (char *read, unsigned int m, char *genome, unsigned int n, gene_t **model, Uint a, Uint b, Uint l, Uint r);

