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
 *	queryalign.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 29.07.2013 17:57:49 CEST  
 *
 */

#ifndef QUERYALIGN_H
#define QUERYALIGN_H

#include "mapfrag.h"
#include "multicharseq.h"

 mapseedlist_t*
bl_getGoodSeeds (matchstem_t **items, unsigned m, unsigned n, karlin_t *stats, 
    segemehl_t *nfo);

mappingset_t*
bl_seedAlign(Suffixarray *s, mappingset_t* set, MultiCharSeq *mseq, 
    char **seqs, char **qual, Uint len, char* qname, mapseedlist_t *l, 
    segemehl_t *nfo, Uint *enctab, bitvector* D, char ismate);


#endif

