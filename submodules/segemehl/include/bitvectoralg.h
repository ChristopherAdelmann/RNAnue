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
 *	bitvectoralg.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 05/26/2008 12:26:03 PM CEST  
 *
 */
#include "basic-types.h"
#include "alignment.h"
#include "bitVector.h"

PairSint myersbitvector( void *space,
    char *query, 
    Uint qlen, 
    char *subject, 
    Uint slen, 
    char *alphabet, 
    Uint asize,
    Uint *enctab,
    Uint k,
    bitvector *peq);

bitvector*
myersbitmatrix( void *space,
    char *query, 
    Uint qlen, 
    char *subject, 
    Uint slen, 
    char *alphabet, 
    Uint asize,
    Uint *enctab,
    Uint k,
    bitvector *peq,
    PairSint *res,
    bitvector *D,
    Uint dim);

bitvector*
myersblockbitmatrix(
    void *space,
    char *query, 
    Uint qlen, 
    char *subject, 
    Uint slen, 
    char *alphabet, 
    Uint asize,
    Uint *enctab,
    Uint k,
    bitvector *peq,
    PairSint *res,
    bitvector *D,
    Uint dim);

bitvector*
getpeq(void *space,
    char *query, 
    Uint qlen,
    char *alphabet,
    Uint asize,
    Uint *enctab);

Uint*
encodetab(char *alphabet, Uint asize) ;


char*
getstringalphabet (void *space, char *string, Uint len, Uint *asize);

Alignment*
bitvectorbacktrack(Alignment *al, bitvector *D, Uint dim, Uint k, 
    Uint l, 
    char *subject, Uint *enctab, bitvector *peq); //patch

