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


#ifndef SEQCLIP_H
#define SEQCLIP_H

/*
 *
 *	seqclip.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 24.04.2010 22:32:24 CEST  
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "biofiles.h"
#include "basic-types.h"


Uint
bl_seqclipPolyA(void *space, char *sequence, Uint len, char* clp, Uint clen);

Uint
bl_seqclipSoft3Prime(void *space, char *sequence, Uint len, 
    char *toClip, Uint clipsize, Uint minclipacc, Uint pAlen);

Uint
bl_seqclipSoft5Prime(void *space, char *s, Uint len, 
    char *C, Uint clen, Uint minclipscr);

Uint
bl_seqclipHard3Prime(Uint len, Uint clipsize);

char*
bl_seqclipHard5Prime(char *s, Uint len, Uint clen);

char*
bl_seqclipFind3Prime (void *space, fasta_t *set, Uint samplesize, Uint fs, int ws);


#endif
