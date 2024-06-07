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


#ifndef _BEDFILES_
#define _BEDFILES_

/*
 *
 *	bedfiles.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 10/30/16 17:41:47 CET  
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stringutils.h"
#include "basic-types.h"
#include "charsequence.h"
#include "gzidx.h"
#include "biofiles.h"

annotationtrack_t* bl_BEDread (void *space, char *filename);
void bl_annotationitemDump(FILE *dev, annotationitem_t *item);
void bl_BEDwrite (annotationitem_t *items, Uint n, FILE *dev);

#endif
