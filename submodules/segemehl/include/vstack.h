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


/**
 * vstack.h
 * implementation of a simple stack for objects of defined size
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Fri Oct 10 11:37:36 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 69 $
 * Author: $Author: steve $
 * Date: $Date: 2008-10-16 15:10:07 +0200 (Thu, 16 Oct 2008) $
 * Id: $Id$
 * Url: $URL$
 */

#ifndef VSTACK_H
#define VSTACK_H

#include <stdio.h>
#include <stdlib.h>
#include "basic-types.h"

#define VSTACKINC 10000

#ifndef BASEINC
#define BASEINC VSTACKINC
#endif

typedef struct{
  void* stackspace;
  Lint allocelem;
  Lint top;
  size_t sizeofelem;
} VStack;

void bl_vstackInit(VStack *s, Lint allocelem, size_t sizeofelem);
void bl_vstackDestruct(VStack *s, void (*rmv)(void*));
BOOL bl_vstackIsEmpty(VStack *s);
void bl_vstackPush(VStack *s, void *elem);
void* bl_vstackTop(VStack *s);
void* bl_vstackTopN(VStack *s, Lint n);
void* bl_vstackPop(VStack *s, void (*rmv)(void*));
Lint bl_vstackSize(VStack *s);

#endif /* VSTACK_H */
