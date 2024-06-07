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
 * container.h
 * implementation of a simple container for objects of defined size
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Tue Oct 14 16:31:33 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 73 $
 * Author: $Author: steve $
 * Date: $Date: 2008-10-29 10:03:28 +0100 (Wed, 29 Oct 2008) $
 * Id: $Id$
 * Url: $URL$
 */

#ifndef CONTAINER_H
#define CONTAINER_H

#include <stdlib.h>
#include "basic-types.h"

#define CONTINC 100
#ifndef BASEINC
#define BASEINC CONTINC
#endif

typedef struct {
  void *contspace;
  int nextfree;
  int allocelem;
  size_t sizeofelem;
} Container;

void bl_containerInit(Container *c, int allocelem, size_t sizeofelem);
void bl_containerDestruct(Container *c, void (*rmv) (void*));
BOOL bl_containerIsEmpty(Container *c);
void bl_containerResize(Container *c, int inc);
void bl_containerAdd(Container *c, void *elem);
void* bl_containerGet(Container *c, int n);
void bl_containerMerge(Container *c, Container *s);
Uint bl_containerSize(Container *c);

#endif /* CONTAINER_H */
