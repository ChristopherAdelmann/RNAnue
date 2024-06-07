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

#ifndef MEMORY_H
#define MEMORY_H
/*
#include "memmac.h"
#include "memman.h"
*/
#include <assert.h>
#include <stdlib.h>

#define ALLOCMEMORY(X, PTR, TYPE, SIZE) bl_realloc((PTR), sizeof(TYPE) * (SIZE))
#define CALLOCMEMORY(X, TYPE, SIZE) bl_realloc((PTR), (SIZE), sizeof(TYPE))

#define FREEMEMORY(X, PTR) \
    free(PTR);             \
    PTR = NULL

void* bl_realloc(void*, size_t);
void* bl_calloc(void*, size_t, size_t);

#endif
