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


#ifndef RADIXSORT_H
#define RADIXSORT_H

/*
 *  radixsort.h
 *  segemehl
 *
 *  Created by Steve Hoffmann on 08.02.10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <basic-types.h>

void
bl_radixSort(void *space, void *toSort, 
			 size_t size, size_t nelem, 
			 Uint (*keyaccess)(void *), 
			 Uint bits);

void
bl_radixSortKeyFirst(void *space, void *toSort, 
					 size_t size, size_t nelem, 
					 Uint bits);

void
bl_radixSortUint(void *space, Uint *toSort, 
				 size_t nelem, 
				 Uint bits);

#endif
