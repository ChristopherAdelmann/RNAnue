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
 *  radixsort.c
 *  a radix sort implementation
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07.02.2010 23:26:52 CET
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <limits.h>
#include <string.h>
#include <basic-types.h>
#include <radixsort.h>



void
bl_radixSort(void *space, void *tosrt, 
			 size_t size, size_t nelem, 
			 Uint (*keyaccess)(void *), 
			 Uint bits) {
	
	char *p, *b, *src, *toSort;
	
	Uint mask, offset=0, i, key;
	Uint cntsize;
	Uint *cnt;

    toSort = (char*) tosrt;
	cntsize = 1 << bits;
	cnt = ALLOCMEMORY(space, NULL, Uint, cntsize);
	
	printf("alloc'd %d bins\n", cntsize);
	
	memset(cnt, 0, sizeof(Uint)*cntsize);
	b = src = malloc(size*nelem);
	
	mask =~ (UINT_MAX<<bits);
	
	for(; mask; mask <<= bits, offset+=bits) {
		for(p=toSort; p < toSort+(nelem*size); p+=size) {
			key = (keyaccess(p) & mask) >> offset;
			++cnt[key];
		}
		
		for(i=1; i < cntsize; ++i) {
			cnt[i]+=cnt[i-1];
		}
		
		for(p=toSort+((nelem-1)*size); p >= toSort; p-=size) {
			key = (keyaccess(p) & mask) >> offset;
			memmove(b+((cnt[key]-1)*size), p, size);
			--cnt[key];
		}
		
		p=b; b=toSort; toSort=p;
		memset(cnt, 0, sizeof(Uint)*cntsize);
	}
	
	if(toSort == src) memcpy(b, toSort, size*nelem);
	FREEMEMORY(space, src);
	FREEMEMORY(space, cnt);
	
	return;
}


void
bl_radixSortKeyFirst(void *space, void *tosrt, 
			 size_t size, size_t nelem, 
			 Uint bits) {
	
	char *p, *b, *src, *toSort;
	Uint *cast;
	
	Uint mask, offset=0, i, key;
	Uint cntsize;
	Uint *cnt;
	
    toSort = (char*) tosrt;
	cntsize = 1 << bits;
	cnt = ALLOCMEMORY(space, NULL, Uint, cntsize);

	memset(cnt, 0, sizeof(Uint)*cntsize);
	b = src = malloc(size*nelem);
	
	mask =~ (UINT_MAX<<bits);
	
	for(; mask; mask <<= bits, offset+=bits) {
		for(p=toSort; p < toSort+(nelem*size); p+=size) {
			cast = (Uint*)p;
			key = (*cast & mask) >> offset;
			++cnt[key];
		}
		
		for(i=1; i < cntsize; ++i) {
			cnt[i]+=cnt[i-1];
		}
		
		for(p=toSort+((nelem-1)*size); p >= toSort; p-=size) {
			cast = (Uint*)p;
			key = (*cast & mask) >> offset;
			memmove(b+((cnt[key]-1)*size), p, size);
			--cnt[key];
		}
		
		p=b; b=toSort; toSort=p;
		memset(cnt, 0, sizeof(Uint)*cntsize);
	}
	
	if(toSort == src) memcpy(b, toSort, size*nelem);
	FREEMEMORY(space, src);
	FREEMEMORY(space, cnt);
	
	return;
}


void
bl_radixSortUint(void *space, Uint *toSort, 
					 size_t nelem, 
					 Uint bits) {
	
	Uint *p, *b, *src;

	Uint mask, offset=0, i, key;
	Uint cntsize;
	Uint *cnt;
	
	cntsize = 1 << bits;
	cnt = ALLOCMEMORY(space, NULL, Uint, cntsize);
	
	memset(cnt, 0, sizeof(Uint)*cntsize);
	b = src = malloc(sizeof(Uint)*nelem);
	
	mask =~ (UINT_MAX<<bits);
	
	for(; mask; mask <<= bits, offset+=bits) {
		for(p=toSort; p < toSort+nelem; ++p) {
			key = (*p & mask) >> offset;
			++cnt[key];
		}
		
		for(i=1; i < cntsize; ++i) {
			cnt[i]+=cnt[i-1];
		}
		
		for(p=toSort+((nelem-1)); p >= toSort; --p) {			
			key = (*p & mask) >> offset;
			b[cnt[key]-1] = *p;
			--cnt[key];
		}
		
		p=b; b=toSort; toSort=p;
		memset(cnt, 0, sizeof(Uint)*cntsize);
	}
	
	if(toSort == src) memcpy(b, toSort, sizeof(Uint)*nelem);
	FREEMEMORY(space, src);
	FREEMEMORY(space, cnt);
	
	return;
}
