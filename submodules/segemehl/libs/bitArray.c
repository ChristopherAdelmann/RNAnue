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
 *  bitArray.c
 *  implementations
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/14/2007 04:15:14 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 33 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-27 10:00:08 +0200 (Tue, 27 May 2008) $
 *
 *  Id: $Id: bitArray.c 33 2008-05-27 08:00:08Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/bitArray.c $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "memory.h"
#include "bitArray.h"

bitarray
initbitarray(void *space, Uint len) {
  Uint n;

  bitarray a;

  n = len/8;
  n += (len % 8 > 0) ? 1 : 0;
  fprintf(stderr, "init bit array of %u\n", n);
  a = ALLOCMEMORY(space, NULL, unsigned char, n);

  return a;
}

bitarray
resizebitarray(void *space, bitarray a, Uint len) {
  Uint n;

  n = len/8;
  n += (len % 8 > 0) ? 1 : 0;
  a = ALLOCMEMORY(space, a, unsigned char, n);

  return a;
}

void
setbitarray(bitarray a, Uint len, unsigned char val) {
  Uint n;

  n = len/8;
  n += (len % 8 > 0) ? 1 : 0;
  memset(a,((val) ? 255 :0), n);
}

unsigned char
valbitarray(bitarray a, Uint len, unsigned char val) {
    Uint i;
    char *array;

    array = (char*) a;

    for(i=0; i < (len/8); i++) {  
        if (array[i] != (char) 255) 
          return 0;
    }
    for(i=0; i < (len%8); i++){
        if (getbit(a, len-i-1)!= val) 
          return 0;
    }
    
    return 1;
}

void
dumpbitarray(bitarray a, Uint len) {
  Uint i;
 
  for(i=0; i < len; i++) {
    printf("%d ", getbit(a, i));    
  }
  printf("\n");
}

