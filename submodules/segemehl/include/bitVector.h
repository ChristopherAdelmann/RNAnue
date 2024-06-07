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


#ifndef BITVECTOR_H
#define BITVECTOR_H

/*
 *
 *	bitVector.h
 *  declarations for bit arrays
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/14/2007 04:15:27 PM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 93 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-12-07 16:58:47 +0100 (Sun, 07 Dec 2008) $
 *
 *  Id: $Id: bitArray.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/bitArray.h $
 */

#include "basic-types.h"

#define BITVECTOR_WORDSIZE  (sizeof(unsigned long long int)*8)


typedef unsigned long long int bitvector_t;
typedef bitvector_t* bitvector;

extern bitvector initbitvector(void *, Uint length);
void dumpbitvector(bitvector a, Uint len);
unsigned char valbitvector(bitvector a, Uint len, unsigned char val);
extern void setbitvector(bitvector a, Uint len, unsigned char val);
bitvector resizebitvector(void *space, bitvector, Uint len);
void wrapBitmatrix(void *space, bitvector *, Uint m);

static inline void
bitvector_setbit(bitvector a, Uint pos, unsigned char val) {
  int  byte,
       bits;
  bitvector_t mask=0;

  byte  = pos/BITVECTOR_WORDSIZE;
  bits = pos & (BITVECTOR_WORDSIZE - 1);
  mask = 1;
  mask <<= bits;

  a[byte] ^= ((bitvector_t)-val ^ a[byte]) & mask;
}

static inline unsigned char
bitvector_getbit(bitvector a, Uint pos) {
  int byte;
  int bits;
  bitvector_t mask=0;


  byte = pos/BITVECTOR_WORDSIZE;
  bits  = pos & (BITVECTOR_WORDSIZE - 1);
  mask = 1;
  mask <<= bits;

  return ((bitvector_t)a[byte] & mask)? 1 : 0;
}


#endif
