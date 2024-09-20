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

#ifndef BITARRAY_H
#define BITARRAY_H

/*
 *
 *	bitArray.h
 *  declarations for bit arrays
 *
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig
 *  @date 07/14/2007 04:15:27 PM CEST
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: bitArray.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/bitArray.h $
 */

#include "basic-types.h"
typedef unsigned char *bitarray;

bitarray initbitarray(void *, Uint length);
void dumpbitarray(bitarray a, Uint len);
unsigned char valbitarray(bitarray a, Uint len, unsigned char val);
void setbitarray(bitarray a, Uint len, unsigned char val);
bitarray resizebitarray(void *space, bitarray, Uint len);
static inline void setbit(bitarray a, Uint pos, unsigned char val) {
    int bytes, bits, shift;
    unsigned char byte, resc;

    bytes = pos >> 3;
    bits = pos & 7;
    resc = a[bytes];

    byte = (unsigned char)val;
    shift = 7 ^ bits;

    byte = byte << shift;
    byte = byte ^ resc;
    byte = byte & (1 << shift);

    a[bytes] = byte ^ resc;
}

static inline unsigned char getbit(bitarray a, Uint pos) {
    int bytes, bits;
    unsigned char byte;

    bytes = pos >> 3;
    bits = pos & 7;

    byte = a[bytes];
    byte = byte >> (7 ^ bits);
    byte = byte << 7;
    byte = byte >> 7;
    byte = byte & 1;

    return byte;
}

#endif
