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

#include "portable_endian.h"
#include <stdio.h>
#include <inttypes.h>
#include <stdint.h>

inline uint16_t bl_endianReverse16(uint16_t value)
{
    return (((value & 0xFF) << 8) |
            ((value & 0xFF00) >> 8));
}

inline uint32_t bl_endianReverse32(uint32_t value) 
{
    return (((value & 0x000000FF) << 24) |
            ((value & 0x0000FF00) <<  8) |
            ((value & 0x00FF0000) >>  8) |
            ((value & 0xFF000000) >> 24));
}

inline uint16_t bl_letoh16 (uint16_t value) {
  #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return value;
 #elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return bl_endianReverse16(value);
 #else
    fprintf(stderr, "Unknown endianness. Does your compiler support __BYTE_ORDER__ ?");
    exit(EXIT_FAILURE);
 #endif
}



inline uint32_t bl_letoh32(uint32_t value) {
  #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return value;
 #elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return bl_endianReverse32(value);
 #else
    fprintf(stderr, "Unknown endianness. Does your compiler support __BYTE_ORDER__ ?");
    exit(EXIT_FAILURE);
 #endif
}

void bl_printEndianness(){

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  fprintf(stderr, "this cpu has LITTLE ENDIAN\n"); 
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
  fprintf(stderr, "this cpu has BIG ENDIAN\n");  
#else
  fprintf(stderr, "ENDIANNESS unknown.\n");
#endif
  return;

}

uint16_t bl_readLittleEndian16(uint16_t value){
  unsigned char *bytes;
  uint16_t val = 0;
  bytes = (unsigned char*) &value;
  val |= (bytes[0]);
  val |= (bytes[1] << 8);

  return val;
}

uint32_t bl_readLittleEndian32(uint32_t value){
  unsigned char *bytes;
  uint32_t val = 0;
  bytes = (unsigned char*) &value;
  val |= (bytes[0]); 
  val |= (bytes[1] << 8);
  val |= (bytes[2] << 16);
  val |= (bytes[3] << 24);

  return val;
}

