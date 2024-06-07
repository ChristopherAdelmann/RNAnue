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


#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

#include <stdint.h>

#ifdef __CYGWIN__ 

#define CRLF '\r'

#else
#define CRLF ' '
#endif

#define MAXBUFFERSIZE 10000
#define BASEINC 10000
#define MAX_INT_LENGTH 50
#define MAGIC_INIT 0x0BADF00D

typedef unsigned char Uchar;
typedef unsigned int Uint;
typedef signed long long Lint;
typedef signed long long int LLint;
typedef signed int Sint;
typedef unsigned char BOOL;

#define True 1
#define False 0

#ifndef TRUE
#define TRUE True
#endif

#ifndef FALSE
#define FALSE False
#endif

typedef struct {
  int  a, 
       b;
} PairSint; 

typedef struct{
  Lint a,
       b;
} PairLSint;

typedef struct {
  Uint  a, 
       b;
} PairUint; 


typedef struct {
  Uint a,
      b,
      c;
} TripleUint;


typedef struct {
  int a,
      b,
      c;
} TripleSint;

typedef struct {
  int a,
      b,
      c,
      d;
} QuadSint;


typedef struct {
  Uint a,
      b,
      c,
      d;
} QuadUint;


#endif

