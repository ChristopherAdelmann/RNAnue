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

/*   This source code is derived  from zran.c from the zlib library 
 *   Namely, the structs have been taken from the zran file and wrapped in
 *   typedefs 
 */

/*  This file has been created on the basis of zran.c form the
    zlib library that is under the following licence:

 
   (C) 1995-2017 Jean-loup Gailly and Mark Adler

   This software is provided 'as-is', without any express or implied
   warranty.  In no event will the authors be held liable for any damages
   arising from the use of this software.
 
   Permission is granted to anyone to use this software for any purpose,
   including commercial applications, and to alter it and redistribute it
   freely, subject to the following restrictions:
 
   1. The origin of this software must not be misrepresented; you must not
      claim that you wrote the original software. If you use this software
      in a product, an acknowledgment in the product documentation would be
      appreciated but is not required.
   2. Altered source versions must be plainly marked as such, and must not be
      misrepresented as being the original software.
   3. This notice may not be removed or altered from any source distribution.

   Jean-loup Gailly        Mark Adler
   jloup@gzip.org          madler@alumni.caltech.edu



   If you use the zlib library in a product, we would appreciate *not* receiving
   lengthy legal documents to sign.  The sources are provided for free but without
   warranty of any kind.  The library has been entirely written by Jean-loup
   Gailly and Mark Adler; it does not include third-party code.

   If you redistribute modified sources, we would appreciate that you include in
   the file ChangeLog history information documenting your changes.  Please read
   the FAQ for more information on the distribution of modified source versions.


   Illustrate the use of Z_BLOCK, inflatePrime(), and inflateSetDictionary()
   for random access of a compressed file.  A file containing a zlib or gzip
   stream is provided on the command line.  The compressed stream is decoded in
   its entirety, and an index built with access points about every SPAN bytes
   in the uncompressed output.  The compressed file is left open, and can then
   be read randomly, having to decompress on the average SPAN/2 uncompressed
   bytes before getting to the desired block of data.

   An access point can be created at the start of any deflate block, by saving
   the starting file offset and bit of that block, and the 32K bytes of
   uncompressed data that precede that block.  Also the uncompressed offset of
   that block is saved to provide a referece for locating a desired starting
   point in the uncompressed stream.  build_index() works by decompressing the
   input zlib or gzip stream a block at a time, and at the end of each block
   deciding if enough uncompressed data has gone by to justify the creation of
   a new access point.  If so, that point is saved in a data structure that
   grows as needed to accommodate the points.

   To use the index, an offset in the uncompressed data is provided, for which
   the latest accees point at or preceding that offset is located in the index.
   The input file is positioned to the specified location in the index, and if
   necessary the first few bits of the compressed data is read from the file.
   inflate is initialized with those bits and the 32K of uncompressed data, and
   the decompression then proceeds until the desired offset in the file is
   reached.  Then the decompression continues to read the desired uncompressed
   data from the file.

   Another approach would be to generate the index on demand.  In that case,
   requests for random access reads from the compressed data would try to use
   the index, but if a read far enough past the end of the index is required,
   then further index entries would be generated and added.

   There is some fair bit of overhead to starting inflation for the random
   access, mainly copying the 32K byte dictionary.  So if small pieces of the
   file are being accessed, it would make sense to implement a cache to hold
   some lookahead and avoid many calls to extract() for small lengths.

   Another way to build an index would be to use inflateCopy().  That would
   not be constrained to have access points at block boundaries, but requires
   more memory per access point, and also cannot be saved to file due to the
   use of pointers in the state.  The approach here allows for storage of the
   index in a file.
 */


#ifndef GZIDX_H
#define GZIDX_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sys/types.h>
#include "zlib.h"

#define SPAN 1048576L       /* desired distance between access points */
#define WINSIZE 32768U      /* sliding window size */
#define CHUNK 16384         /* file input buffer size */
#define LARGECHUNK 1638400000
#define MEDIUMCHUNK 1638400
//#define CHUNK 0x10000         /* file input buffer size */
//#define LARGECHUNK 0x10000
//#define MEDIUMCHUNK 0x10000
//#define SPAN 0x10000



/* access point entry */
typedef struct point {
    off_t out;          /* corresponding offset in uncompressed data */
    off_t in;           /* offset in input file of first full byte */
    int bits;           /* number of bits (1-7) from byte at in - 1, or 0 */
    unsigned char window[WINSIZE];  /* preceding 32K of uncompressed data */
} gzpoint_t;

/* access point list */
typedef struct access {
    char type;                  /* 0 - for gzip/zlib, 1 - for bgzip */
    int have;           /* number of list entries filled in */
    int size;           /* number of list entries allocated */
    struct point *list; /* allocated list */
    off_t end;
} gzaccess_t;

typedef struct gzidxfile {
  FILE *fp;
  struct access *index;
  off_t curap;
  int mychunk;
  unsigned char *buf;
  unsigned char *pos;
  int len;
} gzidx_t;

/* Deallocate an index built by build_index() */
void free_index(struct access *index);

//int extract(FILE *in, struct access *index, off_t offset,
//                  unsigned char *buf, int len);
//int extract2(FILE *in, struct access *index, off_t offset,
//                  unsigned char *buf, int len);


//gzaccess_t* bl_zranGetIndex(char *filename, int *len, off_t myoffset);

void bl_destructgzidxfile(struct gzidxfile *file);

off_t bl_ftellgzidx(struct gzidxfile *f);

//int build_index(FILE *in, off_t span, struct access **built);

struct gzidxfile* bl_initgzidxfile(FILE *fp, struct access *index, off_t offset, int len);


int bl_getgzidxc (struct gzidxfile *f);

//struct access *addpoint(struct access *index, int bits,
//    off_t in, off_t out, unsigned left, unsigned char *window, unsigned char isbgzip);
//int bl_zranParseBlocked2(FILE *fp, struct buildidx *idx);
// int bl_zranDecomp(FILE *fp, z_stream *strm, unsigned char* input, unsigned char *window, 
//                  unsigned char *how, struct buildidx *idx,  char *eof); 
//int bl_zranFetch(FILE *fp, z_stream *strm, unsigned char *input,
//                 unsigned char *window, int windowsz, char *eof,
//                 unsigned char *how, struct buildidx *idx, unsigned char rawinflate);
//int bl_zranLook(FILE *fp, z_stream *strm, unsigned char *input, char *eof, 
//    unsigned char *how, char rawinflate);

gzaccess_t* bl_gzGetIndex(char *filename);
uint64_t bl_zranExtract(FILE *fp, gzaccess_t *idx, off_t off, unsigned char *buf, uint32_t len);
#endif
