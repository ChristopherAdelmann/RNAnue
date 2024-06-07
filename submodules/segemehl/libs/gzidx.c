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
 *   Namely, the main functions build_index and extract and addpoint 
 *   have been taken from zran 
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include "basic-types.h"
#include "memory.h"
#include "zlib.h"
#include "gzidx.h"
#include "bgzip.h"
#include "info.h"


void show_window(FILE *dev, unsigned char *window, uint32_t left) {

  unsigned char output[WINSIZE+1];

  memset(output, '.', WINSIZE+1);

  if (left)
    memcpy(output, window + WINSIZE - left, left);
  if (left < WINSIZE)
    memcpy(output + left, window, WINSIZE - left);

  fprintf(dev, "---------- window follows ----------- \n%s\n", output);

  return;
}

/* Deallocate an index built by build_index() */
void free_index(struct access *index)
{
  if (index != NULL) {
    free(index->list);
    free(index);
  }
}

/* Add an entry to the access point list.  If out of memory, deallocate the
   existing list and return NULL. */

gzaccess_t* bl_gzAddpoint(gzaccess_t *index, int bits,
    off_t in, off_t out, unsigned left, unsigned char *window)
{
  struct point *next;
  Uint i, oldsize;

  /* if list is empty, create it (start with eight points) */
  if (index == NULL) {
    index = (struct access *)malloc(sizeof(struct access));
    if (index == NULL) return NULL;

    //gzip type
    index->type = 0;
    index->list = (struct point*)malloc(sizeof(struct point) << 3);

    for(i=0; i < 8; i++) {
      index->list[i].out = 0ul;
      index->list[i].in = 0ul;
      index->list[i].bits= 0;
      //either change point struct or allocate the memory
      //for the dictonary dynamically. This is necessary to
      //avoid allocs for bgzip files and still use the same
      //index structure.
      memset(index->list[i].window, 0, WINSIZE);
    }
    if (index->list == NULL) {
      free(index);
      return NULL;
    }
    index->size = 8;
    index->have = 0;
  }

  /* if list is full, make it bigger */
  else if (index->have == index->size) {
    oldsize = index->size;
    index->size <<= 1;
    next = (struct point*)realloc(index->list, sizeof(struct point) * index->size);
    if (next == NULL) {
      free_index(index);
      return NULL;
    }

    for(i=oldsize; i < index->size; i++) {
      next[i].out = 0;
      next[i].in = 0;
      next[i].bits= 0;
      //either change point struct or allocate the memory
      //for the dictonary dynamically. This is necessary to
      //avoid allocs for bgzip files and still use the same
      //index structure.
      memset(next[i].window, 0, WINSIZE);
    }
    index->list = next;
  }

  /* fill in entry and increment how many we have */
  next = index->list + index->have;
  next->bits = bits;
  next->in = in;
  next->out = out;

  if (left)
    memcpy(next->window, window + WINSIZE - left, left);
  if (left < WINSIZE)
    memcpy(next->window + left, window, WINSIZE - left);
  index->have++;


  /* return list, possibly reallocated */
  return index;
}

gzaccess_t * bl_bgzAddpoint(gzaccess_t *index, off_t in, off_t out) {

  uint64_t i, oldsize;
  gzpoint_t *next;

  /* if list is empty, create it (start with eight points) */
  if (index == NULL) {
    index = malloc(sizeof(gzaccess_t));

    if (index == NULL) {
      fprintf(stderr, "error during allocation of bgz index structure.\n");
      exit(EXIT_FAILURE);
    }

    //bgzip type
    index->type = 1;

    //acclocation of 8 points at once
    index->list = (gzpoint_t*) malloc(sizeof(gzpoint_t) << 3);

    for(i=0; i < 8; i++) {
      index->list[i].out = 0ul;
      index->list[i].in = 0ul;
      index->list[i].bits= 0;
      memset(index->list[i].window, 0, WINSIZE);
    }
    if (index->list == NULL) {
      free(index);
      return NULL;
    }
    index->size = 8;
    index->have = 0;
  }

  /* if list is full, make it bigger */
  else if (index->have == index->size) {

    oldsize = index->size;
    index->size <<= 1;
    next = (gzpoint_t*) realloc(index->list, sizeof(gzpoint_t) * index->size);

    if (next == NULL) {
      fprintf(stderr, "error during re-allocation of bgz index structure.\n");
      exit(EXIT_FAILURE);
    }

    for(i=oldsize; i < index->size; i++) {
      next[i].out = 0;
      next[i].in = 0;
      next[i].bits= 0;
      //either change point struct or allocate the memory
      //for the dictonary dynamically. This is necessary to
      //avoid allocs for bgzip files and still use the same
      //index structure.
      memset(next[i].window, 0, WINSIZE);
    }
    index->list = next;
  }

  /* fill in entry and increment how many we have */
  next = index->list + index->have;
  next->in = in;
  next->out = out;
  index->have++;

  return index;
}


gzaccess_t* bl_gzBuildIndex(char *filename, off_t span) {
  FILE *fp;
  int ret;
  off_t totin, totout;        /* our own total counters to avoid 4GB limit */
  off_t last;                 /* totout value of last access point */
  struct access *index;       /* access points being generated */
  z_stream strm;
  unsigned char input[CHUNK];
  unsigned char window[WINSIZE];

  fp = fopen(filename, "rb");

  if (fp == NULL) {
    fprintf(stderr, "zran: could not open %s for reading\n", filename);
    exit(-1);
  }

  memset(window, 0, WINSIZE);
  /* initialize inflate */
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = 0;
  strm.next_in = Z_NULL;
  ret = inflateInit2(&strm, 47);      /* automatic zlib or gzip decoding */
  if (ret != Z_OK) {
    fprintf(stderr, "initialization of inflate failed.\n");
    exit(EXIT_FAILURE);
  }

  /* inflate the input, maintain a sliding window, and build an index -- this
     also validates the integrity of the compressed data using the check
     information at the end of the gzip or zlib stream */
  totin = totout = last = 0;
  index = NULL;               /* will be allocated by first addpoint() */
  strm.avail_out = 0;
  do {
    /* get some compressed data from input file */
    strm.avail_in = fread(input, 1, CHUNK, fp);
    if (ferror(fp)) {
      ret = Z_ERRNO;
      goto build_index_error;
    }
    if (strm.avail_in == 0) {
      ret = Z_DATA_ERROR;
      goto build_index_error;
    }
    strm.next_in = input;

    /* process all of that, or until end of stream */
    do {
      /* reset sliding window if necessary */
      if (strm.avail_out == 0) {
        strm.avail_out = WINSIZE;
        strm.next_out = window;
      }

      /* inflate until out of input, output, or at end of block --
         update the total input and output counters */
      totin += strm.avail_in;
      totout += strm.avail_out;
      ret = inflate(&strm, Z_BLOCK);      /* return at end of block */
      totin -= strm.avail_in;
      totout -= strm.avail_out;
      if (ret == Z_NEED_DICT)
        ret = Z_DATA_ERROR;
      if (ret == Z_MEM_ERROR || ret == Z_DATA_ERROR)
        goto build_index_error;
      if (ret == Z_STREAM_END)
        break;

      /* if at end of block, consider adding an index entry (note that if
         data_type indicates an end-of-block, then all of the
         uncompressed data from that block has been delivered, and none
         of the compressed data after that block has been consumed,
         except for up to seven bits) -- the totout == 0 provides an
         entry point after the zlib or gzip header, and assures that the
         index always has at least one access point; we avoid creating an
         access point after the last block by checking bit 6 of data_type
         */
      if ((strm.data_type & 128) && !(strm.data_type & 64) &&
          (totout == 0 || totout - last > span)) {
        index = bl_gzAddpoint(index, strm.data_type & 7, totin,
            totout, strm.avail_out, window);
        if (index == NULL) {
          ret = Z_MEM_ERROR;
          goto build_index_error;
        }
        last = totout;
      }
    } while (strm.avail_in != 0);
  } while (ret != Z_STREAM_END);

  /* clean up and return index (release unused entries in list) */
  (void)inflateEnd(&strm);

  index = realloc(index, sizeof(struct point) * index->have);
  index->size = index->have;

  fclose(fp);

  return index;

  /* return error */
build_index_error:
  (void)inflateEnd(&strm);
  if (index != NULL)
    free_index(index);
  fclose(fp);

  return NULL;
}

gzaccess_t* bl_bgzBuildIndex(char *filename) {

  FILE *fp;
  off_t total=0, end;
  uint64_t blockcount = 0;
  gzaccess_t *idx = NULL;

  fp = fopen(filename, "rb");
  if (fp == NULL) {
    fprintf(stderr, "zran: could not open %s for reading\n", filename);
    exit(-1);
  }

  end = fseek(fp, 0, SEEK_END);
  if(end != 0) {
    perror("bl_bgzBuildIndex: error seeking end of file");
    exit(EXIT_FAILURE);
  }

  end = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  gzip_Header gzipheader = gzip_Header_default();

  if(gzip_readHeader(fp, &gzipheader) != 0){
    fprintf(stderr, "Failed to read gzip header\n");
    exit(EXIT_FAILURE);
  }

  bgzip_Header bgzipheader = bgzip_Header_default();

  if(bgzip_extractBgzHeader(&gzipheader, &bgzipheader) != 0){
    fprintf(stderr, "Failed to read bgzip header.\n");
    exit(EXIT_FAILURE);
  }

  while (1) {
 
    off_t in = gzipheader.offsetInFile;
    off_t lenUncompressed =  bgzip_findLenUncompressedData(fp,
        &gzipheader, &bgzipheader);

    if(lenUncompressed == -1){
      fprintf(stderr, "error readling length of bgzip header.\n");
      exit(EXIT_FAILURE);
    }

    /*empty blocks are allowed at the end of the file, rasing EOF*/

    if(lenUncompressed == 0){
      // EOF block
      break;
    }

    //fprintf(stderr, "adding access point %ld -> %ld\n", in, total);
    idx = bl_bgzAddpoint(idx, in, total);

    total += lenUncompressed;
    blockcount++;
   
    /*the following code fixes a EOF violation produced by bcl2fastq2*/
    if(ftell(fp) == end) {
      NFO("EOF block in file '%s' missing. This is a violation of the bgzf convention.\n", filename);
      if(!feof(fp)) {
        NFO("EOF byte at the end of file missing.\n", NULL);
      }
      NFO("This might happen because you are using Illumina data processed by a new bcl2fastq.", NULL);
      NFO("I am continuing anyways.\n", NULL);
      break;
    }

    assert(gzip_readHeader(fp, &gzipheader) == 0);
    bgzip_extractBgzHeader(&gzipheader, &bgzipheader);
  }

  idx->end = end;

  fclose(fp);
  return idx;
}


int bl_bgzWindow2Buffer(unsigned char *buf, uint32_t pos, uint32_t buflen,
    unsigned char *window, uint32_t srcsz, char where) {

  uint32_t dstsz, cpysz;


  if(pos+1 < buflen) {
    dstsz = buflen - pos;
    cpysz = (srcsz > dstsz) ? dstsz : srcsz;
   /*
    uint32_t hdrsz=0;
    char hdr[]= "\n--------- new window follows ---------\n";
    hdrsz = strlen(hdr);
    char *outstring = NULL;
    outstring = ALLOCMEMORY(NULL, NULL, char, hdrsz+cpysz+1);
    memcpy(outstring, hdr, hdrsz);
    memcpy(&outstring[hdrsz], window, cpysz);
    outstring[hdrsz+cpysz] = 0;
    fprintf(stderr, "\n srcsz:%d %d: %s", cpysz==srcsz, where, outstring);
    FREEMEMORY(space, outstring);
  */
    memcpy(&buf[pos], window, cpysz);

    return cpysz;
  }


  return 0;
}

int bl_bgzFeof(FILE *fp, off_t end) {
  if(ftell(fp) == end) return 1;

  return feof(fp);
}

void bl_bgzInitStream(z_stream *strm) {

  strm->zalloc = Z_NULL;
  strm->zfree = Z_NULL;
  strm->opaque = Z_NULL;
  strm->avail_in = 0;
  strm->next_in = Z_NULL;

  return;
}

void bl_bgzResetOutWindow(unsigned char *window, uint32_t wsz, z_stream *strm) {

  memset(window, 0, wsz);
  strm->avail_out = wsz;
  strm->next_out = window;

  return;
}

int bl_bgzFillStream(FILE *fp, unsigned char *input, z_stream *strm) {

  int n;

  if(strm->avail_in) {
    unsigned char *p = input;
    unsigned const char *q = strm->next_in;
    unsigned n = strm->avail_in;
    do {
      *p++ = *q++;
    } while (--n);
  }

  n = fread(&input[strm->avail_in], 1, CHUNK-strm->avail_in, fp);
  if (ferror(fp)) {
    fprintf(stderr, "error reading bgz file.\n");
    perror("The following error occurred:");
    exit(EXIT_FAILURE);
  }

  strm->next_in = input;
  strm->avail_in += n;

  return n;
}

int bl_gzExtract(FILE *fp, gzaccess_t *index, off_t off, unsigned char *buf, uint32_t len) {
  int ret, skip;
  z_stream strm;
  struct point *here;
  unsigned char input[CHUNK];
  unsigned char discard[WINSIZE];

  /* proceed only if something reasonable to do */
  if (len == 0)
    return 0;

  /* find where in stream to start */
  here = index->list;
  ret = index->have;
  while (--ret && here[1].out <= off)
    here++;

  /* initialize file and inflate state to start there */
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = 0;
  strm.avail_out = 0;
  strm.next_in = Z_NULL;
  ret = inflateInit2(&strm, -15);         /* raw inflate */
  if (ret != Z_OK)
    return ret;
  ret = fseeko(fp, here->in - (here->bits ? 1 : 0), SEEK_SET);
  if (ret == -1)
    goto extract_ret;
  if (here->bits) {
    ret = getc(fp);
    if (ret == -1) {
      ret = ferror(fp) ? Z_ERRNO : Z_DATA_ERROR;
      goto extract_ret;
    }
    (void)inflatePrime(&strm, here->bits, ret >> (8 - here->bits));
  }
  (void)inflateSetDictionary(&strm, here->window, WINSIZE);

  /* skip uncompressed bytes until offset reached, then satisfy request */
  off -= here->out;
  strm.avail_in = 0;
  skip = 1;                               /* while skipping to offset */
  do {
    /* define where to put uncompressed data, and how much */
    if (off == 0 && skip) {          /* at offset now */
      strm.avail_out = len;
      strm.next_out = buf;
      skip = 0;                       /* only do this once */
    }
    if (off > WINSIZE) {             /* skip WINSIZE bytes */
      strm.avail_out = WINSIZE;
      strm.next_out = discard;
      off -= WINSIZE;
    }
    else if (off != 0) {             /* last skip */
      strm.avail_out = (unsigned)off;
      strm.next_out = discard;
      off = 0;
    }

    /* uncompress until avail_out filled, or end of stream */
    do {

      if (strm.avail_in == 0) {
        strm.avail_in = fread(input, 1, CHUNK, fp);
        if (ferror(fp)) {
          ret = Z_ERRNO;
          goto extract_ret;
        }
        if (strm.avail_in == 0) {
          ret = Z_DATA_ERROR;
          goto extract_ret;
        }
        strm.next_in = input;
      }
      ret = inflate(&strm, Z_NO_FLUSH);       /* normal inflate */
      if (ret == Z_NEED_DICT)
        ret = Z_DATA_ERROR;
      if (ret == Z_MEM_ERROR || ret == Z_DATA_ERROR)
        goto extract_ret;
      if (ret == Z_STREAM_END)
        break;
    } while (strm.avail_out != 0);

    /* if reach end of stream, then don't keep trying to get more */
    if (ret == Z_STREAM_END){
      //fprintf(stderr, "Z_STREAM_END break.\n");
      break;
    }

    /* do until offset reached and requested data read, or stream ends */
  } while (skip);

  /* compute number of uncompressed bytes read after offset */
  ret = skip ? 0 : len - strm.avail_out;

  /* clean up and return bytes read or error */
extract_ret:
  (void)inflateEnd(&strm);

  return ret;
}

int bl_bgzExtract(FILE *fp, gzaccess_t *index, off_t off, unsigned char *buf, uint32_t len) {

  int64_t ret; 
  uint32_t ptr =0;
  z_stream strm;
  struct point *here;
  unsigned char input[CHUNK];
  unsigned char wndw[WINSIZE];

  if(len <= 0) return 0;

  

  /* find where in stream to start, in bgz we dont need to take
   * care of any bits stored in predeceeding byte since we jump
   * into a complete, decompressable block     */
  here = index->list;
  ret = index->have;
  while (--ret && here[1].out <= off)
    here++;

  //fprintf(stderr, "requested offset:%ld len:%d\n", off, len);
  //fprintf(stderr, "selected block %ld\n", index->have-ret);

  /* initialize inflate */
  bl_bgzInitStream(&strm);
  ret = inflateInit2(&strm, 47);      /* automatic zlib or gzip decoding */

  if (ret != Z_OK)
    return ret;

  ret = fseeko(fp, here->in, SEEK_SET);
  off -= here->out;

  //fprintf(stderr, "remaining offset %ld, here->out %ld\n", off, here->out);
  //fprintf(stderr, "should be skipping %ld windows of size %d\n", off/WINSIZE, WINSIZE);
  
  /* inflate the input, maintain a sliding wndw, and build an index -- this
     also validates the integrity of the compressed data using the check
     information at the end of the gzip or zlib stream */

  strm.avail_out = 0;
  do {

    /* fill input stream with compressed data from input file */
    bl_bgzFillStream(fp, input, &strm);
    /* reset the sliding wndw */
    bl_bgzResetOutWindow(wndw, WINSIZE, &strm);
    /* process all of that, or until end of stream */

    do {

      /* inflate until out of input, output, or at end of block --
         update the total input and output counters */

      ret = inflate(&strm, Z_NO_FLUSH);      /* return at end of block */

      if (ret == Z_NEED_DICT) {
        fprintf(stderr, "dictionary error\n");
        exit(EXIT_FAILURE);
      }
      if (ret == Z_MEM_ERROR || ret == Z_DATA_ERROR) {
        fprintf(stderr, "data or dictionary error\n");
        exit(EXIT_FAILURE);
      }
      if (ret == Z_STREAM_END) {
        break;
      }

      /* reset sliding wndw if necessary */
      if (strm.avail_out == 0) {

        if(ptr+1 < len) {
          if (off > WINSIZE) {
            off -= WINSIZE;
          } else if(off != 0) {
            if(off < (WINSIZE-strm.avail_out)) {
              //fprintf(stderr, "\nentered 1: off:%ld avail:%u ptr:%d len:%d\n", off, strm.avail_out, ptr, len);
              ptr +=  bl_bgzWindow2Buffer(buf, ptr, len, &wndw[off], (WINSIZE-strm.avail_out)-off, 1);
              //fprintf(stderr, "after leaving pointer is %d\n", ptr);
              off = 0;
            } else {
             //fprintf(stderr, "skipped 1: off:%ld avail:%u\n", off, strm.avail_out);
             off -= (WINSIZE-strm.avail_out);
            }
          } else {
            ptr += bl_bgzWindow2Buffer(buf, ptr, len, wndw, (WINSIZE-strm.avail_out), 2);
          }
        }
        /* reset the sliding wndw */
        bl_bgzResetOutWindow(wndw, WINSIZE, &strm);
      }

    } while (strm.avail_in != 0 && ptr+1 < len);

    //show_window(stderr, wndw, strm.avail_out);
    if(ptr+1 < len) {

      if (off > WINSIZE) {
        off -= WINSIZE;
      } else if(off != 0) {
        if(off < (WINSIZE-strm.avail_out)) {
          ptr +=  bl_bgzWindow2Buffer(buf, ptr, len, &wndw[off], (WINSIZE-strm.avail_out)-off, 3);
          off = 0;
        } else {
          //fprintf(stderr, "skipped 3: off:%ld avail:%u\n", off, strm.avail_out);
          off -= (WINSIZE-strm.avail_out);
        }
      } else {
        //fprintf(stderr, "\nentered off:%ld avail:%u ptr:%d, len:%d\n", off, strm.avail_out, ptr, len);
        ptr += bl_bgzWindow2Buffer(buf, ptr, len, wndw, (WINSIZE-strm.avail_out),4);
        //fprintf(stderr, "after leaving pointer is %d\n", ptr);
      }

      /* reset the sliding wndw */
      bl_bgzResetOutWindow(wndw, WINSIZE, &strm);

      if(ret == Z_STREAM_END && !bl_bgzFeof(fp, index->end)) {
        ret = inflateReset(&strm);
        if(ret != Z_OK) {
          fprintf(stderr, "resetting of the bgz inflate stream failed.\n");
          exit(EXIT_FAILURE);
        }
      }
    }

  } while (!bl_bgzFeof(fp, index->end)  && ptr+1 < len);

//  fprintf(stderr, "at off: %ld buffer filled with %ld bytes (requested:%d), ptr says:%d EOF:%d\n", oldoff, strlen(buf), len, ptr, feof(fp));
  (void)inflateEnd(&strm);

  return ptr;
}


int
bl_gzGetType(char *filename) {
  FILE *fp;

  fp = fopen(filename, "rb");
  if (fp == NULL) {
    fprintf(stderr, "zran: could not open %s for reading\n", filename);
    exit(-1);
  }

  gzip_Header gzipheader = gzip_Header_default();

  if(gzip_readHeader(fp, &gzipheader) != 0){
    fprintf(stderr, "Failed to read gzip header\n");
    exit(1);
  }

  bgzip_Header bgzipHeader = bgzip_Header_default();
  if(bgzip_extractBgzHeader(&gzipheader, &bgzipHeader) == 0){
    NFO("bgzip format detected, compressed size: %i\n", bgzipHeader.lenCompressedData);
    return 1;
  }

  fclose(fp);

  return 0;
}

gzaccess_t*
bl_gzGetIndex(char *filename) {


  int type;
  gzaccess_t *idx;

  type = bl_gzGetType(filename);

  if(type) {
    /* bgzip */
    idx = bl_bgzBuildIndex(filename);
  } else {
    /*gzip/zlib*/
    idx = bl_gzBuildIndex(filename, SPAN);
  }

  //fprintf(stderr, "index is build.\n");
  idx->type = type;
  return idx;
}

uint64_t
bl_zranExtract(FILE *fp, gzaccess_t *idx, off_t off, unsigned char *buf, uint32_t len) {

  uint64_t n=0;

  if(idx->type) {
   n = bl_bgzExtract(fp, idx, off, buf, len);
  } else {
   n = bl_gzExtract(fp, idx, off, buf, len);
  }

  return n;
}

/*------------------------------ bl_initgzfile -------------------------------
 *
 * @brief init gzfile
 * @author Steve Hoffmann
 *
 */

  struct gzidxfile*
bl_initgzidxfile(FILE *fp, struct access *index, off_t offset, int mychunk)
{
  struct gzidxfile *file;

  if(mychunk == 0) mychunk = LARGECHUNK;
  file =(struct gzidxfile *) calloc(1, sizeof(struct gzidxfile));
  file->fp = fp;
  file->index = index;
  file->mychunk = mychunk;
  file->buf =(unsigned char*) calloc(file->mychunk+1, sizeof(unsigned char));
  file->buf[file->mychunk] = 0;
  file->curap = offset;
  file->len = 0;
  file->pos = file->buf;

  return file;
}
/*-------------------------- bl_destructgzidxfile ----------------------------
 *
 * @brief destruct idx file
 * @author Steve Hoffmann
 *
 */

void
bl_destructgzidxfile(struct gzidxfile *file) {
  free(file->buf);
}
/*-------------------------------- bl_getgzc ---------------------------------
 *
 * @brief get char from gz stream
 * @author Steve Hoffmann
 *
 */

int
bl_getgzidxc (struct gzidxfile *f)
{


  if(f->len == 0 || f->pos - f->buf >= f->len) {

    memset(f->buf, 0, f->mychunk+1);
    f->curap += f->pos - f->buf;
    //fprintf(stderr, "curap: %ld, flen:%d\n", f->curap, f->len);
    f->len = bl_zranExtract(f->fp, f->index, f->curap, (unsigned char*) f->buf, f->mychunk);
   // fprintf(stderr, "|");

    if (f->len == 0) {

      //fprintf(stderr, "returning EOF\n");
      return EOF;
    }
    if (f->len  < 0) {
      fprintf (stderr, "zran: extraction failed: %s error \n",
        f->len == Z_MEM_ERROR ? "out of memory" : "input corrupted");
      exit(EXIT_FAILURE);
      return EOF;
    }  else {
 //       fprintf(stderr, "zran: extracted %d bytes at %llu\n", f->len, f->curap);
    }

    f->pos = f->buf;
  }

  return *f->pos++;
}

off_t
bl_ftellgzidx(struct gzidxfile *f) {
  return f->curap + (f->pos - f->buf);
}


