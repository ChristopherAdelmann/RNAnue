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
 *  filebuffer.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 22.11.2016 16:45:13 CET
 *  
 */

#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE 1
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <pthread.h>
#include "filebuffer.h"
#include "stringutils.h"
#include "basic-types.h"
#include "memory.h"
#include "mathematics.h"




char 
bl_circBufferRead(circbuffer_t *cb) {
  char mychar;

  if(cb->beg == cb->end) return 0;

  mychar = cb->buffer[cb->beg];
  cb->beg = (cb->beg + 1) % cb->size;

  return mychar;
}


void 
bl_circBufferAddByte(circbuffer_t *cb, char c) {
  cb->buffer[cb->end] = c;
  cb->end = ( ((cb->end) + 1) % cb->size);
  if (cb->end == cb->beg)
  {
    cb->beg = (cb->beg + 1) % cb->size; 
  }
}

void 
bl_circBufferSaveAddByte(circbuffer_t *cb, char c) {
  char mychar;
  if(bl_circBufferIsFull(cb)) {

    if(cb->mtx) pthread_mutex_lock(cb->mtx);
    
      while(!bl_circBufferIsEmpty(cb)) {
      mychar = bl_circBufferRead(cb);
      fprintf(cb->dev, "%c", mychar);
    }

    if(cb->mtx) pthread_mutex_unlock(cb->mtx); 
  }

  bl_circBufferAddByte(cb, c);
}

void
bl_circBufferEmpty(circbuffer_t *cb) {
   size_t right, left=0;

  if(cb->beg <= cb->end) { 
    right = cb->end - cb->beg;
  } else {
    //cb->beg > cb->end
    right = cb->size - cb->beg;
    left = cb->end;
  }

   if(cb->mtx) pthread_mutex_lock(cb->mtx);

    if(cb->beg <= cb->end) {
      fwrite(&cb->buffer[cb->beg], right*sizeof(char), 1, cb->dev);
    } else {
      fwrite(&cb->buffer[cb->beg], right*sizeof(char), 1, cb->dev);
      fwrite(cb->buffer, left*sizeof(char), 1, cb->dev);
    }
    
    cb->beg = 0;
    cb->end = 0;

   if(cb->mtx) pthread_mutex_unlock(cb->mtx);
}

int
bl_circBufferAddSave(circbuffer_t *cb, char *data, size_t len) {

  size_t fspace=0, right = 0;


  //if data does not fit into buffer - empty the buffer and
  //write straight into file
  if(len > cb->size) {

    if(!bl_circBufferIsEmpty(cb)) {
      bl_circBufferEmpty(cb);
    }
    
    if(cb->mtx) pthread_mutex_lock(cb->mtx);
    fwrite(data, len*sizeof(char), 1, cb->dev);
    if(cb->mtx) pthread_mutex_unlock(cb->mtx);

    return 0;
  } 

  if(cb->beg <= cb->end) { 
    fspace = cb->size - (cb->end - cb->beg);
    right = cb->size - cb->end;
  } else {
    fspace = cb->beg - cb->end;
  }

  //if remaining space does not suffice, dump buffer
  if (len > fspace) { 
    bl_circBufferEmpty(cb);
    right = cb->size;
  } 

  if(cb->beg <= cb->end) {
    if(len > right) {
      memmove(&cb->buffer[cb->end], data, right);
      memmove(cb->buffer, &data[right], len-right);
      cb->end = len-right;
    } else {
      memmove(&cb->buffer[cb->end], data, len);
      cb->end += len;
    }
  } else {
    memmove(&cb->buffer[cb->end], data, len);
    cb->end += len;
  }

  return 0;
}

void 
bl_circBufferInit(circbuffer_t *cb, size_t size, FILE *dev, pthread_mutex_t *mtx) {
    cb->buffer = calloc(size + 1, sizeof(char));
    cb->size  = size + 1;
    cb->beg = 0;
    cb->end = 0;
    cb->dev = dev;
    cb->mtx = mtx;
    cb->feof = 0;
    cb->mode = 0;
}

size_t 
bl_circBufferSize(circbuffer_t *cb) {
  return cb->size-1;
}


//distance to the end of the buffer, the end of data
//or the pointer of the begin - whatever comes first

size_t
bl_circBufferScanLimitLength(circbuffer_t *cb, int ptr) {
  
  size_t m;
  
  if (ptr >= cb->beg) {
    if(cb->end <= ptr) {
      m = cb->size - ptr;
    } else {
      m = cb->end - ptr;
    }
  } else {
    m = MIN(cb->beg,cb->end) - ptr;
  }

  return m;
}

//circular distance from pointer to beginning
size_t
bl_circBufferDist(circbuffer_t *cb, int ptr) {
  if(ptr >= cb->beg) {
    return ptr - cb->beg; 
  } else {
    return cb->size - cb->beg + ptr;
  }
  return 0;
}

void
printbuffer(char *buffer, int len) {
  int i;

  fprintf(stderr, "'");
  for(i=0; i < len; i++) {
    fprintf(stderr, "%c",buffer[i]);
  }
  fprintf(stderr, "'");
}


void
setcircbuffer(char *buffer, int len) {
  int i;

  buffer[0] = '^';
  for(i=1; i < len; i++) {
   buffer[i] = '_';
  }

}

char*
bl_circBufferReadLine(circbuffer_t *cb, uint32_t *len) {
  size_t n,m,newsz;
  int32_t ptr, rr, ell;
  char *crlf = NULL, *str = NULL, noeol=0;
  
  *len = 0;
 
  if(feof(cb->dev) && cb->beg >= cb->end+1) {
    cb->feof=1;
    return NULL;
  }

  ptr = cb->beg;

  while(crlf == NULL) {
 
    m = cb->end - ptr + 1;   

    //scan remainder of circ buffer
    crlf = (char*)memchr(&cb->buffer[ptr], '\n', m); 

    if(crlf == NULL) {
      ptr += m-1;
      
      if(!feof(cb->dev)) {

        if((ptr+1) % cb->size == cb->beg) {  
          
          newsz = cb->size * 2;
          cb->buffer = ALLOCMEMORY(NULL, cb->buffer, char, newsz);
          
          if(ptr != cb->size-1) {    
            //shift the block to the right 
            rr = cb->size - cb->beg;
            memmove(&cb->buffer[newsz-rr], &cb->buffer[cb->beg], rr);  
            //set a new beginning
            cb->beg = newsz-rr;
          }

          //set the increased size
          cb->size = newsz;
        } else {
          //the pointer has not reached the beginning but
          //is at the end of the buffer, set to 0 for circle

          ptr = -1;
        }

        //determine the limit (order: new size needs to be set)        
        if (ptr == -1) { 
          m = cb->beg;
        } else if (ptr < cb->beg) {
          m = cb->beg - (ptr+1);
        } else {
          m = cb->size - (ptr+1);
        }
        ptr++; 

        //fill the new buffer space
        n = fread(&cb->buffer[ptr], sizeof(char), m, cb->dev); 
        cb->end = ptr+n-1;       

      } else {
        //last scan round unsuccessful
        noeol = 1;
        break;
      }
    }
  }

  if(crlf) {
    //get postion of crlf pointer in buffer
    ell = crlf-cb->buffer;
  } else {
    //can only occur at the end of file
    assert(noeol);
    assert(feof(cb->dev));
    fprintf(stderr, "no end of line\n");
    ell = cb->end+1;
    cb->feof = 1;
  }

  //copy the string
  m = bl_circBufferDist(cb, ell); 
  str = (char*) calloc(m+1, sizeof(char));

  if(ell > cb->beg) {
    memcpy(str, &cb->buffer[cb->beg], m); 
  } else if (ell < cb->beg) {
    memcpy(str, &cb->buffer[cb->beg], cb->size-cb->beg);
    memcpy(&str[cb->size - cb->beg], cb->buffer, ell);
  }

  *len = m;

  //fill up the buffer for next iteration
  if(!feof(cb->dev) && ell < cb->beg) { 
    m = cb->size - cb->beg;
    n = fread(&cb->buffer[cb->beg], sizeof(char), m, cb->dev); 
    cb->end = cb->beg +n-1;
  }

  //the new line starts after crlf at position ell
  cb->beg = ell+1;


  return str;
}

    
  int 
bl_circBufferIsEmpty(circbuffer_t *cb) {
    return cb->end == cb->beg;
}

void 
bl_circBufferDestruct(circbuffer_t *cb) {
    free(cb->buffer); /* OK if null */ 
}



int 
bl_circBufferIsFull(circbuffer_t *cb) {
    return (cb->end + 1) % cb->size == cb->beg; 
}

circbuffer_t* 
bl_circBufferInitArray(int n, size_t size, FILE *dev, pthread_mutex_t *mtx) {
  unsigned int i;
  circbuffer_t *bufarr;

  bufarr = ALLOCMEMORY(NULL, NULL, circbuffer_t, n);
  for(i=0; i < n; i++) {
    bl_circBufferInit(&bufarr[i], size, dev, mtx);
  }

  return bufarr;
}

void 
bl_circBufferDestructArray(circbuffer_t *bufarr, int n) {
  unsigned int i;

  for(i=0; i < n; i++) {
    bl_circBufferDestruct(&bufarr[i]);
  }
}

void 
bl_circBufferEmptyArray(circbuffer_t *bufarr, int n) {
  unsigned int i;

  for(i=0; i < n; i++) {
    bl_circBufferEmpty(&bufarr[i]);
  }
}

