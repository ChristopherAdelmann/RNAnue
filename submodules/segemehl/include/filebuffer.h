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
//https://codereview.stackexchange.com/questions/75339/read-file-line-by-line-with-circular-buffer

#ifndef FILEBUFFER_H
#define FILEBUFFER_H

/*
 *
 *	filebuffer.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 22.11.2016 16:46:12 CET  
 *
 */


#define _DEFAULT_SOURCE 1
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <inttypes.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <pthread.h>


typedef struct {
    int beg; // Index of first element added to buffer.
    int end; // Index of most recent element added to buffer.
    size_t size; // Number of elements in circular buffer.
    char *buffer; 
    FILE *dev;
    pthread_mutex_t *mtx;
    int line_beg;
    char feof;
    char mode;
} circbuffer_t;


char  bl_circBufferRead(circbuffer_t *cb);
void  bl_circBufferAddByte(circbuffer_t *cb, char c);
void bl_circBufferSaveAddByte(circbuffer_t *cb, char c);
int bl_circBufferAddSave(circbuffer_t *cb, char *data, size_t len);
void bl_circBufferInit(circbuffer_t *cb, size_t size, FILE *dev, pthread_mutex_t *mtx);
int   bl_circBufferIsEmpty(circbuffer_t *cb);
void  bl_circBufferDestruct(circbuffer_t *cb);
int  bl_circBufferIsFull(circbuffer_t *cb);
circbuffer_t* bl_circBufferInitArray(int n, size_t size, FILE *dev, pthread_mutex_t *mtx);
void  bl_circBufferDestructArray(circbuffer_t *bufarr, int n);
void bl_circBufferEmptyArray(circbuffer_t *bufarr, int n);
char* bl_circBufferReadLine(circbuffer_t *cb, uint32_t *len);
void bl_circBufferEmpty(circbuffer_t *cb) ;

#endif
