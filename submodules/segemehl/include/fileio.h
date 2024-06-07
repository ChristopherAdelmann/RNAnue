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


#ifndef FILEIO_H
#define FILEIO_H

/*
 * fileio.h
 * declarations for file io
 *
 * @author Steve Hoffmann
 * @date Sat 25 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: fileio.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/fileio.h $
 */

#ifndef ALLOCMEMORY
	#include "memory.h"
#endif
#include <math.h>
#include <sys/types.h>
#include <pthread.h>
#include "stringutils.h"
char * bl_getTempFile(char *path, char *tmp);
int bl_UnixSort(void *space, char *filename, const char *fieldstring, const char delim);
char* readfile(void *, char *, size_t*);
stringset_t **readcsv(void *, char *, char*, Uint *);
void writeY(char *, double  *, Uint, Uint, Uint);
void writeXYUint(char *filename, Uint *X, Uint *Y, Uint len);
int bl_fgets(void *space, FILE *fp, char **str);
char * bl_basename (const char *name);
char* bl_replacealphanum(char *s, Uint len);
char*  bl_replacenonalphanum(char *s, Uint len);
int bl_fileprefixlen(char *filename);
int bl_UnixSortMerge(void *space, char **filenames, Uint nooffiles, const char *fieldstring, const char delim, char *outfile);
void bl_writeFileHeader(char *filename, char *header) ;
void bl_freplace(char *filename, char oldchar, char newchar, char stop);
void bl_freplacearr(char *filename, char* oldchars, char *newchars, Uint len, char stop);
void bl_freplacestr(char *filename, char* str, Uint len, char stop);
void writeXYZ(char *filename, double *X, double *Y, double *Z, Uint len);
int bl_rm(void *space, char *filename);
void writeYUint(char *filename, Uint *Y, Uint len, Uint xoff, Uint yoff);
void writeYUintNorm(char *filename, Uint *Y, Uint len, Uint yoff);
double* readX(void *space, char *filename, Uint *nvals);
char *bl_changefilesuffix(char *base, char *suffix);
char* bl_dirname(char *origpath);
const char *bl_fsuffix(char *fn);
#endif
