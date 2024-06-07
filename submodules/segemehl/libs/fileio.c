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
 * fileio.c
 * functions to manipulate and read files
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: fileio.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/fileio.c $
 *  
 */


#ifndef _DEFAULT_SOURCE
 #define _DEFAULT_SOURCE 1
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <pthread.h>
#include "stringutils.h"
#include "basic-types.h"
#include "fileio.h"
#include "info.h"

#ifndef DIR_SEPARATOR
#define DIR_SEPARATOR '/'
#endif

#if defined (_WIN32) || defined (__MSDOS__) || defined (__DJGPP__) || \
  defined (__OS2__)
#define HAVE_DOS_BASED_FILE_SYSTEM
#ifndef DIR_SEPARATOR_2 
#define DIR_SEPARATOR_2 '\\'
#endif
#endif

/* Define IS_DIR_SEPARATOR.  */
#ifndef DIR_SEPARATOR_2
# define IS_DIR_SEPARATOR(ch) ((ch) == DIR_SEPARATOR)
#else /* DIR_SEPARATOR_2 */
# define IS_DIR_SEPARATOR(ch) \
  (((ch) == DIR_SEPARATOR) || ((ch) == DIR_SEPARATOR_2))
#endif /* DIR_SEPARATOR_2 */



void
bl_writeFileHeader(char *filename, char *header) {
  char *tmpfilename, *buffer;
  FILE *out, *in;
  size_t buffersize= 1024, len;

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  tmpfilename = bl_getTempFile(".","headerwrite");

  out = fopen(tmpfilename, "w");

  if(!out) {
    fprintf(stderr, "Couldnt open file %s for writing. Exit forced.", tmpfilename);
    exit(-1);
  }

  fprintf(out, "%s\n", header);
  fclose(out);

  out = fopen(tmpfilename, "a");
  in = fopen(filename, "r");

  if(!in) {
    fprintf(stderr, "Couldnt open file %s for reading. Exit forced.", filename);
    exit(-1);
  }

  while((len = fread(buffer, 1, buffersize, in)) > 0) {
    fwrite(buffer, 1, len, out);
  }

  fclose(out);
  fclose(in);
  int ret = rename(tmpfilename, filename);

  if(ret == 0) {
    NFO("renamed file '%s' successfully to '%s'\n.", tmpfilename, filename);
  } else {
    char* errstring = strerror(errno);
    NFO("renaming file '%s' to '%s' unsuccessful: %s\n.", tmpfilename, filename, errstring);
  }

  return;
}


char* 
bl_replacenonalphanum(char *s, Uint len) {
  Uint i, u=0, lastalphanum=0;
  int ch;
  char *new = calloc(len+1, sizeof(char));

  for(i=0; i < len; i++) {
    ch = s[i];
    if(((Uint)((ch | 0x20) - 'a') < 26u)|| ((Uint)(ch-'0') < 10u)) {
      lastalphanum = u;
      new[u++] = ch;
    } else {
      if (u > 0)
        new[u++] = '_';
    }
  }

  new[lastalphanum+1] ='\0';
  return new;
}



/*------------------------------ bl_getTempFile ------------------------------
 *    
 * @brief get a temporary file
 * @author Steve Hoffmann 
 *   
 */


  char *
bl_getTempFile(char *path, char *tmp)
{

  int res;
  char *fname=NULL;

  fname = ALLOCMEMORY(NULL, NULL, char, strlen(path) + strlen(tmp)+11);

  if(strlen(tmp) > 0)
    sprintf(fname, "%s/%sXXXXXX", path, tmp);
  else
    sprintf(fname,"%s/XXXXXX", path);

  if ((res = mkstemp(fname)) == -1) {
    fprintf(stderr, "Error in creating temporary file '%s'. Exit forced.\n",
        fname);
    exit(-1);
  }
  if (close(res) == -1){
    fprintf(stderr, "Error in closing temporary file '%s'. Exit forced.\n",
        fname);
    exit(-1);   
  }
  return fname;
}

int
bl_UnixSortMerge(void *space, char **filenames, Uint nooffiles, 
    const char *fieldstring, const char delim, char *outfile) {
  int ret;
  Uint i, filenamestringpos;
  char *prg = "LC_COLLATE=C sort";
  char *cmd;
  char *filenamestring = NULL;

  filenamestringpos = 0;
  for(i = 0; i < nooffiles; i++) {

    filenamestring = 
      ALLOCMEMORY(space, filenamestring, char, filenamestringpos+strlen(filenames[i])+2);

    memmove(&filenamestring[filenamestringpos], filenames[i], strlen(filenames[i]));
    filenamestringpos += strlen(filenames[i]);
    filenamestring[filenamestringpos] = ' ';
    filenamestringpos++;
    filenamestring[filenamestringpos] = 0;
  }


  cmd = ALLOCMEMORY(space, NULL, char, strlen(prg) + strlen(fieldstring) 
      + strlen(filenamestring) + strlen(outfile) + 15);
  sprintf(cmd, "%s -m -t '%c' %s %s > %s", prg, delim, fieldstring, filenamestring, outfile);
  ret = system(cmd);

  return ret;
}

int
bl_rm(void *space, char *filename) {
  int ret=0;
  char *prg = "rm";
  char *cmd;

  cmd = calloc(strlen(prg) + strlen(filename) + 10, 1);

  sprintf(cmd, "%s -f %s", prg, filename);
  system(cmd);

  free(cmd);
  return ret;
}

char* 
bl_memrchr(const void* s, char c, size_t n) {
  char *beg = (char *)s;
  char *p = NULL;

  if (s && n) { 

    for(p = beg+n-1;  p > beg; --p) {
      if(*p == c) break;
    }

    if(*p == c) return p; 
  }

  return NULL;
}

char* bl_dirname(char *origpath) {

  char *s=NULL, *p=NULL;
  static const char dot[]=".";
  char *path;

  path = bl_strdup(origpath);

  if(path == NULL) {
    s = NULL;
  } else {
    //set last slash
    s = strrchr(path, DIR_SEPARATOR);
  }

  //strrchr has found a value, at the end
  if(s && s != path && s[1]=='\0') { 
    //scan back to remove multple slashes, e.g. '//'
    for(p=s; p != path; --p) {
      if(!IS_DIR_SEPARATOR(p[-1])) break;
    }

    if(p != path) {
      s= bl_memrchr(path, DIR_SEPARATOR, p - path);
    }
  }

  if(s != NULL) {

    for(p=s; p != path; --p) {
      if(!IS_DIR_SEPARATOR(p[-1])) break;
    }

    if(p == path) {

      if(s == path + 1) {
        s++;
      } else {
        s = path + 1;
      }

    } else {
      s = p;
    }

    *s = '\0';
  } else {

    path = bl_strdup((char *) dot);

  }

  return path;
}

const char *bl_fsuffix(char *fn) {
  const char *start = strrchr(fn,'.');
  if(start == NULL || start == fn) 
    return "";

  return start;
}

int
bl_fileprefixlen(char *filename) {
  Uint i, suf;

  suf = strlen(filename);
  for(i=1; i < strlen(filename); i++) {
    if(filename[i] == '.') {
      suf = i;
    }
  }
  return suf;
}

  char *
bl_basename (const char *name)
{
  const char *base;

#if defined (HAVE_DOS_BASED_FILE_SYSTEM)
  /* Skip over the disk name in MSDOS pathnames. */
  if (ISALPHA (name[0]) && name[1] == ':') 
    name += 2;
#endif

  for (base = name; *name; name++)
  {
    if (IS_DIR_SEPARATOR (*name))
    {
      base = name + 1;
    }
  }
  return (char *) base;
}

int
bl_UnixSort(void *space, char *fn, const char *fieldstring, const char delim) {
  int ret=0;
  char *prg = "LC_COLLATE=C sort";
  char *mv = "mv";
  char *cmd=NULL, *cmd2=NULL;
  char *dn;
  char *tmpfn;

  dn = bl_dirname(fn); 
  tmpfn = bl_getTempFile(dn,"sort");

  bl_bsprintf(&cmd, "%s -o %s -t '%c' %s %s", prg, tmpfn, delim, fieldstring, fn);
  NFO("sorting to '%s'\n", tmpfn);
  NFO("%s.\n", cmd);

  ret = system(cmd);

  if(ret == -1) {
    char* errstring = strerror(errno);
    NFO("sorting to '%s' with '%s' failed:\n'%s'\n", tmpfn, fieldstring, errstring);
    FREEMEMORY(NULL, cmd);
    FREEMEMORY(NULL, dn);
    return ret;
  }

  bl_rm(space, fn);
  ret = rename(tmpfn, fn);

  if(ret == 0) {
    NFO("renamed '%s'\n", tmpfn);
  } else {
    //char* errstring = strerror(errno);
    //NFO("renaming of '%s' failed:\n'%s'\n", tmpfn, errstring);
    NFO("moving '%s' file instead.\n", tmpfn);
    bl_bsprintf(&cmd2, "%s %s %s", mv, tmpfn, fn);
    ret = system(cmd2);
    if(ret == -1) {
      char* errstring = strerror(errno);
      NFO("renaming of '%s' failed: '%s'.\n", tmpfn, errstring);
    } else {
      NFO("renaming of '%s' successful.\n", tmpfn);
    }

    FREEMEMORY(NULL, cmd2); 
  }

  FREEMEMORY(NULL, cmd);
  FREEMEMORY(NULL, dn);
  FREEMEMORY(NULL, tmpfn);

  return ret;
}


char*
bl_changefilesuffix(char *filename, char *suffix) {
  Uint prefixlen;
  Uint suffixlen;
  char *new;

  prefixlen = bl_fileprefixlen(filename);
  suffixlen = strlen(suffix);
  new = ALLOCMEMORY(NULL, NULL, char, prefixlen+1+suffixlen+1);
  memmove(new, filename, prefixlen);
  new[prefixlen] = '.';
  memmove(&new[prefixlen+1], suffix, suffixlen);
  new[prefixlen+suffixlen+1]=0;

  return new;
}


void
bl_fnreplace(char *filename, char oldchar, char newchar, Uint nreplace) {
  int ch, n=0;
  FILE *fp;

  fp = fopen(filename, "rb+");
  if(!fp) {
    fprintf(stderr,"Couldnt open file '%s'. Exit forced!\n", filename);
    exit(-1);
  }

  while((ch = fgetc(fp)) != EOF) {
    if(ch == oldchar) {
      fseek(fp, -1, SEEK_CUR);
      fputc(newchar, fp);
      n++;
    }
    if(nreplace == n) break;
  }


  fclose(fp);
  return;
}

void
bl_freplacearr(char *filename, char* oldchars, char* newchars, Uint len, char stop) {
  int ch, i;
  char oldchar;
  FILE *fp;

  fp = fopen(filename, "rb+");
  if(!fp) {
    fprintf(stderr,"Couldnt open file '%s'. Exit forced!\n", filename);
    exit(-1);
  }

  while((ch = fgetc(fp)) != EOF) {
    for(i=0; i < len; i++) {
      oldchar = oldchars[i];
      if(ch == oldchar) {
        fseek(fp, -1, SEEK_CUR);
        fputc(newchars[i], fp);
        break;
      }
    }

    if(ch == stop) {
      fseek(fp, -1, SEEK_CUR);
      fputc(' ', fp);
      break;
    }
  }

  fclose(fp);
  return;
}

void
bl_freplace(char *filename, char oldchar, char newchar, char stop) {
  int ch;
  FILE *fp;

  fp = fopen(filename, "rb+");
  if(!fp) {
    fprintf(stderr,"Couldnt open file '%s'. Exit forced!\n", filename);
    exit(-1);
  }

  while((ch = fgetc(fp)) != EOF) {
    if(ch == oldchar) {
      fseek(fp, -1, SEEK_CUR);
      fputc(newchar, fp);
    }
    if(ch == stop) {
      fseek(fp, -1, SEEK_CUR);
      fputc('\n', fp);
      break;
    }
  }

  fclose(fp);
  return;
}

void
bl_freplacestr(char *filename, char *str, Uint len, char stop){
  int i = 0;
  char ch;
  FILE *fp;

  fp = fopen(filename, "rb+");  
  if (!fp) {
    fprintf(stderr, "Couldn't open file '%s'. Exit forced.\n", filename);
    assert(0);
    //exit(EXIT_FAILURE);
  }

  while((ch = fgetc(fp)) != EOF){
    if (ch == stop){
      break;
    }
    fseek(fp, -1, SEEK_CUR);
    fputc(str[i%len], fp);
    i++;
  }

  fclose(fp);
  return;
}

int
bl_fgets(void *space, FILE *fp, char **str) {
  char ch, *buffer;
  size_t buffersize = MAXBUFFERSIZE;
  size_t len = 0;

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  while((ch=getc(fp)) != EOF && ch != '\n') {
    if(len == buffersize - 1) {
      buffersize += MAXBUFFERSIZE + 1;
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
    }
    buffer[len++] = (char) ch;
  }

  if(ch == EOF) return EOF;

  buffer[len] = '\0';
  *str = buffer;

  return len;
}


char* 
readfile(void* space, char* filename, size_t* strlen) {

  char ch;
  char *buffer;
  FILE *fp;
  size_t buffersize = MAXBUFFERSIZE;
  size_t len=0;

  fp = fopen(filename, "r");
  if (fp == NULL){
    fprintf(stderr, "Opening of file %s failed. Exit forced.\n", filename);
    exit(EXIT_FAILURE);
  }

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  while((ch=getc(fp)) != EOF) {
    if(len == buffersize-1) {
      buffersize += MAXBUFFERSIZE+1;
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
    }
    len++;
    buffer[len-1]=(char)ch;	
  }
  buffer[len]='\0';
  fclose(fp);

  *strlen = len;
  return buffer;
}



stringset_t **
readcsv(void *space, 
    char* filename, 
    char *delim, 
    Uint *linecount) {

  size_t i, contentlen;
  char *content;
  stringset_t *lines, **csv;

  content = readfile(space, filename, &contentlen);
#ifndef _CRLF_
  lines = tokensToStringset(space, "\n", content, contentlen);
#else
  lines = tokensToStringset(space, "\r\n", content, contentlen);
#endif
  FREEMEMORY(space, content);
  *linecount=lines->noofstrings;
  csv=ALLOCMEMORY(space, NULL, stringset_t *, lines->noofstrings);

  for(i=0; i < lines->noofstrings; i++) {
    csv[i] = tokensToStringset(space, delim, lines->strings[i].str, lines->strings[i].len);
  }

  destructStringset(space, lines);	
  return csv;
}

double*
readX(void *space, char *filename, Uint *nvals) {

  double *X, r;
  Uint n, i, j = 0;
  stringset_t **csv;

  csv = readcsv(space, filename, "\t ", &n);
  X = ALLOCMEMORY(space, NULL, double, n);

  for(i=0; i < n; i++) {
    if(csv[i]->noofstrings){
      //fprintf(stderr, "%s\n", csv[i]->strings[0].str);
      r= atof(csv[i]->strings[0].str);
      if(!isinf(r)) {
        X[j] = r;
        j++;    
        //  fprintf(stderr, "%d\t%f\n",j, r);
      }
    }
    destructStringset(space, csv[i]);
  }


  FREEMEMORY(space, csv);
  X = ALLOCMEMORY(space, X, double, j);
  *nvals = j;
  return X;
}

void
writeY(char *filename, double *Y, Uint len, Uint xoff, Uint yoff) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=yoff; i < len; i++) {
    fprintf(file,"%d\t%f\n", i+xoff, Y[i]);
  }

  fclose(file);
  return;
}

void
writeYUint(char *filename, Uint *Y, Uint len, Uint xoff, Uint yoff) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=yoff; i < len; i++) {
    fprintf(file,"%d\t%d\n", i+xoff, Y[i]);
  }

  fclose(file);
  return;
}

void
writeYUintNorm(char *filename, Uint *Y, Uint len, Uint off) {
  FILE *file;
  Uint i, norm=0;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }
  for(i=0; i < len; i++) {
    norm += Y[i];
  }

  for(i=off; i < len; i++) {
    fprintf(file,"%d\t%f\n", i, (double)Y[i]/norm);
  }

  fclose(file);
  return;
}


void 
writeXYUint(char *filename, Uint *X, Uint *Y, Uint len) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=0; i < len; i++) {
    fprintf(file,"%d\t%d\t%d\n", i, X[i], Y[i]);
  }

  fclose(file);
}

void
writeXYZ(char *filename, double *X, double *Y, double *Z, Uint len) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=0; i < len; i++) {
    fprintf(file,"%f\t%f\t%f\n", X[i], Y[i], Z[i]);
  }

  fclose(file);
}
