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
 *  debug.c
 *  debug messages
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 08/26/2007 06:49:02 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: debug.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/debug.c $
 *  
 */
 
 #include <stdarg.h>
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include "debug.h"

 FILE *dbgdevice = NULL;

 int
 debugmsg( const char *file, 
           const int line, 
           const char *fmt, ...) {

   if (dbgdevice == NULL) {
     dbgdevice = DBGDEVICE;
   }
   
   int ret;
   va_list ap;
   va_start(ap, fmt);
#ifdef DBGNFO   
   fprintf(dbgdevice, "[%s] file: %s, line: %d: ", "SEGEMEHL", file, line);
#endif
   ret = vfprintf(dbgdevice, fmt, ap);
   va_end(ap);

   return ret; 
 }


void 
setdebugdevice(char *filename) {
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    DBG("Couldn't open file '%s'. Exit forced.\n", filename);
    exit(-1);
  }

  dbgdevice = fp;
}

int
debuglevel( const char *file,
    const int line,
    int level,
    const char *fmt, ... ) {

  int ret=0;
  va_list ap;

   if (dbgdevice == NULL) {
     dbgdevice = DBGDEVICE;
   }
  
   if (DBGLEVEL >= level) {

    va_start(ap, fmt);
#ifdef DBGNFO
    fprintf(dbgdevice, "[%s] file: %s, line: %d: ", "segemehl", file, line);
#endif
    ret = vfprintf(dbgdevice, fmt, ap);
    va_end(ap);
  }

  return ret;
}




