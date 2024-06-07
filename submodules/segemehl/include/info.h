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


 #ifndef INFO_H
 #define INFO_H

/*
 *
 *	info.h
 *  nfo messages
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 08/26/2007 07:17:44 PM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: info.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/info.h $
 */

 #include <stdarg.h>
 #include <stdio.h>
 #include <string.h>

#ifndef NFOLEVEL
#define NFOLEVEL 0
#endif

#ifndef NFODEVICE
#define NFODEVICE stderr
#endif

#define NFOL(L, X, ... ) debuglevel (__FILE__, __LINE__, L, X, __VA_ARGS__)
#define NFO(X, ...) infomsg(__FILE__, __LINE__, X, __VA_ARGS__)
#define INFO(X, ...) infomsg(__FILE__, __LINE__, X, __VA_ARGS__)
#define MSG(X) infomsg(__FILE__, __LINE__, X)

#define handle_error_en(en, msg) \
               do { errno = en; perror(msg); exit(EXIT_FAILURE); } while (0)

#define handle_error(msg) \
               do { perror(msg); exit(EXIT_FAILURE); } while (0)

extern unsigned char mute;

int infomsg(char *, int, const char *fmt, ...);
int infolevel(char *, int, int, const char *fmt, ...);

 #endif
