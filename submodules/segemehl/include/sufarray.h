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


#ifndef SUF_ARRAY_H
#define SUF_ARRAY_H

/*
 *
 *	sufarray.h
 *  declarations for enhanced suffix arrays
 *  for char alphabets
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 12/10/06 22:01:33 CET  
 *
 *  SVN
 *  Revision of last commit: $Rev: 72 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-10-28 18:14:42 +0100 (Tue, 28 Oct 2008) $
 *
 *  Id: $Id: sufarray.h 72 2008-10-28 17:14:42Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sufarray/sufarray.h $
 */

#include <sys/types.h>
#include "basic-types.h"
#include "falphabet.h"
#include "charsequence.h"
#include "container.h"
#include "multicharseq.h"


#define LCP_TAB_STORED      ((unsigned char) (1 << 0))
#define CHLD_TAB_STORED     ((unsigned char) (1 << 1))
#define SUFLINK_TAB_STORED  ((unsigned char) (1 << 2))
#define SUFLINK_COMPRESSED  ((unsigned char) (1 << 3))
#define MD5_STORED          ((unsigned char) (1 << 4))
#define LINT_SUFLINKS       ((unsigned char) (1 << 5)) 

typedef struct {
    Uint *suf;
    Uint pos;
} suffix_t;

typedef struct {
 /* int up;
    int down;
    int nextlIndex;
 */

  Uint val;
} childtab;


typedef struct {
  MultiCharSeq  *seq;

  Uint		numofsuffixes;
  Uint 		*suftab;
  Uint		*inv_suftab;
  Uint      *suflink;
  Uint      *suflink_l;
  Uint      *suflink_r;

  Uint 		*bwttab; 	     /* burrows-wheeler array*/
  Uint		*lcptab;         /* alternative: Abouelhoda et al.*/

  unsigned char *lcpctab;    /* nB to store lcp values < 255*/
  PairUint	    *llvtab;     /* array of 8B to store lcp val >=255*/
  Uint           llvcnt;
  Uint           maxlcp;

  signed char   *id;
  PairLSint     *idvtab;
  Uint          idvcnt;

  childtab      *chldtab;    /* a child table*/
  Uint		    *bcktab;     /* the bucket container*/

  unsigned char *mdfive;
  unsigned char llint;
#ifdef SUFLINK_MMAP 
  int pagediff_id;
  int pagediff_sl;
#endif
  int fd;
  off_t off_sl;
  off_t off_id;

} Suffixarray;


 
Uint getMultiCharSeqIndex (MultiCharSeq *mseq, char *ptr);
Suffixarray* readSuffixarray(void *, char *, CharSequence **, Uint, unsigned char silent);
void writeSuffixarray(Suffixarray *s, char *filename);
Suffixarray* constructSufArr(void *, CharSequence **, Uint, FAlphabet *, unsigned char silent);
void constructchildtab(void *, Suffixarray *);
void constructsuflinks(void *, Suffixarray *, Uint *);
void constructLcp (void *, Suffixarray *);
void computeId(void*, Suffixarray *);
Uint* getsufsucc(void *, Suffixarray *);
void checksuflinks(Suffixarray *s, Uint i, Uint j);
void destructSufArr (void *, Suffixarray *); 
extern Container* getChildintervals(void *, Suffixarray *, Uint, Uint, BOOL); 
extern Lint* getChildintervalsArr(void *, Suffixarray *, Uint, Uint, Uint *, BOOL); 
extern PairUint getCharInterval(void *, Suffixarray *, Uint, Uint, Uint, char);
extern PairUint getCharIntervalArr(void *, Suffixarray *, Uint, Uint, Uint, char);
extern PairUint getSuflink(Suffixarray *, Uint, Uint);
extern PairUint jumpkSuflinks(Suffixarray *s, Uint i, Uint j, Uint k);
extern Uint getlcpval(Suffixarray *, Uint, Uint);
extern Uint getfirstlindex(Suffixarray *, Uint, Uint);
extern Lint id (Suffixarray *, Uint);
extern Uint lcp(Suffixarray *s, Uint i);
extern Uint maxlcp(Suffixarray *s);
extern unsigned char isnextlIndex(Suffixarray *s, Uint i);
extern unsigned char isdownIndex(Suffixarray *s, Uint i);
extern unsigned char isupIndex(Suffixarray *s, Uint i);
extern PairUint jumpkSuflinks(Suffixarray *s, Uint i, Uint j, Uint k);
extern void addinterval(void *space, Container *c, Uint a, Uint b);
void destructinterval(void *space, void *data);
void dumplcps (Suffixarray *);
PairUint searchSuffix (void *space, Suffixarray *arr, char *p, Uint len);
void dumpSufArr (Suffixarray *arr);
Uint* searchSuffixList (void *space, Suffixarray *arr, char *p, Uint m, Uint maxsize, Uint *size);
PairUint* searchSuffixPath (void *space, Suffixarray *arr, char *p, Uint len, Uint *last);
PairUint* suffixArrayCreateHash(void* space, Suffixarray *arr, Uint size);
PairUint searchSuffixOff (void *space, Suffixarray *arr, Uint i, Uint j, char *p, Uint pstart, Uint len);
PairUint searchSuffixArrayHash (void *space, Suffixarray *arr, PairUint *hash, Uint size, char *p, Uint len);
unsigned int translateSufArrPos(Suffixarray *esa, Uint pos, Uint *idxi, Uint *l, Uint *r);
void constructInvSufTab(Suffixarray *esa); 
#endif

