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
 *  alignment.c
 *  implementation to handle alignments
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 02/03/2009 02:50:06 PM CET
 *  
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include "mathematics.h"
#include "basic-types.h"
#include "iupac.h"
#include "debug.h"
#include "charsequence.h"
#include "biofiles.h"
#include "locus.h"

const char ntdecode[]  = {'A', 'C', 'G', 'T', '-', 'N' };
const char decodeEop[] = {'R','D','I','N','S','=','X'};

char*
getNTcodekey(void *space) {
  char *code;
  code = ALLOCMEMORY(space, NULL, char, 256);
  memset(code, 5, sizeof(char)*256);

  code['A'] = 0;
  code['a'] = 0;
  code['C'] = 1;
  code['c'] = 1;
  code['G'] = 2;
  code['g'] = 2;
  code['T'] = 3;
  code['t'] = 3;
  code['-'] = 4;
  
  return code;
}

void 
initAlignment(Alignment *al, 
    char *u, Uint ulen, Uint uoff, 
    char *v, Uint vlen, Uint voff) {
  
  assert(uoff < ulen && voff < vlen);
  
  al->u = u;
  al->v = v;
  al->ulen = ulen;
  al->vlen = vlen;
  al->uoff = uoff;
  al->voff = voff;
  al->numofmeops = 0;
  al->meops = malloc(sizeof(Multieop)*(ulen+vlen));
  memset(al->meops, 0, sizeof(Multieop)*(ulen+vlen));
  al->rmvseq = 0;
  al->rmuseq = 0;
}

void 
wrapAlignment(Alignment *al) {
  free(al->meops);
  al->meops = NULL;
  al->numofmeops = 0;
  al->vlen = 0;
  al->ulen = 0;
  al->uoff = 0;
  al->voff = 0;
  if(al->rmvseq) FREEMEMORY(NULL, al->v);
  if(al->rmuseq) FREEMEMORY(NULL, al->u);
  al->u = NULL;
  al->v = NULL;
}

void
copyAlignment(Alignment *to, Alignment *from) {
  //to->meops = malloc(sizeof(Multieop)*(from->ulen+from->vlen));

  to->meops = malloc(sizeof(Multieop)*(from->numofmeops));
  memmove(to->meops, from->meops, sizeof(Multieop)*(from->numofmeops));
  to->numofmeops = from->numofmeops; 
  to->uoff = from->uoff;
  to->voff = from->voff;
  to->rmvseq = from->rmvseq;
  to->rmuseq = from->rmuseq;

  if(from->rmvseq) {
    to->v = malloc(sizeof(char)*(from->vlen+1));
    memmove(to->v, from->v, sizeof(char)*from->vlen);
    to->v[from->vlen] = 0;
  } else {
    to->v = from->v; 
  }
 
  if(from->rmuseq) {
    to->u = malloc(sizeof(char)*(from->ulen+1));
    memmove(to->u, from->u, sizeof(char)*from->ulen);
    to->u[from->ulen] = 0;
  } else {
    to->u = from->u; 
  }
  to->vlen = from->vlen;
  to->ulen = from->ulen;
}


/*--------------------------------- isMatch ----------------------------------
 *    
 * @brief check if alignment has u[i] == v[j]
 * @author Steve Hoffmann 
 *   
 */
 
char
isMatch (Alignment *al, Uint i, Uint j)
{
  if(!matchIUPAC(al->u[i+al->uoff], al->v[j+al->voff])) 
    return 0;

  return 1;
}


void
countEops(Alignment *al, Uint *mat, Uint *mis, Uint *ins, Uint *del, Uint *lmat) {
  Uint i,j,k=0,l=0,matchrun=0;

  *mat = 0;
  *mis = 0;
  *ins = 0;
  *del = 0;
  *lmat = 0; //longest match

  for(i=0; i < al->numofmeops; i++) {
    if(al->meops[i].eop != Replacement && al->meops[i].eop != Match 
        && al->meops[i].eop != Mismatch) {

      if (al->meops[i].eop == Deletion) {
        *del += al->meops[i].steps;
        l += al->meops[i].steps;
        matchrun = 0;
      } 

      if(al->meops[i].eop == Insertion){
        *ins += al->meops[i].steps;
        k += al->meops[i].steps;
        matchrun = 0;
      }

      if(al->meops[i].eop == Softclip){
        k += al->meops[i].steps;
        matchrun = 0;
      }

      if(al->meops[i].eop == Skipped) {
        //do nothing, count nothing
      }

    } else {
      if(al->meops[i].eop == Replacement) { 
        for(j=0; j < al->meops[i].steps; j++) {
          if(!matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff])) { 
            matchrun =0;
            *mis += 1; 
          } else { 
            if(*lmat < matchrun+1) {
              *lmat = matchrun+1;
            }
            matchrun += 1;
            *mat += 1;
          }

          k++; l++;
        }
      }  
      
      if(al->meops[i].eop == Mismatch) {
        matchrun = 0;
        *mis += al->meops[i].steps;
        k+= al->meops[i].steps;
        l+= al->meops[i].steps;
      }

      if (al->meops[i].eop == Match){
        
        if (*lmat < matchrun + al->meops[i].steps) {
          *lmat = matchrun + al->meops[i].steps; 
        }

        matchrun += al->meops[i].steps;
        *mat += al->meops[i].steps;
        k+= al->meops[i].steps;
        l+= al->meops[i].steps;
      }
    }
  }

  return;
}

Uint
getEdist(Alignment *al) {

  Uint mat;
  Uint mis;
  Uint ins;
  Uint del;
  Uint lmat;

  countEops(al, &mat, &mis, &ins, &del, &lmat) ;
  return mis+ins+del;
}

Uint*
getSplitEdist(Alignment *al, Uint *noofsplits) {
  Uint i,j,k=0,l=0,u=0;
  Uint *edist=NULL, n=1;
  edist = ALLOCMEMORY(NULL, NULL, Uint, n);
  edist[u]=0;

  for(i=0; i < al->numofmeops; i++) {
    if(al->meops[i].eop != Replacement && al->meops[i].eop != Match 
        && al->meops[i].eop != Mismatch) {

      if (al->meops[i].eop == Deletion) {
        l+= al->meops[i].steps;
        edist[u] += al->meops[i].steps;
      } 

      if(al->meops[i].eop == Insertion) {
        k+= al->meops[i].steps;
        edist[u] += al->meops[i].steps;
      }

      if(al->meops[i].eop == Softclip) {
        k+= al->meops[i].steps;
        //edist[u] += al->meops[i].steps;
      }

      if(al->meops[i].eop == Skipped) {
        edist = ALLOCMEMORY(NULL, edist, Uint, n+1);
        edist[n] = 0;
        u+=1;
        n+=1;
      }

    } else {
      if(al->meops[i].eop == Replacement) { 
        for(j=0; j < al->meops[i].steps; j++) {
          if(!matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff])) edist[u]++;
          k++; l++;
        }
      }

      if(al->meops[i].eop == Mismatch) {
        edist[u] += al->meops[i].steps;
        k+= al->meops[i].steps;
        l+= al->meops[i].steps;
      }

      if (al->meops[i].eop == Match){
        k+= al->meops[i].steps;
        l+= al->meops[i].steps;
      }
    }
  }

  *noofsplits = n;
  return edist;
}


Uint
getBisulfiteMismatches(Alignment *al, Uint bisulfite){
  Uint i,j,k=0,l=0;
  Uint mis=0;
  char *seq;

  seq = malloc(al->ulen);
  memmove(seq, al->u, al->ulen);
  bl_reconvertBisulfite(seq, al->ulen, bisulfite);

  for(i=0; i < al->numofmeops; i++) {
    if(al->meops[i].eop != Replacement && al->meops[i].eop != Match 
       && al->meops[i].eop != Mismatch) {

      if (al->meops[i].eop == Deletion) {
        l+= al->meops[i].steps;
      } 
      
      if(al->meops[i].eop == Insertion){
        k+= al->meops[i].steps;
      }

      if(al->meops[i].eop == Softclip){
        k+= al->meops[i].steps;
      }
      if(al->meops[i].eop == Skipped){
        //do nothing
      }

    } else {
      for(j=0; j < al->meops[i].steps; j++) {
        /* bisulfite mismatch: IUPAC match on converted query but mismatch on recoverted one */
        if(matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff]) &&
           seq[k+al->uoff] != al->v[l+al->voff]) mis++;
        k++; l++;
      }
    }
  }
  free(seq);
  return mis;
}

Uint getWrongStrandBisulfiteMismatches(Alignment *al, Uint bisulfite){
  Uint i,j,k=0,l=0;
  Uint mis=0;
  char *seq;

  seq = malloc(al->ulen);
  memmove(seq, al->u, al->ulen);
  /* get other bisulfite run by altering bisulfite (e.g. 1=>2, 2=>1) */
  bl_convertBisulfite(seq, al->ulen, (bisulfite % 2) ? bisulfite + 1 : bisulfite - 1, 0);

  for(i=0; i < al->numofmeops; i++) {
    if(al->meops[i].eop != Replacement && al->meops[i].eop != Match 
        && al->meops[i].eop != Mismatch) {
      if (al->meops[i].eop == Deletion) {
        l+= al->meops[i].steps;
      } 
      if (al->meops[i].eop == Insertion){
        k+= al->meops[i].steps;
      }

      if (al->meops[i].eop == Softclip){
        k+= al->meops[i].steps;
      }

      if(al->meops[i].eop == Skipped) {
        //do nothing
      }

    } else {
      for(j=0; j < al->meops[i].steps; j++) {
        /* 
         * wrong strand bisulfite mismatch: IUPAC mismatch on converted query
         * but match on coverted one in other bisulfite run
         */
        if(!matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff]) &&
           matchIUPAC(seq[k+al->uoff], al->v[l+al->voff])) mis++;
        k++; l++;
      }
    }
  }
  free(seq);
  return mis;
}

void
getSoftClipScores(Alignment *al, int polyAlen, int *scores, int indel, 
    int *pAscr, int *adscr, int *adlen) {

  Uint i,j,k=0,l=0;
  int polyAscore=0, adapterScore=0, adapterLen=0;

  for(i=0; i < al->numofmeops; i++) {

    if(k+al->uoff < polyAlen) { 
      if(al->meops[i].eop != Replacement && al->meops[i].eop != Match 
        && al->meops[i].eop != Mismatch) {
        polyAscore += indel * al->meops[i].steps;
        if (al->meops[i].eop == Deletion) {
          l+= al->meops[i].steps;
        } 
        
        if(al->meops[i].eop == Insertion){
          k+= al->meops[i].steps;
        }
        
        if(al->meops[i].eop == Softclip){
          k+= al->meops[i].steps;
        }

        if(al->meops[i].eop == Skipped) {
        //do nothing
        }

      } else {
        for(j=0; j < al->meops[i].steps; j++) {
          if(al->u[k+al->uoff] != 'N' && al->v[k+al->voff] != 'N' 
              && !matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff])) {
            polyAscore += scores[1];
          } else {
            polyAscore += scores[0];
          }

          k++; l++;
        }
      }
    } else {
      adapterLen++;
      if(al->meops[i].eop != Replacement && al->meops[i].eop != Match 
        && al->meops[i].eop != Mismatch) {
        adapterScore += indel * al->meops[i].steps;
        if (al->meops[i].eop == Deletion) {
          l+= al->meops[i].steps;
        } 
        if (al->meops[i].eop == Insertion){
          k+= al->meops[i].steps;
        }
  
        if (al->meops[i].eop == Softclip){
          k+= al->meops[i].steps;
        }

        if(al->meops[i].eop == Skipped) {
        //do nothing
        }

      } else {
        for(j=0; j < al->meops[i].steps; j++) {
          if(!matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff])) {
            adapterScore += scores[1];
          } else {
            adapterScore += scores[0];
          }

          k++; l++;
        }
      }
    }
  }

  *pAscr = polyAscore;
  *adscr = adapterScore;
  *adlen = adapterLen;

  return;
}

int
getAlignScore(Alignment *al, int *scores, int indel) {
  Uint i,j,k=0,l=0;
  int score=0;

  for(i=0; i < al->numofmeops; i++) {

    if(al->meops[i].eop != Replacement && al->meops[i].eop != Match 
        && al->meops[i].eop != Mismatch) {

      if (al->meops[i].eop == Deletion) {
        score += indel * al->meops[i].steps;
        l+= al->meops[i].steps;
      } 

      if(al->meops[i].eop == Insertion){
        score += indel * al->meops[i].steps;
        k+= al->meops[i].steps;
      }

      if(al->meops[i].eop == Softclip){
        k+= al->meops[i].steps;
      }

      if(al->meops[i].eop == Skipped){
        //no moves, no penalties
      }

    } else {
      if(al->meops[i].eop == Replacement || al->meops[i].eop == Mismatch) { 
        for(j=0; j < al->meops[i].steps; j++) {
          if(!matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff])) {
            score += scores[1];
          } else {
            score += scores[0];
          }

          k++; l++;
        }
      }
     
/*      if(al->meops[i].eop == Mismatch) {
        score += al->meops[i].steps*scores[1];
        k+= al->meops[i].steps;
        l+= al->meops[i].steps;
      }
*/
      if (al->meops[i].eop == Match){
        score +=  al->meops[i].steps*scores[0];
        k+= al->meops[i].steps;
        l+= al->meops[i].steps;
      }
    }
  }
  return score;
}

/*int
getAlignScore(Alignment *al, int *scores, int indel) {
  Uint i,j,k=0,l=0;
  int score=0;

  for(i=0; i < al->numofmeops; i++) {

    if(al->meops[i].eop != Replacement) {
      score += indel * al->meops[i].steps;
      if (al->meops[i].eop == Deletion) {
        l+= al->meops[i].steps;
      } 
      if (al->meops[i].eop == Insertion){
        k+= al->meops[i].steps;
      }

      if (al->meops[i].eop == Softclip){
        k+= al->meops[i].steps;
      }

      if(al->meops[i].eop == Skipped) {
        // no moves, no penalties
      }

    } else {
      for(j=0; j < al->meops[i].steps; j++) {
        if(!matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff])) {
          score += scores[1];
        } else {
          score += scores[0];
        }

        k++; l++;
      }
    }
  }
  return score;
}*/


/*---------------------------- getSubstringEdist -----------------------------
 *    
 * @brief get the edist of a substring of an aligned pair of strings
 * @author Steve Hoffmann 
 *   
 */
int
getSubstringEdist(Alignment *al, Uint u, Uint v) {

  Uint i, j, k=0, l=0, edist=0, scr=0;

  for(i=0; i < al->numofmeops; i++) {

    if(al->meops[i].eop != Replacement && al->meops[i].eop != Match 
        && al->meops[i].eop != Mismatch) {

      if(al->meops[i].eop != Skipped) { 
        if(k >= u && k < v) { 
          edist += al->meops[i].steps;
        }

        if (al->meops[i].eop == Deletion) {
          l+= al->meops[i].steps;
        } 

        if (al->meops[i].eop == Insertion){
          k+= al->meops[i].steps;
        }

        if (al->meops[i].eop == Softclip){
          k+= al->meops[i].steps;
        }
      } else {
        //        k+= al->meops[i].steps; no penalties, no moves
      }

    } else {

      Uint check1 =0;
      Uint check2 =0;
      Uint lprime =l;
      Uint kprime =k;
      
      if(al->meops[i].eop == Replacement || al->meops[i].eop == Mismatch) { 
        for(j=0; j < al->meops[i].steps; j++) {
          if(k >=u && k < v && !matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff])) { 
            check1++;
          } else if (k >=u && k<v) scr += 1;
 
          k++; l++;
        }
      }

      if(al->meops[i].eop == Mismatch) {
        k = kprime;
        l = lprime;
        if(k + al->meops[i].steps -1 >= u && k < v) {
          
          //right overlap
          if(k >= u && k + al->meops[i].steps-1 >= v) {
            check2 += v - k;
          } 
          //complete overlap
          if(k < u && k + al->meops[i].steps -1 >= v) {
            check2 += v - u;
          }
          //left overlap
          if(k < u && k + al->meops[i].steps -1 >= u && k + al->meops[i].steps -1 < v) {
            check2 += k + al->meops[i].steps - u; 
          }
          //inclusion
          if(k>=u && k + al->meops[i].steps -1 < v) {
            check2 += al->meops[i].steps;
          }
        }

//        k+= al->meops[i].steps;
//        l+= al->meops[i].steps;
      }

      //check!
      if(al->meops[i].eop == Mismatch) { 
        assert(check1 == check2);
        edist += check1;

        k+= al->meops[i].steps;
        l+= al->meops[i].steps;
      }

      //check!
      if(al->meops[i].eop == Replacement) {
        edist += check1;

        k+= al->meops[i].steps;
        l+= al->meops[i].steps;
      }


      if (al->meops[i].eop == Match){
        k+= al->meops[i].steps;
        l+= al->meops[i].steps;
      }

    }
  }
  return edist;
}


/*----------------------------- getSubAlignment ------------------------------
 *    
 * @brief get partial alignment
 * @author Steve Hoffmann 
 *   
 */
 
Alignment*
getSubAlignment (Alignment *al, Uint l, Uint r)
{
	  
  Uint i, j, p=0, q=0;
  Alignment *new;

  new = ALLOCMEMORY(NULL, NULL, Alignment, 1);
  initAlignment(new, al->u, al->ulen, al->uoff, al->v, al->vlen, al->voff);


  for(i=0; i < al->numofmeops; i++) {
    for(j=0; j < al->meops[i].steps; j++) { 

      if(q+al->voff <= l) {
        new->uoff = al->uoff+p;
        new->voff = al->voff+q;
      }

      if(q+al->voff >= l && q+al->voff <= r) {
        insertEop(new, al->meops[i].eop);
      }
      
      if(al->meops[i].eop == Replacement || al->meops[i].eop == Match 
          || al->meops[i].eop == Mismatch) {
        p++;
        q++;
      }

      if (al->meops[i].eop == Deletion) {
        q++;
      } 

      if (al->meops[i].eop == Insertion){
        p++;
      }
    
      if (al->meops[i].eop == Softclip){
        p++;
      }

    }
  }

  return new;
}

Alignment*
expandAlignment(Alignment *al, Uint* expandpos, Uint *expandlen) { 
  Uint i, j, p=0, q=0, k=0;
  Alignment *new;

  new = ALLOCMEMORY(NULL, NULL, Alignment, 1);
  initAlignment(new, al->u, al->ulen, al->uoff, al->v, al->vlen, al->voff);
 
  fprintf(stdout, "uoff: %d, voff: %d\n%s(%d) -> (%d)\n", al->uoff, al->voff, al->u, al->ulen, al->vlen);

  if(al->voff > expandpos[k]) {
    new->voff += expandlen[k];
    fprintf(stdout, "new expandpos %d, len:%d\n", expandpos[k+1], expandlen[k+1]);
    k++;
  }


  for(i=0; i < al->numofmeops; i++) {
    for(j=0; j < al->meops[i].steps; j++) { 

      if(q+al->voff > expandpos[k]) {
        insertMeop(new, Deletion, expandlen[k]);
        fprintf(stdout, "expandlen[%d] %d nucleotides at meop:%d step:%d; expandpos:%d\n", k, expandlen[k], i, j, expandpos[k]);
     //   q += expandlen[k];
        k++;
      }
         
      insertEop(new, al->meops[i].eop);     
      
      if(al->meops[i].eop == Replacement || al->meops[i].eop == Match 
          || al->meops[i].eop== Mismatch) {
        p++;
        q++;
      }

      if (al->meops[i].eop == Deletion) {
        q++;
      } 

      if (al->meops[i].eop == Insertion){
        p++;
      }

      if (al->meops[i].eop == Softclip){
        p++;
      }


    }
  }

  fprintf(stdout, "expanding alignment ended\n");
  return new;

}

 
/*----------------------------- showmultieoplist -----------------------------
 *    
 * @brief show the multi edit operation list
 * @author Steve Hoffmann 
 * dumps visual representation of alignments and should be shorter!
 */

void 
showmultieoplist(FILE *dev, Alignment *al) {

  Uint i=0;
  fprintf(dev, "[");
  if(al->numofmeops) {
    for(i=0; i < al->numofmeops; i++) {
      fprintf(dev, "%c %d, ", decodeEop[al->meops[i].eop], al->meops[i].steps);
    }
  fprintf(dev, "%c %d",decodeEop[al->meops[i].eop], al->meops[i].steps);    
  }
  fprintf(dev, "]\n");
}

//dumps visual representation of alignments and should be shorter!
char *
multieopstring(Alignment *al, Uint leftclip, Uint rightclip, unsigned char rev) {
  Uint i, j, k, q=0, p=0, cur=0, strsize, steps, msteps, ssteps;
  char *meopstr;
  char eopc=0;

  meopstr = (char*) malloc(sizeof(char)*(3*(al->vlen+al->ulen+leftclip+rightclip)+1));

  if(leftclip || (rightclip && rev)) {
    steps = (rev) ? rightclip : leftclip;
    //strsize = floor(log(steps)/log(10))+3;
    strsize = snprintf(NULL, 0, "%d", steps)+2;
    meopstr[cur] = 'C';
    sprintf(&meopstr[cur+1], "%d;", steps);
    cur+=strsize;
  }

  for(k=0; k < al->numofmeops; k++) {
    i = (rev) ? al->numofmeops - k -1 : k;
    //if Replacement occured
    steps=0;
    if (al->meops[i].eop == Replacement) {
      //iter over all steps
      msteps=0;
      ssteps=0;
      for (j=0; j < al->meops[i].steps; j++) {
        if (!matchIUPAC(al->u[j+p+al->uoff], al->v[j+q+al->voff])) {
          if (j==0 || eopc == 'S') {
            ssteps++;
          } else {
            //strsize = floor(log(msteps)/log(10))+3;
            strsize = snprintf(NULL, 0, "%d", msteps)+2;
            meopstr[cur] = eopc;
            sprintf(&meopstr[cur+1], "%d;", msteps);
            cur+=strsize;
            msteps=0;
            ssteps=1;
          }
          eopc = 'S';
        } else {
          if (j==0 || eopc == 'M') {
            msteps++;
          } else {
            //strsize = floor(log(ssteps)/log(10))+3;
            strsize = snprintf(NULL, 0, "%d", ssteps)+2;
            meopstr[cur] = eopc;
            sprintf(&meopstr[cur+1], "%d;", ssteps);
            cur+=strsize;
            msteps=1;
            ssteps=0;
          }
          eopc = 'M';
        }
      }
      steps = msteps + ssteps;
      assert(msteps == 0 || ssteps == 0);
      //set string ptrs
      p+=j;
      q+=j;
    }
    //if a match occured
    if(al->meops[i].eop == Match) {
      eopc = 'M';
      steps = al->meops[i].steps;
      p+=steps;
      q+=steps;
    }
    //if a mismatch occured 
    if(al->meops[i].eop == Mismatch) {
      eopc = 'S';
      steps = al->meops[i].steps;
      p+=steps;
      q+=steps;
    }
    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      eopc = 'D';
      //set ptrs
      steps = al->meops[i].steps;
      q+=steps;
    }

    //if insertions occured
    if(al->meops[i].eop == Insertion)  {
      eopc = 'I';  
      steps = al->meops[i].steps;
      p+=steps;
    }

    //Softclipping
    if(al->meops[i].eop == Softclip)  {
      eopc = '^';  
      steps = al->meops[i].steps;
      p+=steps;
    }

    //if skipping occured
    if(al->meops[i].eop == Skipped){
      eopc = 'N';
      steps = al->meops[i].steps;
    //  p +=steps; no moves
    }

    //strsize = floor(log(steps)/log(10))+3;
    strsize = snprintf(NULL, 0, "%d", steps)+2;
    meopstr[cur] = eopc;
    sprintf(&meopstr[cur+1], "%d;", steps);
    cur+=strsize;
  }

  if(rightclip || (leftclip && rev)) {
    steps = (rev) ? leftclip : rightclip;
    //strsize = floor(log(steps)/log(10))+3;
    strsize = snprintf(NULL, 0, "%d", steps)+2;
    meopstr[cur] = 'C';
    sprintf(&meopstr[cur+1], "%d;", steps);
    cur+=strsize;
  }
  return meopstr;
}

char*
mdstring(Alignment *al, unsigned char rev) {
  Uint i, j, k, q=0, p=0, cur=0, strsize, steps, msteps=0, ssteps=0;
  char *mdstr;
  char eopc=0;

  mdstr = (char*) malloc(sizeof(char)*(3*(al->vlen+al->ulen)+1));
  memset(mdstr, 0, sizeof(char)*(3*(al->vlen+al->ulen)+1));

    msteps=0;
    ssteps=0;
  for(k=0; k < al->numofmeops; k++) {
    i = (rev) ? al->numofmeops - k - 1 : k;
    //if Replacement occured
    steps=0;
    
    if (al->meops[i].eop == Replacement || al->meops[i].eop == Match || al->meops[i].eop == Mismatch) {
      //iter over all steps

      for (j=0; j < al->meops[i].steps; j++) {  
        if (!matchIUPAC(al->u[j+p+al->uoff], al->v[j+q+al->voff])) {
          if(msteps) {
            //strsize = floor(log(msteps)/log(10))+1;
            strsize = snprintf(NULL, 0, "%d", msteps);
            sprintf(&mdstr[cur], "%d", msteps);
            cur += strsize;
            msteps = 0;
          }
          if(eopc != 'M') {
            sprintf(&mdstr[cur], "0");
            cur += 1;
          }
          sprintf(&mdstr[cur], "%c", al->v[j+q+al->voff]);
          cur += 1;
          ssteps++;
          eopc = 'S';
        } else {
          if (msteps) {
            msteps++;
          } else {
            msteps=1;
            ssteps=0;
          }
          eopc = 'M';
        }
      }
      steps = msteps + ssteps;
      assert(msteps == 0 || ssteps == 0);
      //set string ptrs
      p+=j;
      q+=j;
    } 

    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      if (msteps) {
        //strsize = floor(log(msteps)/log(10))+1;
        strsize = snprintf(NULL, 0, "%d", msteps);
        sprintf(&mdstr[cur], "%d", msteps);
        cur += strsize;
        msteps = 0;
      } else { 
        sprintf(&mdstr[cur], "0");
        cur += 1;
      }
      
      eopc = 'D';
      //set ptrs
      steps = al->meops[i].steps;
      sprintf(&mdstr[cur], "^");
      cur+=1;
      for(j=0; j < steps; j++) {
        sprintf(&mdstr[cur], "%c", al->v[j+q+al->voff]);
        cur+=1;
      }
      q+=steps;
    }

    //if insertions occured
    if(al->meops[i].eop == Insertion)  {
      //eopc = 'I'; 
      steps = al->meops[i].steps;
      p+=steps;
    }

    //if insertions occured
    if(al->meops[i].eop == Softclip)  {
      //eopc = 'I'; 
      steps = al->meops[i].steps;
      p+=steps;
    }


    //if skipping occured
    if(al->meops[i].eop == Skipped) {
      steps = al->meops[i].steps;
     // p+=steps; no moves!
    }
  }

  if(eopc != 'M') {
    sprintf(&mdstr[cur], "0");
    cur += 1;
  }
  
  if (msteps) {
    //strsize = floor(log(msteps)/log(10))+1;
    strsize = snprintf(NULL, 0, "%d", msteps);
    sprintf(&mdstr[cur], "%d", msteps);
    cur += strsize;
    msteps = 0;
  } 

  return mdstr;
}



/*-------------------------- bl_cigarGetAlignString --------------------------
 *    
 * @brief decode a cigar string
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_cigarGetAlignString(char *cigar, uint64_t **usplits, uint64_t **vsplits, 
    Uint *nsplits) {
  Uint i, len, allen=0, cur=0;
  Uint upos=0, vpos=0;
  uint64_t *myusplits = NULL;
  uint64_t *myvsplits = NULL;
  Uint mynsplits = 0;;
  char *buffer, *string = NULL;

  len = strlen(cigar);
  buffer = calloc(len, sizeof(char));

  for(i=0; i < len; i++) {
    switch (cigar[i]) {
      case 'S':
        string = realloc(string, allen+atoi(buffer)+1);
        memset(&string[allen], 'S', atoi(buffer));
        allen += atoi(buffer);
        upos+=atoi(buffer);
        string[allen] = 0;
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'M':
      case 'X':
      case '=':
        string = realloc(string, allen+atoi(buffer)+1);
        memset(&string[allen], 'M', atoi(buffer));
        allen += atoi(buffer); 
        upos+= atoi(buffer);
        vpos+= atoi(buffer);
        string[allen] = 0;
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'D':         
        string = realloc(string, allen+atoi(buffer)+1);
        memset(&string[allen], 'D', atoi(buffer));
        allen += atoi(buffer);
        vpos+=atoi(buffer);
        string[allen] = 0;
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'I':
        string = realloc(string, allen+atoi(buffer)+1);
        memset(&string[allen], 'I', atoi(buffer));
        allen += atoi(buffer);
        upos += atoi(buffer);
        string[allen] = 0;
        memset (buffer, 0, len);
        cur =0;
        break;
      case 'N':

        myvsplits = ALLOCMEMORY(NULL, myvsplits, uint64_t, mynsplits+2);
        myusplits = ALLOCMEMORY(NULL, myusplits, uint64_t, mynsplits+2);
        myvsplits[mynsplits] = vpos;
        vpos += atoi(buffer);
        myvsplits[mynsplits+1] = vpos;
        myusplits[mynsplits] = upos;
        upos+=1;
        myusplits[mynsplits+1] = upos;
        mynsplits += 2;

        memset(buffer, 0, len);
        cur = 0;
        break;  //nomoves on skipping
      default :
        buffer[cur++] = cigar[i];
    }
  }

  free(buffer);
  *usplits = myusplits;
  *vsplits = myvsplits;
  *nsplits = mynsplits;
  return string;
}


/*--------------------------- bl_cigarGetAlignLen ----------------------------
 *    
 * @brief get alignment length from a cigar string
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_cigarGetAlignLen(char *cigar) {
  Uint i, len, allen=0, cur=0;
  char *buffer;

  len = strlen(cigar);
  buffer = calloc(len, sizeof(char));

  for(i=0; i < len; i++) {
    switch (cigar[i]) {
      case 'M':
      case 'X':
      case '=':
        allen += atoi(buffer);
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'I': 
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'D':
        allen += atoi(buffer);
        memset (buffer, 0, len);
        cur =0;
        break;
      case 'N':
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'S':
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'H':
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'P':
        memset (buffer, 0, len);
        cur = 0;
        break;
      default :
        buffer[cur++] = cigar[i];
    }
  }

  free(buffer);
  return allen;
}

char*
bl_mdGetDiffString(char *MD) {
  Uint i, k=0, MDlen, allen=0, nops=0, buffersize=100;
  char chr;
  unsigned char del = 0;
  char *buffer = NULL;
  char *alignstr = NULL;

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);
  memset(buffer, 0, sizeof(char)*buffersize);
  MDlen = strlen(MD); 

  for(i=0; i < MDlen; i++) {
    
    if(!isalpha((int)MD[i]) && MD[i] != '^') {
      if(k >= buffersize-2) {
        buffersize += 100;
        buffer = ALLOCMEMORY(space, buffer, char, buffersize);
      }
      buffer[k] = MD[i];
      buffer[k+1] = 0;
      k++;
      del = 0;
    } else {
      
      nops = atoi(buffer);
      if(nops) {
        alignstr = ALLOCMEMORY(space, alignstr, char, allen+nops+1);
        memset(&alignstr[allen], 'M', nops);
        alignstr[allen+nops] = 0;
        allen+=nops;
       }
       memset(buffer, 0, sizeof(char)*buffersize);
       k=0;

      if(MD[i] == '^') {
        del = 1;
      } else {
        chr = (del) ? 'D' : MD[i];
        alignstr = ALLOCMEMORY(space, alignstr, char, allen+2);
        alignstr[allen] = chr;
        alignstr[allen+1] = 0;
        allen += 1;
      }
    }
  }

  nops = atoi(buffer);
  if(nops) {
    alignstr = ALLOCMEMORY(space, alignstr, char, allen+nops+1);
    memset(&alignstr[allen], 'M', nops);
    alignstr[allen+nops] = 0;
    allen+=nops;
  }

  FREEMEMORY(space, buffer);
  return alignstr;
}

char*
cigarstring(Alignment *al, Uint leftclip, Uint rightclip, char clipch, unsigned char rev, char brief) {
  Uint i, j, k, q=0, p=0, cur=0, strsize, steps=0, msteps;
  char *cigarstr;
  char eopc=0;

  cigarstr = (char*) malloc(sizeof(char)*(3*(al->vlen+al->ulen+leftclip+rightclip)+1));

  if(leftclip || (rightclip && rev)) {
    steps = (rev) ? rightclip : leftclip;
    //    strsize = floor(log(steps)/log(10))+2;
    strsize = snprintf(NULL, 0, "%d", steps)+1;
    sprintf(&cigarstr[cur], "%d%c", steps, clipch);
    cur+=strsize;
  }

  for(steps=0, k=0; k < al->numofmeops; k++) {
    i = (rev) ? al->numofmeops - k - 1 : k;
    //if Replacement occured 
    if (al->meops[i].eop == Replacement) {
      //iter over all steps
      msteps=0;
      for (j=0; j < al->meops[i].steps; j++) {

        if (brief)
          eopc = 'M';
        else { 
          if (!matchIUPAC(al->u[j+p+al->uoff], al->v[j+q+al->voff])) {
            if(eopc == '=') {
              strsize = snprintf(NULL, 0, "%d", msteps)+1;
              sprintf(&cigarstr[cur], "%d%c", msteps, eopc);
              cur+=strsize;
              msteps = 0;
            }
            eopc = 'X';
          } else {
            if(eopc == 'X') {
              strsize = snprintf(NULL, 0, "%d", msteps)+1;
              sprintf(&cigarstr[cur], "%d%c", msteps, eopc);
              cur+=strsize;
              msteps = 0;
            }
            eopc = '=';
          }
        }
        msteps++;
      }

      steps += msteps;
      //set string ptrs
      p+=j;
      q+=j;
    }
    //if match occured
    if (al->meops[i].eop == Match) {
      if(brief)
        eopc = 'M';
      else
        eopc = '=';
      //set ptrs
      //steps = al->meops[i].steps;
      steps += al->meops[i].steps;
      p+=steps;
      q+=steps;
    }
    //if mismatch occured
    if (al->meops[i].eop == Mismatch) {
      if(brief)
        eopc = 'M';
      else
        eopc = 'X';
      //set ptrs
      //steps = al->meops[i].steps;
      steps += al->meops[i].steps;
      p+=steps;
      q+=steps;
    }

    //dump previous matches and mismatches if in brief mode
    //if a different operation is coming
    if(brief && steps && eopc == 'M' && 
        al->meops[i].eop != Match && al->meops[i].eop != Mismatch) {
      
      strsize = snprintf(NULL, 0, "%d", steps)+1;
      sprintf(&cigarstr[cur], "%d%c", steps, eopc);
      cur+=strsize;
      steps = 0;
    }

    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      eopc = 'D';
      //set ptrs
      steps = al->meops[i].steps;
      q+=steps;
    }
    //if insertions occured
    if(al->meops[i].eop == Insertion)  {
      eopc = 'I';  
      steps = al->meops[i].steps;
      p+=steps;
    }

    if(al->meops[i].eop == Softclip)  {
      eopc = 'S';  
      steps = al->meops[i].steps;
      p+=steps;
    }

    //if skipping occured
    if(al->meops[i].eop == Skipped)  {
      eopc = 'N';  
      steps = al->meops[i].steps;
      //     p+=steps; no moves!
    }

    //dump everything except "M"s in the brief mode
    //unless we are at the end
    if(eopc != 'M' || k+1 == al->numofmeops) {
      //strsize = floor(log(steps)/log(10))+2;
      strsize = snprintf(NULL, 0, "%d", steps)+1;
      sprintf(&cigarstr[cur], "%d%c", steps, eopc);
      cur+=strsize;
      steps = 0;
    }
  }

  if(rightclip || (leftclip && rev)) {
    steps = (rev) ? leftclip : rightclip;
    //strsize = floor(log(steps)/log(10))+2;
    strsize = snprintf(NULL, 0, "%d", steps)+1;
    sprintf(&cigarstr[cur], "%d%c", steps, clipch);
    cur+=strsize;
  }
  return cigarstr;
}

//shows muliteoplist of all Alignments in *al
void 
showDynmultieoplist(Alignment* al, int size) {

  int i;
  for (i=0; i < size; i++) {
    showmultieoplist(stdout, &al[i]);
  }
}

//dumps visual representation of alignments and should be shorter!
char*
getAlignString(Alignment* al, char lf) {
  int i, j , k, nlines, len;
  char *mystring = NULL;

  Uint p = 0, q = 0, r = 0;
  //output strings
  char* utemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* vtemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* comp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));

  //iter over all multieops
  for (i=0; i < al->numofmeops; i++) {
    //if Replacement occured
    if (al->meops[i].eop == Replacement || al->meops[i].eop == Match 
        || al->meops[i].eop == Mismatch) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = al->v[j+q+al->voff];
        //real Replacement?
        if (!matchIUPAC(utemp[j+r], vtemp[j+r]))
          comp[j+r] = ' ';
        else
          comp[j+r] = '|';
      }
      //set string ptrs
      p+=j;
      q+=j;
      r+=j;
    }
    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = '-';
        vtemp[j+r] = al->v[j+q+al->voff];
        comp[j+r] = ' ';
      }
      //set ptrs
      r+=j;
      q+=j;
    }
    //if insertions occured
    if (al->meops[i].eop == Insertion) {
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = '-';
        comp[j+r] = ' ';
      }
      r+=j;
      p+=j;
    }
    
    if (al->meops[i].eop == Softclip) {
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = '^';
        comp[j+r] = ' ';
      }
      r+=j;
      p+=j;
    }

    //if skipping occured
    if (al->meops[i].eop == Skipped) {
    //  for (j=0; j < al->meops[i].steps; j++) {
         j=0;
         utemp[j+r] = '*';
         vtemp[j+r] = '*';
         comp[j+r]  = '*';
    //  }
      r+=1;
//      p+=j; no moves!
    }
    
    if(i == al->numofmeops-1) {
      //terminate strings
      utemp[r]='\0';
      vtemp[r]='\0';
      comp[r] ='\0';
      
      nlines = r/60;
      nlines += (r % 60) ? 1 : 0;
      //dump strings
      for(k=0; k < nlines; k++) {
        len = (k*60 > r) ? r % 60 : 60;
        bl_bsprintf(&mystring, "%.*s%c", len, &utemp[k*60], lf);
        bl_bsprintf(&mystring, "%.*s%c", len, &comp[k*60], lf);
        bl_bsprintf(&mystring, "%.*s%c", len, &vtemp[k*60], lf);
      }
      bl_bsprintf(&mystring, "%c", lf);
      memset(utemp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(comp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(vtemp, 0, sizeof(char)*(al->ulen+al->vlen));
    }
  }

  free(utemp);
  free(comp);
  free(vtemp);

  return mystring;
}


//dumps visual representation of alignments and should be shorter!
void 
showAlign(Alignment* al, FILE *dev) {
  int i, j , k, nlines, len;

  Uint p = 0, q = 0, r = 0;
  //output strings
  char* utemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* vtemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* comp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));

  //iter over all multieops
  for (i=0; i < al->numofmeops; i++) {
    //if Replacement occured
    if (al->meops[i].eop == Replacement || al->meops[i].eop == Match 
        || al->meops[i].eop == Mismatch) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = al->v[j+q+al->voff];
        //real Replacement?
        if (!matchIUPAC(utemp[j+r], vtemp[j+r]))
          comp[j+r] = ' ';
        else
          comp[j+r] = '|';
      }
      //set string ptrs
      p+=j;
      q+=j;
      r+=j;
    }
    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = '-';
        vtemp[j+r] = al->v[j+q+al->voff];
        comp[j+r] = ' ';
      }
      //set ptrs
      r+=j;
      q+=j;
    }
    //if insertions occured
    if (al->meops[i].eop == Insertion) {
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = '-';
        comp[j+r] = ' ';
      }
      r+=j;
      p+=j;
    }
    
    if (al->meops[i].eop == Softclip) {
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = '^';
        comp[j+r] = ' ';
      }
      r+=j;
      p+=j;
    }

    //if skipping occured
    if (al->meops[i].eop == Skipped) {
    //  for (j=0; j < al->meops[i].steps; j++) {
         j=0;
         utemp[j+r] = '*';
         vtemp[j+r] = '*';
         comp[j+r]  = '*';
    //  }
      r+=1;
//      p+=j; no moves!
    }
    
    if(i == al->numofmeops-1) {
      //terminate strings
      utemp[r]='\0';
      vtemp[r]='\0';
      comp[r] ='\0';
      
      nlines = r/60;
      nlines += (r % 60) ? 1 : 0;
      //dump strings
      for(k=0; k < nlines; k++) {
        len = (k*60 > r) ? r % 60 : 60;
        fprintf(dev, "%.*s\n", len, &utemp[k*60]);
        fprintf(dev, "%.*s\n", len, &comp[k*60]);
        fprintf(dev, "%.*s\n", len, &vtemp[k*60]);
      }
      fprintf(dev, "\n");
      memset(utemp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(comp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(vtemp, 0, sizeof(char)*(al->ulen+al->vlen));
    }
  }

  free(utemp);
  free(comp);
  free(vtemp);
}

//dumps visual representation of alignments and should be shorter!
void 
showAlignModel(Alignment* al, FILE *dev, gene_t *model) {
  int i, j , k, nlines, len;

  Uint p = 0, q = 0, r = 0;

  //output strings
  char* utemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* vtemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* comp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* pred = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  memset(pred, 0, sizeof(char)*(al->vlen+al->ulen));
  //iter over all multieops
  for (i=0; i < al->numofmeops; i++) {
    //if Replacement occured
    if (al->meops[i].eop == Replacement || al->meops[i].eop == Match 
        || al->meops[i].eop == Mismatch) {

      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {

        if(bl_isExon(model, j+q+al->voff)) {
            pred[j+r]='E';
        } else {
            pred[j+r]='I';
        }

        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = al->v[j+q+al->voff];
        //real Replacement?
        if (!matchIUPAC(utemp[j+r], vtemp[j+r]))
          comp[j+r] = ' ';
        else
          comp[j+r] = '|';
      }
      //set string ptrs
      p+=j;
      q+=j;
      r+=j;
    }
    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        if(bl_isExon(model, j+q+al->voff)) {
            pred[j+r]='E';
        } else {
            pred[j+r]='I';
        }
        utemp[j+r] = '-';
        vtemp[j+r] = al->v[j+q+al->voff];
        comp[j+r] = ' ';
      }
      //set ptrs
      r+=j;
      q+=j;
    }
    //if insertions occured
    if (al->meops[i].eop == Insertion) {
      for (j=0; j < al->meops[i].steps; j++) {

      
        if(bl_isExon(model, j+q+al->voff)) {
            pred[j+r]='E';
        } else {
            pred[j+r]='I';
        }

        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = '-';
        comp[j+r] = ' ';
      }
      r+=j;
      p+=j;
    }

   if (al->meops[i].eop == Softclip) {
      for (j=0; j < al->meops[i].steps; j++) {

      
        if(bl_isExon(model, j+q+al->voff)) {
            pred[j+r]='E';
        } else {
            pred[j+r]='I';
        }

        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = '^';
        comp[j+r] = ' ';
      }
      r+=j;
      p+=j;
    }


    //if skipping occured
    if (al->meops[i].eop == Skipped) {
    //  for (j=0; j < al->meops[i].steps; j++) {
         j=0;
         utemp[j+r] = '*';
         vtemp[j+r] = '*';
         comp[j+r]  = '*';
         pred[j+r] = '*';
    //  }
      r+=1;
//      p+=j; no moves!
    }
    
    if(i == al->numofmeops-1) {
      //terminate strings
      utemp[r]='\0';
      vtemp[r]='\0';
      comp[r] ='\0';
      pred[r] ='\0';
      
      nlines = r/60;
      nlines += (r % 60) ? 1 : 0;
      //dump strings
      for(k=0; k < nlines; k++) {
        len = (k*60 > r) ? r % 60 : 60;
        fprintf(dev, "%.*s\n", len, &utemp[k*60]);
        fprintf(dev, "%.*s\n", len, &comp[k*60]);
        fprintf(dev, "%.*s\n", len, &vtemp[k*60]);
        fprintf(dev, "%.*s\n", len, &pred[k*60]);
      }
      fprintf(dev, "\n");
      memset(utemp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(comp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(vtemp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(pred, 0, sizeof(char)*(al->ulen+al->vlen));
    }
  }

  free(utemp);
  free(comp);
  free(vtemp);
  free(pred);
}



//dumps visual representation of alignments and should be shorter!
void 
showAlignLF(Alignment* al, FILE *dev, char lf) {
  int i, j , k, nlines, len;
  Uint p = 0, q = 0, r = 0;

  if(dev == NULL) return;

  //output strings
  char* utemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* vtemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* comp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));

  //iter over all multieops
  for (i=0; i < al->numofmeops; i++) {
    //if Replacement occured
    if (al->meops[i].eop == Replacement
        || al->meops[i].eop == Match 
        || al->meops[i].eop == Mismatch) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = al->v[j+q+al->voff];
        //real Replacement?
        if (!matchIUPAC(utemp[j+r], vtemp[j+r]))
          comp[j+r] = ' ';
        else
          comp[j+r] = '|';
      }
      //set string ptrs
      p+=j;
      q+=j;
      r+=j;
    }
    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = '-';
        vtemp[j+r] = al->v[j+q+al->voff];
        comp[j+r] = ' ';
      }
      //set ptrs
      r+=j;
      q+=j;
    }
    //if insertions occured
    if (al->meops[i].eop == Insertion) {
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = '-';
        comp[j+r] = ' ';
      }
      r+=j;
      p+=j;
    }

    if (al->meops[i].eop == Softclip) {
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = '^';
        comp[j+r] = ' ';
      }
      r+=j;
      p+=j;
    }


    //if skipping occured
    if (al->meops[i].eop == Skipped) {
    //  for (j=0; j < al->meops[i].steps; j++) {
         j=0;
         utemp[j+r] = '*';
         vtemp[j+r] = '*';
         comp[j+r]  = '*';
    //  }
      r+=1;
//      p+=j; no moves!
    }
    
    if(i == al->numofmeops-1) {
      //terminate strings
      utemp[r]='\0';
      vtemp[r]='\0';
      comp[r] ='\0';
      
      nlines = r/60;
      nlines += (r % 60) ? 1 : 0;
      //dump strings
      for(k=0; k < nlines; k++) {
        len = (k*60 > r) ? r % 60 : 60;
        fprintf(dev, "%.*s%c", len, &utemp[k*60], lf);
        fprintf(dev, "%.*s%c", len, &comp[k*60],  lf);
        fprintf(dev, "%.*s%c", len, &vtemp[k*60], lf);
      }
      fprintf(dev, "%c", lf);
      memset(utemp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(comp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(vtemp, 0, sizeof(char)*(al->ulen+al->vlen));
    }
  }

  free(utemp);
  free(comp);
  free(vtemp);
}


//dumps visual representation of alignments and should be shorter!
/*void 
showAlignLF(Alignment* al, FILE *dev, char lf) {
  int i, j , k, nlines, len;

  Uint p = 0, q = 0, r = 0;
  //output strings
  char* utemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* vtemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* comp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));

  //iter over all multieops
  for (i=0; i < al->numofmeops; i++) {
    //if Replacement occured
    if (al->meops[i].eop == Replacement) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = al->v[j+q+al->voff];
        //real Replacement?
        if (!matchIUPAC(utemp[j+r], vtemp[j+r]))
          comp[j+r] = ' ';
        else
	     comp[j+r] = '|';
      }
      //set string ptrs
      p+=j;
      q+=j;
      r+=j;
    }
    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = '-';
        vtemp[j+r] = al->v[j+q+al->voff];
        comp[j+r] = ' ';
      }
      //set ptrs
      r+=j;
      q+=j;
    }
    //if insertions occured
    if (al->meops[i].eop == Insertion) {
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = '-';
        comp[j+r] = ' ';
      }
      r+=j;
      p+=j;
    }
    //if skipping occured
    if (al->meops[i].eop == Skipped) {
    //  for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = '*';
        vtemp[j+r] = '*';
        comp[j+r] = '*';
    //  }
      r+=1;
//      p+=j; no moves!
    }


    if(i == al->numofmeops-1) {
      //terminate strings
      utemp[r]='\0';
      vtemp[r]='\0';
      comp[r] ='\0';
      
      nlines = r/60;
      nlines += (r % 60) ? 1 : 0;
      //dump strings
      for(k=0; k < nlines; k++) {
        len = (k*60 > r) ? r % 60 : 60;
        fprintf(dev, "%.*s%c", len, &utemp[k*60], lf);
        fprintf(dev, "%.*s%c", len, &comp[k*60], lf);
        if(k < nlines-1)
        fprintf(dev, "%.*s%c", len, &vtemp[k*60], lf);
        else
        fprintf(dev, "%.*s", len, &vtemp[k*60]);
      }
    
      memset(utemp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(comp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(vtemp, 0, sizeof(char)*(al->ulen+al->vlen));
    }
  }

  free(utemp);
  free(comp);
  free(vtemp);
}
*/

void
insertEop(Alignment *al, Eoptype eop) {
  Uint pos = 0;
  //if previous multieops have been set up
  if (al->numofmeops)  {
    pos = al->numofmeops-1;
    //inc steps if curr eop matches prev eops
    if (al->meops[pos].eop == eop) {
      al->meops[pos].steps++;
      //set up new multieop otherwise
    } else {
      al->numofmeops++;
      pos = al->numofmeops-1;
      al->meops[pos].eop =  eop;
      al->meops[pos].steps = 1;
    }
    //set up first multieop
  } else {
    al->numofmeops = 1;
    al->meops[0].eop = eop;
    al->meops[0].steps = 1;
  }
}


Alignment*
reevalAlignment(Alignment *al) {
  Alignment *new;
  Uint i, j, p=0, q=0;

  new = ALLOCMEMORY(space, NULL, Alignment, 1);
  initAlignment(new, al->u, al->ulen, al->uoff, al->v, al->vlen, al->voff);
  
  for(i=0; i < al->numofmeops; i++) {
    
    switch(al->meops[i].eop) {
      case Replacement:
      case Match:
      case Mismatch:
        for(j=0; j < al->meops[i].steps; j++) {
          if(matchIUPAC(al->u[p+al->uoff+j], al->v[q+al->voff+j])) {
            insertEop(new, Match);
          } else { 
            insertEop(new, Mismatch);
          }
        }
        p+=al->meops[i].steps;
        q+=al->meops[i].steps;
        break;
      case Deletion:
        insertMeop(new, Deletion, al->meops[i].steps);
        q+=al->meops[i].steps;
        break;
      case Insertion:
        insertMeop(new, Insertion, al->meops[i].steps);
        p+=al->meops[i].steps;
        break;
      case Softclip:
        insertMeop(new, Softclip, al->meops[i].steps);
        p+=al->meops[i].steps;
        break;
      case Skipped:
        insertMeop(new, al->meops[i].eop, al->meops[i].steps);
        break;
    }
  }
  return new;
}

void
insertMeop(Alignment *al, Eoptype eop, Uint steps) {
  
  if (al->numofmeops > 0)  { 
    //inc steps if curr eop matches prev eops
    if (al->meops[al->numofmeops-1].eop == eop) {
      al->meops[al->numofmeops-1].steps+=steps;
      //set up new multieop otherwise
    } else {
      al->numofmeops++;
      al->meops[al->numofmeops-1].eop =  eop;
      al->meops[al->numofmeops-1].steps = steps;
    }
    //set up first multieop
  } else {
    al->numofmeops = 1;
    al->meops[0].eop = eop;
    al->meops[0].steps = steps;
  }

  return;
}

void
revMeops(Alignment *al) {
  Uint start = 0;
  Uint end = al->numofmeops-1;
  Multieop *meops = al->meops; 

  if (al->numofmeops == 0) return;

  while (start<end) {
    meops[start].eop ^= meops[end].eop;
    meops[start].steps ^= meops[end].steps;
    meops[end].eop ^= meops[start].eop;
    meops[end].steps ^= meops[start].steps;
    meops[start].eop ^= meops[end].eop;
    meops[start].steps ^= meops[end].steps;

    start++;
    end--;
  } 
}

Uint
getValignlenAndSkipped(Alignment *al) {

 Uint i, vallen=0, steps;
  Eoptype eop;

  for(i=0; i < al->numofmeops; i++) {
    eop = al->meops[i].eop;
    steps = al->meops[i].steps;
    switch(eop) {
      case Replacement:
      case Mismatch:
      case Match:
        vallen += steps;
        break;
      case Deletion:
        vallen += steps;
        break;
      case Insertion:
        break;
     case Softclip:
        break;
      case Skipped: //no acutal moves - but on the genome it is
        vallen += steps;
        break;
    }
  }

  return vallen;

}

Uint
getValignlen(Alignment *al) {
  Uint i, vallen=0, steps;
  Eoptype eop;

  for(i=0; i < al->numofmeops; i++) {
    eop = al->meops[i].eop;
    steps = al->meops[i].steps;
    switch(eop) {
      case Replacement:
      case Match:
      case Mismatch:
        vallen += steps;
        break;
      case Deletion:
        vallen += steps;
        break;
      case Insertion:
        break;
     case Softclip:
        break;
      case Skipped: //no acutal moves: skipped are masked from reference
        break;
    }
  }

  return vallen;
}

Uint
getUalignlen(Alignment *al) {
  Uint i, uallen=0, steps;
  Eoptype eop;

  for(i=0; i < al->numofmeops; i++) {
    eop = al->meops[i].eop;
    steps = al->meops[i].steps;
    switch(eop) {
      case Replacement:
      case Match:
      case Mismatch:
        uallen += steps;
        break;
      case Deletion:
        break;
      case Insertion:
        uallen += steps;
        break;
     case Softclip:
        uallen += steps;
        break;
      case Skipped: // no actual moves: skipped are masked from reference
        break;
    }
  }

  return uallen;
}

int
getUalignlenNoClip(Alignment *al) {
  Uint i, uallen=0, steps;
  Eoptype eop;

  for(i=0; i < al->numofmeops; i++) {
    eop = al->meops[i].eop;
    steps = al->meops[i].steps;
    switch(eop) {
      case Replacement:
      case Match:
      case Mismatch:
        uallen += steps;
        break;
      case Deletion:
        break;
      case Insertion:
        uallen += steps;
        break;
     case Softclip:
        break;
      case Skipped: // no actual moves: skipped are masked from reference
        break;
    }
  }

  return uallen;
}



/*------------------------------- getEopString -------------------------------
 *    
 * @brief 
 * @author Steve Hoffmann 
 *   
 */
 
char*
getEopString (Alignment *al)
{
  Uint i, k, len=0;
  char *string = NULL;


  for(i=0; i < al->numofmeops; i++) {
    for(k=0; k < al->meops[i].steps; k++) { 
      string = ALLOCMEMORY(NULL, string, char, len+3);
      len += snprintf(&string[len], 2, "%c", decodeEop[al->meops[i].eop]);
    }
  }
	
  return string;
}

char predictStrand(Alignment *al, char *vseq) {

  Uint i, vallen=0, steps;
  Eoptype eop;
  char skippedalign = 0;
  char strand  =0, *v;

  for(i=0; i < al->numofmeops; i++) {
    eop = al->meops[i].eop;
    steps = al->meops[i].steps;
    switch(eop) {
      case Replacement:
      case Mismatch:
      case Match:
        vallen += steps;
        break;
      case Deletion:
        vallen += steps;
        break;
      case Insertion:
        break;
     case Softclip:
        break;
     case Skipped: //no acutal moves - but on the genome it is
        skippedalign = 1;
        if(steps >= 2) {
          v = &vseq[al->voff+vallen];
          if( v[0] == 'G' && 
              v[1] == 'T' &&
              v[steps-2] == 'A' &&
              v[steps-1] == 'G') { 
            strand =1;
          }
 
          if( v[0] == 'C' && 
              v[1] == 'T' &&
              v[steps-2] == 'A' &&
              v[steps-1] == 'C') { 
            strand =2;
          }

          if( v[0] == 'G' && 
              v[1] == 'C' &&
              v[steps-2] == 'A' &&
              v[steps-1] == 'G') { 
            strand =1;
          }
 
          if( v[0] == 'C' && 
              v[1] == 'T' &&
              v[steps-2] == 'G' &&
              v[steps-1] == 'C') { 
            strand =2;
          }
        
          if( v[0] == 'A' && 
              v[1] == 'T' &&
              v[steps-2] == 'A' &&
              v[steps-1] == 'C') { 
            strand =1;
          }
 
          if( v[0] == 'G' && 
              v[1] == 'T' &&
              v[steps-2] == 'A' &&
              v[steps-1] == 'T') { 
            strand =2;
          }
        }

        vallen += steps;
        break;
    }
  }

  if(!strand && skippedalign) {
    strand = 3;
  }

  return strand;
}


Uint*
getUPartialAlignlen(Alignment *al, Uint *noofparts) {
  Uint i, steps;
  Eoptype eop;
  Uint *uallen, n=1;

  uallen = ALLOCMEMORY(NULL, NULL, Uint, n);
  uallen[n-1] = 0;

  for(i=0; i < al->numofmeops; i++) {
    eop = al->meops[i].eop;
    steps = al->meops[i].steps;
    switch(eop) {
      case Replacement:
      case Match:
      case Mismatch:
        uallen[n-1] += steps;
        break;
      case Deletion:
        break;
      case Insertion:
        uallen[n-1] += steps;
        break;
     case Softclip:
        uallen[n-1] += steps;
        break;
      case Skipped: // no actual moves: skipped are masked from reference
        uallen = ALLOCMEMORY(NULL, uallen, Uint, n+1);
        n++;
        uallen[n-1] = 0;
        break;
    }
  }

  *noofparts = n;
  return uallen;
}

Uint
getPartialAlignNumber(Alignment *al) {
  Uint i;
  Eoptype eop;
  Uint n=1;


  for(i=0; i < al->numofmeops; i++) {
    eop = al->meops[i].eop;
    switch(eop) {
      case Replacement:
      case Match:
      case Mismatch:
        break;
      case Deletion:
        break;
      case Insertion:
        break;
     case Softclip:  
        break;
      case Skipped: // no actual moves: skipped are masked from reference
        n++;
        break;
    }
  }

  return n;
}


/*------------------------- bl_getPartialAlignments --------------------------
 *    
 * @brief extract the partial alignments (separated by N)
 * @author Steve Hoffmann 
 *   
 */
 
Alignment*
bl_getPartialAlignments (Alignment *al, char *vseq, Uint *nalign)
{
 
  Uint i, n, u, v, k=0;
  Alignment* new=NULL;
  Eoptype eop;
  Uint vlen; 
  char *ref;

  n = getPartialAlignNumber(al);
  vlen = getValignlenAndSkipped(al);
  u = al->uoff;
  v = al->voff;

  ref = &vseq[v];
  new = ALLOCMEMORY(NULL, NULL, Alignment, n);  
  initAlignment(&new[0], al->u, al->ulen, u, ref, vlen, v);

  for(i=0; i < al->numofmeops; i++) {
    eop = al->meops[i].eop;
    switch(eop) {
      case Replacement:
        insertMeop(&new[k], Replacement, al->meops[i].steps);
        u+=al->meops[i].steps;
        v+=al->meops[i].steps;
        break;
      case Match:
        insertMeop(&new[k], Match, al->meops[i].steps);
        u+=al->meops[i].steps;
        v+=al->meops[i].steps;
        break;
      case Mismatch:
        insertMeop(&new[k], Mismatch, al->meops[i].steps);
        u+=al->meops[i].steps;
        v+=al->meops[i].steps;
        break;
      case Deletion:         
        insertMeop(&new[k], Deletion, al->meops[i].steps);
        v+=al->meops[i].steps;
        break;
      case Insertion:
        insertMeop(&new[k], Insertion, al->meops[i].steps);
        u+=al->meops[i].steps;
        break;
     case Softclip:  
        insertMeop(&new[k], Insertion, al->meops[i].steps);
        u+=al->meops[i].steps;
        break;
      case Skipped: // no actual moves: skipped are masked from reference
        k++; 
        v+=al->meops[i].steps;
        initAlignment(&new[k], al->u, al->ulen, u, ref, vlen, v);
        break;
    }
  }

  *nalign =k+1;
  return new;
}


/*------------------------------- bl_get5primeU -------------------------------
 *    
 * @brief get 5'-end from alignment on query
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_alignGet5PrimeU (Alignment *al, char rc)
{  
  unsigned int end = al->uoff;
    
  if(rc == 1) end = al->ulen - al->uoff - getUalignlen(al);

  return end;
}




/*------------------------------- bl_get3primeU -------------------------------
 *    
 * @brief get 3'-end from alignment on query
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_alignGet3PrimeU (Alignment *al, char rc)
{  
  unsigned int end = al->uoff + getUalignlen(al) - 1;
    
  if(rc == 1) end = al->ulen - al->uoff - 1;

  return end;
}



/*------------------------------- bl_get5primeV -------------------------------
 *    
 * @brief get 5'-end from alignment on reference
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_alignGet5PrimeV (Alignment *al, char rc)
{
  unsigned int end = al->voff;

    if(rc == 1) end += getValignlen(al) - 1;

	return end;
}


/*--------------------------- bl_compareAlignments ---------------------------
 *    
 * @brief compare alignments
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_compareAlignments (Alignment *a, Alignment *b)
{
  Uint i;

  assert(a->voff == b->voff);
  assert(a->uoff == b->uoff);
	
  assert(a->numofmeops == b->numofmeops);
  for(i=0; i < a->numofmeops; i++) {
    assert(a->meops[i].steps == b->meops[i].steps);
    assert(a->meops[i].eop == b->meops[i].eop);
  }

  return ;
}
