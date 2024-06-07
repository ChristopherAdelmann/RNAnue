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
 *  bitvectoralg.c
 *  implementation of Gene Myers 
 *  bitvector algorithm
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 05/23/2008 06:12:54 PM CEST
 *  
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "basic-types.h"
#include "mathematics.h"
#include "bitVector.h"
#include "memory.h"
#include "alignment.h"
#include "iupac.h"


/*-------------------------------- alphacheck ---------------------------------
 *    
 * @brief alphabet check
 * @author Steve Hoffmann 
 *   
 */

unsigned char 
alphacheck (char c) {
  if (   c == 'A' || c == 'a' || c == 'C' || c == 'c' 
      || c == 'T' || c == 't' || c == 'G' || c == 'g') {
    return 1;
  }
  return 0;
}


/*---------------------------- getstringalphabet -----------------------------
 *    
 * @brief return the alphabet of a string
 * @author Steve Hoffmann 
 *   
 */
 
char*
getstringalphabet (void *space, char *string, Uint len, Uint *asize)
{
  char *alphabet=NULL;
  unsigned char *found;
  Uint i, l = 0;

  found = ALLOCMEMORY(space, NULL, char, 256);
  memset(found, 0, 256);

  for(i=0; i < len; i++) {
    if(!found[(Uint)string[i]]) {
      alphabet = ALLOCMEMORY(space, alphabet, char, l+2);
      alphabet[l] = string[i];
      alphabet[l+1] = 0;
      found[(Uint)string[i]] = 1;
      l++;
    }
  }

  FREEMEMORY(space, found);
  *asize = l;
  return alphabet;
}



/*-------------------------------- encodetab ---------------------------------
 *    
 * @brief brutal ascii encoding
 * @author Steve Hoffmann 
 *   
 */

Uint*
encodetab(char *alphabet, Uint asize) {
    Uint i;
    Uint *tab;

    tab = ALLOCMEMORY(space, NULL, Uint, 255);
    memset(tab, asize, sizeof(Uint)*255);
    for(i=0; i < asize; i++) {
        tab[(Uint)alphabet[i]] = i;
    }
    return tab;
} 


/*---------------------------------- getpeq ----------------------------------
 *    
 * @brief returns pattern mask for each char in alphabet
 * @author Steve Hoffmann 
 *   
 */
 
bitvector*
getpeq(void *space,
    char *query, 
    Uint qlen,
    char *alphabet,
    Uint asize,
    Uint *enctab) {

  bitvector* peq;
  Uint i,j, wordno;

  wordno = qlen/BITVECTOR_WORDSIZE;
  wordno++;
  peq = ALLOCMEMORY(space, NULL, bitvector*, asize);

//  printf("query: %s\n", query);
  for(i=0; i < asize; i++) {
    peq[i] = initbitvector(space, BITVECTOR_WORDSIZE*wordno);
    setbitvector(peq[i], BITVECTOR_WORDSIZE*wordno, 0);
    //peq[i] = 0;
    for(j=0; j < qlen; j++) {
      if (matchIUPAC(query[j], alphabet[i])){
	bitvector_setbit(peq[i], j, 1);	  
      }
    }
  //  printf("char:%c\n", alphabet[i]);
  //  dumpbitvector(peq[i], BITVECTOR_WORDSIZE*wordno);
  }

  return peq;
}


/*------------------------------ myersbitvector ------------------------------
 *    
 * @brief approx. string matching: calculate min edist of query and subject
 * @author Steve Hoffmann 
 *   
 */
 
PairSint
myersbitvector(
    void *space,
    char *query, 
    Uint qlen, 
    char *subject, 
    Uint slen, 
    char *alphabet, 
    Uint asize,
    Uint *enctab,
    Uint k,
    bitvector *peq) {

  bitvector 
  Pv,
  Mv,
  Eq;
  PairSint res;

  Uint score=qlen,
       i,
       j,
       wordno,
       bits;
  bitvector_t check,
              temp,
              carryxh,
              carryph,
              carrymh,
              Ph=0,
              Mh=0,
              Xv=0,
              Xh=0;

  res.a = -1;
  res.b = qlen;
  wordno = qlen/BITVECTOR_WORDSIZE;
  bits   = qlen & (BITVECTOR_WORDSIZE-1);
  wordno++;

  Pv  = initbitvector(space, wordno*BITVECTOR_WORDSIZE);
  Mv  = initbitvector(space, wordno*BITVECTOR_WORDSIZE);

  check = 0;
  bitvector_setbit(&check, bits, 1);
 
  setbitvector(Pv, wordno*BITVECTOR_WORDSIZE, 1);
  setbitvector(Mv, wordno*BITVECTOR_WORDSIZE, 0);


  for(i=0; i < slen; i++) {
   
    Eq = peq[enctab[(Uint)subject[i]]];
    carryxh = carryph = carrymh = 0;

    for(j=0; j < wordno; j++) {
     
      Xv = Eq[j] | Mv[j];
      temp = ((Eq[j] & Pv[j]) + Pv[j] + carryxh); 
      Xh = (temp ^ Pv[j]) | Eq[j];
      
      if (carryxh)
        carryxh = (temp <= (Eq[j] & Pv[j]) || temp <= Pv[j]);
      else
        carryxh = (temp < (Eq[j] & Pv[j]) || temp < Pv[j]);

      Ph = Mv[j] | ~(Xh | Pv[j]);
      Mh = Pv[j] & Xh;

      //check if last word
      if (j == wordno-1) {
        if (Ph & check) 
          score+=1;
        else if(Mh & check) 
          score-=1;
      }

      /*Ph = Ph << 1; with carry*/
      temp = (Ph << 1) | carryph;
      carryph = Ph >> (BITVECTOR_WORDSIZE-1);
      Ph = temp; 

      temp = (Mh << 1) | carrymh;
      carrymh = Mh >> (BITVECTOR_WORDSIZE-1);
      Mh = temp;

      Pv[j] = Mh | ~(Xv | Ph);
      Mv[j] = Ph & Xv;

    }

    if (score <= k && score <= res.b) {
        res.a = i;
        res.b = score;
    } 
  }
 
  FREEMEMORY(space, Pv);
  FREEMEMORY(space, Mv);       

  return res;
}


/*------------------------------ myersbitmatrix ------------------------------
 *    
 * @brief modified bitvector algorithm to return bitmatrix for backtracking
 * @author Steve Hoffmann 
 *   
 */
 
bitvector*
myersbitmatrix(
    void *space,
    char *query, 
    Uint qlen, 
    char *subject, 
    Uint slen, 
    char *alphabet, 
    Uint asize,
    Uint *enctab,
    Uint k,
    bitvector *peq,
    PairSint *res,
    bitvector *D,
    Uint dim) {

  bitvector 
  *Pv,
  *Mv,
  *myMh,
  MvP,
  PvP,
  Eq;

  Uint score=qlen,
       i,
       j,
       wordno,
       bits;
  bitvector_t check,
              temp,
              carryxh,
              carryph,
              carrymh,
              Ph=0,
              Mh=0,
              Xv=0,
              Xh=0;

  res->a = -1;
  res->b = qlen;
  wordno = qlen/BITVECTOR_WORDSIZE;
  bits   = qlen & (BITVECTOR_WORDSIZE-1);
  wordno++;

  Pv = D; 
  Mv = &Pv[dim+1];  
  myMh = &D[2*(dim+1)+1];

 
  memset(Pv[0], 255, wordno*(sizeof(bitvector_t)));
  memset(Mv[0], 0, wordno*(sizeof(bitvector_t)));
  memset(myMh[0], 255, wordno*(sizeof(bitvector_t)));

  check = 0;
  bitvector_setbit(&check, bits, 1);

 
  for(i=0; i < slen; i++) {
   
    Eq = peq[enctab[(Uint)subject[i]]];
    carryxh = carryph = carrymh = 0;

    MvP = Mv[i];
    PvP = Pv[i];

    for(j=0; j < wordno; j++) {
     
      Xv = Eq[j] | MvP[j];
      temp = ((Eq[j] & PvP[j]) + PvP[j] + carryxh); 
      Xh = (temp ^ PvP[j]) | Eq[j];
      
      if (carryxh)
        carryxh = (temp <= (Eq[j] & PvP[j]) || temp <= PvP[j]);
      else
        carryxh = (temp < (Eq[j] & PvP[j]) || temp < PvP[j]); 

      Ph = MvP[j] | ~(Xh | PvP[j]);
      Mh = PvP[j] & Xh;

      //check if last word
      if (j == wordno-1) {
        if (Ph & check) {
          score+=1;
       //   printf("%d,%d:%d hout: %d\n", i, j, score, 1);
        }
        else if(Mh & check) {
          score-=1;
        //  printf("%d,%d:%d hout: %d\n", i, j, score, -1);
        }
      }

      /*Ph = Ph << 1; with carry*/
      temp = (Ph << 1) | carryph;
      carryph = Ph >> (BITVECTOR_WORDSIZE-1);
      Ph = temp; 

      temp = (Mh << 1) | carrymh;
      carrymh = Mh >> (BITVECTOR_WORDSIZE-1);
      Mh = temp;

      Pv[i+1][j] = Mh | ~(Xv | Ph);
      Mv[i+1][j] = Ph & Xv;
      myMh[i+1][j] = Xh;
    }

    if (score <= k && score <= res->b) { // && i < slen - 1) {
        res->a = i;
        res->b = score;
      // fprintf(stderr, "%d: %d wordno:%d\n", i, score, wordno);
    } 
  }

  return Pv;
}




/*---------------------------- myersbitblockmatrix ----------------------------
 *    
 * @brief modified bitvector algorithm to return bitmatrix for backtracking
 * @author Steve Hoffmann 
 *   
 */
 
bitvector*
myersblockbitmatrix(
    void *space,
    char *query, 
    Uint qlen, 
    char *subject, 
    Uint slen, 
    char *alphabet, 
    Uint asize,
    Uint *enctab,
    Uint k,
    bitvector *peq,
    PairSint *res,
    bitvector *D,
    Uint dim) {

  bitvector 
  *Pv,
  *Mv,
  MvP,
  PvP,
  Eq,
  W;

  Uint *score,
       i,
       j,
       y,
       wordno,
       bits;
  int hout=0;
  bitvector_t check,
              last,
              first,
              temp,
              carryxh,
              carryph,
              carrymh,
              Ph=0,
              Mh=0,
              Xv=0,
              Xh=0,
              w;

  res->a = -1;
  res->b = qlen;
  wordno = qlen/BITVECTOR_WORDSIZE;
  bits   = qlen & (BITVECTOR_WORDSIZE-1);
  wordno++;
  
  Pv = D; 
  Mv = &Pv[dim+1];
 
  memset(Pv[0], 255, wordno*(sizeof(bitvector_t)));
  memset(Mv[0], 0, wordno*(sizeof(bitvector_t)));

  check = 0;
  bitvector_setbit(&check, bits, 1);
  last = 0;
  bitvector_setbit(&last, BITVECTOR_WORDSIZE-1, 1);
  first = 0;
  bitvector_setbit(&first, 0, 1);

  score = calloc(wordno+1, sizeof(Uint));
  W = ALLOCMEMORY(space, NULL, bitvector_t, wordno);
  y = MIN((Uint) ceil((double)k/(double)BITVECTOR_WORDSIZE)-1, (wordno-1));

  for(i=0; i <= y ; i++) {
    score[i] = (i+1)*BITVECTOR_WORDSIZE;
    W[i] = last; 
  }
  score[wordno-1] = qlen;
  W[wordno-1] = check; 

  for(i=0; i < slen; i++) {
   
    Eq = peq[enctab[(Uint)subject[i]]];
    carryxh = carryph = carrymh = 0;

    MvP = Mv[i];
    PvP = Pv[i];

    for(j=0; j <= y; j++) {

      Xv = Eq[j] | MvP[j];
      temp = ((Eq[j] & PvP[j]) + PvP[j] + carryxh); 
      Xh = (temp ^ PvP[j]) | Eq[j];

      if (carryxh)
        carryxh = (temp <= (Eq[j] & PvP[j]) || temp <= PvP[j]);
      else
        carryxh = (temp < (Eq[j] & PvP[j]) || temp < PvP[j]); 

      Ph = MvP[j] | ~(Xh | PvP[j]);
      Mh = PvP[j] & Xh;

      w = W[j];

      hout = 0;
      if (Ph & w) { 
        score[j] += 1;
        hout = 1;
      } else if (Mh & w) {
        score[j] -= 1; 
        hout = -1;
      } 
      
      /*Ph = Ph << 1; with carry*/
      temp = (Ph << 1) | carryph;
      carryph = Ph >> (BITVECTOR_WORDSIZE-1);
      Ph = temp; 

      temp = (Mh << 1) | carrymh;
      carrymh = Mh >> (BITVECTOR_WORDSIZE-1);
      Mh = temp;

      Pv[i+1][j] = Mh | ~(Xv | Ph);
      Mv[i+1][j] = Ph & Xv;
    }

  
    if (y < wordno-1 && score[y]-hout <= k 
        && ((Eq[y+1] & first) || hout < 0)) {
 
      y += 1;

      memset(&Pv[i][j], 255,(sizeof(bitvector_t)));
      memset(&Mv[i][j], 0,  (sizeof(bitvector_t)));
  
      MvP = Mv[i];
      PvP = Pv[i];

      //since we open a new zone here we have to update the score
      Xv = Eq[j] | MvP[j];
      temp = ((Eq[j] & PvP[j]) + PvP[j] + carryxh); 
      Xh = (temp ^ PvP[j]) | Eq[j];

      Ph = MvP[j] | ~(Xh | PvP[j]);
      Mh = PvP[j] & Xh;

      //check if last word
      if (j < wordno-1) {
        score[j] = score[j-1] + BITVECTOR_WORDSIZE - hout;
        W[j] = last; 
        w = last;
      } else { 
        score[j] = score[j-1] + bits - hout;
        w = check;
      }

      if (Ph & w) { 
        score[j] += 1;
      } else if (Mh & w) {
        score[j] -= 1; 
      }

      temp = (Ph << 1) | carryph;
      carryph = Ph >> (BITVECTOR_WORDSIZE-1);
      Ph = temp; 

      temp = (Mh << 1) | carrymh;
      carrymh = Mh >> (BITVECTOR_WORDSIZE-1);
      Mh = temp; 

      Pv[i+1][j] = Mh | ~(Xv | Ph);
      Mv[i+1][j] = Ph & Xv;

    } else  {

      while(y > 0 && score[y] >= k + BITVECTOR_WORDSIZE) {
        y -= 1;
      }
    } 

    if (score[wordno-1] <= k && score[wordno-1] <= res->b && i < slen - 1) {
      res->a = i;
      res->b = score[wordno-1];
     // fprintf(stderr, "block score: %d:%d (%d); wordno:%d\n", i, score[wordno-1], score[0], wordno-1);
    } 
  }

  FREEMEMORY(space, W);
  FREEMEMORY(space, score);
  return Pv;
}


/*---------------------------- bitvectorbacktrack ----------------------------
 *    
 * @brief backtracking in bitmatrix
 * @author Steve Hoffmann 
 *   
 */
 
Alignment*
bitvectorbacktrack(Alignment *al, bitvector *D, Uint dim, Uint k, Uint l, 
    char *subject, Uint *enctab, bitvector *peq) { //patch
  int i=k-1, 
  j=l;
  bitvector *Pv = D;
  bitvector *Mv = &D[dim+1];
  bitvector *Xh = &D[2*(dim+1)+1];
//  Uint match1=0, match2=0, match3=0;
//  Uint mismatch1=0, mismatch2=0, mismatch3=0;

  while (i > 0 && j > 0) {
    if (bitvector_getbit(Pv[j], i)) {
   /*   fprintf(stderr, "[%d,%d]: Xh: %d, %d, %d, %d, Pv[j]: %d, Mv[j-1]:%d, equal: %d (%c)\n",            i,   j, 
          bitvector_getbit(Xh[j], i-1), bitvector_getbit(Xh[j-1], i-1), 
          bitvector_getbit(Xh[j-1], i), bitvector_getbit(Xh[j], i),
          bitvector_getbit(Pv[j], i), bitvector_getbit(Pv[j-1], i),     
          bitvector_getbit(peq[enctab[(Uint)subject[j-1]]]    , i)  , subject[j-1]); //patch
     */ 
      if(!bitvector_getbit(Xh[j], i) && !bitvector_getbit(peq[enctab[(Uint)subject[j]]], i)) { //patch2
      //if(bitvector_getbit(Pv[j-1], i) || bitvector_getbit(peq[enctab[(Uint)subject[j-1]]], i)) { //patch
        

      //  fprintf(stderr, "[%d,%d]: following bitvector: %d, equal: %d (%c)\n",i+2, j+2, bitvector_getbit(Pv[j+1],  i+2),  bitvector_getbit(peq[enctab[(Uint)subject[j+1]]], i+2), subject[j+1]); //patch
      //  fprintf(stderr, "[%d,%d]: following bitvector: %d, equal: %d (%c)\n",i+1, j+1, bitvector_getbit(Pv[j],  i+1),    bitvector_getbit(peq[enctab[(Uint)subject[j]]]  , i+1), subject[j]); //patch
      //  fprintf(stderr, "[%d,%d]: bitvector: %d, equal: %d (%c)\n",            i,   j, bitvector_getbit(Pv[j-1], i),     bitvector_getbit(peq[enctab[(Uint)subject[j-1]]]    , i)  , subject[j-1]); //patch
      //  fprintf(stderr, "[%d,%d]: previous bitvector: %d, equal: %d (%c)\n", i-1, j-1, bitvector_getbit(Pv[j-2], i-1),   bitvector_getbit(peq[enctab[(Uint)subject[j-2]]]  , i-1), subject[j-2]); //patch
      //  fprintf(stderr, "[%d,%d]: previous bitvector: %d, equal: %d (%c)\n", i-2, j-2, bitvector_getbit(Pv[j-3], i-2),   bitvector_getbit(peq[enctab[(Uint)subject[j-3]]]  , i-2), subject[j-3]); //patch
        
     //   if(bitvector_getbit(peq[enctab[(Uint)subject[j-1]]], i)) {
     //     match1++;
     //   } else {
     //     mismatch1++;
     //   }
        insertEop(al, Mismatch); //Replacement //patch
        i--; j--; //patch
      } else {  //patch
       insertEop(al, Insertion);
       i--;
      } //patch
    } else {
      if (bitvector_getbit(Mv[j-1], i) 
     //     && !bitvector_getbit(peq[enctab[(Uint)subject[j-1]]], i)   //patch
         ) { 
        insertEop(al, Deletion);
      } else {
        if(bitvector_getbit(peq[enctab[(Uint)subject[j-1]]], i)) {
//          match2++;
          insertEop(al, Match);
        } else {
          insertEop(al, Mismatch);
//          mismatch2++;
        }
       // insertEop(al, Replacement);
        i--;
      }
      j--;
    }
  }

  if (i==0 && j > 0) {
    if(bitvector_getbit(peq[enctab[(Uint)subject[j-1]]], i)) {
//      match3++;
      insertEop(al, Match);
    } else {
//      mismatch3++;
      insertEop(al, Mismatch);
    }
//    insertEop(al, Replacement); 
  } else {
   //  if(i>=0 && j==0 && bitvector_getbit(peq[enctab[(Uint)subject[0]]],i)) { //patch
   //    insertEop(al, Replacement); //patch
   //    i--; //patch
   //  } //patch

     while(i>=0) {
      insertEop(al, Insertion);
      i--;
    }
    /*no insertions in app. string matching at the end*/ 
  }

  /*adjust subject boundaries*/
  if(j>0) {
    al->voff = j-1;
    al->vlen -= (al->vlen > j) ?  j : al->vlen;
  }
  
  revMeops(al);
 // fprintf(stdout, "bitvector align match1:%d/%d match2:%d/%d match3:%d/%d\n",match1, mismatch1, match2, mismatch2, match3, mismatch3);
//  showAlign(al,stdout);
//  fprintf(stdout, "-----\n");

  return al;
}

/*------------------------------ wrapBitmatrix -------------------------------
 *    
 * @brief destruct bit matrix
 * @author Steve Hoffmann 
 *   
 */
 
void
wrapBitmatrix(void *space, bitvector *D, Uint m) {
  Uint i;
  for(i=0; i < m; i++) {
    FREEMEMORY(space, D[i]);
  }

  return;
}

