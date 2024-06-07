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
 *  nw.c
 *  needleman wunsch
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 11/15/2010 12:12:27 AM CET
 *  
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "basic-types.h"
#include "memory.h"
#include "mathematics.h"
#include "alignment.h"
#include <limits.h>
#include "sw.h"
#include "iupac.h"


/*--------------------------------- nwalign ----------------------------------
 *      
 *  needleman-wunsch global similarity alignment
 *  returns a matrix of size (m+1)*(n+1) where m is length of given sequence a
 *  and n the length of sequence b, respectively. Function expects
 *  a function to calculate a substitution score
 *     
 */

  int*
nwmatrix (void *space, symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo)
{
  int i, j, cols, rows, size;
  int *L;

  rows = m+1;
  cols = n+1;

  size = rows*cols;
  L = ALLOCMEMORY(space, NULL, int, size);
  L = memset(L, 0, sizeof(int)*size);

  for(i=1; i < m+1; i++) {
    MATRIX2D(L, cols, i, 0) = i*indel;
    for(j=1; j < n+1; j++) {
      MATRIX2D(L, cols, 0, j) = j*indel;

      MATRIX2D(L, cols, i, j) = 
        MAX3(
            MATRIX2D(L, cols, (i-1), j) + indel ,    
            MATRIX2D(L, cols, i, (j-1)) + indel , 
            MATRIX2D(L, cols, (i-1), (j-1)) + sub(a[i-1], b[j-1], nfo)
            );
    }
  }

  return L;
}


/*------------------------------- nwtraceback --------------------------------
 *      
 *  traceback to find optimal global alignment path
 *   
 */

  void
nwtraceback (void *space, int *M,  
    symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo, Alignment *al)
{
  Uint i, j, ncol, cur; 

  ncol = n+1;
  i = m;
  j = n;

  al->uoff = 0;
  al->voff = 0;

  while(i > 0 && j > 0) {

    cur = MATRIX2D(M, ncol, i, j);
 //   fprintf(stderr, "enter while cur = %d, %c-%c, %d\n", cur, a[i-1], b[j-1],
   //   sub(a[i-1], b[j-1], nfo) );
    if (MATRIX2D(M, ncol, i-1, j) + indel == cur){
      insertEop(al, Insertion);
      i--;   

  //        fprintf(stderr, "insertion\n");
    } else {
      if (MATRIX2D(M, ncol, i, j-1) + indel == cur) {
        insertEop(al, Deletion);
        j--;

    //      fprintf(stderr, "deletion\n");
      } else {
        if (MATRIX2D(M, ncol, i-1, j-1)+sub(a[i-1], b[j-1], nfo) 
            == cur){
          insertEop(al, Replacement);
          i--; j--;
      //    fprintf(stderr, "replacement\n");
        }
        else {
          assert(cur == 0);
        //  fprintf(stderr, "asserting.\n");
          al->uoff = i;
          al->voff = j;

          revMeops(al);
          return;
        }
      }
    }
  }

  al->uoff = i;
  al->voff = j;
  revMeops(al);

  return;
}

/*--------------------------------- sgalign ----------------------------------
 *      
 *  semi-global similarity alignment
 *  returns a matrix of size (m+1)*(n+1) where m is length of given sequence a
 *  and n the length of sequence b, respectively. Function expects
 *  a function to calculate a substitution score
 *     
 */

  int*
sgmatrix (void *space, symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo)
{
  int i, j, cols, rows, size;
  int *L;

  rows = m+1;
  cols = n+1;

  size = rows*cols;
  L = ALLOCMEMORY(space, NULL, int, size);
  L = memset(L, 0, sizeof(int)*size);

  for(i=1; i < m+1; i++) {
    MATRIX2D(L, cols, i, 0) = i*indel;
    for(j=1; j < n+1; j++) {
      MATRIX2D(L, cols, 0, j) = 0;

      MATRIX2D(L, cols, i, j) = 
        MAX3(
            MATRIX2D(L, cols, (i-1), j) + indel ,    
            MATRIX2D(L, cols, i, (j-1)) + indel , 
            MATRIX2D(L, cols, (i-1), (j-1)) + sub(a[i-1], b[j-1], nfo)
            );
    }
  }

  return L;
}

/*------------------------------- sgtraceback --------------------------------
 *      
 *  traceback to find optimal semi global alignment path
 *   
 */

  void
sgtraceback (void *space, int *M,  
    symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo, Alignment *al)
{
  Uint i, j, ncol, cur; 
  Uint maxcol;
  char explicitmatches = 1;
  ncol = n+1;
  i = m;
  maxcol = n;

  for(j=0; j < n; j++) {
    if(MATRIX2D(M, ncol, i, maxcol) < MATRIX2D(M, ncol, i, j))
      maxcol = j;
  }

  j = maxcol;

  al->uoff = 0;
  al->voff = 0;

  while(i > 0 && j > 0) {

    cur = MATRIX2D(M, ncol, i, j);
 //   fprintf(stderr, "enter while cur = %d, %c-%c, %d\n", cur, a[i-1], b[j-1],
   //   sub(a[i-1], b[j-1], nfo) );
    if (MATRIX2D(M, ncol, i-1, j) + indel == cur){
      insertEop(al, Insertion);
      i--;   

  //        fprintf(stderr, "insertion\n");
    } else {
      if (MATRIX2D(M, ncol, i, j-1) + indel == cur) {
        insertEop(al, Deletion);
        j--;

    //      fprintf(stderr, "deletion\n");
      } else {
        if (MATRIX2D(M, ncol, i-1, j-1)+sub(a[i-1], b[j-1], nfo) 
            == cur){
          if(explicitmatches) { 
            if (matchIUPAC(a[i-1], b[j-1])) {
              insertEop(al, Match);
            } else {
              insertEop(al, Mismatch);
            } 
          } else {
            insertEop(al, Replacement);
          }
          i--; j--;
          //    fprintf(stderr, "replacement\n");
        }
        else {
          assert(cur == 0);
        //  fprintf(stderr, "asserting.\n");
          al->uoff = i;
          al->voff = j;

          revMeops(al);
          return;
        }
      }
    }
  }

  al->uoff = i;
  al->voff = j;
  revMeops(al);

  return;
}


/*--------------------------------- sgaffine ---------------------------------
 *    
 * @brief affine gap cost semi global alignment
 * @author Steve Hoffmann 
 *   
 */
 
void
sgaffinematrix (void *space, int **A, int **B, int **S, symtype *a, Uint m, symtype *b,
    Uint n, int open, int ext, int close, Sint (*sub)(symtype, symtype, void*),
    void *nfo)
{
 
  int i, j, cols, rows, size;

  rows = m+1;
  cols = n+1;

  size = rows*cols;
  //A is gaps in a, B is gaps in b and S is substitutions
  *A = ALLOCMEMORY(space, NULL, int, size);
  *B = ALLOCMEMORY(space, NULL, int, size);
  *S = ALLOCMEMORY(space, NULL, int, size);
  
  *A = memset(*A, 0, sizeof(int)*size);
  *B = memset(*B, 0, sizeof(int)*size);
  *S = memset(*S, 0, sizeof(int)*size);

  //init for semi global gap affine alignment
  for(i=0; i < m+1; i++) {
    MATRIX2D(*S, cols, i, 0) = 0;
    MATRIX2D(*B, cols, i, 0) = INT_MIN;
  }

  for(j=0; j < n+1; i++) {
    MATRIX2D(*S, cols, 0, j) = 0;
    MATRIX2D(*B, cols, 0, j) = INT_MIN;
  }

  for(i=1; i < m+1; i++) {
    for(j=1; j < n+1; j++) {
      MATRIX2D(*A, cols, i, j) = 
        MAX(
            MATRIX2D(*A, cols, (i-1), j) + ext,    
            MATRIX2D(*S, cols, (i-1), j) + open + ext
            );

      MATRIX2D(*B, cols, i, j) = 
        MAX(
            MATRIX2D(*B, cols, i, (j-1)) + ext,    
            MATRIX2D(*S, cols, i, (j-1)) + open + ext 
            );

      MATRIX2D(*S, cols, i, j) = 
        MAX3(
            MATRIX2D(*S, cols, (i-1), (j-1)) + sub(a[i-1], b[j-1], nfo),    
            MATRIX2D(*A, cols, i, j),
            MATRIX2D(*B, cols, i, j)
            );
    }
  }


  return ;
}


/*---------------------------- sgaffinetraceback -----------------------------
 *    
 * @brief trace back the affine gap cost alignment
 * @author Steve Hoffmann 
 *   
 */

  void
sgaffinetraceback(void *space, int *A, int *B, int *S, symtype *a, Uint m, symtype *b,
    Uint n, int open, int ext, int close, Sint (*sub)(symtype, symtype, void*),
    void *nfo, Alignment *al)

{

  Uint i, j, ncol, cur; 
  Uint maxcol;
  char state = 'S';

  ncol = n+1;
  i = m;
  maxcol = n;

  for(j=0; j < n; j++) {
    if(MATRIX2D(S, ncol, i, maxcol) < MATRIX2D(S, ncol, i, j))
      maxcol = j;
  }

  j = maxcol;

  al->uoff = 0;
  al->voff = 0;


  while(i > 0 && j > 0) {
    if (state == 'S') { 
      cur = MATRIX2D(S, ncol, i, j);
      if (cur == MATRIX2D(S, ncol, (i-1), (j-1)) + sub(a[i-1], b[j-1], nfo)) {
        //go diagonally dont change the state
        insertEop(al, Replacement);
        i--; j--;
      } else if(cur == MATRIX2D(A, ncol, i, j)) {
        //move up
        state = 'A';
      } else if(cur == MATRIX2D(B, ncol, i, j)) {
        //move left
        state = 'B';
      }
    }

    if(state == 'A') {
      cur = MATRIX2D(A, ncol, i, j);
      if(cur == MATRIX2D(A, ncol, (i-1), j) + open) {
        //if this was an opening - change the state to S
        state = 'S';
      }
      //but move upwards anyways
      insertEop(al, Insertion);
      i--;
    }

    if(state == 'B') {
      cur = MATRIX2D(B, ncol, i, j);
      if(cur == MATRIX2D(B, ncol, i, (j-1)) + open) {
        //if this was an opening - change the state to S
        state = 'S';
      }
      //but move leftwards anyways
      insertEop(al, Deletion);
      j--;
    }
  }

  al->uoff = 0;
  al->voff = 0;
  revMeops(al);

  return ;
}


/*----------------------------- sgaffinecircular -----------------------------
 *    
 * @brief get circular matches with affine gap cost semi global alignment
 * @author Steve Hoffmann 
 *   
 */
 
void
sgaffinecircular (  )
{

  //take reference and double 
	
  
  return ;
}
