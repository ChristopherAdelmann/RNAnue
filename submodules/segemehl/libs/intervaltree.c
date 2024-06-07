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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <inttypes.h>
#include "mathematics.h"
#include "basic-types.h"
#include "iupac.h"
#include "debug.h"
#include "charsequence.h"
#include "biofiles.h"
#include "locus.h"
#include "intervaltree.h"
#include <pthread.h>
#include "segemehl.h"
#include "manopt.h"
#include "version.h"
#include "bedfiles.h"
#include "info.h"



int64_t 
gethigh(void *x, void *nfo) {
  int64_t *elem;
  elem = (int64_t*) x;

  return elem[1];
}

int64_t 
getlow(void *x, void *nfo) {
  int64_t *elem;
  elem = (int64_t*) x;

  return elem[0];
}


int64_t 
gethigh_annotitem(void *x, void *nfo) {
  annotationitem_t *elem;
  annotationoffs_t *offs;
  elem = (annotationitem_t*) x;
  offs = (annotationoffs_t*) nfo;
  uint64_t high = elem->end;
  uint64_t low = elem->start;
  int64_t d = high - low;
  /*
   * TODO concurring offsets!
   *
   */

  if(offs) {
    /*
     * first add any global right offset to end, high
     */
    if(offs->right < 0 && labs(offs->right) >= d) {
      high = low;
    } else {
      high += offs->right;
    }
    /*
     * now add directional offset to end, high
     */
    if(elem->strand == '+') {

      if(offs->dir3prime < 0 && labs(offs->dir3prime) >= d) {
        high = low;
      } else {
        high += offs->dir3prime;
      }
    } else {

      if(offs->dir5prime < 0 && labs(offs->dir5prime) >= d) {
        high = low;
      } else {
        high += offs->dir5prime;
      }

    }
  }

  return high;
}

int64_t 
getlow_annotitem(void *x, void *nfo) {
  annotationitem_t *elem;
  annotationoffs_t *offs;
  elem = (annotationitem_t*) x;
  offs = (annotationoffs_t*) nfo;
  uint64_t low = elem->start;
  uint64_t high = elem->end;
  int64_t d = high - low;

  /*
   * TODO concurring offsets!
   *
   */
  if(offs) {
    /*
     * first add any global left offset to start, low
     */
    
    if(offs->left > 0 && labs(offs->left) >= d) {
      low = high;
    } else {
      low += offs->left;
    }

    /*
     * now add directional offset to end, high
     */
    if(elem->strand == '+') {

     if(offs->dir5prime > 0 && labs(offs->dir5prime) >= d) {
        low = high;
      } else {
        low += offs->dir5prime;
      }

    } else {

     if(offs->dir3prime > 0 && labs(offs->dir3prime) >= d) {
        low = high;
      } else {
        low += offs->dir3prime;
      } 
    }
  }

  return low;
}


char
bl_intervaltreeCompare(uint64_t a, uint64_t b, uint64_t u, uint64_t v) {
  if(a > u) return 1;
  if(u > a) return -1;
  if(b > v) return 1;
  if(v > b) return -1;
  return 0;
}


void 
bl_intervaltreeInit(intervaltree_t *T) {
  T->l = NULL;
  T->r = NULL;
  T->max = 0;
  T->elem = NULL;
  T->ncoll = 0;
  T->collisions = NULL;
  T->h = 1;
}

uint64_t
bl_intervaltreeGetHeight(intervaltree_t *T){
  if(T == NULL) {
    return 0;
  } else {
    return T->h;
  }
}

uint64_t
bl_intervaltreeGetMaxPos(intervaltree_t *T){
  if(T == NULL) {
    return 0;
  } else {
    return T->max;
  }
}

//the left subinterval is the new root, its right
//child is the new left child of T,
//T becomes the rigth child of root
intervaltree_t *
bl_intervaltreeRightRotate(intervaltree_t *T, int64_t (*high)(void*, void*), void* nfo) {
  intervaltree_t *root, *tmp;

  root = T->l;
  tmp = T->l->r;

  root->r =  T;
  T->l =  tmp;

  T->h = MAX(bl_intervaltreeGetHeight(T->l), bl_intervaltreeGetHeight(T->r))+1;
  T->max = MAX3(high(T->elem, nfo), bl_intervaltreeGetMaxPos(T->l), bl_intervaltreeGetMaxPos(T->r));

  root->h = MAX(bl_intervaltreeGetHeight(root->r), bl_intervaltreeGetHeight(root->l))+1; //T->l->l
  root->max = MAX3(high(root->elem, nfo),bl_intervaltreeGetMaxPos(root->r), bl_intervaltreeGetMaxPos(root->l)); //T->l->l



  return root;
}



//the right subinterval is the new root, its left
//child is the new right child of T
intervaltree_t *
bl_intervaltreeLeftRotate(intervaltree_t *T, int64_t (*high)(void*, void*), void* nfo) {
  intervaltree_t *root, *tmp;

  root = T->r;
  tmp = T->r->l;

  root->l = T;
  T->r = tmp;

  T->h = MAX(bl_intervaltreeGetHeight(T->l), bl_intervaltreeGetHeight(T->r))+1;
  T->max = MAX3(high(T->elem, nfo),bl_intervaltreeGetMaxPos(T->l), bl_intervaltreeGetMaxPos(T->r));

  root->h = MAX(bl_intervaltreeGetHeight(root->r), bl_intervaltreeGetHeight(root->l))+1; //T->r->r
  root->max = MAX3(high(root->elem, nfo),bl_intervaltreeGetMaxPos(root->r), bl_intervaltreeGetMaxPos(root->l)); //T->r->r

  return root;
}


intervaltree_t*
bl_intervaltreeInsert(intervaltree_t *T, void *intervals, uint64_t k, 
    size_t elemsz, int64_t (*low)(void*, void*), int64_t (*high)(void *, void*), void *nfo) {

  char *ptr;
  int64_t imbalance, left_h, right_h;
  int64_t a,b,u,v,p,q,r,s;


  ptr = intervals;
  ptr = &ptr[k*elemsz];
  a = low(ptr, nfo);
  b = high(ptr, nfo);
  //fprintf(stderr, "[%"PRId64",%"PRId64"]->T\n", low(ptr), high(ptr));

  if(T == NULL) {      
    T = ALLOCMEMORY(NULL, NULL, intervaltree_t, 1);
    bl_intervaltreeInit(T);
    T->low = low(ptr, nfo);
    T->max = high(ptr, nfo);
    T->elem = ptr;
    return T;
  }

  u = low(T->elem, nfo);
  v = high(T->elem, nfo);

  if(bl_intervaltreeCompare(a, b, u, v) < 0) {
    //if(low(ptr) < T->low) 
    T->max = MAX(T->max, high(ptr, nfo));
    T->l = bl_intervaltreeInsert(T->l, intervals, k, elemsz, low, high, nfo);
  } else if (bl_intervaltreeCompare(a, b, u, v) > 0) {
    // else if (low(ptr) > T->low) 
    T->max = MAX(T->max, high(ptr, nfo));
    T->r = bl_intervaltreeInsert(T->r, intervals, k, elemsz, low, high, nfo);
  } else {
    //handle key collision
    //fprintf(stderr, "adding collisions\n");
    T->collisions= ALLOCMEMORY(NULL, T->collisions, void*, T->ncoll+1);
    T->collisions[T->ncoll] = ptr;
    T->ncoll++;
    return T;
  }

  left_h = bl_intervaltreeGetHeight(T->l);
  right_h = bl_intervaltreeGetHeight(T->r);
  T->h = MAX(left_h,right_h) + 1;

  imbalance = left_h - right_h;

  if(imbalance > 1) {
    p = low(T->l->elem, nfo);
    q = high(T->l->elem, nfo);
  }

  //left tree is larger, left value is larger (left, left)
  if(imbalance > 1 &&  bl_intervaltreeCompare(a, b, p, q) < 0) {
    //if(imbalance > 1 && low(ptr) < T->l->low) 
    //right rotate
    //fprintf(stderr, "ll\n");
    return bl_intervaltreeRightRotate(T, high, nfo);
  }



  //left tree is larger, left value is smaller (left, right) 
  if(imbalance > 1 && bl_intervaltreeCompare(a, b, p, q) > 0) {
    //if(imbalance > 1 && low(ptr) > T->l->low) 
    //left and right rotate
    //fprintf(stderr, "lr\n");
    T->l = bl_intervaltreeLeftRotate(T->l, high, nfo);
    return bl_intervaltreeRightRotate(T, high, nfo);
  }

  if(imbalance < -1) {
    r = low(T->r->elem, nfo);
    s = high(T->r->elem, nfo);
  }

  //right tree is larger, right value is smaller (right, right) 
  if(imbalance < -1 &&  bl_intervaltreeCompare(a, b, r, s) > 0) {
    //if(imbalance < -1 && low(ptr) > T->r->low) 
    //left rotate  
    //fprintf(stderr, "rr \n");
    return bl_intervaltreeLeftRotate(T, high, nfo);
  }


  //right tree is larger, right value is larger (right, left)
  if(imbalance < -1 &&bl_intervaltreeCompare(a, b, r, s) < 0) {
    //if(imbalance < -1 && low(ptr) < T->r->low) 
    //right and left rotate
    //fprintf(stderr, "rl\n");
    T->r = bl_intervaltreeRightRotate(T->r, high,nfo);
    return bl_intervaltreeLeftRotate(T, high, nfo);
  }

  return T;
}


void* 
bl_intervaltreeSearch(intervaltree_t *T, void* interval, 
    int64_t (*low)(void *, void*), int64_t (*high)(void *, void*), void ***results, 
    uint64_t *nelems, void* nfo) {

  char *ptr;
  uint64_t k;

  ptr = interval;

  if(T == NULL) {
    return NULL;
  }

  //query interval
  int64_t q_lo = low(ptr, nfo);
  int64_t q_hi = high(ptr, nfo);
  //current tree interval
  int64_t t_lo = low(T->elem, nfo);
  int64_t t_hi = high(T->elem, nfo);
  int64_t t_max = T->max;

  //check left 
  if(t_max < q_lo) {
    return NULL;
  }

  //search left
  if(T->l) {
    bl_intervaltreeSearch(T->l, interval, low, high, results, nelems, nfo);
  }


  //this interval overlaps with query
  if(t_lo <= q_hi && t_hi >= q_lo) {
    uint64_t n = *nelems;
    void **ptr = *results;

    if(n == 0) { 
      assert(ptr == NULL); 
    } else {
    //  printf("0x%" PRIXPTR "\n", (uintptr_t)results);
    //  printf("0x%" PRIXPTR "\n", (uintptr_t)*results);
    }

    ptr = ALLOCMEMORY(NULL, ptr, void*, n+1);
    //fprintf(stderr, "adding interval [%"PRId64",%"PRId64"] ovl w/ [%"PRId64",%"PRId64"] \n", 
    //   t_lo, t_hi, q_lo, q_hi);
    ptr[n] = T->elem;
    n++;

    for(k=0; k < T->ncoll; k++) {
      ptr = ALLOCMEMORY(NULL, ptr, void*, n+1);
      //fprintf(stderr, "adding collision [%"PRId64",%"PRId64"] ovl w/ [%"PRId64",%"PRId64"] \n", 
      //    t_lo, t_hi, q_lo, q_hi);
      ptr[n] = T->collisions[k];
      n++;
    }

    *results = ptr;
    *nelems = n;
  }

  //check right
  if(t_lo > q_hi) {
    return NULL;
  }

  //search right
  if(T->r) {
    bl_intervaltreeSearch(T->r, interval, low, high, results, nelems, nfo);
  }

  return NULL;
}

void 
preOrder(FILE *dev, intervaltree_t *root, int64_t (*low)(void *, void*), int64_t (*high)(void *, void*), void *nfo) {
  if(root != NULL)
  {
    fprintf(dev, "[%"PRId64",%"PRId64"] (max:%"PRId64")\n", low(root->elem, nfo), high(root->elem,nfo), root->max);
    preOrder(dev, root->l, low, high, nfo);
    preOrder(dev, root->r, low, high, nfo);
  }
}

void
bl_intervaltreeDestruct(intervaltree_t *root) {
  if(root->l != NULL) {
    bl_intervaltreeDestruct(root->l);
  }
  if(root->r != NULL) {
    bl_intervaltreeDestruct(root->r);
  }
  FREEMEMORY(NULL, root->l);
  FREEMEMORY(NULL, root->r);

  if(root->ncoll) FREEMEMORY(NULL, root->collisions);
}

void* 
bl_intervalforestSearch(intervalforest_t *T, char *treename, void* interval, 
    int64_t (*low)(void *, void*), int64_t (*high)(void *, void*), void ***results, 
    uint64_t *nelems, void *nfo) {

  uint32_t j;


  for(j=0; j < T->n; j++) { 
    if(!strcmp(T->treenames[j], treename)) {
      break;
    }
  }

  if(j < T->n) { 
    bl_intervaltreeSearch(T->trees[j], interval, getlow_annotitem, 
        gethigh_annotitem, results, nelems, NULL);
  }

  return NULL;
}


void
bl_intervalforestInsert(intervalforest_t *T, char *treename, void *intervals, 
    uint64_t k, size_t elemsz,int64_t (*low)(void *,void*), 
    int64_t (*high)(void *, void*), void *nfo) {

  uint32_t j;

  for(j=0; j < T->n; j++) {
    if(!strcmp(T->treenames[j], treename)) {
      break;
    }
  }
  if(j == T->n) {
    T->treenames = ALLOCMEMORY(NULL, T->treenames, char*, T->n+1);
    T->trees = ALLOCMEMORY(NULL, T->trees, intervaltree_t, T->n+1);
    T->trees[j] = NULL;
    T->treenames[j] = bl_strdup(treename);
    T->n++;
  }

  T->trees[j] = bl_intervaltreeInsert(T->trees[j], intervals, k, 
      elemsz, low, high, nfo);

}

intervalforest_t*
bl_intervalforestBuildAnnot(annotationmultitrack_t *mannot) {
  uint32_t i;
  intervalforest_t *forest;

  forest = ALLOCMEMORY(NULL, NULL, intervalforest_t, 1);
  forest->n = 0;
  forest->trees = NULL;
  forest->treenames = NULL;

  for(i=0; i < mannot->noofitems; i++) {   
    bl_intervalforestInsert(forest, mannot->items[i].chromname, mannot->items, i, 
            sizeof(annotationitem_t), getlow_annotitem, gethigh_annotitem, NULL);

  }

  return forest;
}

void
bl_intervalforestDestruct(intervalforest_t *forest) {

  uint32_t i; 
  for(i=0; i < forest->n; i++) {
    bl_intervaltreeDestruct(forest->trees[i]);
    FREEMEMORY(NULL, forest->trees[i]); 
    FREEMEMORY(NULL, forest->treenames[i]); 
  }
  FREEMEMORY(NULL, forest->treenames);
  FREEMEMORY(NULL, forest->trees);
}


/* Constructing tree given in the above figure */
/* 

   uint64_t intervals[] = {10,15,20,25,30,35,40,45,50,55,25,30,23,28,19,20,38,52,25,30};
   root = bl_intervaltreeInsert(root, intervals, 8, sizeof(uint64_t), getlow, gethigh);
   root = bl_intervaltreeInsert(root, intervals, 10, sizeof(uint64_t), getlow, gethigh); 
   root = bl_intervaltreeInsert(root, intervals, 4, sizeof(uint64_t), getlow, gethigh);
   root = bl_intervaltreeInsert(root, intervals, 6, sizeof(uint64_t), getlow, gethigh);
   root = bl_intervaltreeInsert(root, intervals, 12, sizeof(uint64_t), getlow, gethigh);
   root = bl_intervaltreeInsert(root, intervals, 14, sizeof(uint64_t), getlow, gethigh);
   root = bl_intervaltreeInsert(root, intervals, 16, sizeof(uint64_t), getlow, gethigh);
   root = bl_intervaltreeInsert(root, intervals, 0, sizeof(uint64_t), getlow, gethigh);
   root = bl_intervaltreeInsert(root, intervals, 2, sizeof(uint64_t), getlow, gethigh);
   root = bl_intervaltreeInsert(root, intervals, 18, sizeof(uint64_t), getlow, gethigh);
   */
/*
   The constructed AVL Tree would be
   30
   /  \
   20   40
   /  \     \
   10  25    50

*/
/* 
   uint64_t myinterval[]={26,31};
//uint64_t myinterval[]={10,25};
*/



