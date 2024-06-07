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


#ifndef INTERVALTREE_H
#define INTERVALTREE_H

#include "basic-types.h"
#include "biofiles.h"


//this is an AVL tree node with an extra field max with the
//rightmost postion in the corresponding subtree

struct intervaltree_s {
    int64_t max; 
    int64_t low; 
    uint64_t k;
    uint64_t h;
    void *elem;
    uint64_t ncoll;
    void **collisions;
    struct intervaltree_s *l;
    struct intervaltree_s *r;
};

typedef struct intervaltree_s intervaltree_t;

typedef struct intervaltrees_s {
  intervaltree_t **trees;
  char **treenames;
  uint32_t n;
} intervalforest_t;


int64_t 
gethigh_annotitem(void *x, void* nfo) ;
int64_t 
getlow_annotitem(void *x, void *nfo) ;
intervaltree_t*
bl_intervaltreeInsert(intervaltree_t *T, void *intervals, uint64_t k, 
    size_t elemsz, int64_t (*low)(void *,void*), int64_t (*high)(void *, void*), void *nfo) ;
void* 
bl_intervaltreeSearch(intervaltree_t *T, void* interval, 
    int64_t (*low)(void *, void*), int64_t (*high)(void *, void*), void ***results, 
    uint64_t *nelems, void *nfo);

void* 
bl_intervalforestSearch(intervalforest_t *T, char *treename, void* interval, 
    int64_t (*low)(void *, void*), int64_t (*high)(void *, void*), void ***results, 
    uint64_t *nelems, void *nfo);

void
bl_intervalforestInsert(intervalforest_t *T, char *treename, void *intervals, uint64_t k, size_t elemsz,int64_t (*low)(void *,void*), int64_t (*high)(void *, void*), void *nfo);

intervalforest_t *
bl_intervalforestBuildAnnot(annotationmultitrack_t *mannot);

  void 
preOrder(FILE *dev, intervaltree_t *root, int64_t (*low)(void *, void*), int64_t (*high)(void *, void*), void *nfo);

void
bl_intervaltreeDestruct(intervaltree_t *root);

void
bl_intervalforestDestruct(intervalforest_t *forest); 


void
bl_annotationitemInitOffset(annotationoffs_t *off); 
#endif 
