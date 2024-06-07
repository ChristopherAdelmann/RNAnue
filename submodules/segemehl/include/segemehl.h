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

#ifndef SEGEMEHL_H
#define SEGEMEHL_H

/*
 *
 *	segemehl.h
 *  declarations for threaded segemehl
 *
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig
 *  @date 01/02/2008 10:12:46 AM CET
 *
 *  SVN
 *  Revision of last commit: $Rev: 101 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-12-08 02:29:27 +0100 (Mon, 08 Dec 2008) $
 *
 *  Id: $Id: segemehl.h 101 2008-12-08 01:29:27Z steve $
 *  Url: $URL:
 *http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/src/segemehl.h $
 *
 */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include "alignment.h"
#include "basic-types.h"
#include "biofiles.h"
#include "charsequence.h"
#include "fileBins.h"
#include "filebuffer.h"
#include "fileio.h"
#include "htslib/sam.h"
#include "multicharseq.h"
#include "stringutils.h"
#include "sufarray.h"

#define MAXFILEBINS 50

typedef struct mappingstats_s {
    uint64_t total;
    uint64_t mapped;
    uint64_t uniquemapped;
    uint64_t multiplemapped;
    uint64_t unmapped;
    uint64_t paired;
    uint64_t uniquepaired;
    uint64_t multiplepaired;
    uint64_t singlequerymapped;
    uint64_t singlematemapped;
    uint64_t splitpair;
    uint64_t singlesplit;

} mappingstats_t;

typedef struct {
    Uint beg;
    Uint end;
    Uint n;
    Uint *binids;
    Uint *queuepos;
    Uint chk;
    char processed;
} mergemap_t;

typedef struct {
    mergemap_t *map;
    Uint mapsz;
} merge_t;

typedef struct segemehl_s {
    void *space;
    void *head;
    char *outfn;
    char *queryfilename;
    char *matefilename;
    char *idxfilename;
    char *idx2filename;
    char *dbfilename;
    char *nomatchname;
    char *samnomatchname;
    char *matelessfilename;
    char *softclip3Prime;
    char *softclip5Prime;
    char *filebinbasename;
    char *splitfilebasename;
    FILE *dev;
    FILE *nomatchdev;
    FILE *samnomatchdev;
    FILE *splitdev;
    htsFile *bamdev;
    unsigned char bam;
    bam_hdr_t *bamhdr;
    mappingstats_t *stats;
    bl_fileBinDomains_t *bins;
    bl_fileBinDomains_t *splitbins;
    int (*slct)(void *, void *);
    Suffixarray *sarray;
    MultiCharSeq *seq;
    fasta_t *reads;
    fasta_t *fasta;
    char *cmdline;
    Uint hardclip;
    Uint hardclip3Prime;
    Uint hardclip5Prime;
    Uint softclip3PrimeLen;
    Uint softclip5PrimeLen;
    Uint totallength;
    Uint minsize;
    Uint *mapmatches;
    Uint *counter;
    Uint M;
    Uint maxmyers;
    Uint maxout;
    Uint jump;
    Uint s_ext;
    Uint p_mis;
    Uint Xoff;
    Uint k;
    Uint k_p;
    Uint rep_type;
    Uint check;
    Uint kMis;
    Uint threadno;
    Uint threadid;
    Uint bestonly;
    Uint maxpairinsertsize;
    Uint maxaligninsertsize;
    Uint polyAlen;
    Uint minclipscr3;
    Uint minclipscr5;
    Uint bisulfiterun;
    Uint bisulfiteprotocol;
    Uint bisulfitemerging;
    Uint bisulfite;
    Uint strand;
    Uint minfragmentalignlen;
    int scores[3];
    int indel;
    int transition;
    int localparams[3];
    int minfragmentalignscore;
    Uint minsplicedaligncover;
    int minsplicedalignscore;
    int minscore;
    unsigned char split;
    // unsigned char index;
    unsigned char bining;
    unsigned char matchingstat;
    unsigned char align;
    unsigned char mute;
    unsigned char oldversion;
    unsigned char showMateless;
    unsigned char polyA;
    unsigned char order;
    unsigned char sam;
    unsigned char autoclip;
    unsigned char SAMmeop;
    unsigned char SAMpairstat;
    unsigned char nohead;
    unsigned char nosuflinks;
    double maxevalue;
    double maxsplitevalue;
    int accuracy;
    int clipacc;
    Uint bedist;
    Uint fusion;
    Eoptype cliptype;
    unsigned char mappingqual;
    double chainscorescale;
    pthread_mutex_t *mtx;
    pthread_mutex_t *mtx2;
    pthread_mutex_t *mtx3;
    pthread_mutex_t *mtx4;
    pthread_mutex_t *mtx5;
    pthread_mutex_t *mtx6;
    pthread_mutex_t *mtx7;
    pthread_mutex_t *mtx8;
    pthread_mutex_t *mtx9;
    unsigned int *nextfastaidx;
    PairUint *hash;
    Uint hashsize;
    unsigned char briefcigar;
    char *readgroupid;
    char *readgroupinfo;
    char *readgroupfile;
    char *singlesplitfilename;
    char *multisplitfilename;
    char *transsplitfilename;
    FILE *singlesplitdev;
    FILE *multisplitdev;
    FILE *transsplitdev;
    Uint maxcollinearlen;
    Uint maxcircularlen;
    double minentropy;
    Uint maxfragocc;
    Uint splitfragmargin;
    Uint maxsplitfragmargin;
    Uint fixsplitlen;
    Uint maxfixsplitocc;
    Uint fixsplitdist;
    int fixsplittransition;
    int bufferedwrite;
    /*  circbuffer_t *sambuffer;
      circbuffer_t *snglbuffer;
      circbuffer_t *multbuffer;
      circbuffer_t *trnsbuffer;
      */
    merge_t *merge;
    char orientation;
    char fullname;
} segemehl_t;

typedef struct checkthread_s {
    Uint noofseqs;
    Uint *counter;
} checkthreadinfo_t;

inline static void se_setdefault(segemehl_t *info) {
    info->bins = NULL;
    info->split = 0;
    info->splitbins = NULL;
    info->bining = 0;
    info->slct = bl_fileBinsCClassSelect;
    info->dev = stdout;
    info->idxfilename = NULL;
    info->idx2filename = NULL;
    info->dbfilename = NULL;
    info->queryfilename = NULL;
    info->matefilename = NULL;
    info->splitfilebasename = NULL;
    info->multisplitfilename = NULL;
    info->singlesplitfilename = NULL;
    info->transsplitfilename = NULL;
    info->singlesplitdev = NULL;
    info->multisplitdev = NULL;
    info->transsplitdev = NULL;
    info->sarray = NULL;
    info->maxout = 0;
    info->nohead = 0;
    info->seq = NULL;
    info->fasta = NULL;
    info->polyAlen = 0;
    info->reads = NULL;
    info->outfn = NULL;
    info->totallength = 0;
    info->counter = 0;
    info->minsize = 12;
    info->k_p = 1;
    // info->index = 1;
    info->threadno = 1;
    info->threadid = 0;
    info->accuracy = 90;
    info->clipacc = 70;
    info->maxsplitevalue = 50;
    info->M = 100;
    info->maxmyers = 100;
    info->jump = 0;
    info->s_ext = 2;
    info->p_mis = 4;
    info->Xoff = 8;
    info->rep_type = 12;
    info->kMis = 0;
    info->mute = 1;
    info->maxaligninsertsize = 1000;
    info->maxpairinsertsize = 200000;
    info->maxcollinearlen = 100000;
    info->maxcircularlen = 20000;
    info->matchingstat = 0;
    info->bestonly = 1;
    info->check = 0;
    info->maxevalue = 5;
    info->space = NULL;
    info->nomatchname = NULL;
    info->nomatchdev = NULL;
    info->samnomatchname = NULL;
    info->samnomatchdev = NULL;
    info->splitdev = NULL;
    info->showMateless = 1;
    info->bedist = 0;
    info->autoclip = 0;
    info->align = 0;
    info->order = 0;
    info->sam = 1;
    info->bisulfiterun = 0;
    info->bisulfiteprotocol = 0;
    info->bisulfitemerging = 0;
    info->bisulfite = 0;
    info->strand = 0;
    info->oldversion = 0;
    info->fusion = 0;
    info->hardclip = 0;
    info->hardclip3Prime = 0;
    info->hardclip5Prime = 0;
    info->softclip3Prime = NULL;
    info->softclip5Prime = NULL;
    info->filebinbasename = NULL;
    info->polyA = 0;
    info->softclip3PrimeLen = 0;
    info->softclip5PrimeLen = 0;
    info->minclipscr3 = 5;
    info->minclipscr5 = 5;
    info->SAMmeop = 0;
    info->SAMpairstat = 1;
    info->scores[0] = 1;
    info->scores[1] = -2;
    info->scores[2] = -2;
    info->localparams[0] = 20;
    info->localparams[1] = 50;
    info->indel = -3;
    info->transition = -10;
    info->chainscorescale = 0.9;
    info->minfragmentalignlen = 20;
    info->minfragmentalignscore = 18;
    info->minsplicedaligncover = 80;
    info->cmdline = NULL;
    info->nosuflinks = 0;
    info->minsplicedalignscore = 2 * 18;
    info->minscore = 20;
    info->mappingqual = 0;
    info->cliptype = Softclip;
    info->mtx = (pthread_mutex_t *)ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
    info->mtx2 = (pthread_mutex_t *)ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
    info->mtx3 = (pthread_mutex_t *)ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
    info->mtx4 = (pthread_mutex_t *)ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
    info->mtx5 = (pthread_mutex_t *)ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
    info->mtx6 = (pthread_mutex_t *)ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
    info->mtx7 = (pthread_mutex_t *)ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
    info->mtx8 = (pthread_mutex_t *)ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
    info->mtx9 = (pthread_mutex_t *)ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
    info->nextfastaidx = (unsigned int *)ALLOCMEMORY(space, NULL, Uint, 1);
    info->nextfastaidx[0] = 0;
    pthread_mutex_init(info->mtx, NULL);
    pthread_mutex_init(info->mtx2, NULL);
    pthread_mutex_init(info->mtx3, NULL);
    pthread_mutex_init(info->mtx4, NULL);
    pthread_mutex_init(info->mtx5, NULL);
    pthread_mutex_init(info->mtx6, NULL);
    pthread_mutex_init(info->mtx7, NULL);
    pthread_mutex_init(info->mtx8, NULL);
    pthread_mutex_init(info->mtx9, NULL);
    info->hash = NULL;
    info->hashsize = 10;
    info->stats = (mappingstats_t *)ALLOCMEMORY(space, NULL, mappingstats_t, 1);
    memset(info->stats, 0, sizeof(mappingstats_t));
    info->briefcigar = 0;
    info->readgroupid = NULL;
    info->readgroupinfo = NULL;
    info->readgroupfile = NULL;
    info->minentropy = 1.3;
    info->maxfragocc = 50;
    info->splitfragmargin = 50;
    info->maxsplitfragmargin = 150;
    info->fixsplitlen = 10;
    info->maxfixsplitocc = 30000;
    info->fixsplitdist = 20000;
    info->fixsplittransition = -2;
    info->bufferedwrite = 0;
    info->merge = (merge_t *)ALLOCMEMORY(space, NULL, merge_t, 1);
    info->merge->map = NULL;
    info->merge->mapsz = 0;
    info->bamdev = NULL;
    info->bam = 0;
    info->bamhdr = NULL;
    info->head = NULL;
    info->orientation = 1;
    info->fullname = 0;
}

inline static void se_destructInfo(void *space, segemehl_t *info) {
    bl_fastaDestruct(space, info->fasta);
    FREEMEMORY(space, info->fasta);

    if (info->bisulfiteprotocol) destructMultiCharSeq(space, info->seq);
    if (info->softclip3Prime) {
        FREEMEMORY(space, info->softclip3Prime);
    }

    if (info->queryfilename) {
        bl_fastaDestruct(space, info->reads);
        FREEMEMORY(space, info->reads);
    }

    if (info->singlesplitfilename) {
        FREEMEMORY(space, info->singlesplitfilename);
    }

    if (info->multisplitfilename) {
        FREEMEMORY(space, info->multisplitfilename);
    }

    if (info->transsplitfilename) {
        FREEMEMORY(space, info->transsplitfilename);
    }

    if (info->merge) {
        FREEMEMORY(space, info->merge);
    }

    FREEMEMORY(space, info->mtx);
    FREEMEMORY(space, info->mtx2);
    FREEMEMORY(space, info->mtx3);
    FREEMEMORY(space, info->mtx4);
    FREEMEMORY(space, info->mtx5);
    FREEMEMORY(space, info->mtx6);
    FREEMEMORY(space, info->mtx7);
    FREEMEMORY(space, info->mtx8);
    FREEMEMORY(space, info->mtx9);
    FREEMEMORY(space, info->nextfastaidx);
    FREEMEMORY(space, info->cmdline);
    FREEMEMORY(space, info->stats);
}

int segemehl(int argc, char **argv);

extern const char *SORT[];
extern const char *SORTBIN[];
extern const char SORTDELIM;
// void se_updateProgressBar(Uint k, segemehl_t *nfo);
// fasta_t* se_fastaMaster(void *space, fasta_t *f, Uint size, segemehl_t* nfo);
#endif
