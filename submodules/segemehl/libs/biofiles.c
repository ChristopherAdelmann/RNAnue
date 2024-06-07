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
 *  biofiles.c
 *  helper functions to handle file types
 *  used in bioinformatics
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/10/2007 01:56:15 PM CEST
 *
 *  SVN
 *  Revision of last commit: $Rev: 76 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-11-11 16:34:21 +0100 (Tue, 11 Nov 2008) $
 *
 *  Id: $Id: biofiles.c 76 2008-11-11 15:34:21Z steve $
 *  Url: $URL:
 * http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/biofiles.c $
 *
 */

#include "biofiles.h"

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "annotation.h"
#include "basic-types.h"
#include "bedfiles.h"
#include "charsequence.h"
#include "debug.h"   //(skipped)
#include "fileio.h"  //(skipped)
#include "gzidx.h"
#include "info.h"  //(skipped)
#include "mathematics.h"
#include "memory.h"       //(skipped)
#include "seqclip.h"      //(skipped)
#include "stringutils.h"  //(skipped)

/*------------------------------- bl_fastaInit -------------------------------
 *
 * @brief initialize the fasta struct
 * @author Steve Hoffmann
 *
 */

fasta_t *bl_fastaInit(void *space) {
    fasta_t *f;

    f = ALLOCMEMORY(space, NULL, fasta_t, 1);
    f->seqs = NULL;
    f->quals = NULL;
    f->active_noofseqs = 0;
    f->active_noofmates = 0;
    f->noofseqs = 0;
    f->minlen = 0;
    f->maxlen = 0;
    f->matestart = NULL;
    f->hasIndex = 0;
    f->gzip = 0;
    f->nooffiles = 0;
    f->filenames = NULL;
    f->matefilenames = NULL;
    f->filetotal = 0;
    f->findex = NULL;
    f->gzindex = NULL;
    f->chunkindex = NULL;
    f->matechunkindex = NULL;
    f->curchunk = 0;
    f->chunkIsActive = 0;
    f->hasMates = 0;
    f->checkid = 1;

    return f;
}

/*---------------------------- bl_fastaHasIndex -----------------------------
 *
 * @brief returns 1 if fasta is divided into chunks, 0 otherwise
 * @author Steve Hoffmann
 *
 */

unsigned char bl_fastaHasIndex(fasta_t *f) { return (f->hasIndex); }

/*--------------------------- bl_fastxInitSeqIndex ---------------------------
 *
 * @brief initialize a sequence index
 * @author Steve Hoffmann
 *
 */

fastxfileindex_t *bl_fastxInitFileIndex(void *space, Uint noofaccesspoints) {
    fastxfileindex_t *idx;
    idx = ALLOCMEMORY(space, NULL, fastxfileindex_t, 1);
    idx->ap = ALLOCMEMORY(space, NULL, seqaccesspoint_t, noofaccesspoints);
    memset(idx->ap, 0, sizeof(seqaccesspoint_t) * noofaccesspoints);
    idx->size = 0;
    idx->allocated = noofaccesspoints;

    return idx;
}

/*--------------------------- bl_fastxInitSeqIndex ---------------------------
 *
 * @brief initialize a sequence index
 * @author Steve Hoffmann
 *
 */

fastxseqindex_t *bl_fastxInitSeqIndex(void *space, Uint noofaccesspoints) {
    fastxseqindex_t *idx;
    idx = ALLOCMEMORY(space, NULL, fastxseqindex_t, 1);
    idx->ap = ALLOCMEMORY(space, NULL, seqaccesspoint_t, noofaccesspoints);
    memset(idx->ap, 0, sizeof(seqaccesspoint_t) * noofaccesspoints);
    idx->size = 0;
    idx->allocated = noofaccesspoints;

    return idx;
}

/*----------------------- bl_fastaGetDescriptionLength -----------------------
 *
 * @brief get the length of the fasta description
 * @author Steve Hoffmann
 *
 */

Uint bl_fastaGetDescriptionLength(fasta_t *f, Uint elem) {
    Uint k;
    if (!bl_fastaHasIndex(f)) return f->seqs[elem]->descrlen;

    k = bl_fastxGetChunkElem(NULL, f, elem);
    if (k == -1) {
        DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
        exit(-1);
    }
    return f->seqs[k]->descrlen;
}

/*-------------------------- bl_fastaGetDescription --------------------------
 *
 * @brief returns the descriptions of a fasta sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetDescription(fasta_t *f, Uint elem) {
    Uint k;
    if (!bl_fastaHasIndex(f)) return f->seqs[elem]->description;

    k = bl_fastxGetChunkElem(NULL, f, elem);
    if (k == -1) {
        DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
        exit(-1);
    }

    return f->seqs[k]->description;
}

/*--------------------- bl_fastaGetMateDescriptionLength ---------------------
 *
 * @brief get the length of the fasta mate description
 * @author Steve Hoffmann
 *
 */

Uint bl_fastaGetMateDescriptionLength(fasta_t *f, Uint elem) {
    Uint k;
    if (!bl_fastaHasIndex(f)) return f->seqs[elem]->noofinfo;

    k = bl_fastxGetChunkElem(NULL, f, elem);
    if (k == -1) {
        DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
        exit(-1);
    }
    return f->seqs[k]->noofinfo;
}

/*-------------------------- bl_fastaGetDescription --------------------------
 *
 * @brief returns the descriptions of a fasta sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetMateDescription(fasta_t *f, Uint elem) {
    Uint k;
    if (!bl_fastaHasIndex(f)) return (char *)f->seqs[elem]->info;

    k = bl_fastxGetChunkElem(NULL, f, elem);
    if (k == -1) {
        DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
        exit(-1);
    }

    return (char *)f->seqs[k]->info;
}

/*------------------------ bl_fastaGetSequenceLength -------------------------
 *
 * @brief returns the length of a fasta sequence
 * @author Steve Hoffmann
 *
 */

Uint bl_fastaGetSequenceLength(fasta_t *f, Uint elem) {
    Uint k, clip3 = 0, clip5 = 0;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
        }
    } else {
        k = elem;
    }

    clip3 = f->seqs[k]->clip3[0];
    clip5 = f->seqs[k]->clip5[0];

    if (bl_fastaHasMate(f)) {
        return f->matestart[k] - 1 - clip3 - clip5;
    }

    return f->seqs[k]->length - clip3 - clip5;
}

/*--------------------------- bl_fastaGetSequence ----------------------------
 *
 * @brief return the fasta sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetSequence(fasta_t *f, Uint elem) {
    Uint k, clip5 = 0;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        // fprintf(stdout, "found chunk elem %d\n", k);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    clip5 = f->seqs[k]->clip5[0];
    return &f->seqs[k]->sequence[clip5];
}

/*--------------------------- bl_fastaGetSequenceNoClip ----------------------------
 *
 * @brief return the fasta sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetSequenceNoClip(fasta_t *f, Uint elem) {
    Uint k;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        // fprintf(stdout, "found chunk elem %d\n", k);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    return f->seqs[k]->sequence;
}

/*--------------------------- bl_fastaGetLeftClip ----------------------------
 *
 * @brief return the fasta sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetLeftClip(fasta_t *f, Uint elem) {
    Uint k, clip5 = 0;
    char *seq = NULL;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        // fprintf(stdout, "found chunk elem %d\n", k);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    clip5 = f->seqs[k]->clip5[0];
    if (clip5) {
        seq = ALLOCMEMORY(NULL, NULL, char, clip5 + 1);
        memmove(seq, f->seqs[k]->sequence, sizeof(char) * clip5);
        seq[clip5] = 0;
    }

    return seq;
}

/*--------------------------- bl_fastaGetLeftClipQual ----------------------------
 *
 * @brief return the fasta sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetLeftClipQual(fasta_t *f, Uint elem) {
    Uint k, clip5 = 0;
    char *seq = NULL;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        // fprintf(stdout, "found chunk elem %d\n", k);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    clip5 = f->seqs[k]->clip5[0];
    if (clip5) {
        seq = ALLOCMEMORY(NULL, NULL, char, clip5 + 1);
        memmove(seq, f->quals[k]->sequence, sizeof(char) * clip5);
        seq[clip5] = 0;
    }
    return seq;
}

/*--------------------------- bl_fastaGetRightClip ----------------------------
 *
 * @brief return the fasta sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetRightClip(fasta_t *f, Uint elem) {
    Uint k, clip3 = 0;
    char *seq = NULL;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        // fprintf(stdout, "found chunk elem %d\n", k);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    clip3 = f->seqs[k]->clip3[0];
    if (clip3) {
        seq = ALLOCMEMORY(NULL, NULL, char, clip3 + 1);
        memmove(seq, f->seqs[k]->sequence, sizeof(char) * clip3);
        seq[clip3] = 0;
    }

    return seq;
}
/*--------------------------- bl_fastaGetRightClipQual ----------------------------
 *
 * @brief return the fasta sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetRightClipQual(fasta_t *f, Uint elem) {
    Uint k, clip3 = 0;
    char *seq = NULL;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        // fprintf(stdout, "found chunk elem %d\n", k);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    clip3 = f->seqs[k]->clip3[0];
    if (clip3) {
        seq = ALLOCMEMORY(NULL, NULL, char, clip3 + 1);
        memmove(seq, f->quals[k]->sequence, sizeof(char) * clip3);
        seq[clip3] = 0;
    }

    return seq;
}

/*----------------------------- bl_fastaHasMate ------------------------------
 *
 * @brief returns 1 if fasta has a mate pair, 0 otherwise
 * @author Steve Hoffmann
 *
 */

unsigned char bl_fastaHasMate(fasta_t *f) {
    return (unsigned char)(f->matestart != NULL || (f->hasMates && f->hasIndex));
}

/*--------------------------- bl_fastaGetMateStart ---------------------------
 *
 * @brief get the start pos of mate sequence in string
 * @author Steve Hoffmann
 *
 */

Uint bl_fastaGetMateStart(fasta_t *f, Uint elem) {
    Uint k = 0;

    if (!bl_fastaHasMate(f)) {
        return 0;
    }

    if (!bl_fastaHasIndex(f)) return f->matestart[elem];

    k = bl_fastxGetChunkElem(NULL, f, elem);
    if (k == -1) {
        DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
        exit(-1);
    }

    return f->matestart[k];
}

/*-------------------------- bl_fastaGetMateLength ---------------------------
 *
 * @brief return length of mate sequence
 * @author Steve Hoffmann
 *
 */

Uint bl_fastaGetMateLength(fasta_t *f, Uint elem) {
    Uint k, clip3 = 0, clip5 = 0;

    if (!bl_fastaHasMate(f)) {
        return 0;
    }

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);

        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
        }
    } else {
        k = elem;
    }

    clip3 = f->seqs[k]->clip3[1];
    clip5 = f->seqs[k]->clip5[1];

    return (f->seqs[k]->length - f->matestart[k]) - clip3 - clip5;
}

/*----------------------------- bl_fastaGetMate ------------------------------
 *
 * @brief returns the mate sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetMate(fasta_t *f, Uint elem) {
    char *res, clip5 = 0;
    Uint k;

    if (!bl_fastaHasMate(f)) {
        return NULL;
    }

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);

        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
        }
    } else {
        k = elem;
    }

    clip5 = f->seqs[k]->clip5[1];

    res = f->seqs[k]->sequence;
    return &res[f->matestart[k] + clip5];
}
/*----------------------------- bl_fastaGetMateNoClip ------------------------------
 *
 * @brief returns the mate sequence w/o clipping
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetMateNoClip(fasta_t *f, Uint elem) {
    char *res;
    Uint k;

    if (!bl_fastaHasMate(f)) {
        return NULL;
    }

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);

        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
        }
    } else {
        k = elem;
    }

    res = f->seqs[k]->sequence;
    return &res[f->matestart[k]];
}

/*---------------------------- bl_fastaGetClipPos ----------------------------
 *
 * @brief return the clipping positions
 * @author Steve Hoffmann
 *
 */

void bl_fastaGetClipPos(fasta_t *f, Uint elem, Uint *p5, Uint *p3) {
    Uint k;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);

        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
        }
    } else {
        k = elem;
    }

    *p5 = f->seqs[k]->clip5[0];
    *p3 = f->seqs[k]->clip3[0];

    return;
}

/*-------------------------- bl_fastaGetMateClipPos --------------------------
 *
 * @brief return the clipping positions
 * @author Steve Hoffmann
 *
 */

void bl_fastaGetMateClipPos(fasta_t *f, Uint elem, Uint *p5, Uint *p3) {
    Uint k;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);

        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
        }
    } else {
        k = elem;
    }

    *p5 = f->seqs[k]->clip5[1];
    *p3 = f->seqs[k]->clip3[1];

    return;
}
/*--------------------------- bl_fastaGetMateLeftClip ----------------------------
 *
 * @brief return the fasta sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetMateLeftClip(fasta_t *f, Uint elem) {
    Uint k, clip5 = 0;
    char *seq = NULL;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        // fprintf(stdout, "found chunk elem %d\n", k);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    clip5 = f->seqs[k]->clip5[1];

    if (clip5) {
        seq = ALLOCMEMORY(NULL, NULL, char, clip5 + 1);
        memmove(seq, f->seqs[k]->sequence, sizeof(char) * clip5);
        seq[clip5] = 0;
    }

    return seq;
}

/*--------------------------- bl_fastaGetMateLeftClipQual ----------------------------
 *
 * @brief return the fasta sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetMateLeftClipQual(fasta_t *f, Uint elem) {
    Uint k, clip5 = 0;
    char *seq = NULL;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        // fprintf(stdout, "found chunk elem %d\n", k);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    clip5 = f->seqs[k]->clip5[1];
    if (clip5) {
        seq = ALLOCMEMORY(NULL, NULL, char, clip5 + 1);
        memmove(seq, f->quals[k]->sequence, sizeof(char) * clip5);
        seq[clip5] = 0;
    }

    return seq;
}

/*--------------------------- bl_fastaGetMateRightClip ----------------------------
 *
 * @brief return the fasta sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetMateRightClip(fasta_t *f, Uint elem) {
    Uint k, clip3 = 0;
    char *seq = NULL;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        // fprintf(stdout, "found chunk elem %d\n", k);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    clip3 = f->seqs[k]->clip3[1];

    if (clip3) {
        seq = ALLOCMEMORY(NULL, NULL, char, clip3 + 1);
        memmove(seq, f->seqs[k]->sequence, sizeof(char) * clip3);
        seq[clip3] = 0;
    }

    return seq;
}

/*--------------------------- bl_fastaGetMateRightClipQual ----------------------------
 *
 * @brief return the fasta sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetMateRightClipQual(fasta_t *f, Uint elem) {
    Uint k, clip3 = 0;
    char *seq = NULL;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        // fprintf(stdout, "found chunk elem %d\n", k);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    clip3 = f->seqs[k]->clip3[1];
    if (clip3) {
        seq = ALLOCMEMORY(NULL, NULL, char, clip3 + 1);
        memmove(seq, f->quals[k]->sequence, sizeof(char) * clip3);
        seq[clip3] = 0;
    }
    return seq;
}

/*--------------------------- bl_fastaDestructMate ---------------------------
 *
 * @brief free mate information
 * @author Steve Hoffmann
 *
 */

void bl_fastaDestructMate(void *space, fasta_t *f) {
    assert(f->matestart);
    FREEMEMORY(space, f->matestart);

    return;
}

/*---------------------------- bl_fastaHasQuality ----------------------------
 *
 * @brief returns 1 if fasta has quality information, 0 otherwise
 * @author Steve Hoffmann
 *
 */

unsigned char bl_fastaHasQuality(fasta_t *f) { return (unsigned char)(f->quals != NULL); }

/*---------------------------- bl_fastaGetQuality ----------------------------
 *
 * @brief returns the quality information of the sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetQuality(fasta_t *f, Uint elem) {
    Uint k, clip5 = 0;

    if (!bl_fastaHasQuality(f)) {
        return NULL;
    }

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);

        if (k == -1) {
            DBG("retrieval of quality %d failed. Exit forced.\n", elem);
        }
    } else {
        k = elem;
    }

    clip5 = f->seqs[k]->clip5[0];

    return &f->quals[k]->sequence[clip5];
}

/*---------------------------- bl_fastaGetQualityNoClip ----------------------------
 *
 * @brief returns the quality information of the sequence
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetQualityNoClip(fasta_t *f, Uint elem) {
    Uint k;

    if (!bl_fastaHasQuality(f)) {
        return NULL;
    }

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);

        if (k == -1) {
            DBG("retrieval of quality %d failed. Exit forced.\n", elem);
        }
    } else {
        k = elem;
    }

    return f->quals[k]->sequence;
}

/*-------------------------- bl_fastaGetMateQuality --------------------------
 *
 * @brief returns the quality string of the mate sequence elem
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetMateQuality(fasta_t *f, Uint elem) {
    char *res, clip5 = 0;
    Uint k;

    if (!bl_fastaHasMate(f) || !bl_fastaHasQuality(f)) {
        return NULL;
    }

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);

        if (k == -1) {
            DBG("retrieval of quality %d failed. Exit forced.\n", elem);
        }
    } else {
        k = elem;
    }

    clip5 = f->seqs[k]->clip5[1];

    res = f->quals[k]->sequence;
    return &res[f->matestart[k] + clip5];
}

/*-------------------------- bl_fastaGetMateQualityNoClip --------------------------
 *
 * @brief returns the quality string of the mate sequence elem
 * @author Steve Hoffmann
 *
 */

char *bl_fastaGetMateQualityNoClip(fasta_t *f, Uint elem) {
    char *res;
    Uint k;

    if (!bl_fastaHasMate(f) || !bl_fastaHasQuality(f)) {
        return NULL;
    }

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);

        if (k == -1) {
            DBG("retrieval of quality %d failed. Exit forced.\n", elem);
        }
    } else {
        k = elem;
    }

    res = f->quals[k]->sequence;
    return &res[f->matestart[k]];
}

/*------------------------- bl_fastaDestructQuality --------------------------
 *
 * @brief free quality information
 * @author Steve Hoffmann
 *
 */

void bl_fastaDestructQuality(void *space, fasta_t *f) {
    Uint i;
    assert(f->quals);

    for (i = 0; i < f->active_noofseqs; i++) {
        destructSequence(space, f->quals[i]);
    }

    return;
}

#ifdef HASHING
/*---------------------------- bl_fastaGetQuantity -----------------------------
 *
 * @brief get quantity of fasta sequence in terms of tag count
 *        (number of input sequences with equal nucleotides)
 * @author Christian Otto
 *
 */
Uint bl_fastaGetQuantity(fasta_t *f, Uint elem) {
    Uint k;

    if (!bl_fastaHasIndex(f)) return f->seqs[elem]->quantity;

    k = bl_fastxGetChunkElem(NULL, f, elem);

    if (k == -1) {
        DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
    }

    return f->seqs[k]->quantity;
}

/*------------------------ bl_fastaSetQuantity ------------------------
 *
 *
 * @brief set quantity of fasta sequence in terms of tag count
 *        (number of input sequences with equal nucleotides)
 * @author Christian Otto
 *
 */
void bl_fastaSetQuantity(fasta_t *f, Uint elem, Uint quantity) {
    Uint k;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    f->seqs[k]->quantity = quantity;
}

#endif

/*------------------------ bl_fastaSetMateDescription ------------------------
 *
 * @brief set the description for a fasta sequence elem
 * @author Steve Hoffmann
 *
 */

void bl_fastaSetMateDescription(fasta_t *f, Uint elem, char *descr, Uint len) {
    assert(descr[0] == '@' || descr[0] == '>');
    memmove(descr, &descr[1], len - 1);
    descr[len - 1] = 0;
    f->seqs[elem]->info = descr;
    f->seqs[elem]->noofinfo = len - 1;
    return;
}

/*-------------------------- bl_fastaSetDescription --------------------------
 *
 * @brief set the description for a fasta sequence elem
 * @author Steve Hoffmann
 *
 */

void bl_fastaSetDescription(fasta_t *f, Uint elem, char *descr, Uint len) {
    assert(descr[0] == '@' || descr[0] == '>');
    memmove(descr, &descr[1], len - 1);
    descr[len - 1] = 0;
    f->seqs[elem]->description = descr;
    f->seqs[elem]->descrlen = len - 1;

    return;
}

/*--------------------------- bl_fastaSetSequence ----------------------------
 *
 * @brief set the fasta sequence elem
 * @author Steve Hoffmann
 *
 */

void bl_fastaSetSequence(void *space, fasta_t *f, Uint elem, char *seq, Uint len) {
    Uint oldlen, newlen;
    char *oldseq;

    oldlen = f->seqs[elem]->length;
    oldseq = f->seqs[elem]->sequence;
    newlen = oldlen + len;

    if (oldlen) {
        newlen = len + oldlen + 1;
        seq = ALLOCMEMORY(space, seq, char, newlen + 1);
        f->matestart = ALLOCMEMORY(space, f->matestart, Uint, elem + 1);
        f->matestart[elem] = oldlen + 1;
        memmove(&seq[oldlen + 1], seq, len);
        memmove(seq, oldseq, oldlen);
        seq[oldlen] = '\0';
        seq[newlen] = '\0';
        FREEMEMORY(space, oldseq);
    }

    f->seqs[elem]->length = newlen;
    f->seqs[elem]->sequence = seq;

    return;
}

/*----------------------------- bl_fastaSoftClip -----------------------------
 *
 * @brief clip sequence
 * @author Steve Hoffmann
 *
 */

Uint bl_fastaSoftClip(void *space, fasta_t *f, Uint elem, char *p5, Uint p5len, Uint p5scr,
                      char *p3, Uint p3len, Uint p3acc, Uint pAlen) {
    Uint len, l = 0, r = 0, k;
    char *seq;

    seq = bl_fastaGetSequence(f, elem);
    len = bl_fastaGetSequenceLength(f, elem);

    // fprintf(stderr, "***************** clipping '%s' - '%s' of length %d - %d\n", p5, p3, p5len,
    // p3len);

    if (p3 && p3len) r = bl_seqclipSoft3Prime(space, seq, len, p3, p3len, p3acc, pAlen);
    if (p5 && p5len) l = bl_seqclipSoft5Prime(space, seq, len, p5, p5len, p5scr);

    //  fprintf(stderr, "clipping returned l:%d and r:%d\n", l, r);

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    if (f->seqs[k]->clip5[0] + f->seqs[k]->clip3[0] + l + r >= len) {
        return 0;
    }

    f->seqs[k]->clip5[0] += l;
    f->seqs[k]->clip3[0] += r;

    return l + r;
}

/*--------------------------- bl_fastaMateSoftClip ---------------------------
 *
 * @brief clip sequence
 * @author Steve Hoffmann
 *
 */

Uint bl_fastaMateSoftClip(void *space, fasta_t *f, Uint elem, char *p5, Uint p5len, Uint p5scr,
                          char *p3, Uint p3len, Uint p3acc, Uint pAlen) {
    Uint len, l = 0, r = 0, k;
    char *seq;

    if (!bl_fastaHasMate(f)) {
        return 0;
    }

    seq = bl_fastaGetMate(f, elem);
    len = bl_fastaGetMateLength(f, elem);

    if (p3 && p3len) r = bl_seqclipSoft3Prime(space, seq, len, p3, p3len, p3acc, pAlen);
    if (p5 && p5len) l = bl_seqclipSoft5Prime(space, seq, len, p5, p5len, p5scr);

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    if (f->seqs[k]->clip5[1] + f->seqs[k]->clip3[1] + l + r >= len) {
        return 0;
    }

    f->seqs[k]->clip5[1] += l;
    f->seqs[k]->clip3[1] += r;

    return l + r;
}

/*----------------------------- bl_fastaHardClip -----------------------------
 *
 * @brief hard clip sequence
 * @author Steve Hoffmann
 *
 */

Uint bl_fastaHardClip(void *space, fasta_t *f, Uint elem, Uint p5, Uint p3) {
    Uint k, len;

    len = bl_fastaGetSequenceLength(f, elem);

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    if (f->seqs[k]->clip5[0] + f->seqs[k]->clip3[0] + p3 + p5 >= len) {
        return 0;
    }

    f->seqs[k]->clip5[0] += p5;
    f->seqs[k]->clip3[0] += p3;

    return p5 + p3;
}

/*--------------------------- bl_fastaMateHardClip ---------------------------
 *
 * @brief hard clip sequence
 * @author Steve Hoffmann
 *
 */

Uint bl_fastaMateHardClip(void *space, fasta_t *f, Uint elem, Uint p5, Uint p3) {
    Uint k, len;

    if (!bl_fastaHasMate(f)) return 0;

    len = bl_fastaGetMateLength(f, elem);

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    if (f->seqs[k]->clip5[1] + f->seqs[k]->clip3[1] + p3 + p5 >= len) {
        return 0;
    }

    f->seqs[k]->clip5[1] += p5;
    f->seqs[k]->clip3[1] += p3;

    return p5 + p3;
}

void bl_fastaSetClip(fasta_t *f, Uint elem, Uint p5, Uint p3) {
    Uint k;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    f->seqs[k]->clip5[0] = p5;
    f->seqs[k]->clip3[0] = p3;
}

void bl_fastaSetMateClip(fasta_t *f, Uint elem, Uint p5, Uint p3) {
    Uint k;

    if (bl_fastaHasIndex(f)) {
        k = bl_fastxGetChunkElem(NULL, f, elem);
        if (k == -1) {
            DBG("retrieval of sequence %d failed. Exit forced.\n", elem);
            exit(-1);
        }
    } else {
        k = elem;
    }

    f->seqs[k]->clip5[1] = p5;
    f->seqs[k]->clip3[1] = p3;
}

/*---------------------------- bl_fastaSeqQuality ----------------------------
 *
 * @brief set the quality information for sequence eleme
 * @author Steve Hoffmann
 *
 */

void bl_fastaSetQuality(void *space, fasta_t *f, Uint elem, char *qlty, Uint len) {
    Uint oldlen, newlen;
    char *oldqlty;

    assert(f->quals && f->quals[elem]);

    oldlen = f->quals[elem]->length;
    oldqlty = f->quals[elem]->sequence;
    newlen = len + oldlen;

    if (oldlen) {
        newlen = len + oldlen + 1;
        qlty = ALLOCMEMORY(space, qlty, char, newlen + 1);
        memmove(&qlty[oldlen + 1], qlty, len);
        memmove(qlty, oldqlty, oldlen);
        qlty[oldlen] = '\0';
        qlty[newlen] = '\0';
        FREEMEMORY(space, oldqlty);
    }

    f->quals[elem]->length = newlen;
    f->quals[elem]->sequence = qlty;

    return;
}

/*---------------------------- bl_fastaAddQuality ----------------------------
 *
 * @brief allocate memory for a new quality information
 * @author Steve Hoffmann
 *
 */

void bl_fastaAddQuality(void *space, fasta_t *f) {
    Uint n = f->active_noofseqs;

    f->quals = ALLOCMEMORY(space, (f->quals), CharSequence *, n + 1);
    assert(f->quals != NULL);
    f->quals[n] = initSequence(space);

    return;
}

/*--------------------------- bl_fastaAddSequence ----------------------------
 *
 * @brief alloc memory for a new fasta sequence
 * @author Steve Hoffmann
 *
 */

void bl_fastaAddSequence(void *space, fasta_t *f) {
    Uint n = f->active_noofseqs;

    f->seqs = ALLOCMEMORY(space, (f->seqs), CharSequence *, n + 1);
    assert(f->seqs != NULL);
    f->seqs[n] = initSequence(space);
    f->seqs[n]->clip3[0] = 0;
    f->seqs[n]->clip3[1] = 0;
    f->seqs[n]->clip5[0] = 0;
    f->seqs[n]->clip5[1] = 0;
    return;
}

/*------------------------------- bl_fastaAdd --------------------------------
 *
 * @brief add a new fasta structure
 * @author Steve Hoffmann
 *
 */

void bl_fastaAdd(void *space, fasta_t *f, char *descr, Uint descrlen, char *sequence, Uint seqlen,
                 Uint n) {
    assert(n == f->active_noofseqs);

    bl_fastaAddSequence(space, f);
    bl_fastaSetDescription(f, n, descr, descrlen);
    bl_fastaSetSequence(space, f, n, sequence, seqlen);

    f->minlen = (seqlen < f->minlen) ? seqlen : f->minlen;
    f->maxlen = (seqlen > f->maxlen) ? seqlen : f->maxlen;
    f->active_noofseqs++;

    if (!f->hasIndex) {
        f->noofseqs++;
    }

    return;
}

/*------------------------------- bl_fastaxAdd -------------------------------
 *
 * @brief add a new fasta structure for fastq or fasta information
 * @author Steve Hoffmann
 *
 */

void bl_fastxAdd(void *space, fasta_t *f, char *descr, Uint descrlen, char *sequence, char *quality,
                 Uint seqlen, Uint n) {
    assert(n == f->active_noofseqs);

    bl_fastaAddSequence(space, f);
    bl_fastaSetDescription(f, f->active_noofseqs, descr, descrlen);
    bl_fastaSetSequence(space, f, f->active_noofseqs, sequence, seqlen);

    f->minlen = (seqlen < f->minlen) ? seqlen : f->minlen;
    f->maxlen = (seqlen > f->maxlen) ? seqlen : f->maxlen;

    if (quality) {
        assert(n == 0 || bl_fastaHasQuality(f));
        bl_fastaAddQuality(space, f);
        bl_fastaSetQuality(space, f, f->active_noofseqs, quality, seqlen);
    } else {
        assert(!bl_fastaHasQuality(f));
    }

    f->active_noofseqs++;

    if (!f->hasIndex) {
        f->noofseqs++;
    }

    return;
}

/*--------------------------- bl_fastaCheckMateID ----------------------------
 *
 * @brief check if fasta description matches with mate pair description
 * @author Steve Hoffmann
 *
 */

unsigned char bl_fastaCheckMateID(fasta_t *f, Uint elem, char *mateid, Uint matelen) {
    char *id, *id2, *tok1, *tok2, *desc;
    Uint descrlen;
    unsigned char res;

    descrlen = f->seqs[elem]->descrlen;
    desc = f->seqs[elem]->description;

    id = ALLOCMEMORY(space, NULL, char, descrlen + 2);
    id2 = ALLOCMEMORY(space, NULL, char, matelen + 2);

    strcpy(id, desc);
    strcpy(id2, mateid);

    tok1 = strtok(id, "/");
    tok2 = strtok(id2, "/");
    res = (strcmp(tok1, tok2) == 0);

    if (!res) {
        FREEMEMORY(space, id);
        FREEMEMORY(space, id2);

        id = ALLOCMEMORY(space, NULL, char, descrlen + 2);
        id2 = ALLOCMEMORY(space, NULL, char, matelen + 2);

        strcpy(id, desc);
        strcpy(id2, mateid);

        tok1 = strtok(id, " ");
        tok2 = strtok(id2, " ");
        res = (strcmp(tok1, tok2) == 0);
    }

    FREEMEMORY(space, id);
    FREEMEMORY(space, id2);
    return res;
}

/*----------------------------- bl_fastaAddMate ------------------------------
 *
 * @brief add a new mate pair
 * @author Steve Hoffmann
 *
 */

void bl_fastaAddMate(void *space, fasta_t *f, char *descr, Uint descrlen, char *sequence,
                     Uint seqlen, Uint n) {
    bl_fastaSetMateDescription(f, n, descr, descrlen);
    if (f->checkid) {
        assert(bl_fastaCheckMateID(f, n, descr, descrlen));
    }
    bl_fastaSetSequence(space, f, n, sequence, seqlen);

    f->minlen = (seqlen < f->minlen) ? seqlen : f->minlen;
    f->maxlen = (seqlen > f->maxlen) ? seqlen : f->maxlen;

    f->active_noofmates++;
    return;
}

/*----------------------------- bl_fastxAddMate ------------------------------
 *
 * @brief add a new mate pair with or without quality information
 * @author Steve Hoffmann
 *
 */

void bl_fastxAddMate(void *space, fasta_t *f, char *descr, Uint descrlen, char *sequence,
                     char *quality, Uint seqlen, Uint n) {
    bl_fastaSetMateDescription(f, n, descr, descrlen);
    if (f->checkid) {
        if (bl_fastaCheckMateID(f, n, descr, descrlen) == 0) {
            NFO("The fasta/fastq IDs in both mate files do not match.\n", NULL);
            NFO("The first mismatch occurred at fastq number %u\n", n);
            NFO("Exiting.\n", NULL);
            exit(EXIT_FAILURE);
        }
    }
    bl_fastaSetSequence(space, f, n, sequence, seqlen);

    if (quality) {
        assert(bl_fastaHasQuality(f));
        bl_fastaSetQuality(space, f, n, quality, seqlen);
    } else {
        assert(!bl_fastaHasQuality(f));
    }

    f->minlen = (seqlen < f->minlen) ? seqlen : f->minlen;
    f->maxlen = (seqlen > f->maxlen) ? seqlen : f->maxlen;

    f->active_noofmates++;
    return;
}

/*------------------------------- bl_fastaChop -------------------------------
 *
 * @brief chop the fasta structure into pieces
 * @author Steve Hoffmann
 *
 */

fasta_t **bl_fastaChop(void *space, fasta_t *f, Uint pieces) {
    Uint size, r, i, j, offset = 0;
    fasta_t **chops;

    size = f->active_noofseqs / pieces;
    r = f->active_noofseqs - (size * pieces);
    assert((pieces * size) + r == f->active_noofseqs);

    chops = ALLOCMEMORY(space, NULL, fasta_t *, pieces);
    for (i = 0; i < pieces; i++) {
        chops[i] = bl_fastaInit(NULL);
        if (i < pieces - 1) {
            chops[i]->active_noofseqs = size;
            chops[i]->noofseqs = size;
        } else {
            chops[i]->active_noofseqs = size + r;
            chops[i]->noofseqs = size + r;
        }

        chops[i]->seqs = ALLOCMEMORY(space, NULL, CharSequence *, chops[i]->active_noofseqs);
        if (bl_fastaHasQuality(f)) {
            chops[i]->quals = ALLOCMEMORY(space, NULL, CharSequence *, chops[i]->active_noofseqs);
        }

        if (bl_fastaHasMate(f)) {
            chops[i]->matestart = &f->matestart[offset];
        }

        for (j = 0; j < chops[i]->active_noofseqs; j++) {
            chops[i]->seqs[j] = f->seqs[j + offset];

            if (bl_fastaHasQuality(f)) {
                chops[i]->quals[j] = f->quals[j + offset];
            }

            chops[i]->minlen = (f->seqs[j + offset]->length < chops[i]->minlen)
                                   ? f->seqs[j + offset]->length
                                   : chops[i]->minlen;
            chops[i]->maxlen = (f->seqs[j + offset]->length > chops[i]->maxlen)
                                   ? f->seqs[j + offset]->length
                                   : chops[i]->maxlen;
        }

        offset += chops[i]->active_noofseqs;
    }
    return chops;
}

/*------------------------------- bl_fastxCopy -------------------------------
 *
 * @brief copy a fasta set
 * @author Steve Hoffmann
 *
 */

fasta_t *bl_fastxCopy(void *space, fasta_t *f, Uint offset, Uint size) {
    Uint newsize, j;
    fasta_t *copy;

    assert(offset < f->noofseqs);

    copy = bl_fastaInit(NULL);

    if (offset + size > f->noofseqs)
        newsize = f->noofseqs - offset;
    else
        newsize = size;

    copy->active_noofseqs = newsize;
    copy->noofseqs = newsize;

    copy->seqs = ALLOCMEMORY(space, NULL, CharSequence *, copy->active_noofseqs);
    if (bl_fastaHasQuality(f)) {
        copy->quals = ALLOCMEMORY(space, NULL, CharSequence *, copy->active_noofseqs);
    }

    if (bl_fastaHasMate(f)) {
        copy->matestart = &f->matestart[offset];
    }

    for (j = 0; j < copy->active_noofseqs; j++) {
        copy->seqs[j] = f->seqs[j + offset];

        if (bl_fastaHasQuality(f)) {
            copy->quals[j] = f->quals[j + offset];
        }

        copy->minlen = (f->seqs[j + offset]->length < copy->minlen) ? f->seqs[j + offset]->length
                                                                    : copy->minlen;
        copy->maxlen = (f->seqs[j + offset]->length > copy->maxlen) ? f->seqs[j + offset]->length
                                                                    : copy->maxlen;
    }

    return copy;
}

/*---------------------------- bl_fastxChopIndex -----------------------------
 *
 * @brief chop the index to pieces
 * @author Steve Hoffmann
 *
 */

fasta_t **bl_fastxChopIndex(void *space, fasta_t *f, Uint pieces) {
    Uint i, j, offset, noofchunks, chunksperpiece, rest = 0, chunks, chunkoff = 0;

    fasta_t **chops;
    fastxseqindex_t *chunkindex = NULL;
    fastxseqindex_t *matechunkindex = NULL;

    assert(f->hasIndex);
    assert(pieces <= f->chunkindex->size);

    chops = ALLOCMEMORY(space, NULL, fasta_t *, pieces);
    noofchunks = f->chunkindex->size;
    chunksperpiece = noofchunks / pieces;
    rest = noofchunks - (pieces * chunksperpiece);

    for (i = 0; i < pieces; i++) {
        chops[i] = bl_fastaInit(space);
        memmove(chops[i], f, sizeof(fasta_t));
        chops[i]->chunkIsActive = 0;
        chops[i]->active_noofseqs = 0;
        chops[i]->active_noofmates = 0;
        chops[i]->seqs = NULL;
        chops[i]->quals = NULL;
        chops[i]->matestart = NULL;
        chops[i]->chunkIsActive = 0;
        chunks = chunksperpiece;

        if (rest > 0) {
            rest--;
            chunks++;
        }

        chunkindex = ALLOCMEMORY(space, NULL, fastxseqindex_t, 1);
        chunkindex->ap = ALLOCMEMORY(space, NULL, seqaccesspoint_t, chunks);
        memmove(chunkindex->ap, &f->chunkindex->ap[chunkoff], chunks * sizeof(seqaccesspoint_t));
        chunkindex->size = chunks;
        chunkindex->allocated = chunks;

        if (chunkoff) {
            offset = f->chunkindex->ap[chunkoff - 1].cumnoofseqs;
        } else {
            offset = 0;
        }

        for (j = 0; j < chunks; j++) {
            chunkindex->ap[j].cumnoofseqs -= offset;
        }

        if (f->matechunkindex) {
            matechunkindex = ALLOCMEMORY(space, NULL, fastxseqindex_t, 1);
            matechunkindex->ap = ALLOCMEMORY(space, NULL, seqaccesspoint_t, chunks);
            memmove(matechunkindex->ap, &f->matechunkindex->ap[chunkoff],
                    chunks * sizeof(seqaccesspoint_t));
            matechunkindex->size = chunks;
            matechunkindex->allocated = chunks;

            if (chunkoff) {
                offset = f->matechunkindex->ap[chunkoff - 1].cumnoofseqs;
            } else {
                offset = 0;
            }

            for (j = 0; j < chunks; j++) {
                matechunkindex->ap[j].cumnoofseqs -= offset;
            }
        }

        chops[i]->chunkindex = chunkindex;
        chops[i]->matechunkindex = matechunkindex;
        if (matechunkindex)
            assert(chunkindex->ap[chunks - 1].cumnoofseqs ==
                   matechunkindex->ap[chunks - 1].cumnoofseqs);
        chops[i]->noofseqs = chunkindex->ap[chunks - 1].cumnoofseqs;
        chunkoff += chunks;
    }

    return chops;
}

/*---------------------------- bl_fastxCopyIndex -----------------------------
 *
 * @brief copy with index
 * @author Steve Hoffmann
 *
 */

fasta_t *bl_fastxCopyIndex(void *space, fasta_t *f, Uint offset, Uint size) {
    Uint j, chunks, localoffset;
    fasta_t *copy;
    fastxseqindex_t *chunkindex = NULL;
    fastxseqindex_t *matechunkindex = NULL;

    assert(f->hasIndex);
    assert(offset < f->chunkindex->size);

    copy = bl_fastaInit(space);

    memmove(copy, f, sizeof(fasta_t));
    copy->chunkIsActive = 0;
    copy->active_noofseqs = 0;
    copy->active_noofmates = 0;
    copy->seqs = NULL;
    copy->quals = NULL;
    copy->matestart = NULL;
    copy->chunkIsActive = 0;

    if (offset + size > f->chunkindex->size)
        chunks = f->chunkindex->size - offset;
    else
        chunks = size;

    chunkindex = ALLOCMEMORY(space, NULL, fastxseqindex_t, 1);
    chunkindex->ap = ALLOCMEMORY(space, NULL, seqaccesspoint_t, chunks);
    memmove(chunkindex->ap, &f->chunkindex->ap[offset], chunks * sizeof(seqaccesspoint_t));
    chunkindex->size = chunks;
    chunkindex->allocated = chunks;

    if (offset) {
        localoffset = f->chunkindex->ap[offset - 1].cumnoofseqs;
    } else {
        localoffset = 0;
    }

    for (j = 0; j < chunks; j++) {
        chunkindex->ap[j].cumnoofseqs -= localoffset;
    }

    if (f->matechunkindex) {
        matechunkindex = ALLOCMEMORY(space, NULL, fastxseqindex_t, 1);
        matechunkindex->ap = ALLOCMEMORY(space, NULL, seqaccesspoint_t, chunks);
        memmove(matechunkindex->ap, &f->matechunkindex->ap[offset],
                chunks * sizeof(seqaccesspoint_t));
        matechunkindex->size = chunks;
        matechunkindex->allocated = chunks;

        if (offset) {
            localoffset = f->matechunkindex->ap[offset - 1].cumnoofseqs;
        } else {
            localoffset = 0;
        }

        for (j = 0; j < chunks; j++) {
            matechunkindex->ap[j].cumnoofseqs -= localoffset;
        }
    }

    copy->chunkindex = chunkindex;
    copy->matechunkindex = matechunkindex;
    if (matechunkindex)
        assert(chunkindex->ap[chunks - 1].cumnoofseqs ==
               matechunkindex->ap[chunks - 1].cumnoofseqs);
    copy->noofseqs = chunkindex->ap[chunks - 1].cumnoofseqs;

    return copy;
}

/*------------------------------- bl_fastxScan -------------------------------
 *
 * @brief scan fasta or fastq format
 * @author Steve Hoffmann
 *
 */

Uint bl_fastxScan(void *space, char *filename, struct access *index, off_t offset,
                  fastxfileindex_t *findex, Uint max, Uint *minlen, Uint *maxlen,
                  unsigned char *minq, unsigned char *maxq) {
    char ch;
    char idchar = 0;
    int ret = 0;
    off_t curseqoffset, lastindexoffset = 0;
    unsigned char desc = 0;
    unsigned char fastq = 0;
    unsigned char qualdesc = 0;
    unsigned char qual = 0;
    unsigned char seq = 0;
    unsigned char minqual = 255, maxqual = 0;
    FILE *fp;
    Uint seqlen = 0;
    Uint n = 0;
    Uint len = 0;

    struct gzidxfile *gzf = NULL;

    if (index) {
        fp = fopen(filename, "rb");

        if (fp == NULL) {
            fprintf(stderr, "Couldnt open %s for reading. Exit forced.\n", filename);
            exit(-1);
        }

        gzf = bl_initgzidxfile(fp, index, offset, MEDIUMCHUNK);

    } else {
        fp = fopen(filename, "r");
        if (fp == NULL) {
            fprintf(stderr, "Couldnt open %s for reading. Exit forced.\n", filename);
            exit(-1);
        }

        ret = fseeko(fp, offset, SEEK_SET);

        if (ret == -1) {
            MSG("fseeko failed. Exit forced.\n");
            exit(-1);
        }
    }

    if (findex->size + 2 >= findex->allocated) {
        findex->ap = ALLOCMEMORY(space, findex->ap, seqaccesspoint_t, findex->allocated + 11);
        findex->allocated += 11;
    }

    findex->ap[findex->size].offset = offset;
    findex->ap[findex->size].noofseqs = 0;
    findex->size++;

    lastindexoffset = offset;

    while ((ch = (index) ? bl_getgzidxc(gzf) : getc(fp)) != EOF) {
        // fprintf(stderr, "################## scanning with %ld\n", offset);

        if ((ch == '@' || ch == '>') && !idchar) {
            desc = 1;
            fastq = (ch == '@');
            idchar = ch;
        }

        if (ch == idchar && ((!fastq && len > 0) || (qual && len > 0 && len == seqlen))) {
            seq = 0;
            qual = 0;
            desc = 1;

            if (n == 0 || *minlen > len) {
                *minlen = len;
            }

            if (n == 0 || *maxlen < len) {
                *maxlen = len;
            }

            n++;
            seqlen = 0;
            len = 0;

            if (index) {
                curseqoffset = bl_ftellgzidx(gzf);
            } else {
                curseqoffset = ftello(fp);
                if (curseqoffset == -1) {
                    MSG("ftello failed. Exit forced.\n");
                    exit(-1);
                }
            }

            if (curseqoffset > lastindexoffset + SPAN || (max && n == max)) {
                if (findex->size + 2 >= findex->allocated) {
                    findex->ap =
                        ALLOCMEMORY(space, findex->ap, seqaccesspoint_t, findex->allocated + 11);
                    findex->allocated += 11;
                }

                findex->ap[findex->size].offset = curseqoffset - 1;
                findex->ap[findex->size].noofseqs = n - 1;
                findex->size++;

                lastindexoffset = curseqoffset;
            }

            if (max && n == max) {
                break;
            }
        }

        if (qual && len > seqlen) {
            NFO("fastq error: qual string > nt string: %d\n", n);
            exit(-1);
        }

        if (qual) {
            if (ch < minqual) minqual = ch;
            if (ch > maxqual) maxqual = ch;
        }

        if (fastq && ch == '+' && seq && len > 0) {
            seq = 0;
            qualdesc = 1;
            seqlen = len;
            len = 0;
        }

        if (!desc && !qualdesc && ch == '\n') {
            /*do nothing.*/
        } else {
            if (desc && ch == '\n') {
                len = 0;
                desc = 0;
                seq = 1;
            } else if (fastq && qualdesc && ch == '\n') {
                len = 0;
                qualdesc = 0;
                qual = 1;
            } else {
                if (ch == '\r') continue;
                len++;
            }
        }
    }

    if ((!fastq && seq && len > 0) || (fastq && qual && len > 0 && len == seqlen)) {
        if (n == 0 || *minlen > len) {
            *minlen = len;
        }
        if (n == 0 || *maxlen < len) {
            *maxlen = len;
        }
        n++;
    }

    *minq = minqual;
    *maxq = maxqual;

    fclose(fp);
    if (index) bl_destructgzidxfile(gzf);
    FREEMEMORY(space, gzf);

    return n;
}

/*--------------------------- bl_fastxChunkIndex ----------------------------
 *
 * @brief adjust offset-index for n seqs to chunks with approx k seqs
 * @author Steve Hoffmann
 *
 */

fastxseqindex_t *bl_fastxChunkIndex(void *space, char **filenames, struct access **gzindex,
                                    fastxfileindex_t **findex, Uint *n, Uint nooffiles, Uint total,
                                    Uint k) {
    Uint i = 0, j, chunks, cur = 0, skip, chunksize, curchunksize = 0, rest = 0, nseqs = 0,
         nseqsfile = 0, minlen = 0, maxlen = 0, off = 0;
    unsigned char minqual = 255, maxqual = 0;
    // Uint l;
    // char ch;
    //  FILE *fp=NULL;

    fastxseqindex_t *idx;
    fastxfileindex_t *tmp;

    chunks = (k > total) ? 1 : total / k;
    idx = bl_fastxInitSeqIndex(space, chunks + 1000);

    // fprintf(stderr, "start fastxchunkindex\n");

    for (j = 0; j < nooffiles; j++) {
        chunks = (k > n[j]) ? 1 : n[j] / k;
        chunksize = (k > n[j]) ? n[j] : k + ((n[j] % k) / chunks);
        rest = (n[j] % k) % chunks;

        curchunksize = chunksize;

        if (rest > 0) {
            rest--;
            curchunksize++;
        }

        nseqs += curchunksize;
        idx->ap[idx->size].noofseqs = curchunksize;
        idx->ap[idx->size].cumnoofseqs = nseqs;
        idx->ap[idx->size].offset = 0;
        idx->ap[idx->size].fileid = j;
        idx->size++;

        tmp = bl_fastxInitFileIndex(space, chunks);
        nseqsfile = curchunksize;

        for (i = 1; i < chunks; i++) {
            cur = 0;
            while (cur + 1 < findex[j]->size && findex[j]->ap[cur + 1].noofseqs < nseqsfile) cur++;

            if (findex[j]->ap[cur].noofseqs) {
                off = findex[j]->ap[cur].noofseqs + 1;
            } else {
                off = 0;
            }

            skip = nseqsfile - off;

            //   fprintf(stderr, "chunk %d in file %d (curchunksize %d): skiping %d (= nseqsfile:%d
            //   - findex.noofseqs:%d)\n",
            //       i, j, curchunksize, skip, nseqsfile, findex[j]->ap[cur].noofseqs);
            //   fprintf(stderr, "at offset: %lu\n", findex[j]->ap[cur].offset);
            //   fprintf(stderr, "@%d\n", nseqsfile);

            if (skip) {
                if (gzindex) {
                    bl_fastxScan(space, filenames[j], gzindex[j], findex[j]->ap[cur].offset, tmp,
                                 skip, &minlen, &maxlen, &maxqual, &minqual);
                } else {
                    bl_fastxScan(space, filenames[j], NULL, findex[j]->ap[cur].offset, tmp, skip,
                                 &minlen, &maxlen, &maxqual, &minqual);
                }
            } else {
                if (tmp->size + 2 >= tmp->allocated) {
                    tmp->ap = ALLOCMEMORY(space, tmp->ap, seqaccesspoint_t, tmp->allocated + 11);
                    tmp->allocated += 11;
                }

                tmp->ap[tmp->size].offset = findex[j]->ap[cur].offset;
                tmp->ap[tmp->size].noofseqs = curchunksize;
                tmp->size++;
            }

            /*===================

            fp = fopen(filenames[j],"r");
            if (fp == NULL) {
              fprintf(stderr, "Couldnt open %s for reading. Exit forced.\n", filenames[j]);
              exit(-1);
            }

            fseeko(fp, tmp->ap[tmp->size-1].offset, SEEK_SET);

            for(l=0; l < 50; l++) {
              ch = getc(fp);
              fprintf(stderr, "%c", ch);
            }

            fprintf(stderr, "\n");
            fclose(fp);
            ===================== */

            curchunksize = chunksize;

            if (rest > 0) {
                rest--;
                curchunksize++;
            }

            nseqs += curchunksize;
            nseqsfile += curchunksize;
            idx->ap[idx->size].noofseqs = curchunksize;
            idx->ap[idx->size].cumnoofseqs = nseqs;
            idx->ap[idx->size].offset = tmp->ap[tmp->size - 1].offset;
            idx->ap[idx->size].fileid = j;
            idx->size++;
        }

        FREEMEMORY(space, tmp->ap);
        FREEMEMORY(space, tmp);
    }

    return idx;
}

/*----------------------------- bl_fastxGetChunk -----------------------------
 *
 * @brief return chunknumber for sequence k
 * @author Steve Hoffmann
 *
 */

int bl_fastxGetChunk(fasta_t *fasta, Uint k) {
    Uint chunksize, i, chunks;

    if (k >= fasta->noofseqs || !fasta->chunkindex) return -1;

    chunks = fasta->chunkindex->size;
    chunksize = fasta->chunkindex->ap[0].noofseqs;
    i = k / chunksize;
    if (i >= chunks) i = 0;

    while (i < chunks && fasta->chunkindex->ap[i].cumnoofseqs <= k) i++;
    while (i > 0 && fasta->chunkindex->ap[i - 1].cumnoofseqs > k) i--;

    if (fasta->chunkindex->ap[i].cumnoofseqs <= k ||
        (i > 0 && fasta->chunkindex->ap[i - 1].cumnoofseqs > k)) {
        DBG("chunk not found: chunks:%d, i:%d, idx[i]:%d, idx[i-1]:%d, k:%d\n", chunks, i,
            fasta->chunkindex->ap[i].cumnoofseqs, fasta->chunkindex->ap[i - 1].cumnoofseqs, k);

        i = 0;
        DBG("list: chunks:%d, i:%d, idx[i]:%d, idx[i-1]:%d, k:%d, fid:%d\n", chunks, i,
            fasta->chunkindex->ap[i].cumnoofseqs, 0, k, fasta->chunkindex->ap[i].fileid);

        for (i = 1; i < chunks; i++) {
            DBG("list: chunks:%d, i:%d, idx[i]:%d, idx[i-1]:%d, k:%d, fid:%d\n", chunks, i,
                fasta->chunkindex->ap[i].cumnoofseqs, fasta->chunkindex->ap[i - 1].cumnoofseqs, k,
                fasta->chunkindex->ap[i].fileid);
        }

        exit(-1);
        return -1;
    }

    return i;
}

/*------------------------- bl_fastxGetChunkElem -------------------------
 *
 * @brief load chunk if necessary and return position of sequence k
 * @author Steve Hoffmann
 *
 */

Uint bl_fastxGetChunkElem(void *space, fasta_t *f, Uint k) {
    Uint off = 0, fid;
    int cur;

    cur = f->curchunk;

    if (f->chunkIsActive) {
        if (f->chunkindex->ap[cur].cumnoofseqs > k &&
            (cur == 0 || f->chunkindex->ap[cur - 1].cumnoofseqs <= k)) {
            if (cur > 0) off = f->chunkindex->ap[cur - 1].cumnoofseqs;
            assert(k >= off);
            return k - off;
        }
        bl_fastxDestructSequence(space, f);
        f->active_noofseqs = 0;
    }

    cur = bl_fastxGetChunk(f, k);

    assert(cur > -1);
    fid = f->chunkindex->ap[cur].fileid;

    if (f->gzip) {
        f = bl_fastxgzRead(space, f, f->filenames[fid], f->gzindex[fid], f->upper, f->lower,
                           f->chunkindex->ap[cur].offset, 0, f->chunkindex->ap[cur].noofseqs,
                           bl_fastxAdd);

        if (f->hasMates) {
            f = bl_fastxgzRead(space, f, f->matefilenames[fid], f->mategzindex[fid], f->upper,
                               f->lower, f->matechunkindex->ap[cur].offset, 0,
                               f->matechunkindex->ap[cur].noofseqs, bl_fastxAddMate);
        }
    } else {
        f = bl_fastxRead(space, f, f->filenames[fid], f->upper, f->lower,
                         f->chunkindex->ap[cur].offset, 0, f->chunkindex->ap[cur].noofseqs,
                         // f->gzip, f->index[fid],
                         bl_fastxAdd);
        if (f->hasMates) {
            f = bl_fastxRead(space, f, f->matefilenames[fid], f->upper, f->lower,
                             f->matechunkindex->ap[cur].offset, 0,
                             f->matechunkindex->ap[cur].noofseqs,
                             // f->gzip, f->index[fid],
                             bl_fastxAddMate);
        }
    }

    if (cur > 0) off = f->chunkindex->ap[cur - 1].cumnoofseqs;

    f->chunkIsActive = 1;
    f->curchunk = cur;

    assert(k >= off);
    return k - off;
}

/*------------------------------ bl_fastxIndex ------------------------------
 *
 * @brief  build index of fasta or fastq file
 * @author Steve Hoffmann
 *
 */

fasta_t *bl_fastxIndex(void *space, fasta_t *set, char **filenames, Uint nooffiles,
                       unsigned char isMate, unsigned char gzip, Uint pieces) {
    int i, len = 0;
    Uint total = 0, *ftotal = NULL;
    Uint chunksize = 10000, minlen = 0, maxlen = 0, seqsperpiece, noofchunks, rest;
    unsigned char minqual = 0, maxqual = 0;

    struct access **gzindex = NULL;
    fastxfileindex_t **findex = NULL;
    fastxseqindex_t *chunkindex;

    if (gzip) {
        gzindex = ALLOCMEMORY(space, NULL, struct access *, nooffiles);
        for (i = 0; i < nooffiles; i++) {
            gzindex[i] = bl_gzGetIndex(filenames[i]);
        }
    }

    findex = ALLOCMEMORY(space, NULL, fastxfileindex_t *, nooffiles);
    ftotal = ALLOCMEMORY(space, NULL, Uint, nooffiles);

    for (i = 0; i < nooffiles; i++) {
        if (gzindex) {
            len = gzindex[i]->size;
        } else {
            len = 0;
        }

        findex[i] = bl_fastxInitFileIndex(space, len);

        if (gzip)
            ftotal[i] = bl_fastxScan(space, filenames[i], gzindex[i], 0, findex[i], 0, &minlen,
                                     &maxlen, &minqual, &maxqual);
        else
            ftotal[i] = bl_fastxScan(space, filenames[i], NULL, 0, findex[i], 0, &minlen, &maxlen,
                                     &minqual, &maxqual);

        total += ftotal[i];
    }

    seqsperpiece = MAX(1, total / pieces);

    if (seqsperpiece < chunksize) {
        chunksize = seqsperpiece;
    } else {
        noofchunks = seqsperpiece / chunksize;
        rest = seqsperpiece - (noofchunks * chunksize);
        chunksize += rest / noofchunks;
    }

    chunkindex =
        bl_fastxChunkIndex(space, filenames, gzindex, findex, ftotal, nooffiles, total, chunksize);

    if (!isMate) {
        set = bl_fastaInit(space);
        set->hasIndex = 1;
        set->active_noofseqs = 0;

        set->chunkIsActive = 0;
        set->chunkindex = chunkindex;
        set->findex = findex;
        set->gzindex = gzindex;
        set->filenames = filenames;
        set->gzip = gzip;
        set->noofseqs = total;
        set->upper = 1;
        set->lower = 0;
        set->nooffiles = nooffiles;
        set->filetotal = ftotal;
        set->minlen = minlen;
        set->maxlen = maxlen;
        set->minqual = minqual;
        set->maxqual = maxqual;

    } else {
        if (set == NULL || set->nooffiles != nooffiles || set->noofseqs != total ||
            set->chunkindex == NULL || set->chunkindex->size != chunkindex->size) {
            MSG("1: Reading mates failed: mate and query files differ in size!\n");
            NFO("set->nooffiles %d = %d nooffiles\n", set->nooffiles, nooffiles);
            NFO("set->noofseqs %d = %d noofseqs\n", set->noofseqs, total);
            NFO("set->chunkindex->size %d = %d chunkindex->size", set->chunkindex->size,
                chunkindex->size);
            exit(-1);
        }

        for (i = 0; i < nooffiles; i++) {
            if (ftotal[i] != set->filetotal[i]) {
                MSG("2: Reading mates failed: mate and query files differ in size!\n");
                exit(-1);
            }
        }

        for (i = 0; i < chunkindex->size; i++) {
            if (chunkindex->ap[i].noofseqs != set->chunkindex->ap[i].noofseqs) {
                MSG("3: Reading mates failed: mate and query files differ in size!\n");
                exit(-1);
            }
        }

        FREEMEMORY(space, ftotal);

        set->matefilenames = filenames;
        set->matechunkindex = chunkindex;
        set->matefindex = findex;
        set->mategzindex = gzindex;
        set->hasMates = 1;

        if (set->maxlen < maxlen) {
            set->maxlen = maxlen;
        }

        if (set->minlen > minlen) {
            set->minlen = minlen;
        }

        if (set->minqual > minqual) {
            set->minqual = minqual;
        }

        if (set->maxqual < maxqual) {
            set->maxqual = maxqual;
        }
    }

    set->active_noofseqs = 0;
    set->chunkIsActive = 0;

    return set;
}

/*------------------------------- bl_fastxRead -------------------------------
 *
 * @brief read fasta or fastq format
 * @author Steve Hoffmann
 *
 */

fasta_t *bl_fastxRead(void *space, fasta_t *fasta, char *filename, unsigned char upper,
                      unsigned char lower, off_t offset, Uint startseq, Uint lastseq,
                      // unsigned char gzip, struct access *index,
                      void (*handler)(void *, fasta_t *, char *, Uint, char *, char *, Uint,
                                      Uint)) {
    FILE *fp;
    char ch;
    char *buffer;
    char *descrbuffer = NULL;
    char *seqbuffer = NULL;
    char *qualbuffer = NULL;
    char idchar = 0;
    int ret = 0;
    unsigned char desc = 0;
    unsigned char fastq = 0;
    unsigned char qualdesc = 0;
    unsigned char qual = 0;
    unsigned char seq = 0;
    unsigned char gzip = 0;
    struct gzidxfile *gzf = NULL;
    struct access *index = NULL;

    Uint descrlength = 0;
    Uint seqlen = 0;
    Uint buffersize = MAXBUFFERSIZE;
    Uint n = startseq;
    Uint len = 0;

    NFO("Reading file %s\n", filename);

    buffer = ALLOCMEMORY(space, NULL, char, buffersize);
    if (fasta == NULL) fasta = bl_fastaInit(space);

    if (gzip) {
        fp = fopen(filename, "rb");
        gzf = bl_initgzidxfile(fp, index, offset, MEDIUMCHUNK);
    } else {
        fp = fopen(filename, "r");
    }

    if (fp == NULL) {
        NFO("fastxRead: Couldn't open file '%s': %d. Exit forced.\n", filename, errno);
        exit(-1);
    }

    if (offset > 0) {
        ret = fseeko(fp, offset, SEEK_SET);
        if (ret == -1) {
            NFO("fastxRead: fseeko failed for file %s. Exit forced.\n", filename);
            exit(-1);
        }
    }

    while ((ch = (gzip) ? bl_getgzidxc(gzf) : getc(fp)) != EOF) {
        if (len == buffersize - 1) {
            buffersize = 2 * buffersize + 1;
            buffer = ALLOCMEMORY(space, buffer, char, buffersize);
        }

        if ((ch == '@' || ch == '>') && !idchar) {
            desc = 1;
            fastq = (ch == '@');
            idchar = ch;
        }

        if (fastq && ch == '+' && seq && len > 0) {
            buffer = ALLOCMEMORY(space, buffer, char, len + 1);
            buffer[len] = '\0';
            seq = 0;
            qualdesc = 1;

            seqbuffer = buffer;
            seqlen = len;
            len = 0;

            buffersize = MAXBUFFERSIZE;
            buffer = ALLOCMEMORY(space, NULL, char, buffersize);
        }

        if (qual && len > seqlen) {
            NFO("fastq format error: quality string longer than nt string: %s\n", descrbuffer);
            exit(-1);
        }

        assert(!qual || len <= seqlen);  //  v-- && seqlen > 0 produces segfault!!
        if (ch == idchar && ((!fastq && len > 0) || (qual && len > 0 && len == seqlen))) {
            if (lastseq && n >= lastseq - 1) {
                break;
            }

            seq = 0;
            qual = 0;
            desc = 1;

            buffer = ALLOCMEMORY(space, buffer, char, len + 1);
            buffer[len] = '\0';

            assert(!fastq || seqbuffer);

            if (!seqbuffer) {
                seqbuffer = buffer;
            } else {
                qualbuffer = buffer;
            }

            handler(space, fasta, descrbuffer, descrlength, seqbuffer, qualbuffer, len, n);
            n++;

            descrlength = 0;
            descrbuffer = NULL;
            seqlen = 0;
            seqbuffer = NULL;
            qualbuffer = NULL;

            len = 0;
            buffersize = MAXBUFFERSIZE;
            buffer = ALLOCMEMORY(space, NULL, char, buffersize);
        }

        if (!desc && !qualdesc && ch == '\n') {
            /*do nothing.*/
        } else {
            if (desc && ch == '\n') {
                buffer = ALLOCMEMORY(space, buffer, char, len + 1);
                buffer[len] = '\0';

                descrbuffer = buffer;
                descrlength = len;

                len = 0;
                buffersize = MAXBUFFERSIZE;
                buffer = ALLOCMEMORY(space, NULL, char, buffersize);
                desc = 0;
                seq = 1;
            } else if (fastq && qualdesc && ch == '\n') {
                FREEMEMORY(space, buffer);
                len = 0;
                buffersize = MAXBUFFERSIZE;
                buffer = ALLOCMEMORY(space, NULL, char, buffersize);
                qualdesc = 0;
                qual = 1;
            } else {
                if (ch == '\r') continue;
                len++;
                if (upper && !desc && !qualdesc && !qual) {
                    buffer[len - 1] = (char)toupper((int)ch);
                } else if (lower && !desc && !qualdesc && !qual) {
                    buffer[len - 1] = (char)tolower((int)ch);
                } else {
                    buffer[len - 1] = (char)ch;
                }
            }
        }
    }

    if ((!fastq && len > 0) || (qual && len > 0 && len == seqlen)) {
        buffer = ALLOCMEMORY(space, buffer, char, len + 1);
        buffer[len] = '\0';

        assert(!fastq || seqbuffer);

        if (!seqbuffer) {
            seqbuffer = buffer;
        } else {
            qualbuffer = buffer;
        }

        if (descrbuffer == NULL) DBG("empty descr buffer after loop n=%d\n", n);

        handler(space, fasta, descrbuffer, descrlength, seqbuffer, qualbuffer, len, n);
    }

    if (gzip) {
        bl_destructgzidxfile(gzf);
        FREEMEMORY(space, gzf);
    }

    fclose(fp);
    return fasta;
}

/*------------------------------ bl_fastgzRead -------------------------------
 *
 * @brief randomly access zipped fasta gz using the index
 * to populate fasta struct
 * @author Steve Hoffmann
 *
 */

fasta_t *bl_fastxgzRead(void *space, fasta_t *fasta, char *filename, struct access *idx,
                        unsigned char upper, unsigned char lower, off_t offset, Uint startseq,
                        Uint lastseq,
                        void (*handler)(void *, fasta_t *, char *, Uint, char *, char *, Uint,
                                        Uint)) {
    FILE *fp = NULL;
    char ch = 0;
    char *buffer;
    // char *check;
    unsigned char *gzbuffer;
    char *descrbuffer = NULL;
    char *seqbuffer = NULL;
    char *qualbuffer = NULL;
    char idchar = 0;
    int nb;
    int gzbufpos = 0;
    unsigned char desc = 0;
    unsigned char fastq = 0;
    unsigned char qualdesc = 0;
    unsigned char qual = 0;
    unsigned char seq = 0;
    unsigned char finalize = 0;
    int ret;
    Uint descrlength = 0;
    Uint seqlen = 0;
    Uint buffersize = MAXBUFFERSIZE;
    Uint n = startseq;
    Uint len = 0;

    if (fasta == NULL) fasta = bl_fastaInit(space);
    fp = fopen(filename, "rb");

    if (fp == NULL) {
        DBG("fastxgzRead: Couldn't open file '%s': %s. Exit forced.\n", filename, strerror(errno));
        exit(-1);
    }

    /* access block and read it entirely */

    buffer = ALLOCMEMORY(space, NULL, char, buffersize);
    gzbuffer = ALLOCMEMORY(space, NULL, char, SPAN);
    //  check = ALLOCMEMORY(space, NULL, char, 1001);

    nb = bl_zranExtract(fp, idx, offset, gzbuffer, SPAN);

    /*
    memset(check, 0, 1001);
    memmove(check, gzbuffer, MIN(nb,1000));
    fprintf(stderr, "%s\n-------\n", check);
    */

    if (nb < 0) {
        DBG("extraction failed (%s)\n", nb == Z_MEM_ERROR ? "out of memory" : "input corrupted");
        fclose(fp);
        exit(-1);
    }

    do {
        finalize = 0;
        while (gzbufpos < nb) {
            ch = gzbuffer[gzbufpos++];

            if (len == buffersize - 1) {
                buffersize = 2 * buffersize + 1;
                buffer = ALLOCMEMORY(space, buffer, char, buffersize);
            }

            if ((ch == '@' || ch == '>') && !idchar) {
                desc = 1;
                fastq = (ch == '@');
                idchar = ch;
                len = 0;
                // fprintf(stderr, "found idchar: %c\n", idchar);
            }

            /*fastq only: store the sequence and read the quality string*/
            if (fastq && ch == '+' && seq && len > 0) {
                buffer = ALLOCMEMORY(space, buffer, char, len + 1);
                buffer[len] = '\0';
                seq = 0;
                qualdesc = 1;

                seqbuffer = buffer;
                seqlen = len;
                len = 0;

                buffersize = MAXBUFFERSIZE;
                buffer = ALLOCMEMORY(space, NULL, char, buffersize);
            }

            if (qual && len > seqlen) {
                DBG("%s: qual longer than nt string (n=%d). Exit.\n", descrbuffer, n);
                exit(-1);
            }

            /*reading a new id*/
            if (ch == idchar && ((!fastq && len > 0) || (qual && len > 0 && len == seqlen))) {
                seq = 0;
                qual = 0;
                desc = 1;

                if (lastseq && n >= lastseq - 1) {
                    finalize = 1;
                    break;
                }

                buffer = ALLOCMEMORY(space, buffer, char, len + 1);
                buffer[len] = '\0';
                assert(!fastq || seqbuffer);

                /*switch fasta:fastq*/
                if (!seqbuffer) {
                    seqbuffer = buffer;
                } else {
                    qualbuffer = buffer;
                }

                if (descrbuffer == NULL) DBG("empty descr buffer in loop n=%d\n", n);

                handler(space, fasta, descrbuffer, descrlength, seqbuffer, qualbuffer, len, n);

                n++;

                descrlength = 0;
                descrbuffer = NULL;
                seqlen = 0;
                seqbuffer = NULL;
                qualbuffer = NULL;

                len = 0;
                buffersize = MAXBUFFERSIZE;
                buffer = ALLOCMEMORY(space, NULL, char, buffersize);
            }

            if (!desc && !qualdesc && ch == '\n') {
                /*do nothing.*/
            } else {
                if (desc && ch == '\n') {
                    /*store description*/
                    buffer = ALLOCMEMORY(space, buffer, char, len + 1);
                    buffer[len] = '\0';

                    descrbuffer = buffer;
                    descrlength = len;

                    len = 0;
                    buffersize = MAXBUFFERSIZE;
                    buffer = ALLOCMEMORY(space, NULL, char, buffersize);
                    desc = 0;
                    seq = 1;
                } else if (fastq && qualdesc && ch == '\n') {
                    /*toss quality string description*/
                    FREEMEMORY(space, buffer);
                    len = 0;
                    buffersize = MAXBUFFERSIZE;
                    buffer = ALLOCMEMORY(space, NULL, char, buffersize);
                    qualdesc = 0;
                    qual = 1;
                } else {
                    if (ch == '\r') continue;
                    len++;
                    if (upper && !desc && !qualdesc && !qual) {
                        buffer[len - 1] = (char)toupper((int)ch);
                    } else if (lower && !desc && !qualdesc && !qual) {
                        buffer[len - 1] = (char)tolower((int)ch);
                    } else {
                        buffer[len - 1] = (char)ch;
                    }
                }
            }
        }

        FREEMEMORY(space, gzbuffer);

        if (!finalize) {
            /*read next chunk*/
            offset += nb;
            gzbuffer = ALLOCMEMORY(space, NULL, char, SPAN);
            nb = bl_zranExtract(fp, idx, offset, gzbuffer, SPAN);

            /*
            memmove(check, gzbuffer, MIN(nb,1000));
            fprintf(stderr, "%s\n-------\n", check);
            */

            if (nb < 0) {
                DBG("extraction failed (%s)\n",
                    nb == Z_MEM_ERROR ? "out of memory" : "input corrupted");
                exit(-1);
            } else if (nb == 0) {
                FREEMEMORY(space, gzbuffer);
                break;
            }

            gzbufpos = 0;
        } else {
            break;
        }

    } while (1);

    buffer = ALLOCMEMORY(space, buffer, char, len + 1);
    buffer[len] = '\0';

    assert(!fastq || seqbuffer);

    if (!seqbuffer) {
        seqbuffer = buffer;
    } else {
        qualbuffer = buffer;
    }

    if (descrbuffer == NULL) DBG("empty descr buffer after loop n=%d\n", n);

    handler(space, fasta, descrbuffer, descrlength, seqbuffer, qualbuffer, len, n);

    ret = fclose(fp);
    if (ret == EOF) {
        fprintf(stderr, "Couldnt close file!\n");
        exit(-1);
    }
    return fasta;
}

/*------------------------------ bl_fastxGetSet ------------------------------
 *
 * @brief read a set of fasta or fastq files
 * @author Steve Hoffmann
 *
 */

fasta_t *bl_fastxGetSet(void *space, char **filenames, unsigned int nooffiles, unsigned char upper,
                        unsigned char lower, unsigned char index, Uint pieces) {
    Uint i, prefixlen, n = 0;
    unsigned char gzip = 2;
    struct access *gzindex;
    fasta_t *set;

    for (i = 0; i < nooffiles; i++) {
        prefixlen = bl_fileprefixlen(filenames[i]);

        if (strncmp(&filenames[i][prefixlen], ".gz", 3) == 0 ||
            strncmp(&filenames[i][prefixlen], ".gzip", 5) == 0 ||
            strncmp(&filenames[i][prefixlen], ".bgz", 4) == 0 ||
            strncmp(&filenames[i][prefixlen], ".bgzip", 6) == 0) {
            if (gzip == 2 || gzip == 1) {
                gzip = 1;
            } else {
                MSG("Provide fastx files either gzipped xor plain. Exit forced.\n");
                exit(-1);
            }
        } else {
            if (gzip == 2 || gzip == 0) {
                gzip = 0;
            } else {
                MSG("Provide fastx files either gzip'd xor plain. Exit forced.\n");
                exit(-1);
            }
        }
    }

    if (!index) {
        set = bl_fastaInit(space);
        set->filenames = filenames;
        set->nooffiles = nooffiles;
        set->gzip = gzip;
        set->upper = upper;
        set->lower = lower;

        for (i = 0; i < nooffiles; i++) {
            if (!gzip) {
                set =
                    bl_fastxRead(space, set, filenames[i], upper, lower, 0, n, 0,  // gzip, gzindex,
                                 bl_fastxAdd);

                n = set->active_noofseqs;

            } else {
                gzindex = bl_gzGetIndex(filenames[i]);
                set = bl_fastxgzRead(space, set, filenames[i], gzindex, upper, lower, 0, n, 0,
                                     bl_fastxAdd);

                n = set->active_noofseqs;

                FREEMEMORY(space, gzindex->list);
                FREEMEMORY(space, gzindex);
            }
        }
    } else {
        set = bl_fastxIndex(space, NULL, filenames, nooffiles, 0, gzip, pieces);
    }

    return set;
}

/*---------------------------- bl_fastxGetMateSet ----------------------------
 *
 * @brief get a set of fastq or fasta mate sequence files
 * @author Steve Hoffmann
 *
 */

fasta_t *bl_fastxGetMateSet(void *space, fasta_t *set, char **filenames, unsigned int nooffiles,
                            unsigned char upper, unsigned lower, unsigned char index, Uint pieces) {
    Uint i, prefixlen, n = 0;
    unsigned char gzip = 2;
    struct access *gzindex;

    assert(set != NULL);

    for (i = 0; i < nooffiles; i++) {
        prefixlen = bl_fileprefixlen(filenames[i]);

        if (strncmp(&filenames[i][prefixlen], ".gz", 3) == 0 ||
            strncmp(&filenames[i][prefixlen], ".gzip", 5) == 0 ||
            strncmp(&filenames[i][prefixlen], ".bgz", 4) == 0 ||
            strncmp(&filenames[i][prefixlen], ".bgzip", 6) == 0) {
            if (gzip == 2 || gzip == 1) {
                gzip = 1;
            } else {
                MSG("Provide fastx files either gzipped xor txt. Exit forced.\n");
                exit(-1);
            }
        } else {
            if (gzip == 2 || gzip == 0) {
                gzip = 0;
            } else {
                MSG("Provide fastx files either gzip'd xor txt. Exit forced.\n");
                exit(-1);
            }
        }
    }

    if (!index) {
        set = bl_fastaInit(space);
        for (i = 0; i < nooffiles; i++) {
            if (!gzip) {
                set =
                    bl_fastxRead(space, set, filenames[i], upper, lower, 0, n, 0, bl_fastxAddMate);

                n = set->active_noofmates;

            } else {
                gzindex = bl_gzGetIndex(filenames[i]);
                set = bl_fastxgzRead(space, set, filenames[i], gzindex, upper, lower, 0, n, 0,
                                     bl_fastxAddMate);

                n = set->active_noofmates;

                FREEMEMORY(space, gzindex->list);
                FREEMEMORY(space, gzindex);
            }
        }
    } else {
        set = bl_fastxIndex(space, set, filenames, nooffiles, 1, gzip, pieces);
    }

    return set;
}

/*------------------------------ bl_fastxIDcmp -------------------------------
 *
 * @brief compare fasta ids
 * @author Steve Hoffmann
 *
 */

int bl_fastxIDcmp(char *a, char *b) {
    char *desc;
    char *temp = "chr";
    int i, j, alen, blen, tlen;
    alen = strlen(a);
    blen = strlen(b);
    tlen = strlen(temp);

    if (strcmp(a, b) == 0) return 0;

    for (i = 0; i < alen; i++) {
        if (isspace((int)a[i])) break;
    }

    if (blen <= i && strncmp(a, b, i) == 0) {
        return 0;
    }

    for (j = 0; j < blen; j++) {
        if (isspace((int)b[j])) break;
    }

    desc = calloc(i + tlen + 1, sizeof(char));
    memmove(desc, temp, tlen);
    memmove(&desc[tlen], a, i);
    desc[tlen + i] = 0;

    if (blen >= tlen + i && strncmp(desc, b, tlen + i) == 0) {
        free(desc);
        return 0;
    }

    free(desc);

    desc = calloc(j + tlen + 1, sizeof(char));
    memmove(desc, temp, tlen);
    memmove(&desc[tlen], b, j);
    desc[tlen + j] = 0;

    if (alen >= tlen + j && strncmp(a, desc, tlen + j) == 0) {
        free(desc);
        return 0;
    }

    free(desc);

    return -1;
}

/*----------------------------- bl_fastxFindIdx ------------------------------
 *
 * @brief return the index for a fasta with name id or description
 * @author Steve Hoffmann
 *
 */

Uint bl_fastxFindIDIdx(char *id, fasta_t *set) {
    Uint i, j;
    char *desc;

    for (i = 0; i < set->noofseqs; i++) {
        desc = bl_fastaGetDescription(set, i);
        // fprintf(stderr, "desc:%s, id:%s\n", desc, id);

        if (strcmp(id, desc) == 0) {
            break;
        }
        // fprintf(stderr, "not break\n");

        for (j = 0; j < bl_fastaGetDescriptionLength(set, i); j++) {
            if (isspace((int)desc[j])) break;
        }

        if (strlen(id) <= j && strncmp(id, desc, j) == 0) {
            break;
        }

        // fprintf(stderr, "not break 2\n");
    }

    return (i < set->noofseqs) ? i : -1;
}

/*------------------------ bl_fastxDestructChunkIndex ------------------------
 *
 * @brief destruct the chunk index
 * @author Steve Hoffmann
 *
 */

void bl_fastxDestructChunkIndex(void *space, fasta_t *f) {
    if (f->chunkindex) {
        FREEMEMORY(space, f->chunkindex->ap);
        FREEMEMORY(space, f->chunkindex);

        if (f->hasMates) {
            FREEMEMORY(space, f->matechunkindex->ap);
            FREEMEMORY(space, f->matechunkindex);
        }

        f->matechunkindex = NULL;
        f->chunkindex = NULL;
    }

    return;
}

/*------------------------- bl_fastaDestructIndex -------------------------
 *
 * @brief free sequence index
 * @author Steve Hoffmann
 *
 */

void bl_fastxDestructIndex(void *space, fasta_t *f) {
    Uint i;

    bl_fastxDestructChunkIndex(space, f);

    if (f->filetotal) {
        FREEMEMORY(space, f->filetotal);
    }

    if (f->findex) {
        for (i = 0; i < f->nooffiles; i++) {
            FREEMEMORY(space, f->findex[i]->ap);
            FREEMEMORY(space, f->findex[i]);
            if (f->hasMates) {
                FREEMEMORY(space, f->matefindex[i]->ap);
                FREEMEMORY(space, f->matefindex[i]);
            }
        }
        FREEMEMORY(space, f->findex);
        if (f->hasMates) {
            FREEMEMORY(space, f->matefindex);
        }
        f->matefindex = NULL;
        f->findex = NULL;
    }

    f->curchunk = 0;
    f->noofseqs = 0;
    f->chunkIsActive = 0;

    if (f->gzip) {
        for (i = 0; i < f->nooffiles; i++) {
            if (f->gzindex) {
                FREEMEMORY(space, f->gzindex[i]->list);
                FREEMEMORY(space, f->gzindex[i]);
            }
            if (f->hasMates && f->mategzindex) {
                FREEMEMORY(space, f->mategzindex[i]->list);
                FREEMEMORY(space, f->mategzindex[i]);
            }
        }
        if (f->gzindex) {
            FREEMEMORY(space, f->gzindex);
        }

        if (f->hasMates && f->mategzindex) {
            FREEMEMORY(space, f->mategzindex);
        }
        f->gzip = 0;
    }

    return;
}

/*------------------------- bl_fastaDestructSequence -------------------------
 *
 * @brief descruct the fasta struct
 * @author Steve Hoffmann
 *
 */

void bl_fastxDestructSequence(void *space, fasta_t *f) {
    Uint i;

    for (i = 0; i < f->active_noofseqs; i++) {
        destructSequence(space, f->seqs[i]);
    }
    FREEMEMORY(space, f->seqs);

    if (bl_fastaHasMate(f)) {
        FREEMEMORY(space, f->matestart);
    }

    if (bl_fastaHasQuality(f)) {
        bl_fastaDestructQuality(space, f);
        FREEMEMORY(space, f->quals);
    }

    f->seqs = NULL;
    f->matestart = NULL;
    f->quals = NULL;
    f->active_noofseqs = 0;
    f->minlen = 0;
    f->maxlen = 0;
    f->chunkIsActive = 0;
}

/*----------------------------- bl_fastaDestruct -----------------------------
 *
 * @brief destruct fasta
 * @author Steve Hoffmann
 *
 */

void bl_fastaDestruct(void *space, fasta_t *f) {
    bl_fastxDestructSequence(space, f);
    bl_fastxDestructIndex(space, f);
    return;
}

/*------------------------------- bl_fastxDump -------------------------------
 *
 * @brief dump the fastx
 * @author Steve Hoffmann
 *
 */

void bl_fastxDump(void *space, fasta_t *fasta, char *desc, Uint desclen, char *sequence,
                  char *quality, Uint quallen, Uint n) {
    fprintf(stderr, "%s\n", desc);
    fprintf(stderr, "%s\n", sequence);
    fprintf(stderr, "+%s\n", &desc[1]);
    fprintf(stderr, "%s\n", quality);
}

/*------------------------------- bl_GFFwrite --------------------------------
 *
 * @brief write GFF to file filename
 * @author Steve Hoffmann
 *
 */

void bl_GFFwrite(char *filename, annotationtrack_t *set) {
    Uint i, j;
    FILE *fp;
    annotationitem_t *g;

    fp = fopen(filename, "w");
    if (fp == NULL) {
        fprintf(stderr, "couldn't open %s - exit forced", filename);
        exit(-1);
    }

    for (i = 0; i < set->noofitems; i++) {
        g = &set->items[i];
        fprintf(fp, "%s\t%s\t%s\t", g->chromname, g->source, g->name);
        // GFF file has a 1-offset as opposed to internal representation w/ 0-offset
        fprintf(fp, "%" PRIu64 "\t%" PRIu64 "\t%c\t", g->start + 1, g->end + 1, g->strand);
        fprintf(fp, "%d", g->frame);

        if (g->noofattributes) {
            fprintf(fp, "\t");
        }

        for (j = 0; j < g->noofattributes; j++) {
            fprintf(fp, "%s", g->attributes[j]);
            if (j < g->noofattributes - 1) {
                fprintf(fp, ";");
            }
        }

        fprintf(fp, "\n");
    }

    fclose(fp);
    return;
}

/*--------------------------- bl_blGFFAddAttribute ---------------------------
 *
 * @brief add an attribute to annotation item
 * @author Steve Hoffmann
 *
 */

void bl_GFFAddAttribute(void *space, annotationitem_t *item, char *attr, Uint len) {
    item->attributes = ALLOCMEMORY(space, item->attributes, char *, item->noofattributes + 1);

    item->attributelen = ALLOCMEMORY(space, item->attributelen, Uint, item->noofattributes + 1);

    item->attributes[item->noofattributes] = ALLOCMEMORY(space, NULL, char, len + 1);

    item->attributelen[item->noofattributes] = len;
    memmove(item->attributes[item->noofattributes], attr, len);

    item->attributes[item->noofattributes][item->attributelen[item->noofattributes]] = 0;

    item->noofattributes++;

    return;
}
/*-------------------------------- bl_annotationRead --------------------------------
 *
 * @brief read a gff or bed file
 * @author Steve Hoffmann
 *
 */

annotationtrack_t *bl_annotationRead(void *space, char *fn) {
    annotationtrack_t *annot = NULL;
    const char *suf = bl_fsuffix(fn);

    if (!strcmp(suf, ".bed")) {
        annot = bl_BEDread(NULL, fn);
    } else if (!strcmp(suf, ".gff") || !strcmp(suf, ".gff3")) {
        annot = bl_GFFread(NULL, fn);
    } else {
        NFO("please provide a bed or gff file with the appropriate extension.\n", NULL);
        exit(EXIT_FAILURE);
    }
    return annot;
}

/*-------------------------------- bl_GFFread --------------------------------
 *
 * @brief read a gff file
 * @author Steve Hoffmann
 *
 */

annotationtrack_t *bl_GFFread(void *space, char *filename) {
    stringset_t **set;
    annotationtrack_t *track;
    annotationitem_t *item;
    char *str, *pch;
    Uint linecount, i, j, p, noofstr, len, pchlen;

    track = ALLOCMEMORY(space, NULL, annotationtrack_t, 1);
    bl_annotationtrackInit(track);
    set = readcsv(space, filename, "\t", &linecount);
    track->items = ALLOCMEMORY(space, NULL, annotationitem_t, linecount);

    for (i = 0; i < linecount; i++) {
        noofstr = set[i]->noofstrings;

        /* fprintf(stderr, "input:'"); */
        /* for(uint64_t y=0; y < noofstr; y++){ */
        /*     fprintf(stderr, "%s\t", set[i]->strings[y].str); */
        /* } */
        /* fprintf(stderr, "'\n"); */

        if (noofstr) {
            str = set[i]->strings[0].str;
            len = strlen(str);

            // comment line
            if (strncmp(str, "#", 1) == 0) {
                destructStringset(space, set[i]);
                continue;
            }

            // track description
            if (len >= 5 && !strncmp(str, "track", 5)) {
                for (j = 1; j < noofstr; j++) {
                    str = set[i]->strings[j].str;
                    len = strlen(str);

                    if (len > 5 && !strncmp(str, "name=", 5)) {
                        track->tracknamelen = len - 5;
                        track->trackname = ALLOCMEMORY(space, NULL, char, len - 4);
                        memmove(track->trackname, &str[5], len - 5);
                        track->trackname[len - 5] = '\0';
                    }

                    if (len > 12 && !strncmp(str, "description=", 12)) {
                        track->descriptionlen = len - 12;
                        track->description = ALLOCMEMORY(space, NULL, char, len - 11);
                        memmove(track->description, &str[5], len - 12);
                        track->description[len - 12] = '\0';
                    }
                }
                destructStringset(space, set[i]);
                continue;
            }

            // real data
            if (noofstr >= 3) {
                item = &track->items[track->noofitems];
                bl_annotationitemInit(item, GFFITEM);

                for (j = 0; j < noofstr; j++) {
                    str = set[i]->strings[j].str;
                    len = strlen(str);

                    switch (j) {
                        case 0:
                            item->chromnamelen = len;
                            item->chromname = ALLOCMEMORY(space, NULL, char, len + 1);
                            memmove(item->chromname, str, len);
                            item->chromname[len] = '\0';
                            break;
                        case 1:
                            item->sourcelen = len;
                            item->source = ALLOCMEMORY(space, NULL, char, len + 1);
                            memmove(item->source, str, len);
                            item->source[len] = '\0';
                            break;
                        case 2:
                            item->namelen = len;
                            item->name = ALLOCMEMORY(space, NULL, char, len + 1);
                            memmove(item->name, str, len);
                            item->name[len] = '\0';
                            break;
                        case 3:
                            item->start = atoi(str);
                            if (!item->start && str[0] != '0') {
                                DBG("GFF '%s' %d:%d: atoi failed", filename, i, j);
                                exit(-1);
                            }
                            // convert GFF 1-offset to internal 0-offset
                            item->start -= 1;
                            break;
                        case 4:
                            item->end = atoi(str);
                            if (!item->end && str[0] != '0') {
                                DBG("GFF '%s' %d:%d: atoi failed", filename, i, j);
                                exit(-1);
                            }
                            // convert GFF 1-offset to internal 0-offset
                            item->end -= 1;
                            break;
                        case 5:
                            item->score = atof(str);
                            if (item->score == 0.0 && str[0] != '0' && str[0] != '.') {
                                DBG("GFF '%s' %d:%d: %f(%s) :atof failed", filename, i, j,
                                    item->score, str);
                                exit(-1);
                            }
                            break;
                        case 6:
                            if (str[0] != '-' && str[0] != '+' && str[0] != '.') {
                                DBG("GFF '%s' %d:%d: strand failed", filename, i, j);
                                exit(-1);
                            }
                            item->strand = str[0];
                            break;
                        case 7:
                            item->frame = atoi(str);
                            if ((!item->frame && str[0] != '.' && str[0] != '0') ||
                                item->frame > 2) {
                                DBG("GFF '%s' %d:%d: %s:atoi frame failed\n", filename, i, j, str);
                                exit(-1);
                            }
                            break;
                        case 8:
                            pch = strtok(str, ";");

                            while (pch != NULL) {
                                pchlen = strlen(pch);

                                for (p = 0; isspace((int)pch[p]) && p < pchlen; p++);
                                if (p < pchlen) {
                                    bl_GFFAddAttribute(space, item, &pch[p], strlen(&pch[p]));
                                }
                                pch = strtok(NULL, ";");
                            }

                            break;
                        default:
                            DBG("'%s' not in GFF format\n", filename);
                            exit(-1);
                            break;
                    }
                }

                track->noofitems++;
            }
        }
        destructStringset(space, set[i]);
    }

    qsort(track->items, track->noofitems, sizeof(annotationitem_t), bl_annotationitem_cmp);

    track->sorted = 1;

    /* bl_annotationtrackAssignTrackLevel(track); */

    FREEMEMORY(space, set);
    return track;
}

/*-------------------------------- bl_GTFread --------------------------------
 *
 * @brief read a GTF file an construct a list of gene models
 * @author Steve Hoffmann
 *
 */

geneset_t *bl_GTFread(char *filename) { return NULL; }

/*------------------------------- bl_initGene --------------------------------
 *
 * @brief init a gene_t record
 * @author Steve Hoffmann
 *
 */

gene_t *bl_initGene(char *identifier, char strand) {
    gene_t *mygene;

    mygene = ALLOCMEMORY(NULL, NULL, gene_t, 1);
    mygene->id = identifier;
    mygene->direction = strand;
    mygene->exons = NULL;
    mygene->noofexons = 0;

    return mygene;
}

/*-------------------------------- bl_addExon --------------------------------
 *
 * @brief add an exon to gene_t
 * @author Steve Hoffmann
 *
 */

void bl_addExon(gene_t *gene, Uint start, Uint end, char strand, char *refchr) {
    Uint k = gene->noofexons;
    gene->exons = ALLOCMEMORY(NULL, gene->exons, exon_t, k + 1);
    gene->exons[k].start = start;
    gene->exons[k].end = end;
    gene->exons[k].strand = strand;
    gene->exons[k].refchr = refchr;
    gene->exons[k].cds = NULL;
    gene->exons[k].noofcds = 0;
    gene->noofexons++;

    return;
}

void bl_dumpGene(gene_t *gene) {
    Uint i;

    for (i = 0; i < gene->noofexons; i++) {
        fprintf(stdout, "exon %d\t%d\t%d\n", i, gene->exons[i].start, gene->exons[i].end);
    }
}

void bl_wrapGene(gene_t *gene) {
    Uint i;
    for (i = 0; i < gene->noofexons; i++) {
        FREEMEMORY(NULL, gene->exons[i].cds);
    }
    FREEMEMORY(NULL, gene->exons);
}

char bl_isExon(gene_t *gene, Uint pos) {
    Uint i;

    for (i = 0; i < gene->noofexons; i++) {
        if (pos >= gene->exons[i].start && pos <= gene->exons[i].end) {
            return 1;
        }
    }

    return 0;
}

/*----------------------------- bl_getExonNumber -----------------------------
 *
 * @brief get the number of exons
 * @author Steve Hoffmann
 *
 */

Uint bl_getExonNumber(gene_t *model) { return model->noofexons; }

/*----------------------------- bl_getExonStart ------------------------------
 *
 * @brief get start of exon with id
 * @author Steve Hoffmann
 *
 */

Uint bl_getExonStart(gene_t *model, Uint i) { return model->exons[i].start; }

/*----------------------------- bl_getExonEnd -------------------------------
 *
 * @brief get start of exon with id
 * @author Steve Hoffmann
 *
 */

Uint bl_getExonEnd(gene_t *model, Uint i) { return model->exons[i].end; }
