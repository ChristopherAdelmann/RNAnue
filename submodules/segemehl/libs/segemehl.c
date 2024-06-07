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
 *  segemehl.c
 *  segemehl
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/10/2007 02:50:57 PM CEST
 *
 *  Revision of last commit:
 *  $Rev: 103 $
 *  $Author: steve $
 *  $Date: 2008-12-10 15:18:18 +0100 (Wed, 10 Dec 2008) $
 *
 *
 *  $Id: segemehl.c 103 2008-12-10 14:18:18Z steve $
 *  $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/src/segemehl.c $
 *
 */

#include "segemehl.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "basic-types.h"
#include "biofiles.h"
#include "citation.h"
#include "fileBins.h"
#include "filebuffer.h"
#include "fileio.h"
#include "info.h"
#include "iupac.h"
#include "manopt.h"
#include "manout.h"
#include "match.h"
#include "mathematics.h"
#include "memory.h"
#include "merge.h"
#include "multicharseq.h"
#include "samio.h"
#include "seqclip.h"
#include "stringutils.h"
#include "sufarray.h"
#include "version.h"
#include "vtprogressbar.h"

#ifdef HASHING
#include "hash.h"
#endif

pthread_mutex_t updatemtx;
pthread_mutex_t fastamtx;
unsigned char mute = 0;
char *ntcode;

/*------------------------------- se_progressBar--------------------------------
 *
 * @brief the pogress bar
 * @author Steve Hoffmann
 *
 */

void *se_progressBar(void *args) {
    checkthreadinfo_t *t;

    sleep(2);
    cursorVisible();
    t = (checkthreadinfo_t *)args;
    initProgressBarVT();

    while (pthread_mutex_trylock(&updatemtx) != 0) {
        progressBarVT("reads matched.", t->noofseqs, (*t->counter), 25);
    }

    cursorVisible();
    fprintf(stderr, "\n");
    return NULL;
}

/*-------------------------------- matchSlave --------------------------------
 *
 * @brief the slave for threaded matching
 * @author Steve Hoffmann
 *
 */

void *matchSlave(void *args) {
    segemehl_t *t;

    t = (segemehl_t *)args;
    match(t->space, t->sarray, t->reads, t);
    return NULL;
}

/*----------------------------- se_readQuery-----------------------------
 *
 * @brief read the queries and mates if available, perpare clipping
 * @author Steve Hoffmann
 *
 */
fasta_t *se_readQuery(segemehl_t *nfo, unsigned char checkid) {
    char ch, *adapter;

    if (nfo->queryfilename) {
        NFO("reading queries in '%s'.\n", nfo->queryfilename);
        /*
         * chunksize is triggered with last 1
         */
        nfo->reads = bl_fastxGetSet(NULL, &nfo->queryfilename, 1, 1, 0, 1, nfo->threadno);

        NFO("%d query sequences found.\n", nfo->reads->noofseqs);

        nfo->reads->checkid = checkid;
#ifdef HASHING
        if (!nfo->index) {
            MSG("Hashing without fasta index\n");
            bl_fastxGetTags(NULL, nfo->reads);
        } else {
            MSG("Hashing with fasta index\n");
            bl_fastxGetTags3(NULL, nfo->reads);
        }
        exit(-1);
#endif

        if (nfo->threadno > nfo->reads->noofseqs) {
            NFO("more threads than queries. Exit forced\n", NULL);
            exit(EXIT_FAILURE);
        }

        if (nfo->reads->noofseqs < 50 && nfo->autoclip) {
            NFO("A minimum of 50 queries is reccommended for autoclip.\n", nfo->reads->noofseqs);
            MSG("Do you want to proceed with autoclip? (y/n): ");
            while ((ch = getchar()) != 'n' && ch != 'y');
            if (ch == 'n') {
                MSG("Do you want to proceed without clipping? (y/n): ");
                while ((ch = getchar()) != 'n' && ch != 'y');
                if (ch == 'n')
                    exit(EXIT_SUCCESS);
                else
                    nfo->autoclip = 0;
            }
        }

        if (nfo->autoclip) {
            adapter = bl_seqclipFind3Prime(NULL, nfo->reads, 100000, 40, 10);
            NFO("found adapter sequence: '%s'\n", adapter);
            MSG("Do you want to clip? (y/n): ");
            while ((ch = getchar()) != 'n' && ch != 'y');
            if (ch == 'n') {
                MSG("Do you want to proceed without clipping? (y/n): ");
                while ((ch = getchar()) != 'n' && ch != 'y');
                if (ch == 'n') exit(EXIT_SUCCESS);
            } else {
                nfo->softclip3Prime = adapter;
            }
        }
    }

    if (nfo->softclip3Prime) {
        nfo->softclip3PrimeLen = strlen(nfo->softclip3Prime);
    }

    if (nfo->softclip5Prime) {
        nfo->softclip5PrimeLen = strlen(nfo->softclip5Prime);
    }

    if (nfo->queryfilename && nfo->matefilename) {
        NFO("reading mates in '%s'.\n", nfo->matefilename);

        nfo->reads =
            bl_fastxGetMateSet(NULL, nfo->reads, &nfo->matefilename, 1, 1, 0, 1, nfo->threadno);
        NFO("%d mate sequences found.\n", nfo->reads->noofseqs);
    }

    return nfo->reads;
}
/*-------------------------- se_prepareBisulfiteRun --------------------------
 *
 * @brief prepare the bisulfite run
 * @author Christian Otto
 *
 */

int se_prepareBisulfiteRun(segemehl_t *nfo, manopt_optionset *optset, Uint k, char *old, char *new,
                           unsigned char *ndx) {
    char oldch = ' ';
    char newch = ' ';
    unsigned char index = 0;
    Uint qfbaselen = 0, filebinbasenamelen;
    char *qfbase, *filebinbasename;

    /* initialize bisulfite matching run */
    if (k == 0) {
        initIUPAC(2, 1);

        if (manopt_isset(optset, 'i', NULL) ^ manopt_isset(optset, 'x', NULL)) {
            nfo->bisulfiterun = 1;
        } else {
            nfo->bisulfiterun = 2;
        }

        /* bisulfite binning in case of two matching runs */
        if (nfo->bisulfitemerging) {
            /* create domain for each matching run with bins for each thread */
            if (!nfo->filebinbasename) {
                qfbase = bl_basename(nfo->queryfilename);
                qfbaselen = strlen(qfbase);
                filebinbasename = ALLOCMEMORY(space, NULL, char, qfbaselen);
                memmove(filebinbasename, qfbase, bl_fileprefixlen(qfbase));
                filebinbasename[bl_fileprefixlen(qfbase)] = 0;
                nfo->filebinbasename = filebinbasename;
            }
            filebinbasenamelen = strlen(nfo->filebinbasename);

            NFO("creating bisulfite bins.\n", NULL);
            nfo->bins = se_createBisulfiteBins(NULL, 2, nfo->threadno, nfo->filebinbasename,
                                               filebinbasenamelen);

            if (nfo->bins == NULL) {
                NFO("Could not create bisulfite bins %s*! Exit forced.\n", nfo->filebinbasename);
                exit(-1);
            }
        }
    } else {
        if (manopt_isset(optset, 'i', NULL) ^ manopt_isset(optset, 'x', NULL) &&
            manopt_isset(optset, 'j', NULL) ^ manopt_isset(optset, 'y', NULL)) {
            nfo->bisulfiterun = 2;
            /* cleanup */
            destructMultiCharSeq(NULL, nfo->seq);
            bl_fastaDestruct(NULL, nfo->fasta);
            FREEMEMORY(NULL, nfo->fasta);
        } else {
            return -1;
        }
    }

    if (nfo->bisulfiterun == 1) {
        oldch = 'C';
        newch = 'T';
        index = manopt_isset(optset, 'x', NULL);
    } else if (nfo->bisulfiterun == 2) {
        // reset fastaMaster pointer
        nfo->nextfastaidx[0] = 0;
        oldch = 'G';
        newch = 'A';
        nfo->idxfilename = nfo->idx2filename;
        index = manopt_isset(optset, 'y', NULL);
    }
    /*
     * set conversion accordingly to run
     * in bisulfite and PARCLIP with 4SG
     * nfo->bisulfite = 1 in run 1
     * nfo->bisulfite = 2 in run 2
     */
    nfo->bisulfite = nfo->bisulfiterun;
    /*
     * adjustment of conversion in PAR-CLIP with 4SU:
     * nfo->bisulfite = 3 in run 1
     * nfo->bisulfite = 4 in run 2
     */
    if (nfo->bisulfiteprotocol == 3) {
        nfo->bisulfite = nfo->bisulfiterun + 2;
    }

    /*
     * set strand accordingly in bisulfite with
     * Lister et al.'s protocol and PARCLIP with 4SU
     * nfo->strand == 1 (plus) in run 1
     * nfo->strand == 2 (minus) in run 2
     */
    if (nfo->bisulfiteprotocol == 1 || nfo->bisulfiteprotocol == 3) {
        nfo->strand = nfo->bisulfiterun;
    }
    /*
     * adjustment to PAR-CLIP with 4SG:
     * nfo->strand == 2 (minus) in run 1
     * nfo->strand == 1 (plus) in run 2
     */
    if (nfo->bisulfiteprotocol == 4) {
        nfo->strand = 3 - nfo->bisulfiterun;
    }

    //      NFO("nfo->bisulfiteprotocol=%d\tnfo->bisulfite=%d\tnfo->strand=%d\tseedconv:%c->%c\n",
    //          nfo->bisulfiteprotocol, nfo->bisulfite, nfo->strand, oldch, newch);
    //      NFO("bisulfite/parclip mapping run %d\n", nfo->bisulfiterun, oldch, newch);
    //

    *old = oldch;
    *new = newch;
    *ndx = index;

    return 1;
}
/*-----------------------------------se_getESA----------------------------------
 *
 * @brief the main function
 * @author Steve Hoffmann
 *
 */

void se_getESA(unsigned char buildindex, segemehl_t *nfo) {
    uint64_t i;
    time_t startsuf, endsuf;
    double difsuf;
    unsigned int *suflinktable;

    if (buildindex) {
        time(&startsuf);
        nfo->sarray = constructSufArr(NULL, nfo->fasta->seqs, nfo->fasta->noofseqs, NULL, mute);

        for (i = 0; i < nfo->fasta->noofseqs; i++) {
            FREEMEMORY(NULL, nfo->fasta->seqs[i]->sequence);
            nfo->fasta->seqs[i]->sequence = NULL;
        }

        MSG("constructing lcp.\n");
        constructLcp(NULL, nfo->sarray);
        MSG("deleting inv_suftab\n");
        //     if(!nfo->check) [
        FREEMEMORY(NULL, nfo->sarray->inv_suftab);
        nfo->sarray->inv_suftab = NULL;

        MSG("constructing child tab.\n");
        constructchildtab(NULL, nfo->sarray);
        MSG("constructing suffix links.\n");
        MSG("constructing id.\n");
        computeId(NULL, nfo->sarray);
        MSG("constructing suflinks - bottom up.\n");
        suflinktable = getsufsucc(NULL, nfo->sarray);
        MSG("constructing suflinks - top down.\n");
        constructsuflinks(NULL, nfo->sarray, suflinktable);
        FREEMEMORY(NULL, suflinktable);
        time(&endsuf);
        difsuf = difftime(endsuf, startsuf);
        NFO("building  the suffix array has taken %f seconds.\n", difsuf);
        NFO("total length of suffix array was %u.\n", nfo->totallength);

        if (nfo->idxfilename) {
            NFO("writing suffix array '%s' to disk.\n", nfo->idxfilename);
            writeSuffixarray(nfo->sarray, nfo->idxfilename);
        }

    } else {
        time(&startsuf);
        NFO("reading suffix array '%s' from disk.\n", nfo->idxfilename);
        nfo->sarray =
            readSuffixarray(NULL, nfo->idxfilename, nfo->fasta->seqs, nfo->fasta->noofseqs, mute);

        for (i = 0; i < nfo->fasta->noofseqs; i++) {
            FREEMEMORY(space, nfo->fasta->seqs[i]->sequence);
            nfo->fasta->seqs[i]->sequence = NULL;
        }

        time(&endsuf);
        difsuf = difftime(endsuf, startsuf);
        NFO("reading the suffix array has taken %f seconds.\n", difsuf);
    }

    if (nfo->check) {
        NFO("checking suffixarray %s\n", nfo->idxfilename);
        for (i = 1; i < nfo->sarray->numofsuffixes - 1; i++) {
            if (!mute) {
                progressBarVT("suffixes checked.", nfo->sarray->numofsuffixes, i, 25);
            }
            if (strcmp(&nfo->sarray->seq->sequences[nfo->sarray->suftab[i - 1]],
                       &nfo->sarray->seq->sequences[nfo->sarray->suftab[i]]) > 0) {
                NFO("suffix array '%s' corrupted! Exit forced.\n", nfo->idxfilename);
                exit(-1);
            }
        }
        checksuflinks(nfo->sarray, 0, nfo->sarray->numofsuffixes - 1);
    }

    return;
}

/*----------------------------------- main -----------------------------------
 *
 * @brief the main function
 * @author Steve Hoffmann
 *
 */

int segemehl(int argc, char **argv) {
    segemehl_t info, *th_info;
    manopt_arg *unflagged;
    manopt_optionset optset;
    manopt_intconstraint threadconstraint;
    manopt_intconstraint accuracyconstraint;
    manopt_intconstraint jumpconstraint;
    manopt_intconstraint bisulfiteconstraint;
    unsigned char skipidcheck = 0;
    int *space = NULL;
    int i = 0, k, qfbaselen = 0;
    Uint filebinbasenamelen = 0, clipseqlen3 = 0, desclen;

    char oldch, newch, *qfbase, *splitfilebasename = NULL, *filebinbasename = NULL, *version,
                                *clipseq3 = NULL, *desc, *adapter = NULL;

    unsigned int counter = 0;
    unsigned char index = 0;

    bl_fileBinDomains_t *bins;

    double difmatch;
    time_t startmatch, endmatch;

    pthread_t *threads;
    pthread_t clockthread;
    checkthreadinfo_t ch_info;
    manopt_arg *dbfilenames;
    threadconstraint.max = 3000;
    threadconstraint.min = 1;
    accuracyconstraint.max = 100;
    accuracyconstraint.min = 0;
    jumpconstraint.max = INT_MAX;
    jumpconstraint.min = 0;
    bisulfiteconstraint.min = 1;
    bisulfiteconstraint.max = 6;

    se_setdefault(&info);

    bl_bsprintf(&info.cmdline, "%s", argv[0]);
    for (i = 1; i < argc; i++) {
        bl_bsprintf(&info.cmdline, " %s", argv[i]);
    }

#ifdef LOCAL
    MSG("Segemehl compiled with -DLOCAL option: resetting minscore, mincoverage and accuracy\n");
#endif

    version = getNiceGitVersion(VERSION, REVISION, TIME);
    manopt_initoptionset(&optset, argv[0], NULL, "Heuristic mapping of short sequences\n",
                         "SEGEMEHL is free software under GPL \n  2008 Bioinformatik Leipzig \n  "
                         "2018 Leibniz Institute on Aging (FLI) ",
                         version, "Please report bugs to steve@bioinf.uni-leipzig.de");

    manopt_blockseparator(&optset, "INPUT");
    manopt(&optset, LISTOPT, 1, 'd', "database",
           "list of path/filename(s) of fasta database sequence(s)", "<file> [<file>]", NULL, NULL);
    manopt(&optset, STRINGOPT, 0, 'q', "query", "path/filename of query sequences", "<file>", NULL,
           &info.queryfilename);
    manopt(&optset, STRINGOPT, 0, 'p', "mate", "path/filename of mate pair sequences", "<file>",
           NULL, &info.matefilename);
    manopt(&optset, REQSTRINGOPT, 0, 'i', "index", "path/filename of db index", "<file>", NULL,
           &info.idxfilename);
    manopt(&optset, REQSTRINGOPT, 0, 'j', "index2", "path/filename of second db index", "<file>",
           NULL, &info.idx2filename);
    manopt(&optset, REQSTRINGOPT, 0, 'x', "generate", "generate db index and store to disk",
           "<file>", NULL, &info.idxfilename);
    manopt(&optset, REQSTRINGOPT, 0, 'y', "generate2", "generate second db index and store to disk",
           "<file>", NULL, &info.idx2filename);
    /*
    manopt(&optset, REQUINTOPT, 0, 'F', "bisulfite",
      "bisulfite mapping with methylC-seq/Lister et al. (=1) or bs-seq/Cokus et al. protocol (=2),
    PAR-CLIP with 4SU (=3) or with 6SG (=4)", "<n>", &bisulfiteconstraint, &info.bisulfiteprotocol);
    */

    manopt(&optset, REQSTRINGOPT, 0, 'G', "readgroupfile", "filename to read @RG header", "<file>",
           NULL, &info.readgroupfile);
    manopt(&optset, REQSTRINGOPT, 0, 'g', "readgroupid", "read group id", "<string>", NULL,
           &info.readgroupid);

    manopt(&optset, REQUINTOPT, 0, 't', "threads", "start <n> threads", "<n>", &threadconstraint,
           &info.threadno);
    manopt_blockseparator(&optset, "OUTPUT");

    manopt(&optset, REQSTRINGOPT, 0, 'o', "outfile", "outputfile", "<string>", NULL, &info.outfn);
    manopt(&optset, FLAG, 0, 'b', "bamabafixoida", "generate a bam output (-o <filename> required)",
           NULL, NULL, &info.bam);
    manopt(&optset, REQSTRINGOPT, 0, 'u', "nomatchfilename", "filename for unmatched reads",
           "<file>", NULL, &info.nomatchname);
    manopt(&optset, FLAG, 0, 'e', "briefcigar", "brief cigar string (M vs X and =)", NULL, NULL,
           &info.briefcigar);
    manopt(&optset, FLAG, 0, 's', "progressbar", "show a progress bar", NULL, NULL, &info.mute);
    manopt(&optset, STRINGOPT, 0, 'B', "filebins",
           "file bins with basename <string> for easier data handling", "<string>", NULL,
           &info.filebinbasename);
    manopt(&optset, FLAG, 0, 'V', "MEOP",
           "output MEOP field for easier variance calling in SAM (XE:Z:)", NULL, NULL,
           &info.SAMmeop);
    manopt_blockseparator(&optset, "ALIGNMENT");

    manopt(&optset, REQUINTOPT, 0, 'F', "bisulfite",
           "bisulfite aln with methylC-seq/Lister et al. (=1) or bs-seq/Cokus et al. protocol (=2)",
           "<n>", &bisulfiteconstraint, &info.bisulfiteprotocol);

    manopt(&optset, STRINGOPT, 0, 'S', "splits", "detect split/spliced reads.", "[<basename>]",
           NULL, &info.splitfilebasename);

#ifdef LOCAL
    manopt(&optset, TRIPLEINTOPT, 0, 'k', "scores",
           "set scores for seeded alignments <match> <mismatch> <indel>", "<n><n><n>", NULL,
           &info.scores);
    manopt(&optset, PAIRINTOPT, 0, 'L', "local",
           "use local alignment; absolute min. score and length dependent min. score <absminscore> "
           "<relminscore> ",
           "<n><n>", NULL, &info.localparams);
#endif
    manopt(&optset, REQUINTOPT, 0, 'A', "accuracy",
           "min percentage of matches per read in semi-global alignment", "<n>",
           &accuracyconstraint, &info.accuracy);
    manopt(&optset, REQUINTOPT, 0, 'D', "differences",
           "search seeds initially with <n> differences", "<n>", NULL, &info.k_p);
    manopt(&optset, REQDBLOPT, 0, 'E', "evalue", "max evalue", "<double>", NULL, &info.maxevalue);
    manopt(&optset, REQUINTOPT, 0, 'H', "hitstrategy",
           "report only best scoring hits (=1) or all (=0)", NULL, NULL, &info.bestonly);

    manopt(&optset, REQUINTOPT, 0, 'm', "minsize", "minimum length of queries", "<n>", NULL,
           &info.minsize);

    manopt(&optset, REQUINTOPT, 0, 'Z', "minfraglen", "min length of a spliced fragment", "<n>",
           NULL, &info.minfragmentalignlen);
    manopt(&optset, REQUINTOPT, 0, 'W', "minsplicecover", "min coverage for spliced transcripts",
           "<n>", &accuracyconstraint, &info.minsplicedaligncover);
    manopt(&optset, REQUINTOPT, 0, 'U', "minfragscore", "min score of a spliced fragment", "<n>",
           NULL, &info.minfragmentalignscore);
    manopt(&optset, REQDBLOPT, 0, 'l', "splicescorescale",
           "report spliced alignment with score s only if <f>*s is larger than next best spliced "
           "alignment",
           "<f>", NULL, &info.chainscorescale);
    manopt(&optset, REQDBLOPT, 0, 'w', "maxsplitevalue", "max evalue for splits", "<double>", NULL,
           &info.maxsplitevalue);

    manopt_blockseparator(&optset, "SPECIAL");

    manopt(&optset, REQUINTOPT, 0, 'X', "dropoff", "dropoff parameter for extension", "<n>", NULL,
           &info.Xoff);
    manopt(&optset, REQUINTOPT, 0, 'J', "jump", "search seeds with jump size <n> (0=automatic)",
           "<n>", &jumpconstraint, &info.jump);
    manopt(&optset, FLAG, 0, 'O', "order",
           "sorts the output by chromsome and position (might take a while!)", "<n>", NULL,
           &info.order);
    manopt(&optset, REQINTOPT, 0, 'I', "maxpairinsertsize",
           "maximum size of the inserts (paired end) in case of multiple hits", "<n>", NULL,
           &info.maxpairinsertsize);
    manopt(&optset, REQUINTOPT, 0, 'M', "maxinterval",
           "maximum width of a suffix array interval, i.e. a query seed will be omitted if it "
           "matches more than <n> times",
           "<n>", NULL, &info.M);

    manopt(&optset, FLAG, 0, 'c', "checkidx", "check index", NULL, NULL, &info.check);
    manopt(&optset, REQUINTOPT, 0, 'n', "extensionpenalty",
           "penalty for a mismatch during extension", "<n>", NULL, &info.p_mis);
    manopt(&optset, REQUINTOPT, 0, 'r', "maxout",
           "maximum number of alignments that will be reported. If set to zero, all alignments "
           "will be reported",
           "<n>", NULL, &info.maxout);

    /* hidden options */
    manopt(&optset, FLAG, 0, 0, "skipidcheck",
           "do not check whether the fastq ids of mates / paired ends match. Instead, the first "
           "mate (-q) will be used for output only.",
           NULL, NULL, &skipidcheck);

    manopt(&optset, FLAG, 0, 0, "showalign", "show alignments", NULL, NULL, &info.align);
    manopt(&optset, FLAG, 0, 0, "nohead", "do not output header", NULL, NULL, &info.nohead);

    manopt(&optset, FLAG, 0, 'f', "fullname",
           "write full query name (no trunctation at whitespace)", NULL, NULL, &info.fullname);

    // manopt(&optset, FLAG, 0, 'z', "nosuflinks",
    //     "dont use suflinks (does not affect index construction, for short reads only, increases
    //     runtime!)", NULL, NULL, &info.nosuflinks);
    /*
    manopt(&optset, REQSTRINGOPT, 0, 'P', "prime5",
     "add 5' adapter", "<string>", NULL, &info.softclip5Prime);
    manopt(&optset, REQSTRINGOPT, 0, 'Q', "prime3",
     "add 3' adapter", "<string>", NULL, &info.softclip3Prime);
    manopt(&optset, REQINTOPT, 0, 'R', "clipacc",
     "clipping accuracy", "<n>", NULL, &info.clipacc);
    manopt(&optset, FLAG, 0, 'T', "polyA",
     "clip polyA tail", NULL, NULL, &info.polyA);
    manopt(&optset, FLAG, 0, 'Y', "autoclip",
     "autoclip unknown 3prime adapter", NULL, NULL, &info.autoclip);
    manopt(&optset, FLAG, 0, 'C', "hardclip",
      "enable hard clipping", NULL, NULL, &info.hardclip);

    */

    /*

    manopt(&optset, FLAG, 0, 'Z', "PAIR",
    "output pairstatus flag XA:Z: field in SAM", NULL, NULL, &info.SAMpairstat);

    */

    /*get unflagged options*/
    unflagged = manopt_getopts(&optset, argc, argv);

    /*
     *  option checks
     *
     */

    if (manopt_isset(&optset, 'O', "order") &&
        (!manopt_isset(&optset, 'o', "outfile") && !manopt_isset(&optset, 'B', "filebins"))) {
        manopt_help(&optset,
                    "option -O, --order may not be used when output is dumped to stdout!\n");
    }

    if (info.nomatchname != NULL) {
        info.nomatchdev = fopen(info.nomatchname, "w");
        if (info.nomatchdev == NULL) {
            manopt_help(&optset,
                        "could not open file for unmapped reads. Writing privileges set?\n");
        }
    }

    if (!manopt_isset(&optset, 'o', "outfile") && manopt_isset(&optset, 'b', "bamabafixoida")) {
        NFO("consider writing output to a file (-o <filename>) instead of stdout.\n", NULL);
    }

    if (manopt_isset(&optset, 'O', "order") && manopt_isset(&optset, 'b', "bamabafixoida")) {
        manopt_help(
            &optset,
            "sorting not supported in bam output mode. Please, use samtools sort manually.\n");
    }

    if (manopt_isset(&optset, 'B', "filebins") && manopt_isset(&optset, 'b', "bamabafixoida")) {
        manopt_help(&optset, "file bins not supported in bam output mode\n");
    }

    if (manopt_isset(&optset, 'U', "minfragscore")) {
        info.minsplicedalignscore = 2 * info.minfragmentalignscore;
    }

    if (manopt_isset(&optset, 'S', "splits") && manopt_isset(&optset, 'q', "query")) {
        info.split = 1;
    }

    if (!manopt_isset(&optset, 'F', "bisulfite")) {
        /*
         * option checks for normal mode
         *
         */
        index = manopt_isset(&optset, 'x', NULL);

        if (!(!manopt_isset(&optset, 'i', NULL) ^ !manopt_isset(&optset, 'x', NULL))) {
            manopt_help(&optset, "please give index filename using -i XOR -x option\n");
        } else if (unflagged->noofvalues > 1) {
            manopt_help(&optset, "unknown argument(s)\n");
        }
    } else {
        /*
         * option checks for bisulifte mode
         *
         */
        if (manopt_isset(&optset, 'i', NULL) && manopt_isset(&optset, 'x', NULL)) {
            manopt_help(&optset, "please give C/T index filename using -i XOR -x option\n");
        }
        if (manopt_isset(&optset, 'j', NULL) && manopt_isset(&optset, 'y', NULL)) {
            manopt_help(&optset, "please give G/A index filename using -j XOR -y option\n");
        }
        if (!manopt_isset(&optset, 'i', NULL) && !manopt_isset(&optset, 'x', NULL) &&
            !manopt_isset(&optset, 'j', NULL) && !manopt_isset(&optset, 'y', NULL)) {
            manopt_help(
                &optset,
                "please give C/T and/or G/A idx using (-i XOR -x) AND/OR (-j XOR -y) option\n");
        } else if (unflagged->noofvalues > 1) {
            manopt_help(&optset, "unknown argument(s)\n");
        }

        if (manopt_isset(&optset, 'q', "query")) {
            if ((manopt_isset(&optset, 'i', NULL) ^ manopt_isset(&optset, 'x', NULL)) &&
                (manopt_isset(&optset, 'j', NULL) ^ manopt_isset(&optset, 'y', NULL))) {
                info.bisulfitemerging = 1;
            }
            if (info.nomatchname != NULL) {
                NFO("Warning: file with unmapped reads may contain reads that are mapped \
            in one but not in both matching runs.\n",
                    NULL);
            }
            if (manopt_isset(&optset, 'S', "splits")) {
                manopt_help(&optset, "split alignments not yet supported with bisulfite mapping\n");
            }
        }
    }

#ifdef LOCAL

    NFO("local aln w/ sigma (%d,%d,%d), min abs sigma %d and len score mod %d\n", info.scores[0],
        info.scores[1], info.scores[2], info.localparams[0], info.localparams[1]);

#endif

    /* read queries */
    se_readQuery(&info, !skipidcheck);

    oldch = newch = ' ';
    for (k = 0; k < 2; k++) {
        /* reset counter variables */
        info.counter = 0;
        counter = 0;
        info.totallength = 0;

        /* normal matching run */
        if (!info.bisulfiteprotocol) {
            if (k == 1) {
                break;
            }
            initIUPAC(1, 1);
            info.bisulfite = 0;

            /* bisulfite matching run*/
        } else {
            if (se_prepareBisulfiteRun(&info, &optset, k, &oldch, &newch, &index) == -1) break;
        }

        MSG("reading database sequences.\n");
        dbfilenames = manopt_getarg(&optset, 'd', "database");
        info.fasta =
            bl_fastxGetSet(space, dbfilenames->values, dbfilenames->noofvalues, 1, 0, 0, 1);

        NFO("%d database sequences found.\n", info.fasta->noofseqs);
        for (i = 0; i < info.fasta->noofseqs; i++) {
            info.totallength += bl_fastaGetSequenceLength(info.fasta, i);
        }

        for (i = 0; i < info.fasta->noofseqs; i++) {
            desclen = bl_fastaGetDescriptionLength(info.fasta, i);
            desc = strclip(space, bl_fastaGetDescription(info.fasta, i), &desclen);
            FREEMEMORY(space, info.fasta->seqs[i]->description);
            info.fasta->seqs[i]->description = desc;
            info.fasta->seqs[i]->descrlen = desclen;
        }

        NFO("total length of db sequences: %u\n", info.totallength);

        if (k == 0) {
            /* compile header only once, necessary because bisulifte mode loops twice */
            info.head = se_header(&info, -1);
            NFO("compiled sam header.\n", NULL);
        }

        if (info.bisulfiteprotocol) {
            info.seq = concatCharSequences(space, info.fasta->seqs, info.fasta->noofseqs, (char)126,
                                           (char)127);

            /* character conversion */
            for (i = 0; i < info.fasta->noofseqs; i++) {
                strconvert(bl_fastaGetSequence(info.fasta, i),
                           bl_fastaGetSequenceLength(info.fasta, i), oldch, newch);
            }
        }

        if (!info.bisulfitemerging && !info.bins && manopt_isset(&optset, 'B', "filebins") &&
            manopt_isset(&optset, 'q', "query")) {
            if (!info.filebinbasename) {
                qfbase = bl_basename(info.queryfilename);
                qfbaselen = strlen(qfbase);
                filebinbasename = ALLOCMEMORY(space, NULL, char, qfbaselen);
                memmove(filebinbasename, qfbase, bl_fileprefixlen(qfbase));
                filebinbasename[bl_fileprefixlen(qfbase)] = 0;
                info.filebinbasename = filebinbasename;
            }

            filebinbasenamelen = strlen(info.filebinbasename);

            info.bins = se_createChromDomains(space, info.fasta, 500, 500, info.filebinbasename,
                                              filebinbasenamelen);

            if (info.bins == NULL) {
                NFO("Could not create bins %s*! Try w/o binning! Exit forced.\n",
                    info.filebinbasename);
                exit(-1);
            }
        }

        se_getESA(index, &info);

        if (info.queryfilename) {
            if (!info.bisulfiteprotocol) {
                info.seq = info.sarray->seq;
            }

            if (!info.bins && k == 0) {
                se_openOutputDevices(space, &info);
            }

            if (info.polyA) {
                info.polyAlen = MIN(80, info.reads->maxlen);
                clipseq3 = ALLOCMEMORY(space, NULL, char, info.polyAlen + 1);
                memset(&clipseq3[0], 'A', info.polyAlen);
                clipseqlen3 = info.polyAlen;
                clipseq3[info.polyAlen] = 0;
                info.minclipscr3 = 5;
            }

            if (info.softclip3Prime) {
                clipseqlen3 += info.softclip3PrimeLen;
                clipseq3 = ALLOCMEMORY(space, clipseq3, char, clipseqlen3 + 1);
                memmove(&clipseq3[info.polyAlen], info.softclip3Prime, info.softclip3PrimeLen);
                clipseq3[clipseqlen3] = 0;
                // info.minclipscr3 = floor((((float)info.softclip3PrimeLen)*info.clipacc)/100.);
                info.minclipscr3 = 5;
            }

            info.softclip3Prime = clipseq3;
            info.softclip3PrimeLen = clipseqlen3;

            if (info.softclip5Prime) {
                info.minclipscr5 = floor((((float)info.softclip5PrimeLen) * info.clipacc) / 100.);
            }

#ifdef SEEDHASH
            NFO("creating seed hash with size.\n", info.hashsize);
            info.hash = suffixArrayCreateHash(space, info.sarray, info.hashsize);
#endif

            if (info.threadno > 1) {
                info.counter = &counter;
                NFO("starting %d threads.\n", info.threadno);

                th_info = ALLOCMEMORY(space, NULL, segemehl_t, info.threadno);
                threads = ALLOCMEMORY(space, NULL, pthread_t, info.threadno);
                ch_info.noofseqs = info.reads->noofseqs;
                ch_info.counter = &counter;

                if (!mute && !info.mute) {
                    pthread_mutex_init(&updatemtx, NULL);
                    pthread_mutex_lock(&updatemtx);
                    pthread_create(&clockthread, NULL, se_progressBar, &ch_info);
                }

                time(&startmatch);

                for (i = 0; i < info.threadno; i++) {
                    memmove(&th_info[i], &info, sizeof(segemehl_t));
                    th_info[i].reads = info.reads;
                    th_info[i].threadid = i;
                    pthread_create(&threads[i], NULL, matchSlave, &th_info[i]);
                }

                for (i = 0; i < info.threadno; i++) {
                    pthread_join(threads[i], NULL);
                }

                if (!mute && !info.mute) {
                    /*notification via mutex - why use signals?*/
                    pthread_mutex_unlock(&updatemtx);
                    pthread_join(clockthread, NULL);
                }

                fflush(info.dev);
                time(&endmatch);
                difmatch = difftime(endmatch, startmatch);
                NFO("threaded matching w/ suffixarray has taken %f seconds.\n", difmatch);

                FREEMEMORY(space, th_info);
                FREEMEMORY(space, threads);

            } else {
                time(&startmatch);
                initProgressBarVT();
                match(info.space, info.sarray, info.reads, &info);  // match

                time(&endmatch);
                difmatch = difftime(endmatch, startmatch);
                NFO("matching w/ suffixarray has taken %f seconds.\n", difmatch);
            }
        }
        destructSufArr(space, info.sarray);
    } /* END OF for (k = 0; k < 2; k++) */

    /* merge thread-bins */
    if (info.bisulfitemerging) {
        bl_fileBinDomainsCloseAll(info.bins);
        bins = NULL;

        /* if no chromosomal binning --> write to output file */
        if (!manopt_isset(&optset, 'B', "filebins")) {
            se_openOutputDevices(space, &info);
        }
        /* otherwise initialize and write to chromosome bins */
        else {
            bins = se_createChromDomains(space, info.fasta, 500, 500, info.filebinbasename,
                                         filebinbasenamelen);
            if (bins == NULL) {
                NFO("Could not create bins %s*! Try w/o binning! Exit forced.\n",
                    info.filebinbasename);
                exit(-1);
            }
        }

        /* reset mapping stat */
        memset(info.stats, 0, sizeof(mappingstats_t));
        /* do bisulfite merging and cleanup */
        se_mergeBisulfiteBinsNew(info.bins, info.reads, info.head, info.dev, bins, 1, info.bestonly,
                                 info.stats, &info);

        /* destruct bins */
        bl_fileBinDomainsDestruct(space, info.bins);
        FREEMEMORY(space, info.bins);

        info.bins = bins;
    }

    se_closeOutputDevices(NULL, &info);

    MSG("Mapping stats:\n");
    se_printStats(stderr, &info);

    if (adapter) {
        FREEMEMORY(space, adapter);
    }

    if (info.filebinbasename && !manopt_isset(&optset, 'B', "filebins")) {
        FREEMEMORY(space, info.filebinbasename);
    }

    if (splitfilebasename) {
        FREEMEMORY(space, splitfilebasename);
    }

    se_destructInfo(space, &info);
    manopt_destructoptionset(&optset);
    // manopt_destructarg(unflagged);
    FREEMEMORY(NULL, unflagged);
    bl_samdestructHeader(info.head);
    FREEMEMORY(space, info.head);
    FREEMEMORY(space, version);

    NFO("\nGoodbye.\n %s\n", citerand());
    return EXIT_SUCCESS;
}
