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
 *  haarz.c
 *  
 *
 *  @author 
 *  @email 
 *  @date 04/01/2018 20:04:30 CEST
 *  
 */



#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <pthread.h>

#include "info.h"
#include "basic-types.h"
#include "mathematics.h"
#include "manout.h"
#include "bamio.h"
#include "manopt.h"
#include "annotation.h"
#include "bedfiles.h"
#include "intervaltree.h"
#include "version.h"

unsigned char mute=0;
char *ntcode;

/*-------------------------- bl_splitReadOffsetString ----------------------
 *    
 * @brief read the offset string
 * @author Steve Hoffmann 
 *   
 */

int
bl_splitReadOffsetString(char *str, annotationoffs_t *offs) {
  /*
   * the format is [>]\[d,d\]
   *
   */

  Uint n,i;
  int64_t a=0,b=0;
  int begptr=0;
  char ch, start=1, directionality=0, open=0, sep=0, closed=0, select=0;
  n = strlen(str);

  for(i=0; i < n; i++) {
    ch = str[i];
    switch(ch) {
      case 's':
      case 'e':
      case '>':
        if(!start) {
          NFO("1: malformed expression '%s'. Exit forced.\n", str);
          return -1;
        }
        directionality = 1;
        start=0;
        if(ch == 's' && !select) select = 1;
        if(ch == 'e' && !select) select = 2;
        break;
      case '[':
        if((!directionality && !start) || open || sep) {
          NFO("2: malformed expression '%s'. Exit forced.\n", str);
          return -1;
        }
        begptr = i+1;
        start=0;
        closed=0;
        open=1;
        break;
      case ',':
        if(!open || start || sep || closed) {
          NFO("malformed expression '%s'. Exit forced.\n", str);
          return -1;
        }
        str[i]=0;

        if(i == begptr) {
          if(directionality) {
            NFO("blank offset instead of 0: selected 3prime end only\n", NULL);
          } else {
            NFO("blank offset instead of 0: selected right end only\n", NULL);
          }
          a = LONG_MIN;
        } else {
          a= strtol(&str[begptr], NULL, 10);
          if (errno == ERANGE){
            NFO("range error for number '%s' Exit forced.\n", &str[begptr]);
            return -1;
          }
        }
        str[i] = ',';
        begptr = i+1;
        sep=1;
        break;
      case ']':
        if(!open || start || !sep || closed) {
          NFO("3: malformed expression '%s'. Exit forced.\n", str);
          return -1;
        }

        str[i]=0;
        if(i == begptr) {
          if(directionality) {
            NFO("blank offset instead of 0: selected 5prime end only\n", NULL);
          } else {
            NFO("blank offset instead of 0: selected left end only\n", NULL);
          }
          b = LONG_MIN;
        } else {
          b= strtol(&str[begptr], NULL, 10);
          if (errno == ERANGE){
            NFO("range error for number '%s' Exit forced.\n", &str[begptr]);
            return -1;
          } 
        }
        str[i] = ']';

        if(directionality) {
          offs->dir5prime = a;
          offs->dir3prime = b;
          offs->select = select;
        } else {
          offs->left=a;
          offs->right=b;
        } 

        closed = 1;
        sep=0;
        start=1;
        open=0;
        directionality = 0;

        break;
      default:
        break;
    }
  }


  return 1;
}


/*-------------------------- bl_splitWriteTableHeader ------------------------
 *    
 * @brief write the header of the split table
 * @author Steve Hoffmann 
 *   
 */

void bl_splitWriteTableHeader(FILE *dev, annotationmultitrack_t *bed) {
  Uint i;

  fprintf(dev, "chr\tleft\tright\tn\tmedian_qual\t");
  for(i=0; i < bed->nooftracks-1; i++) {
    fprintf(dev, "%s\t", bed->filename[i]);
  }

  fprintf(dev, "%s\n", bed->filename[i]);
}

/*--------------------------- bl_splitWriteTableRow --------------------------
 *    
 * @brief write the row of the split table
 * @author Steve Hoffmann 
 *   
 */

void bl_splitWriteTableRow(FILE *dev, annotationitem_t *item, uint64_t *C, 
    Uint n, char *Q, Uint m, double minmedianqual, char nostrand) {  

  Uint i,j; 

  double m_hat = median_char(Q, m);
  if(m_hat >= minmedianqual) {

    fprintf(dev, "%s\t%"PRIu64"\t%"PRIu64"\t%d\t%f\t", item->chromname,
        item->start, item->end, m, m_hat); 

    if(nostrand) {
      fprintf(dev, ".");
    } else {
      fprintf(dev, "%c", item->strand);
    }

    for(i=0; i < n; i++) {
      fprintf(dev, "\t%"PRIu64"", C[i]); 
    }
    
    if(item->noofattributes) {
      for(j =0 ; j < item->noofattributes; j++) {
        fprintf(dev, "\t%s", item->attributes[j]);
      }
    } 

    fprintf(dev, "\n");
  }

  return;
}

/*----------------------------- bl_splitSummarize ----------------------------
 *    
 * @brief summarize the split table 
 * @author Steve Hoffmann 
 *   
 */

void bl_splitSummarize(FILE *fp, annotationmultitrack_t *trk, Uint mintotalsplit, 
    double minmedianqual, char nostrand) {

  uint64_t i=0,j,n,total=0; 
  char *Q;
  uint64_t *C;
  Uint nooftracks=0;

  //in practice Q is much smaller -> change to a more
  //flexible size allocation
  n= trk->noofitems;
  Q=ALLOCMEMORY(NULL, NULL, char, n);
  memset(Q, 0, sizeof(char)*n);

  //count the trackids
  nooftracks = trk->nooftracks;
  C=ALLOCMEMORY(NULL, NULL, uint64_t, trk->nooftracks);
  memset(C, 0, sizeof(uint64_t)*trk->nooftracks);

  bl_splitWriteTableHeader(fp, trk);

  for(j=1; j < n; j++) {

    Q[j-i-1] = (int) trk->items[j-1].score;
    C[trk->items[j-1].trackid] += 1;

    if((nostrand && bl_annotationitem_nostrand_cmp(&trk->items[i], &trk->items[j])) || 
        (!nostrand && bl_annotationitem_cmp(&trk->items[i], &trk->items[j]))) {

      if(j-i >= mintotalsplit) {
        total += j-i;
        bl_splitWriteTableRow(fp, &trk->items[i], C, nooftracks, 
            Q, j-i, minmedianqual, nostrand); 
      }
      memset(C, 0, sizeof(uint64_t)*trk->nooftracks);
      i=j;
    } 
  }
  FREEMEMORY(NULL, Q);
  FREEMEMORY(NULL, C);
}

/*----------------------------- bl_splitOverlap ----------------------------
 *    
 * @brief summarize the split table 
 * @author Steve Hoffmann 
 *   
 */

void bl_splitOverlap(annotationmultitrack_t *mspl, intervalforest_t *forest, 
    char mode,  manopt_arg *attributes){

  uint64_t i, u, k, q, nelems=0;
  void **results=NULL;

  for(i=0; i < mspl->noofitems; i++) {

    bl_intervalforestSearch(forest, mspl->items[i].chromname, &mspl->items[i], 
        getlow_annotitem, gethigh_annotitem, &results, &nelems, NULL);

    Uint minmax=0;

    if(nelems > 0) {
      /*
       * output all overlapping elems
       *
       */
      if(mode == 0) {

        for(u=0; u < nelems; u++) {
          annotationitem_t *item = (annotationitem_t*) results[u];  

          for(k=0; k < item->noofattributes; k++) {
            if(attributes){
              for(q=0; q < attributes->noofvalues; q++) {
                char * attribute = attributes->values[q];
                if (strstr(item->attributes[k], attribute)) {
                  bl_GFFAddAttribute(NULL, &mspl->items[i], item->attributes[k], 
                      item->attributelen[k]);
                }
              }
            } else {

              bl_GFFAddAttribute(NULL, &mspl->items[i], item->attributes[k], 
                  item->attributelen[k]);
            }
          }
        }
      }
      /*
       * output largest or smallest intervals
       *
       */
      if(mode) { 
        for(u=0; u < nelems; u++) {
          annotationitem_t *item = (annotationitem_t*) results[u];
          Uint ell = item->end - item->start;
          annotationitem_t *pre = (annotationitem_t*) results[minmax];
          if ((mode == 1 && ell > pre->end - pre->start) ||
              (mode == 2 && ell < pre->end - pre->start)) {

            if(attributes){
              for(k=0; k < item->noofattributes; k++) {
                for(q=0; q < attributes->noofvalues; q++) {
                  char * attribute = attributes->values[q];
                  if (strstr(item->attributes[k], attribute)) {
                    minmax=u;
                  }
                }
              }
            } else { 
              minmax = u;
            }
          }
        }
        annotationitem_t *item = (annotationitem_t*) results[minmax];

        for(k=0; k < item->noofattributes; k++) {
          if(attributes){
            for(q=0; q < attributes->noofvalues; q++) {
              char * attribute = attributes->values[q];
              if (strstr(item->attributes[k], attribute)) {
                bl_GFFAddAttribute(NULL, &mspl->items[i], item->attributes[k], 
                    item->attributelen[k]);
              }
            }
          } else {
            bl_GFFAddAttribute(NULL, &mspl->items[i], item->attributes[k], 
                item->attributelen[k]);
          }

        }
      }
    }

    nelems = 0;
    FREEMEMORY(NULL, results); 
  }

  FREEMEMORY(NULL, results);
}

int main(int argc,char** argv) {

  manopt_optionset prgset; 
  manopt_arg *unflagged;
  manopt_arg *prg;
  manopt_optionset optset;

  uint32_t i=0;
  uint8_t uniqueonly=0;
  int ret;
  bam_info_t *bam;
  pthread_t *threads;
  pthread_mutex_t devmtx, msmtx,bammtx;
  uint32_t nthreads=10;
  uint8_t isthreaded = 0;
  bl_bam_methyl_worker_t* wk;
  bl_bam_methyl_master_t ms;
  circbuffer_t *devbuffer;
  FILE *outfp = stdout;
  char *bamfn = NULL;
  char *fafn = NULL;
  char *outfn = NULL;
  //  char *gfffn = NULL;
  char *version;

  version = getNiceGitVersion(VERSION, REVISION, TIME);

  manopt_initoptionset(&prgset, argv[0], NULL, 
      "Heuristic mapping of short sequences\n\n  available programs:\n\n  callmethyl \t generate methylation vcf from bam\n  methylstring \t get SAM file with methylation string annotation\n  split summarize and annotate segemehl split info\n",
      "SEGEMEHL is free software under GPL \n  2008 Bioinformatik Leipzig \n  2018 Leibniz Institute on Aging (FLI) ",
      version,
      "Please report bugs to steve@bioinf.uni-leipzig.de");


  prg = manopt_getopts(&prgset, MIN(argc,2), argv);

  if(prg->noofvalues == 1) { 
    manopt_help(&prgset, "program needs to be selected\n");
  }

  manopt_initoptionset(&optset, argv[0], NULL, 
      "Heuristic mapping of short sequences\n",
      "SEGEMEHL is free software under GPL \n  2008 Bioinformatik Leipzig \n  2018 Computational Biology, Leibniz Institute on Aging (FLI) ",
      version,
      "Please report bugs to steve@bioinf.uni-leipzig.de");
  

  /*
   *
   *   call methylation string   
   *
   *
   */

  if(!strcmp(prg->values[1], "methylstring")) {
    int32_t i=0; 
    annotationoffs_t offs;

    intervalforest_t *forest=NULL;

    manopt_listconstraint ovlconstraint;
    manopt_arg *annotfns;  
    annotationtrack_t *annot = NULL; 
    annotationmultitrack_t mannot; 


    ovlconstraint.items = ALLOCMEMORY(NULL, NULL, char*, 3);
    ovlconstraint.items[0] = "ALL";
    ovlconstraint.items[1] = "LARGEST";
    ovlconstraint.items[2] = "SMALLEST";
    ovlconstraint.noofitems = 3;
    ovlconstraint.maxlength = 1;
    ovlconstraint.minlength = 1;

    manopt_blockseparator(&optset, "methylation string");
    manopt_blockseparator(&optset, "INPUT");

    manopt(&optset, STRINGOPT, 1, 'd', "database", 
        "list of path/filename(s) of fasta database sequence(s)", "<file> [<file>]", 
        NULL, &fafn);    
    manopt(&optset, STRINGOPT, 1, 'b', "bam", 
        "path/filename of sorted and indexed (!) bamfile", "<file>", NULL, &bamfn);
    manopt(&optset, LISTOPT, 0, 'a', "annotationfiles", 
        "list of path/filename(s) of BED or GFF file (s))", "<file> [<file>]", 
        NULL, NULL);
     manopt(&optset, LISTOPT, 0, 'A', "attributes", 
        "attributes that shall be selected for overlap annotation","<string>", NULL, NULL);
    manopt(&optset, LISTOPT, 0, 'O', "offsets", 
        "offsets for the annotation","<string>", NULL, NULL);
    manopt(&optset, STRINGOPT, 0, 'o', "output", 
        "path/filename of output file (will be sorted)", "<file>", NULL, &outfn);
    manopt(&optset, FLAG, 0, 'u', "uniqueonly", 
      "generate a bam output (-o <filename> required)", NULL, NULL, &uniqueonly); 

    manopt(&optset, REQUINTOPT, 0, 't', "threads", 
        "start <n> threads", "<n>", NULL, &nthreads);



    unflagged = manopt_getopts(&optset, argc, argv);
    pthread_mutex_init(&devmtx, NULL);

    if(nthreads>1) { 
      pthread_mutex_init(&msmtx, NULL);
      pthread_mutex_init(&bammtx, NULL);
      isthreaded = 1;
    }

    if(outfn) {
      outfp = fopen(outfn, "w");
      if(!outfp) {
        NFO("Couldnt open file '%s'. Exit forced!\n", outfn);
        exit(-1);
      }
    }   
    
    //OFFSETS    
    if(manopt_isset(&optset, 'O', NULL)) {
      bl_annotationitemInitOffset(&offs);
      manopt_arg *offsets = manopt_getarg(&optset, 'O', "offsets");
      bl_splitReadOffsetString(offsets->values[0], &offs); 
    }


    //READ ANNOTATION 
    
    if(manopt_isset(&optset, 'a', NULL)) {

      annotfns = manopt_getarg(&optset, 'a', "annotationfiles");
      bl_annotationmultitrackInit(&mannot);

      for(i=0; i < annotfns->noofvalues; i++) {
        NFO("reading annotation '%s'.\n", annotfns->values[i]);
        annot = bl_annotationRead(NULL, annotfns->values[i]); 
        bl_annotationtrackJoin(NULL, &mannot, annot);
      }
  
      if(manopt_isset(&optset, 'O', NULL)) {
        NFO("applying offsets: %lld, %lld, %lld, %lld\n", offs.dir5prime, 
          offs.dir3prime, offs.left, offs.right); 
        bl_annotationmultitrackApplyOffset(&mannot, &offs); 
      }

      forest = bl_intervalforestBuildAnnot(&mannot);
    
      NFO("created intervaltrees for %d chromosomes.\n", forest->n);
    }   
  
    devbuffer = bl_circBufferInitArray(nthreads, BAM_CS_DEFAULT_CIRCBUF_SZ, 
        outfp, &devmtx);

    bam = ALLOCMEMORY(space, NULL, bam_info_t, nthreads);
    wk = ALLOCMEMORY(space, NULL, bl_bam_methyl_worker_t, nthreads);
    threads = ALLOCMEMORY(space, NULL, pthread_t, nthreads);

    for(i=0; i < nthreads; i++) {
      bl_bamLoadInfo(&bam[i], bamfn, fafn); 
      if(i==0) {
        bl_bamInitMaster(&ms, bam[i].hdr, nthreads>1, &msmtx); 
      }
      bl_bamInitWorker(&wk[i], &ms, stdout, &bam[i], uniqueonly, isthreaded, &devmtx, 
          &bammtx, &devbuffer[i], forest);
    }

    if(nthreads > 1) {
      //create multiple threads
      for(i=0; i < nthreads; i++){
        pthread_create(&threads[i], NULL, bl_bamMethylStringWorker, &wk[i]);
      }
      //spawn!
      for(i=0; i < nthreads; i++) {
        pthread_join(threads[i], NULL); 
      }
    } else {
      //single thread
       bl_bamMethylStringWorker(&wk[0]);
    }

    bl_circBufferEmptyArray(devbuffer, nthreads);  
    bl_circBufferDestructArray(devbuffer, nthreads);

    for(i=0; i < nthreads; i++){
      bl_bamDestructInfo(&bam[i]);
    }

    FREEMEMORY(NULL, bam);
    FREEMEMORY(NULL, wk);
    FREEMEMORY(NULL, threads); 

    pthread_mutex_destroy(&devmtx);

    if(nthreads>1) {  
      pthread_mutex_destroy(&msmtx);
      pthread_mutex_destroy(&bammtx);
    }

    if(outfn) {
      fclose(outfp);
    }

    if(forest) {
      bl_annotationmultitrackDestruct (NULL, &mannot);
      bl_annotationtrackDestruct(NULL, annot);
      FREEMEMORY(NULL, annot);
     
      bl_intervalforestDestruct(forest);  
      FREEMEMORY(NULL, forest);
    }

    manopt_destructarg(unflagged);
    FREEMEMORY(NULL, unflagged);
    FREEMEMORY(NULL, devbuffer);
    FREEMEMORY(NULL, ovlconstraint.items);

    /*
     *
     *   call methyl vcf     
     *
     *
     */

  } else if(!strcmp(prg->values[1], "callmethyl")) {


    manopt_blockseparator(&optset, "methylation caller");
    manopt_blockseparator(&optset, "INPUT");

    manopt(&optset, STRINGOPT, 1, 'd', "database", 
        "list of path/filename(s) of fasta database sequence(s)", "<file> [<file>]", 
        NULL, &fafn);    
    manopt(&optset, STRINGOPT, 1, 'b', "bam", 
        "path/filename of sorted and indexed (!) bamfile", "<file>", NULL, &bamfn);
    manopt(&optset, REQUINTOPT, 0, 't', "threads", 
        "start <n> threads", "<n>", NULL, &nthreads);
    manopt(&optset, STRINGOPT, 0, 'o', "output", 
        "path/filename of output file (will be sorted)", "<file>", NULL, &outfn);
    manopt(&optset, FLAG, 0, 'u', "uniqueonly", 
      "only use uniquely mapped reads", NULL, NULL, &uniqueonly); 
    

    unflagged = manopt_getopts(&optset, argc, argv);
    pthread_mutex_init(&devmtx, NULL);

    if(nthreads>1) { 
      pthread_mutex_init(&msmtx, NULL);
      pthread_mutex_init(&bammtx, NULL);
      isthreaded = 1;
    }

    if(outfn) {
      outfp = fopen(outfn, "w");
      if(!outfp) {
        NFO("Couldnt open file '%s'. Exit forced!\n", outfn);
        exit(-1);
      }
    }   

    devbuffer = bl_circBufferInitArray(nthreads, BAM_CS_DEFAULT_CIRCBUF_SZ, 
        outfp, &devmtx);

    bam = ALLOCMEMORY(space, NULL, bam_info_t, nthreads);
    wk = ALLOCMEMORY(space, NULL, bl_bam_methyl_worker_t, nthreads);
    threads = ALLOCMEMORY(space, NULL, pthread_t, nthreads);

    for(i=0; i < nthreads; i++) {
      bl_bamLoadInfo(&bam[i], bamfn, fafn); 
      if(i==0) {
        bl_bamInitMaster(&ms, bam[i].hdr, nthreads>1, &msmtx); 
      }
      bl_bamInitWorker(&wk[i], &ms, stdout, &bam[i], uniqueonly, isthreaded, &devmtx, 
          &bammtx, &devbuffer[i], NULL);
    }

    if(nthreads > 1) {
      //create multiple threads
      for(i=0; i < nthreads; i++){
        pthread_create(&threads[i], NULL, bl_bamCrossSectioMethylWorker, &wk[i]);
      }
      //spawn!
      for(i=0; i < nthreads; i++) {
        pthread_join(threads[i], NULL); 
      }
    } else {
      //single thread
      bl_bamCrossSectioMethylWorker(&wk[0]);
    }

    bl_circBufferEmptyArray(devbuffer, nthreads);  
    bl_circBufferDestructArray(devbuffer, nthreads);

    for(i=0; i < nthreads; i++){
      bl_bamDestructInfo(&bam[i]);
    }

    FREEMEMORY(NULL, bam);
    FREEMEMORY(NULL, wk);
    FREEMEMORY(NULL, threads); 

    pthread_mutex_destroy(&devmtx);

    if(nthreads>1) {  
      pthread_mutex_destroy(&msmtx);
      pthread_mutex_destroy(&bammtx);
    }

    manopt_destructarg(unflagged);
    FREEMEMORY(NULL, unflagged);

    if(outfn) {
      fclose(outfp);
      ret = bl_UnixSort(NULL, outfn, "-k1,1V -k2,2n --parallel=10", '\t');
      if(ret==-1) {
        MSG("sort failed. Try to sort vcf w/ 'sort -k1,1V -k2,2n'.\n");
      }
    }   

    FREEMEMORY(NULL, devbuffer);
  /*
   *
   *  summarize split program
   *
   *
   */
  
  } else if(!strcmp(prg->values[1], "split")) {
 
    intervalforest_t *forest;
    annotationoffs_t offs;
   
    char *ovldefault="ALL";
    uint32_t i;

    Uint mintotalsplit = 5;
    double minmedianqual = 25.0;
    char nostrand = 1;

    annotationtrack_t *annot = NULL;
    annotationtrack_t **spl = NULL;
    annotationmultitrack_t mspl;
    annotationmultitrack_t mannot; 
    manopt_listconstraint ovlconstraint;
    manopt_arg *annotfns;
    manopt_arg *splfns;
    manopt_arg *ovlmode;
    manopt_arg *atribte = NULL;

    ovlconstraint.items = ALLOCMEMORY(NULL, NULL, char*, 3);
    ovlconstraint.items[0] = "ALL";
    ovlconstraint.items[1] = "LARGEST";
    ovlconstraint.items[2] = "SMALLEST";
    ovlconstraint.noofitems = 3;
    ovlconstraint.maxlength = 1;
    ovlconstraint.minlength = 1;


    manopt_blockseparator(&optset, "split");
    manopt_blockseparator(&optset, "INPUT");
    manopt(&optset, LISTOPT, 1, 'f', "files", 
        "list of path/filename(s) of bed files with split info (s)", "<file> [<file>]", 
        NULL, NULL);
    manopt(&optset, REQUINTOPT, 0, 'm', "minsplit", 
        "minimum total split number (all samples) of junction", "<n>", NULL, &mintotalsplit);
    manopt(&optset, REQDBLOPT, 0, 'q', "minqual", 
        "minimum median quality of junction", "<f>", NULL, &minmedianqual);
    manopt(&optset, LISTOPT, 0, 'a', "annotationfiles", 
        "list of path/filename(s) of GFF file (s))", "<file> [<file>]", 
        NULL, NULL);
    manopt(&optset, SELECTOPT, 0, 'M', "ovlmode", 
        "annotation mode LARGEST, SMALLEST or ALL","<string>", &ovlconstraint, &ovldefault);
    manopt(&optset, LISTOPT, 0, 'A', "attributes", 
        "attributes that shall be selected for overlap annotation","<string>", NULL, NULL);
    manopt(&optset, LISTOPT, 0, 'O', "offsets", 
        "offsets for the annotation","<string>", NULL, NULL);
    
    unflagged = manopt_getopts(&optset, argc, argv);

    if(!manopt_isset(&optset, 'f', NULL)) {
      manopt_help(&optset, "input bedfiles are missing\n");
      //the selected program is an unflagged value: fix it
    } else if(unflagged->noofvalues > 2) { 
      manopt_help(&optset, "unknown argument(s)\n");
    }

    //READ SPLITS
    splfns = manopt_getarg(&optset, 'f', "files");
    ovlmode = manopt_getarg(&optset, 'M', "ovlmode");
    atribte = manopt_getarg(&optset, 'A', "attributes");
    
    bl_annotationmultitrackInit(&mspl);

    NFO("reading %d files.\n", splfns->noofvalues);
    spl = ALLOCMEMORY(NULL, NULL, annotationtrack_t*, splfns->noofvalues);
    for(i=0; i < splfns->noofvalues; i++) {
      spl[i] = bl_BEDread(NULL, splfns->values[i]); 
      bl_annotationtrackJoin(NULL, &mspl, spl[i]);  
    }

    NFO("sorting %d items.\n", mspl.noofitems);
    annotationitem_t *ptr = &mspl.items[0];
    qsort(ptr, mspl.noofitems, sizeof(annotationitem_t), bl_annotationitem_cmp);

    //OFFSETS    
    if(manopt_isset(&optset, 'O', NULL)) {
      bl_annotationitemInitOffset(&offs);
      manopt_arg *offsets = manopt_getarg(&optset, 'O', "offsets");
      bl_splitReadOffsetString(offsets->values[0], &offs); 
    }

    //READ ANNOTATION
    if(manopt_isset(&optset, 'a', NULL)) {

      annotfns = manopt_getarg(&optset, 'a', "annotationfiles");
      bl_annotationmultitrackInit(&mannot);

      for(i=0; i < annotfns->noofvalues; i++) {
        NFO("reading annotation in GFF'%s'.\n", annotfns->values[i]);
        annot = bl_GFFread(NULL, annotfns->values[i]);  
        bl_annotationtrackJoin(NULL, &mannot, annot); 
      }

      if(manopt_isset(&optset, 'O', NULL)) {
        NFO("applying offsets: %lld, %lld, %lld, %lld\n", offs.dir5prime, 
          offs.dir3prime, offs.left, offs.right); 
        bl_annotationmultitrackApplyOffset(&mannot, &offs); 
      }
       
      forest = bl_intervalforestBuildAnnot(&mannot);
    
      NFO("created intervaltrees for %d chromosomes.\n", forest->n);
      NFO("searching %d items.\n", mspl.noofitems);
      
      Uint mode = 0;
      if(manopt_isset(&optset, 'M', NULL)) {
        ovldefault = ovlmode->values[0];
        if(!strcmp(ovldefault, "LARGEST")) {
          mode = 1;
          NFO("selecting largest element.\n", NULL);
        } else if (!strcmp(ovldefault, "SMALLEST")) {
          mode = 2;
          NFO("selecting smallest element.\n", NULL);
        } else {
          NFO("selecting all elements.\n", NULL);
          mode = 0;
        }
      }
    
      bl_splitOverlap(&mspl, forest, mode, atribte);

      bl_annotationmultitrackDestruct (NULL, &mannot);
      bl_annotationtrackDestruct(NULL, annot);
      FREEMEMORY(NULL, annot);
     
      bl_intervalforestDestruct(forest);  
      FREEMEMORY(NULL, forest);
    }


    //bl_BEDwrite(ptr, mspl.noofitems, stdout);
    NFO("summarizing %d splits.\n", mspl.noofitems);
    bl_splitSummarize(stdout, &mspl, mintotalsplit, minmedianqual, nostrand); 

    for(i=0; i < splfns->noofvalues; i++) {
      bl_annotationtrackDestruct(NULL, spl[i]);
      FREEMEMORY(NULL, spl[i]);
    }
    FREEMEMORY(NULL, spl);
    FREEMEMORY(NULL, ovlconstraint.items);
    bl_annotationmultitrackDestruct(NULL, &mspl);
    manopt_destructarg(unflagged);
    FREEMEMORY(NULL, unflagged);

  } else {
    manopt_help(&prgset, "unknown program selected\n");
  }

  manopt_destructoptionset(&optset);
  manopt_destructoptionset(&prgset);
  manopt_destructarg(prg);
  FREEMEMORY(NULL, prg);
  FREEMEMORY(NULL, version);

  return 0;
}


