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
 *
 *   manout.c
 *   the new output manager. simply takes a mapping set, determines the
 *   mapping policy and selects the right output functions
 *
 *   @author Steve Hoffmann
 *   @email steve@bioinf.uni-leipzig.de
 *   @date 26.12.2012 02:20:31 CET
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <pthread.h>

#include "basic-types.h"
#include "segemehl.h"
#include "manoutformats.h"
#include "mappingqual.h"
#include "junctions.h"
#include "version.h"
#include "multicharseq.h"
#include "mapfrag.h"	//(skipped)
#include "fileBins.h"	//(skipped)
#include "fileio.h"	//(skipped)
#include "biofiles.h"	//(skipped)
#include "info.h"	//(skipped)
#include "memory.h"	//(skipped)
#include "samio.h"	//(skipped)
#include "mathematics.h"	//(skipped)
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "bamio.h"
/*-------------------------- se_printUnmappedFastx ---------------------------
 *    
 * @brief print unmapped fastx
 * @author Steve Hoffmann 
 *   
 */

  void
se_printUnmappedFastx(FILE* dev, char *desc, char *seq, char *qual, mapseed_t *seed)
{

  if(qual) {
    if(seed)
      fprintf(dev, "@%s ef:%d;if:%d %" PRIu64 ":%" PRIu64 " %" PRIu64 ":%" PRIu64 ":%d\n%s\n+%s\n%s\n",		
          desc, 
          seed->maxevalue, seed->maxinterval,
          seed->u, seed->seedlen,
          seed->refidx, seed->refpos, seed->rc, 
          seq, desc, qual);
    else

      fprintf(dev, "@%s ef:0;if:0 0:0 0:0:0\n%s\n+%s\n%s\n",		
          desc, 
          seq, desc, qual);
  } else {
    if(seed)
      fprintf(dev, ">%s ef:%d;if:%d %" PRIu64 ":%" PRIu64 " %" PRIu64 ":%" PRIu64 ":%d\n%s\n", 
          desc, 
          seed->maxevalue, seed->maxinterval,
          seed->u, seed->seedlen,
          seed->refidx, seed->refpos, seed->rc,
          seq); 
    else

      fprintf(dev, ">%s ef:0;if:0 0:0 0:0:0\n%s\n",		
          desc, 
          seq);
  }

  return ;
}
/*------------------------------ se_printUnmappedReads -------------------------------
 *    
 * @brief  print the unmapped reads to a fasta device for further analysis
 * @author Steve Hoffmann 
 *   
 */

void
se_printUnmappedReads(mappingset_t *s, fasta_t *reads, Uint k, 
    mapseed_t *seed, mapseed_t *mateseed, segemehl_t *nfo) {

  if(nfo->nomatchdev){ 
    //now dump to additional no-match-device if requested
    if(!bl_hasQueryMapping(s)) {
      char *seq  = bl_fastaGetSequenceNoClip(reads, k);
      char *qual = bl_fastaGetQualityNoClip(reads, k);
      char *desc = bl_fastaGetDescription(reads, k);


      if (nfo->threadno > 1) pthread_mutex_lock(nfo->mtx2);
      se_printUnmappedFastx(nfo->nomatchdev, desc, seq, qual, seed);     
      fflush(nfo->nomatchdev);
      if (nfo->threadno > 1) pthread_mutex_unlock(nfo->mtx2);
    } 

    if(bl_fastaHasMate(reads) && !bl_hasMateMapping(s)) {
      char *seq  = bl_fastaGetMateNoClip(reads, k);
      char *qual = bl_fastaGetMateQualityNoClip(reads, k);
      char *desc = bl_fastaGetMateDescription(reads, k);

      if (nfo->threadno > 1) pthread_mutex_lock(nfo->mtx2);
      se_printUnmappedFastx(nfo->nomatchdev, desc, seq, qual, mateseed);  
      fflush(nfo->nomatchdev);
      if (nfo->threadno > 1) pthread_mutex_unlock(nfo->mtx2);
    }
  }
}


/*-------------------------------- se_output ---------------------------------
 *    
 * @brief manage the output
 * @author Steve Hoffmann 
 *   
 */
void
se_updateStats(mappingset_t *s, fasta_t *reads, segemehl_t *nfo) {
  Uint i;

  if(nfo->threadno > 1) pthread_mutex_lock(nfo->mtx5);
  //STATS
  nfo->stats->total +=1;
  if(bl_fastaHasMate(reads)) {
    nfo->stats->total +=1;
  }

  if(bl_mappingsetHasPaired(s)) {
 
    //mapped in pair (+2 mapped reads)
    nfo->stats->mapped+=2;
    nfo->stats->paired+=1;

    if(bl_hasMultiplePairedMappings(s)) {
      //multiple pair (+2 uniquely mapped reads) 
      nfo->stats->multiplemapped+=2;
      nfo->stats->multiplepaired+=1;
    } else {
      //uniq pair (+2 multiple mapped reads)
      nfo->stats->uniquemapped+=2;
      nfo->stats->uniquepaired+=1;
    }
  } else {
   
    //unmapped in pair
    if(bl_hasQueryMapping(s)) {
      //query mapped (+1 mapped reads)
      nfo->stats->mapped+=1;
      nfo->stats->singlequerymapped+=1;
      if(bl_hasMultipleQueryMappings(s)) {
        nfo->stats->multiplemapped+=1;
      } else {
        nfo->stats->uniquemapped+=1;
      }
    } else {
      nfo->stats->unmapped+=1;
    }
    if(bl_hasMateMapping(s)) {
      //mate mapped (+1 mapped reads)
      nfo->stats->mapped+=1;
      nfo->stats->singlematemapped+=1;
      if(bl_hasMultipleMateMappings(s)) {
        nfo->stats->multiplemapped+=1;
      } else {
        nfo->stats->uniquemapped+=1;
      }
    } else if(bl_fastaHasMate(reads)) {
      nfo->stats->unmapped+=1;
    }
  }


  for(i=0; i < s->n; i++) {
    char isquerysplit = bl_isSplitMappingQM(&s->elem[i], 0);
    char ismatesplit = bl_isSplitMappingQM(&s->elem[i], 1);
    if(isquerysplit && ismatesplit) {
      //split pair
      nfo->stats->splitpair +=1;
    } else if(isquerysplit || ismatesplit) {
      //split single
      nfo->stats->singlesplit +=1;
    }
  }

  if(nfo->threadno > 1) pthread_mutex_unlock(nfo->mtx5);

}


/*-------------------------------- se_output ---------------------------------
 *    
 * @brief manage the output
 * @author Steve Hoffmann 
 *   
 */


void
se_output(mappingset_t *s, fasta_t * reads, unsigned int k, 
    mapseed_t *seed, mapseed_t *mateseed, Suffixarray *arr, karlin_t* stats, 
    segemehl_t *nfo) {

  unsigned int i, noofqueries=0, noofmates=0;
  unsigned querycnt=0, matecnt=0;

  samlist_t *samlist;
  MultiCharSeq *mseq=nfo->seq;

  /* remove reads with a bad accuracy */
#ifndef LOCAL
  bl_removeBadMatesAcc(s, bl_fastaGetSequenceLength(reads, k), 
      bl_fastaGetMateLength(reads, k), nfo->accuracy);

  /* remove reads that are heavily softclipped */
  bl_removeBadMatesCov(s, bl_fastaGetSequenceLength(reads, k),
      bl_fastaGetMateLength(reads, k), nfo->minsplicedaligncover);

  /*
   fprintf(stderr, "afer removing bad mates:\n"); 
   bl_dumpMappingSet(stderr, s);
  */

#else


  int minscore = MAX(nfo->localparams[0], 
      ((nfo->localparams[1]*bl_fastaGetSequenceLength(reads, k))/100)*nfo->scores[0]);
  int mateminscore = MAX(nfo->localparams[0], 
      ((nfo->localparams[1]*bl_fastaGetMateLength(reads, k))/100)*nfo->scores[0]);

  bl_removeBadMatesScr(s, bl_fastaGetSequenceLength(reads, k), 
      bl_fastaGetMateLength(reads, k), nfo->scores, nfo->scores[2], minscore, 
      mateminscore);
   
  //remove reads that are heavily softclipped
  bl_removeBadMatesCov(s, bl_fastaGetSequenceLength(reads, k),
      bl_fastaGetMateLength(reads, k), nfo->minsplicedaligncover);
#endif

  if(bl_fastaHasMate(reads)) {
    
    if(bl_mappingsetHasPaired(s)) { 
      //proper pairs available
      if(nfo->bestonly) { 
        bl_removeSuboptimalMapping(s, nfo->scores, nfo->scores[2]);
        /*
        fprintf(stderr, "afer removing suboptimal:\n"); 
        bl_dumpMappingSet(stderr, s);
        */
      }
      //still paired?
      if(bl_mappingsetHasPaired(s)) { 
        bl_removeUnpairedMapping(s);
      }

      se_updateStats(s, reads, nfo); 
    } else { 
      //no proper pairs available
       
      if(nfo->bestonly) { 
        //eliminate mates and queries individually
        //fprintf(stdout, "removing mates individually\n");
        bl_removeSuboptimalMappingQM(s, nfo->scores, nfo->scores[2], 0);
        bl_removeSuboptimalMappingQM(s, nfo->scores, nfo->scores[2], 1);
      }

      se_updateStats(s, reads, nfo); 
      //if there are only two mappings left and it is one query and one mate
      //we are allowed to merge them into one
      bl_mergeMappings(s);
    }
  } else {

    se_updateStats(s, reads, nfo); 
    //remove suboptimal if best only is chosen
    if(nfo->bestonly) { 
      bl_removeSuboptimalMapping(s, nfo->scores, nfo->scores[2]);
    }
  }

  bl_countMultipleMappings(s, &noofqueries, &noofmates);

  /*CALCULATE MAPPING QUAL IF NO VALUE IS SET*/
  if(!nfo->mappingqual) { 
    if(bl_fastaHasQuality(reads)) { 
      type3mappingset(s,bl_fastaGetDescription(reads, k));
    }
    longestmatchqual(s, arr->numofsuffixes, stats);
  }

  //Write BED files for splits
  if(nfo->split) { 
    bl_getSpliceJunctionsFromMappingSet(s, mseq, bl_fastaGetDescription(reads,k), nfo);
  }

  //in case there is one uniq alignment of one mate and the other is not mapped
  //RNEXT is filled with the artificial position of the unmapped mate
#ifdef SORTEDUNMAPPED
  char unmappedtornext = (bl_fastaHasMate(reads) && s->n==1 
      && (!bl_hasQueryMapping(s) || !bl_hasMateMapping(s)));
#else
  char unmappedtornext = 0; 
#endif

  for(i=0; i < s->n; i++) {

    samlist = bl_samgetSamList(reads, k, &s->elem[i], mseq, (s->n > 1), 
        noofqueries, noofmates, querycnt, matecnt, unmappedtornext, nfo);
    /*
     * output to device
     *
     */
    bl_samprintSamlist(samlist, &s->elem[i], nfo);

    querycnt += (bl_isQueryMapping(&s->elem[i])) ? 1 : 0;
    matecnt += (bl_isMateMapping(&s->elem[i])) ? 1 : 0;
   
    //use samlist information to store mate info in empty aligns - if there is a uniq hit
    if(s->n == 1) { 
      if(bl_fastaHasMate(reads)&& !bl_hasQueryMapping(s)) {
        char *seq  = bl_fastaGetSequenceNoClip(reads, k);
        char *qual = bl_fastaGetQualityNoClip(reads, k);
        char *desc = bl_fastaGetDescription(reads, k);
        char *nextrname = samlist->recs[0].rname;
        //should always be > 0, because the alignment exists uniqly
        uint64_t nextrpos = samlist->recs[0].pos;
        char nextrc = samlist->recs[0].flag & 0x10;
        char unmappedmate = !bl_hasMateMapping(s);
        /*
         * output to device
         *
         */
        bl_samprintEmptyAlign (desc, seq, qual, 1, 1, unmappedmate, nextrname, 
            nextrpos, nextrc, 0, 0, seed, nfo);
      } 

      if(bl_fastaHasMate(reads) && !bl_hasMateMapping(s)) {
        char *seq  = bl_fastaGetMateNoClip(reads, k);
        char *qual = bl_fastaGetMateQualityNoClip(reads, k);
        char *desc;
        if(reads->checkid) {
          desc = bl_fastaGetMateDescription(reads, k);
        } else {
          desc = bl_fastaGetDescription(reads, k);
        }

        char *nextrname = samlist->recs[0].rname;
        //should always be > 0, because the alignment exists uniqly
        uint64_t nextrpos = samlist->recs[0].pos;
        char nextrc = samlist->recs[0].flag & 0x10; 
        char unmappedmate = !bl_hasQueryMapping(s);
        /*
         * output to device
         *
         */
        bl_samprintEmptyAlign (desc, seq, qual, 1, 0, unmappedmate, nextrname, 
            nextrpos, nextrc, 0, 0, mateseed, nfo);
      }
    }

    bl_samdestructSamList (samlist);
    FREEMEMORY(NULL, samlist);
  }

  /*
   * dump unaligned reads w/o (uniqely) aligned mate
   *
   */
  if(s->n != 1) {

    if(!bl_hasQueryMapping(s)) {
      char *seq  = bl_fastaGetSequenceNoClip(reads, k);
      char *qual = bl_fastaGetQualityNoClip(reads, k);
      char *desc = bl_fastaGetDescription(reads, k);
      char hasPaired = bl_fastaHasMate(reads);
      char *nextrname = NULL;
      int64_t nextrpos = 0;
      char nextrc = 0;
      char unmappedmate = !bl_hasMateMapping(s);
      /*
       * output to device
       *
       */
      bl_samprintEmptyAlign (desc, seq, qual, hasPaired, 1, unmappedmate, nextrname, 
          nextrpos, nextrc, 0, 0, seed, nfo);
    }

    if(bl_fastaHasMate(reads) && !bl_hasMateMapping(s)) {
      char *seq  = bl_fastaGetMateNoClip(reads, k);
      char *qual = bl_fastaGetMateQualityNoClip(reads, k); 
      char *desc;
        if(reads->checkid) {
          desc = bl_fastaGetMateDescription(reads, k);
        } else {
          desc = bl_fastaGetDescription(reads, k);
        }

      char hasPaired = bl_fastaHasMate(reads);
      char *nextrname = NULL;
      int64_t nextrpos = 0;
      char nextrc = 0;
      char unmappedmate = !bl_hasQueryMapping(s);
      /*
       * output to device
       *
       */
      bl_samprintEmptyAlign (desc, seq, qual, hasPaired, 0, unmappedmate, nextrname, 
          nextrpos, nextrc, 0, 0, mateseed, nfo);
    }
  }

  se_printUnmappedReads(s, reads, k, seed, mateseed, nfo);
 
}

/*----------------------------- se_initChromBins -----------------------------
 *    
 * @brief set up bins for chromosomes
 * @author Steve Hoffmann 
 *   
 */

  bl_fileBinDomains_t*
se_createChromDomains (void *space, fasta_t *f, Uint avgbins, Uint maxbins, 
    char *filetemplate, Uint tmplen)
{
  bl_fileBinDomains_t* domains;
  char **desc;
  Uint *size;
  Uint i, no, total=0;

  no = f->noofseqs+1;
  if(no > maxbins) return NULL;

  desc = ALLOCMEMORY(space, NULL, char*, no);
  size = ALLOCMEMORY(space, NULL, Uint, no);

  /*
   * a null domain '*' is created for empty
   * alignments, this domain needs to be taken
   * care of in bl_samHeader
   *
   */
  desc[0] = "*";
  size[0] = 1;
  total = 1;

  for(i=1; i < no; i++) {
    desc[i] = bl_fastaGetDescription(f, i-1);
    size[i] = bl_fastaGetSequenceLength(f, i-1);
    total += size[i];
  } 

  domains = bl_fileBinsDomainsInit(space, desc, size, no, total,  
      avgbins, maxbins, filetemplate, tmplen);

  FREEMEMORY(space, desc);
  FREEMEMORY(space, size);

  return domains;
}

/*----------------------------- se_initChromBins -----------------------------
 *    
 * @brief set up bins for chromosomes
 * @author Steve Hoffmann 
 *   
 */

  bl_fileBins_t*
se_createChromBins (void *space, fasta_t *f, int maxbins, char *template, 
    Uint tmplen)
{
  bl_fileBins_t* bins;
  char **desc;
  Uint i, no;

  no = f->noofseqs;
  if(no > maxbins) return NULL;

  bins = ALLOCMEMORY(space, NULL, bl_fileBins_t, 1);
  desc = ALLOCMEMORY(space, NULL, char*, no);
  bl_fileBinsInit(space, bins);

  for(i=0; i < no; i++) {
    desc[i] = bl_fastaGetDescription(f, i);
  }

  bl_fileBinsAdd(space, bins, no, bl_fileBinCClassAssign, desc, NULL, 
      template, tmplen);

  FREEMEMORY(space, desc);
  return bins;
}

/*-------------------------------- se_createBisulifteBins ---------------------
 *
 * @brief set up bin domains for matching runs and threads,
 *        domain names are simply 0...(noofdomains-1) as strings
 * @author Christian Otto
 *
 */

bl_fileBinDomains_t*
se_createBisulfiteBins (void *space, Uint noofdomains,
    Uint threadno, char *filetemplate, Uint tmplen){
  Uint i, j;
  bl_fileBinDomains_t *domains;
  bl_fileBinClass_t *class;

  domains = ALLOCMEMORY(space, NULL, bl_fileBinDomains_t, 1);
  domains->noofdomains = noofdomains;
  domains->exp = 0;
  domains->domain = ALLOCMEMORY(space, NULL, bl_fileBinDomain_t, noofdomains);

  for (i = 0; i < noofdomains; i++){    
    domains->domain[i].domainsize = threadno;    
    domains->domain[i].domainname = ALLOCMEMORY(space, NULL, char, log10(i+1) + 3);
    snprintf(domains->domain[i].domainname, log10(i+1)+2, "%d", i);

    bl_fileBinsInit(space, &domains->domain[i].bins);
    bl_fileBinsAdd(space, &domains->domain[i].bins, threadno, NULL, NULL, NULL,
        filetemplate, tmplen);

    domains->domain[i].bins.noofbins = threadno;

    for (j = 0; j < domains->domain[i].bins.noofbins; j++){
      class = ALLOCMEMORY(space, NULL, bl_fileBinClass_t, 1);
      class->start = j;
      class->end = j;
      class->classname = NULL;
      domains->domain[i].bins.b[j].id = class;
    }
  }
  /*
     DBG("domains: noofdomains=%u\texp=%u\n", domains->noofdomains, domains->exp);
     for (i = 0; i < domains->noofdomains; i++){
     DBG("domain %u: domainname=%s\tdomainsize=%u\tnoofbins=%u\n", i,
     domains->domain[i].domainname, domains->domain[i].domainsize,
     domains->domain[i].bins.noofbins);
     for (j = 0; j < domains->domain[i].bins.noofbins; j++){
     DBG("bin %u: filename=%s\tstart=%llu\tend=%llu\n", j, domains->domain[i].bins.b[j].fname,
     domains->domain[i].bins.b[j].id->start, domains->domain[i].bins.b[j].id->end);
     }
     }*/
  return domains;

}
/*--------------------------------- se_header ----------------------------------
 *    
 * @brief  compile the segemehl header
 * @author Steve Hoffmann 
 *   
 */

samheader_t*
se_header(segemehl_t *nfo, Uint binno) {
  const char *rginfodefault = "\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina\0";
  samheader_t *head;
  char *tempgroupdefaults = NULL;

  /*
   *  update read groups first
   *
   */
  if(nfo->readgroupfile) {
    head = bl_samreadHeader(nfo->readgroupfile);
  } else  {
    //no readgroup file provided
    head = ALLOCMEMORY(NULL, NULL, samheader_t, 1);
    bl_saminitHeader(head);
    char *tempid = NULL;

    if(nfo->readgroupid) {
      bl_bsprintf(&tempid, "%s", nfo->readgroupid); 
    } else {
      MSG("assigning all reads to default read group 'A1'.\n");
      bl_bsprintf(&tempid, "A1"); 
    }

    NFO("additional read group default values '%s'\n", rginfodefault);
    bl_bsprintf(&tempgroupdefaults, "%s", rginfodefault);
    bl_samaddReadGroup(head, tempid, tempgroupdefaults);
  }

  if(head->nrgroups != 1) {
    NFO("read group file has %d IDs. Exactly 1 required.\n", head->nrgroups);
    exit(EXIT_FAILURE);
  } else if(!head->rgroups[0]) {
    NFO("bad read group id '%s' ('%s').\n", head->rgroups[0], nfo->readgroupfile);
    exit(EXIT_FAILURE);
  } else { 
    nfo->readgroupid = head->rgroups[0];
    nfo->readgroupinfo = head->rgroupsinfo[0];
  }
  /*
   * update chromosomes 
   *
   */
  char **seq;
  Uint *seqlen, size,i;

  if(nfo->bins && binno != -1) {
    size = bl_fileBinsDomainsGetList(NULL, nfo->bins, &seq, &seqlen);
  } else { 
    if(nfo->fasta) {
      size = nfo->fasta->noofseqs;
      seq = ALLOCMEMORY(NULL, NULL, char**, size);
      seqlen = ALLOCMEMORY(NULL, NULL, Uint*, size);
    } else {
      size = 0;
      seq = NULL;
      seqlen = NULL;
    }
  }
  for(i=0; i < size; i++) {
    char *tempname = NULL;
    seq[i] = bl_fastaGetDescription(nfo->fasta, i); 
    seqlen[i] = bl_fastaGetSequenceLength(nfo->fasta, i); 
    bl_bsprintf(&tempname, "%s", seq[i]);
    bl_samaddReferenceSequence (head, tempname, seqlen[i]);
  }

  FREEMEMORY(NULL, seq);
  FREEMEMORY(NULL, seqlen);

  bl_bsprintf(&head->version,"%s",VERSION);
  bl_bsprintf(&head->cmd,"%s",nfo->cmdline);

  NFO("reads assigned to read group '%s'\n", nfo->readgroupid);
  return head;
}
/*------------------------------ openOutputDevices -----------------------------
 *    
 * @brief select the output device for matches and print header
 * @author Steve Hoffmann 
 *   
 */
void se_openOutputDevices(void *space, segemehl_t *info) {

  char *buffer = NULL;
  char *qfbase;

  if(!info->bam) {
    if(info->outfn) {
      info->dev=fopen(info->outfn, "wb"); //w -> wb
      setvbuf(info->dev, NULL, _IOFBF, 524288);
      NFO("opening sam file '%s'.\n", info->outfn);
      if (info->dev == NULL) {
        NFO("Couldn't open file '%s'. Exit forced.\n", info->outfn);
        exit(EXIT_FAILURE);
      }
    }

    if(!info->order) {
      buffer = bl_samwriteHeader(info->head, -1, 0, '\t', '\n');
      fprintf(info->dev, "%s\n", buffer);
    } else {
      if(!info->outfn) {
        NFO("For sorting, an output filename is needed.\n", NULL);
        exit(EXIT_FAILURE);
      }
      /*
       * effectively this is a place holder
       */
      buffer = bl_samwriteHeader(info->head, -1, 1, 8, 7);
      buffer[strlen(buffer)-1]=29;
    }

    FREEMEMORY(space, buffer);

  } else {

    if(info->outfn) {
      NFO("opening bam file '%s'.\n", info->outfn);
      info->bamdev = bl_bamOpenFile(info->outfn); 
      info->bamhdr = bl_bamGetHeader(info->head, -1);
      int ret = sam_hdr_write(info->bamdev, info->bamhdr);
      if(ret < 0) {
        NFO("error writing header to bam.\n", NULL);
      }
    } else {
      NFO("BAM will be written to stdout.\n", NULL);
      info->bamdev = bl_bamOpenFile("-"); 
      info->bamhdr = bl_bamGetHeader(info->head, -1);
      int ret = sam_hdr_write(info->bamdev, info->bamhdr);
      if(ret < 0) {
        NFO("error writing header to bam.\n", NULL);
      }

    }
  }

  if(info->split) {

    if (info->splitfilebasename) {
      qfbase = info->splitfilebasename;
    } else if (info->outfn){ 
      qfbase = info->outfn;
    } else {
      qfbase = bl_basename(info->queryfilename); 
    }

    info->multisplitfilename = bl_changefilesuffix(qfbase, "mult.bed");
    info->singlesplitfilename = bl_changefilesuffix(qfbase, "sngl.bed");
    info->transsplitfilename = bl_changefilesuffix(qfbase, "trns.txt");

    NFO("writing multi splits to '%s'\n", info->multisplitfilename);
    NFO("writing sngle splits to '%s'\n", info->singlesplitfilename);
    NFO("writing trans splits to '%s'\n", info->transsplitfilename);
    /*
     * set file buffers to a larger size
     *
     */
    info->multisplitdev = fopen(info->multisplitfilename, "wb");
    setvbuf(info->multisplitdev, NULL, _IOFBF, 524288);
    info->singlesplitdev = fopen(info->singlesplitfilename, "wb");
    setvbuf(info->singlesplitdev, NULL, _IOFBF, 524288);
    info->transsplitdev = fopen(info->transsplitfilename, "wb");
    setvbuf(info->transsplitdev, NULL, _IOFBF, 524288);

    fprintf(info->multisplitdev, 
        "track name=\"MultiSplit:%s\" description=\"segemehl multi pred for %s\" visibility=2 itemRgb=\"On\"\n", 
        info->readgroupid, info->queryfilename);

    fprintf(info->singlesplitdev, 
        "track name=\"SingleSplit:%s\" description=\"segemehl sngl pred for %s\" visibility=2 itemRgb=\"On\"\n", 
        info->readgroupid, info->queryfilename);
  }

  return;
}

/*----------------------------- closeOutputDevices ------------------------------
 *    
 * @brief select the output device for matches and print header
 * @author Steve Hoffmann 
 *   
 */

void
se_closeOutputDevices(void *space, segemehl_t *info) {
  uint64_t i,j;
  uint64_t dmno;
  char **header;

  if(!info->bam) {
    if(info->outfn && !info->bins) {
      NFO("closing output file '%s'.\n", info->outfn);
      fclose(info->dev);
      /*
       * to ensure proper sorting with header at the top
       * we simply convert header s.t. it will be always top
       */
      if(info->order) {
        NFO("sorting output file.\n", NULL);
        /* get header */
        char *buffer = bl_samwriteHeader(info->head, -1, 0, '\t', '\n');
        /* make sure heade placeholder will be sorted first */
        bl_freplacestr(info->outfn, "\000\t", 2, 29);
        /* sort file */
        NFO("starting sort.\n", NULL);
        bl_UnixSort(space, info->outfn, SORT[info->rep_type], SORTDELIM);
        /* write header back*/
        NFO("re-writing header to '%s'.\n", info->outfn);
        bl_freplacestr(info->outfn, buffer, strlen(buffer), '\n');
        /* clean up */
        FREEMEMORY(NULL, buffer);
      }
      /*
       * expand alignments in ordered and merged data
       */
      if (info->align && (info->order || info->bisulfitemerging)){
        NFO("Expanding alignments in '%s'.\n", info->outfn);
        bl_freplacearr(info->outfn, "\007","\n", 1, EOF);
      }
      //bl_freplacearr(info.outfile, "\007\010","\n\t", 2, 29);

    } else if(info->bins) {
      NFO("closing output file bins.\n", NULL);
      bl_fileBinDomainsCloseAll(info->bins);

      if(info->order) {
        NFO("sorting output file bins.\n", NULL);
        bl_fileBinDomainsUnixSort(NULL, info->bins, SORTBIN[info->rep_type], SORTDELIM);
      }
      if(info->align && (info->order || info->bisulfitemerging)) {
        MSG("Expanding alignments in all bins.\n");
        for(i=0; i < info->bins->noofdomains; i++) {
          for(j=0; j < info->bins->domain[i].bins.noofbins; j++) {
            bl_freplacearr(info->bins->domain[i].bins.b[j].fname, "\007",
                "\n", 1, EOF);
          }
        }
      }

      dmno = bl_fileBinsDomainsGetNoOfNames(NULL, info->bins); 
      header = ALLOCMEMORY(space, NULL, char**, dmno);
      for(i=0; i < dmno; i++) {
        header[i]= bl_samwriteHeader(info->head, i, info->order, '\t', '\n'); 
      }

      bl_fileBinDomainsMerge(NULL, info->bins, info->filebinbasename, 
          strlen(info->filebinbasename), "sam", 3, header, 1);

      for(i=0; i < dmno; i++) {
        FREEMEMORY(space, header[i]);
      }

      FREEMEMORY(NULL, header);

      bl_fileBinDomainsDestruct(NULL, info->bins);
      FREEMEMORY(space, info->bins);
    }

  } else {
    NFO("closing bam file.\n", NULL);
    int ret = hts_close(info->bamdev);
    if(ret < 0) {
      NFO("error closing bam file.\n", NULL);
      exit(EXIT_FAILURE);
    }

    bam_hdr_destroy(info->bamhdr);
  }

  if(info->nomatchname != NULL) {
    fclose(info->nomatchdev);
  }

  if(info->split) {
    fclose(info->singlesplitdev);
    fclose(info->multisplitdev);
    fclose(info->transsplitdev);
  }

  return;
}
/*-------------------------------- se_printStats ---------------------------------
 *    
 * @brief print mapping stats
 * @author Steve Hoffmann 
 *   
 */

void se_printStats(FILE *dev, segemehl_t *nfo) {
  mappingstats_t *s = nfo->stats;
  double p;
   
  fprintf(dev,"\ttotal\tmapped\t(%%)\t");
  fprintf(dev,"unique\t(%%)\tmulti\t(%%)\tsplit\t(%%)\n"); 

  fprintf(dev, "all\t");

  p = ((double)s->unmapped/s->total)*100.0;
  //fprintf(dev, "%"PRIu64"\t%"PRIu64"\t%.2f%%", s->total, s->unmapped, p);
  fprintf(dev, "%"PRIu64"\t", s->total);

  p = ((double)s->mapped/s->total)*100.0;
  fprintf(dev, "%"PRIu64"\t%.2f%%\t", s->mapped, p);

  p = ((double)s->uniquemapped/s->total)*100.0;
  fprintf(dev, "%"PRIu64"\t%.2f%%\t", s->uniquemapped, p);

  p = ((double)s->multiplemapped/s->total)*100.0;
  fprintf(dev, "%"PRIu64"\t%.2f%%\t", s->multiplemapped, p);

  uint64_t totalsplits = s->splitpair*2+s->singlesplit;
  p = ((double)totalsplits/s->total)*100.0;
  fprintf(dev, "%"PRIu64"\t%.2f%%\n", totalsplits, p);

  if(s->paired) {
    fprintf(dev, "pair\t");

    uint64_t totalpairs = s->total/2;
    uint64_t unmappedpairs = s->singlequerymapped+s->singlematemapped;
    p = ((double)unmappedpairs/totalpairs)*100.0;
    //fprintf(dev, "%"PRIu64"\t%"PRIu64"\t%.2f%%", totalpairs, unmappedpairs, p);
    fprintf(dev, "%"PRIu64"\t", totalpairs);

    p = ((double)s->paired/totalpairs)*100.0;
    fprintf(dev, "%"PRIu64"\t%.2f%%\t", s->paired, p);

    p = ((double)s->uniquepaired/totalpairs)*100.0;
    fprintf(dev, "%"PRIu64"\t%.2f%%\t", s->uniquepaired, p);

    p = ((double)s->multiplepaired/totalpairs)*100.0;
    fprintf(dev, "%"PRIu64"\t%.2f%%\t", s->multiplepaired, p);
 
    p = ((double)s->splitpair/totalpairs)*100.0;
    fprintf(dev, "%"PRIu64"\t%.2f%%\n", s->splitpair, p);
  }

  return;
}
