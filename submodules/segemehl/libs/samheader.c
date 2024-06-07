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
#include <stdint.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "version.h"
#include "info.h"
#include "samheader.h"
#include "samio.h"


/*------------------------------ bl_samwriteHeader -----------------------------
 *    
 * @brief read the header from a sam file
 * @author Steve Hoffmann 
 *   
 */

char*
bl_samwriteHeader(samheader_t *hdr, Uint binno, char order, char sep, char lf) {

  char *s=NULL;
  uint64_t nseqs, i, ngrps;

  bl_bsprintf(&s, "@HD%cVN:1.0",sep);
 
  if(order) {
    bl_bsprintf(&s, "%cSO:coordinate", sep);
  }

  bl_bsprintf(&s, "%c", lf);


  //templates, their ids and infos
  nseqs = hdr->nrnames;
  if(binno != -1) {
    bl_bsprintf(&s, "@SQ%cSN:%s%cLN:%d%c", 
        sep, hdr->rnames[binno], sep, hdr->rlens[binno], lf);
  } else { 
    for(i=0; i < nseqs; i++) {
      bl_bsprintf(&s,"@SQ%cSN:%s%cLN:%d%c", 
          sep, hdr->rnames[i], sep, hdr->rlens[i], lf);
    }
  }
  //read groups and infos 
  ngrps = hdr->nrgroups;
  for(i=0; i < ngrps; i++) {
    if(hdr->rgroupsinfo && hdr->rgroupsinfo[i]) {
     bl_bsprintf(&s, "@RG%cID:%s", sep, hdr->rgroups[i], sep);
     bl_bsprintf(&s, "%s%c", hdr->rgroupsinfo[i], lf);
    } else {
     bl_bsprintf(&s, "@RG%cID:%s%c", sep, hdr->rgroups[i], lf);
    }
  }

  
  bl_bsprintf(&s,"@PG%cID:segemehl%cVN:%s%cCL:%s", sep, sep, 
      hdr->version, sep, hdr->cmd); 

  return s;
}
/*------------------------------- bl_samHeader -------------------------------
 *    
 * @brief SAM header
 * @author Steve Hoffmann 
 *   
 */

  char*
bl_samHeader (void *space, segemehl_t *info, Uint binno)
{

  Uint i, size=0;
  Uint *seqlen;
  char sep, lf;
  char **seq;
  char *hdr=NULL;


  if(info->order){ 
    sep = 8;
    lf = 7;
  } else {
    sep = '\t';
    lf = '\n';
  }

  if(info->bins && binno != -1) {
    size = bl_fileBinsDomainsGetList(space, info->bins, &seq, &seqlen);
  } else { 
    if(info->fasta) {
      size = info->fasta->noofseqs;
      seq = ALLOCMEMORY(space, NULL, char**, size);
      seqlen = ALLOCMEMORY(space, NULL, Uint*, size);

      for(i=0; i < size; i++) {
        seq[i] = bl_fastaGetDescription(info->fasta, i); 
        seqlen[i] = bl_fastaGetSequenceLength(info->fasta, i); 
      }
    } else {
      size = 0;
      seq = NULL;
      seqlen = NULL;
    }
  }

  bl_bsprintf(&hdr,"@HD%cVN:1.0",sep);

  if(info->order) {
    bl_bsprintf(&hdr, "%cSO:coordinate", sep);
  }

  bl_bsprintf(&hdr,"%c",lf);

  for(i=0; i < size; i++) {
    /*
     * in case of binning, a null domain '*'
     * has been created by se_createChromDomains
     * for empty alignments
     *
     */
    if(seq[i][0] != '*') {
      bl_bsprintf(&hdr,"@SQ%cSN:%s%cLN:%d%c", sep, seq[i], sep, seqlen[i], lf);
    }
  }

  bl_bsprintf(&hdr,"@RG%cID:%s", sep, info->readgroupid);

  if(info->readgroupinfo) {
    bl_bsprintf(&hdr,"%s%c", info->readgroupinfo, lf);
  } else {
    bl_bsprintf(&hdr,"%c",lf);
  }

  bl_bsprintf(&hdr,"@PG%cID:segemehl", sep);
  bl_bsprintf(&hdr,"%cVN:%s", sep, VERSION);

  if(info->cmdline) {
    bl_bsprintf(&hdr,"%cCL:%s", sep, info->cmdline);
  }

  bl_bsprintf(&hdr,"%c",lf);

  if(info->order) { 
    sprintf(&hdr[strlen(hdr)-1], "%c", 29);
    bl_bsprintf(&hdr, "%c", '\n'); 
  } 

  FREEMEMORY(space, seq);
  FREEMEMORY(space, seqlen);

  return hdr;
}
/*------------------------------ bl_saminitHeader ------------------------------
 *    
 * @brief initalize header structure
 * @author Steve Hoffmann 
 *   
 */
  void
bl_saminitHeader (samheader_t *head)
{
  head->version = NULL;
  head->rnames = NULL;
  head->rlens = NULL;
  head->nrnames = 0;
  head->rgroups = NULL;
  head->rgroupsinfo = NULL;
  head->nrgroups = 0;
  head->cmd = NULL;

  return ;
}
/*---------------------------- bl_samdestructHeader ----------------------------
 *    
 * @brief destroy the header
 * @author Steve Hoffmann 
 *   
 */

  void
bl_samdestructHeader (samheader_t *head)
{
  Uint i;
  FREEMEMORY(NULL, head->version);
  FREEMEMORY(NULL, head->rlens);
  FREEMEMORY(NULL, head->cmd);
  for(i=0; i< head->nrnames; i++){
    FREEMEMORY(NULL, head->rnames[i]);
  }

  for(i=0; i< head->nrgroups; i++){
    FREEMEMORY(NULL, head->rgroups[i]);
    FREEMEMORY(NULL, head->rgroupsinfo[i]);
  }

  FREEMEMORY(NULL, head->rnames);
  FREEMEMORY(NULL, head->rgroups);
  FREEMEMORY(NULL, head->rgroupsinfo);

  return ;
}
/*----------------------------- bl_samaddReadGroup -----------------------------
 *    
 * @brief add a read group to a header
 * @author Steve Hoffmann 
 *   
 */

  void
bl_samaddReadGroup (samheader_t *head, char *id, char *info)
{

  head->rgroups = ALLOCMEMORY(NULL, head->rgroups, char**, head->nrgroups+1);
  head->rgroupsinfo = ALLOCMEMORY(NULL, head->rgroupsinfo, char**, head->nrgroups+1);
  head->rgroups[head->nrgroups] = id;
  head->rgroupsinfo[head->nrgroups] = info;
  head->nrgroups++;

  return ;
}
/*----------------------------- bl_samaddReadGroup -----------------------------
 *    
 * @brief add a reference sequence to the header
 * @author Steve Hoffmann 
 *   
 */

  void
bl_samaddReferenceSequence (samheader_t *head, char *name, uint64_t len)
{

  head->rnames = ALLOCMEMORY(NULL, head->rnames, char**, head->nrnames+1);
  head->rlens= ALLOCMEMORY(NULL, head->rlens, uint64_t, head->nrnames+1);
  head->rlens[head->nrnames] = len;
  head->rnames[head->nrnames] = name;
  head->nrnames++;

  return ;
}


/*------------------------ bl_samgetReadGroupFromHeader ------------------------
 *    
 * @brief extract the reference sequence name and its length
 * @author Steve Hoffmann 
 *   
 */

  char*
bl_samgetReadGroupFromHeader (char **xval, Uint ntags, char **rgname)
{

  Uint i;
  char checkname=0;
  char *info = NULL;

  for(i=0; i < ntags; i++) {
    if(!strncmp(xval[i],"ID:", 3)) {
      *rgname = bl_strdup(&xval[i][3]);
      checkname++;
    } else {
      bl_bsprintf(&info,"\t%s", xval[i]);
    }
  }
 

  return info;
}


/*------------------- bl_samgetReferenceSequencesFromHeader --------------------
 *    
 * @brief extract the reference sequence name and its length
 * @author Steve Hoffmann 
 *   
 */

  void
bl_samgetReferenceSequencesFromHeader (char **xval, Uint ntags, char**rname, 
    uint64_t *len)
{

  Uint i;
  char checkname=0, checklen=0;

  for(i=0; i < ntags; i++) {
    if(!strncmp(xval[i],"SN:", 3)) {
      *rname = bl_strdup(&xval[i][3]);
      checkname = 1;
    }
    if(!strncmp(xval[i],"LN:", 3)) {
      *len = strtoull(&xval[i][3], NULL, 10);
      checklen = 1;
    }
  }

  assert(checkname && checklen);

  return ;
}


/*------------------------------ bl_samgetHeader -------------------------------
 *    
 * @brief read the header from a sam file
 * @author Steve Hoffmann 
 *   
 */
samheader_t*
bl_samgetHeader(segemehl_t *nfo) {

  samheader_t *head;
  uint64_t nseqs, i, ngrps;

  head = ALLOCMEMORY(NULL, NULL, samheader_t, 1);
  bl_saminitHeader(head);

  bl_bsprintf(&head->cmd, "%s", nfo->cmdline);
  bl_bsprintf(&head->version, "%s",  VERSION);

  nseqs = nfo->fasta->noofseqs;
  head->nrnames = nseqs;
  head->rnames = ALLOCMEMORY(NULL, NULL, char*, nseqs);
  head->rlens = ALLOCMEMORY(NULL, NULL, uint64_t, nseqs);

  for(i=0; i < nseqs; i++) {
    bl_bsprintf(&head->rnames[i], "%s", bl_fastaGetDescription(nfo->fasta, i));
    head->rlens[i] = bl_fastaGetDescriptionLength(nfo->fasta, i);
  }

  ngrps = 1;
  head->nrgroups = ngrps;
  head->rgroups = ALLOCMEMORY(NULL, NULL, char*, ngrps);
  head->rgroupsinfo = ALLOCMEMORY(NULL, NULL, char*, ngrps);
  for(i=0; i < ngrps; i++) {
    bl_bsprintf(&head->rgroups[i], "%s", nfo->readgroupid);
    if(nfo->readgroupinfo) {
      bl_bsprintf(&head->rgroupsinfo[i], "%s", nfo->readgroupinfo);
    } else {
      head->rgroupsinfo[i] = NULL;
    }
  }


  return head;
}

/*--------------------------- bl_samparseHeaderLine------------------------------
 *    
 * @brief read the header from a sam file
 * @author Steve Hoffmann 
 *   
 */

  samheader_t*
bl_samparseHeaderLine (samheader_t *head, char *line)
{

  char *ptr, *saveptr;
  char *mycopy = bl_strdup(line); 
  char **xval=NULL;
  Uint ntags = 0, i=0;
  char *rname = NULL;
  char *rgroup = NULL;
  uint64_t rlen =0;


  if(line[0] != '@')
    return NULL;

  ptr = strtok_bl(mycopy, "\t", &saveptr); 

  while (ptr != NULL) {
    xval = ALLOCMEMORY(NULL, xval, char*, ntags+1); 
    xval[ntags] = bl_strdup(ptr);
    ptr = strtok_bl(NULL, "\t", &saveptr);
    ntags++;
  }


  //header line
  if(!strncmp(xval[0], "@HD", 3)) {
    //don't care
  } else if (!strncmp(xval[0], "@SQ", 3)) { 
    bl_samgetReferenceSequencesFromHeader(&xval[1], ntags-1, &rname, &rlen);
    head->rnames = ALLOCMEMORY(NULL, head->rnames, char*, head->nrnames+1);
    head->rlens = ALLOCMEMORY(NULL, head->rlens, uint64_t, head->nrnames+1);
    head->rnames[head->nrnames] = rname;
    head->rlens[head->nrnames] = rlen;
    head->nrnames++;
  } else if (!strncmp(xval[0], "@RG", 3)) { 
    char *info = bl_samgetReadGroupFromHeader(&xval[1], ntags-1, &rgroup);
    bl_samaddReadGroup(head, rgroup, info);
  } else if (!strncmp(xval[0], "@PG", 3)) { 
    //don't care
  } else if (!strncmp(xval[0], "@CO", 3)) {
    //don't care
  }

  for(i=0; i < ntags; i++) {
    FREEMEMORY(NULL, xval[i]);
  }

  FREEMEMORY(NULL, xval);
  FREEMEMORY(NULL, mycopy);
  return head;
}
/*----------------------------- bl_samdumpHeader -------------------------------
 *    
 * @brief dump the header
 * @author Steve Hoffmann 
 *   
 */

  void
bl_samdumpHeader (samheader_t  *head)
{
  uint64_t i;

  for(i=0; i < head->nrnames; i++) {
    fprintf(stderr, "found rname %s (%" PRIu64 ")\n", head->rnames[i], head->rlens[i]);
  }

  for(i=0; i < head->nrgroups; i++) {
    fprintf(stderr, "found read group %s\n", head->rgroups[i]);
    fprintf(stderr, "additional info %s\n", head->rgroupsinfo[i]);
  }

  return ;
}



samheader_t*
bl_samreadHeader(char* filename)
  //unsigned char gzip, struct access *index,)

{
  FILE *fp;
  off_t offset = 0;
  char ch;
  char *buffer;
  //  char *descrbuffer = NULL;
  //  char *seqbuffer = NULL;
  //  char *qualbuffer = NULL;
  //  char idchar=0;
  int ret=0;
  //  unsigned char desc = 0;

  //  unsigned char qualdesc = 0;
  //  unsigned char qual = 0;
  //  unsigned char seq = 0;
  unsigned char gzip = 0;
  struct gzidxfile *gzf = NULL;
  struct access * index = NULL;

  //  Uint descrlength = 0; 
  //  Uint seqlen = 0;
  Uint buffersize = MAXBUFFERSIZE;
  //  Uint n = startseq;
  Uint len = 0;  

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  if(gzip) {
    fp = fopen(filename, "rb");
    gzf = bl_initgzidxfile(fp, index, offset, MEDIUMCHUNK);
  } else {
    fprintf(stderr, "open normal.\n");
    fp = fopen(filename, "r");
  }

  if (fp == NULL) {
    NFO("Couldn't open file '%s': %d. Exit forced.\n", filename, errno);
    exit(-1);
  }

  if(offset > 0) {
    ret = fseeko(fp, offset, SEEK_SET);
    if (ret == -1) {
      NFO("fseeko failed for file %s. Exit forced.\n", filename);
      exit(-1);
    }
  }

  samheader_t *head;
  head = ALLOCMEMORY(NULL, NULL, samheader_t, 1);
  bl_saminitHeader(head);

  while((ch= (gzip) ? bl_getgzidxc(gzf) : getc(fp)) != EOF) {


    if(ch == '\n' && buffer){
      //finalize buffer
      buffer = ALLOCMEMORY(space, buffer, char, len+1);  
      buffer[len] = '\0'; 
      //process buffer 
      if(buffer[0] == '@') {
        bl_samparseHeaderLine(head, buffer);
      } else { 
        bl_samline2rec(buffer, len, head);
      }
      //free current buffer
      FREEMEMORY(NULL, buffer);
      //and allocate a new one (or erase the old);
      len = 0;
      buffer = ALLOCMEMORY(space, NULL, char, buffersize);
    }

    if(ch != '\n') {
      len++;
      if(len == buffersize-1) {
        buffersize = 2*buffersize+1;
        buffer = ALLOCMEMORY(space, buffer, char, buffersize);
      }
      buffer[len-1] = ch;
    }
  }
  //bl_samdumpHeader(head);

  fclose(fp);
  FREEMEMORY(NULL, buffer);
  return head;
}
