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
 *  bamio.c
 *  
 *
 *  @author 
 *  @email 
 *  @date 04/01/2018 20:04:30 CEST
 *  
 */

#include "bamio.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <pthread.h>
#include "basic-types.h"
#include "memory.h"
#include "mathematics.h"
#include "info.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "stringutils.h"
#include "filebuffer.h"
#include "samio.h"
#include "segemehl.h"

/*
 * not a G -> H
 */


const unsigned char seq_nt16_GHN_table[256] = {
  
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    0,0,1,0, 2,2,2,2, 2,2,2,2, 2,0/*=*/,2,2,
    2,0,0,0, 0,2,2,1, 0,2,2,0, 2,0,2,2,
    2,2,0,0, 0,2,0,0, 2,0,2,2, 2,2,2,2,
    2,0,0,0, 0,2,2,1, 0,2,2,0, 2,0,2,2,
    2,2,0,0, 0,2,0,0, 2,0,2,2, 2,2,2,2,

    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
};

const char seq_nt16_str_comp[] = "=TGKCYSBAWRDMHVN";

const unsigned char seq_nt16_C_context_table[3][3] = { 
    {0,1,2},  //HH, HG, HN
    {3,4,5},  //GH, GG, GN
    {6,7,8}   //NH, NG, NN
 };

const char seq_nt16_C_context_str[9][4] = {"CHH", "CHG", "CHN", "CGH", "CGG", "CGN", "CNH", "CNG", "CNN"}; 

const char seq_nt16_C_context_char[9] =   {'h'  ,   'x',   'u',   'z',  'z',    'z',   'u',   'u',   'u'};

bl_bam_methyl_master_t*
bl_bamInitMaster(bl_bam_methyl_master_t *ms, bam_hdr_t *hdr, 
    uint8_t isthreaded, pthread_mutex_t *mtx) {

  ms->hdr = hdr;
  ms->last_tid = 0;
  ms->last_beg = 0;
  ms->last_len = 0;
  ms->mtx = mtx;
  ms->isthreaded = isthreaded;

  return ms;
}

uint32_t
bl_bamCrossSectionCoverage(bam_cs_data_t *csd, uint8_t allowdel) {
  uint32_t i, cov=0;

  if(allowdel) return csd->n;

  for(i=0; i < csd->n; i++){
    if(BAM_CS_GET_NT(csd->d[i])) cov++;
  }

  return cov;
}

void
bl_bamCrossSectionCount(bam_cs_data_t *csd, uint32_t cnt[]) {

  uint32_t i;
  for(i=0; i < csd->n; i++) {
    cnt[BAM_CS_GET_NT(csd->d[i])]++;
  } 
}

void
bl_bamCrossSectionBSCount(bam_cs_data_t *csd, uint32_t cnt[]) {

  uint32_t i;
  for(i=0; i < csd->n; i++) {
    if ((csd->c == BAM_CS_NT_C && BAM_CS_GET_CP(csd->d[i]) == 1) ||
        (csd->c == BAM_CS_NT_G && BAM_CS_GET_CP(csd->d[i]) == 2)){
      cnt[BAM_CS_GET_NT(csd->d[i])]++;
    }
  } 
}


uint32_t
bl_bamCrossSectionBSCoverage(bam_cs_data_t *csd) {
  uint32_t i,c=0;

  for (i = 0; i < csd->n; i++){ 
    if ((csd->c == BAM_CS_NT_C && BAM_CS_GET_CP(csd->d[i]) == 1) ||
        (csd->c == BAM_CS_NT_G && BAM_CS_GET_CP(csd->d[i]) == 2)){
      if(BAM_CS_GET_NT(csd->d[i]))
        c++;
    }
  }

  return c;
}

uint8_t
bl_bamIsBisulfiteTarget(uint8_t ref) {
  switch(ref){
    case BAM_CS_NT_G:
      return METHYL_REV_STRAND;
      break;
    case BAM_CS_NT_C:
      return METHYL_FWD_STRAND;
      break;
    default:
      return 0;
  }
}

uint8_t 
bl_bamGetBSBase(uint8_t strand, unsigned char conv) {

  if(!strand) {
    return BAM_CS_NT_N;
  }

  if(!conv){
    return (strand == METHYL_FWD_STRAND) ? BAM_CS_NT_C : BAM_CS_NT_G;
  } else {
    return (strand == METHYL_FWD_STRAND) ? BAM_CS_NT_T : BAM_CS_NT_A;
  }

  return BAM_CS_NT_N;
}

uint32_t
bl_bamGetBSCount(bam_cs_data_t *csd, uint8_t strand, unsigned char conv) {
  uint32_t _cnt[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  uint8_t bs=0;

  if(!strand) return 0;
  bl_bamCrossSectionBSCount(csd, _cnt); 
  bs = bl_bamGetBSBase(strand, conv);

  return _cnt[bs];
}



void
bl_bamCrossSectionGetContext2(bam_cs_t *cs, uint32_t k, uint8_t strand, char* ctxt) {  

  if(strand == METHYL_FWD_STRAND) {     
    ctxt[0] = seq_nt16_str[cs->x[k].c];
    if(k+1>= cs->n) {
      ctxt[1] = seq_nt16_str[cs->next];
    } else {
      ctxt[1] = seq_nt16_str[cs->x[k+1].c];
    }
  } else if (strand == METHYL_REV_STRAND) {
    ctxt[0] = seq_nt16_str_comp[cs->x[k].c];
    if(k == 0) {
      ctxt[1]  = seq_nt16_str_comp[cs->prev];
    } else {
      ctxt[1] = seq_nt16_str_comp[cs->x[k-1].c];
    }
  }

  return;
}

unsigned char
bl_bamGetContextCode(char *rseq, uint64_t pos, uint8_t strand) {
  uint8_t i1=2, i2=2; 

  if(strand == METHYL_FWD_STRAND) {     

    if(rseq[pos+1] != 0) {
      i1 = seq_nt16_GHN_table[(int)rseq[pos+1]];
      if(rseq[pos+2] != 0) {
        i2 = seq_nt16_GHN_table[(int)rseq[pos+2]];
      }
    }

  } else if (strand == METHYL_REV_STRAND) {     
    if(pos > 0) {
      i1 = seq_nt16_str_comp[seq_nt16_table[(int)rseq[pos-1]]];
      i1 = seq_nt16_GHN_table[i1];
      if(pos > 1) {
        i2 = seq_nt16_str_comp[seq_nt16_table[(int)rseq[pos-2]]];
        i2 = seq_nt16_GHN_table[i2];
      }
    }
  }

  return seq_nt16_C_context_table[i1][i2];
}

double
bl_bamCrossSectionRateSimple(bam_cs_data_t *csd, uint8_t strand, uint8_t allowdel) { 
  uint32_t i, j, max, m, u, _cnt[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  uint8_t mch, uch;


  double rate = -1.0;
  if(!strand) return rate;

  bl_bamCrossSectionBSCount(csd, _cnt); 
  if(!allowdel) _cnt[0] = 0;
  mch = bl_bamGetBSBase(strand, 0);
  uch = bl_bamGetBSBase(strand, 1);
  m = _cnt[mch];
  u = _cnt[uch];

  if(m+u) {
    i = uarraymax(_cnt, 16);
    max = _cnt[i];
    _cnt[i] = 0;
    j = uarraymax(_cnt, 16);

    /* Two valid cases:
     * case 1: 
     * unique maximal occurring character in cross
     * which is either the methylated or unmethylated
     * character
     * case 2:
     * ambiguous maximal occurring character in cross
     * which must be the methylation and unmethylated
     * character
     */       
    if ((max > _cnt[j] && (i == mch || i == uch)) ||
        (max == _cnt[j] && ((i == mch && j == uch)|| (i == uch && j == mch)))){

      rate = (double) m/((double) m + (double) u);
    }

  }

  return rate;
}

char*
bl_bamCrossSectionMethylSimple(FILE *dev, bam_cs_t *cs, uint32_t k, uint32_t n, char *chrom) {
  char *tmp = NULL;
  char ctxt[3]={0,0,0};
  char strandch[]={'.','+','-'};
  char missing = '.';
  uint32_t cov,bscov,meth,unmeth,other;
  double rate;
  uint8_t ref= cs->x[k].c;

  //get BS strand
  uint8_t strand = bl_bamIsBisulfiteTarget(ref); 
  //return if non-methylation site (<=> st eq 0)
  if(!strand) {
    return NULL;
  }

  //get context
  bl_bamCrossSectionGetContext2(cs, k, strand, ctxt);
  //get coverage: note that deletions are ignored here
  cov = bl_bamCrossSectionCoverage(&cs->x[k], 0);
  bscov = bl_bamCrossSectionBSCoverage(&cs->x[k]);

  /* get base counts in bisulfite cross */
  meth = unmeth = other = 0;
  rate = -1;

  if (bscov > 0){
    assert(strand);

    /* count unconverted/methylated and converted/unmethylated bases */

    meth = bl_bamGetBSCount(&cs->x[k], strand, 0);
    unmeth = bl_bamGetBSCount(&cs->x[k], strand, 1);
    other = bscov - meth - unmeth;

    /* calculate methylation rate */
    rate = bl_bamCrossSectionRateSimple(&cs->x[k], strand, 0);
  }

  if(!strcmp(ctxt, "CG")) {
    if(k > 0 && strand == METHYL_REV_STRAND) {
      //get data from previous cs
      
    } else {
      //should never happen!
    }


    if(k+1 < n) {
      //get data from successive cs
    }
  }
  /* report output*/ 
  //  if (rate != -1){    
  //  /* output first fields */

  //    fprintf(dev, "%s\t%d\t%c\t%c\t%c\t%c\t%c", chrom, cs->beg+k+1, missing, seq_nt16_str[ref], missing, missing, missing);
  //   fprintf(dev, "\t");
  /* info field */
  //   fprintf(dev, "CS=%c;CC=%s;NS=1;MMR=%.2f;DMR=.", strandch[strand], ctxt, rate);
  //   fprintf(dev, "\t");
  /* format field */
  //   fprintf(dev, "DP:MDP:MDP3:MRDP:CM:CU:MR");
  //   fprintf(dev, "\t");
  /* data field */
  //   fprintf(dev, "%d:%d:%d,%d,%d:%d:%d:%d:%.2f", cov, bscov, meth, unmeth, other,
  //           meth + unmeth, meth, unmeth, rate);
  //   fprintf(dev, "\n");
  //  }
  //
  if (rate != -1){ 
    /* output first fields */
    bl_bsprintf(&tmp, "%s\t%d\t%c\t%c\t%c\t%c\t%c", chrom, cs->beg+k+1, missing, seq_nt16_str[ref], missing, missing, missing);
    bl_bsprintf(&tmp, "\t");
    /* info field */
    bl_bsprintf(&tmp, "CS=%c;CC=%s;NS=1;MMR=%.2f;DMR=.", strandch[strand], ctxt, rate);
    bl_bsprintf(&tmp, "\t");
    /* format field */
    bl_bsprintf(&tmp, "DP:MDP:MDP3:MRDP:CM:CU:MR");
    bl_bsprintf(&tmp, "\t");
    /* data field */
    bl_bsprintf(&tmp, "%d:%d:%d,%d,%d:%d:%d:%d:%.2f", cov, bscov, meth, unmeth, other,
        meth + unmeth, meth, unmeth, rate);
    bl_bsprintf(&tmp, "\n");
  }

  return tmp;
}

/*
 * lower -> higher bits
 *  0:  4-bit : char encoding
 *  4:  1-bit : rc 
 *  5:  6-bit : nucleotide qual (64)
 *  11: 6-bit : mapping qual (64)
 *  17: 8-bit : mm/nh
 *  25: 2-bit : paired
 *  27: 2-bit : conversion protocol
 *  29: 17-bit: query position 
 *  46: 16-bit: number mismatches/nm 
 *  62: 2-bit : reserved (indel) 
 */

void
bl_bamCrossSectionAddCoded(bam_cs_t *cs, uint64_t rpos, uint8_t nt, uint8_t rc, uint8_t nq, 
    uint8_t mq, uint8_t mm, uint8_t pp, uint8_t cp, uint32_t qp, uint16_t nm) {

  uint64_t k;
  uint64_t val =0;

  //allocation is done elsewhere
  assert(rpos >= cs->beg);
  k = rpos - cs->beg;
  assert(k < cs->n);

  //if(cs->x[k].n >= BAM_CS_DEFAULT_MAX_COVERAGE) return;

  val |=  (((uint64_t)nt) & BAM_CS_NT_MASK) << BAM_CS_NT_LBIT; 
  val |=  (((uint64_t)rc) & BAM_CS_RC_MASK) << BAM_CS_RC_LBIT; 
  val |=  (((uint64_t)nq) & BAM_CS_NQ_MASK) << BAM_CS_NQ_LBIT; 
  val |=  (((uint64_t)mq) & BAM_CS_MQ_MASK) << BAM_CS_MQ_LBIT; 
  val |=  (((uint64_t)mm) & BAM_CS_MM_MASK) << BAM_CS_MM_LBIT; 
  val |=  (((uint64_t)pp) & BAM_CS_PP_MASK) << BAM_CS_PP_LBIT; 
  val |=  (((uint64_t)cp) & BAM_CS_CP_MASK) << BAM_CS_CP_LBIT; 
  val |=  (((uint64_t)qp) & BAM_CS_QP_MASK) << BAM_CS_QP_LBIT; 
  val |=  (((uint64_t)nm) & BAM_CS_NM_MASK) << BAM_CS_NM_LBIT; 

  uint64_t m = (cs->x[k].n+1) / BAM_CS_DEFAULT_COVERAGE;
  uint64_t r = (cs->x[k].n+1) % BAM_CS_DEFAULT_COVERAGE;

  if(r == 0 || cs->x[k].n == 0) {
    cs->x[k].d = ALLOCMEMORY(NULL, cs->x[k].d, uint64_t, (m+1)*BAM_CS_DEFAULT_COVERAGE);
  }
  cs->x[k].d[cs->x[k].n] = val;
  cs->x[k].n++;

  /*
     assert(BAM_CS_GET_NM(val) == nm);
     assert(BAM_CS_GET_NT(val) == nt); 
     assert(BAM_CS_GET_NQ(val) == nq);
     assert(BAM_CS_GET_MQ(val) == mq);
     assert(BAM_CS_GET_MM(val) == mm);
     assert(BAM_CS_GET_QP(val) == qp);
     assert(BAM_CS_GET_NM(val) == nm);
     assert(BAM_CS_GET_CP(val) == cp);
     assert(BAM_CS_GET_RC(val) == rc);
     assert(BAM_CS_GET_PP(val) == pp);
     assert(BAM_CS_GET_NT(val) == nt);
     */

  return;
}


void
bl_bamCrossSectionInit(bam_cs_t *cs, uint32_t beg, uint32_t len,  uint32_t tid) {

  if(cs==NULL) return;

  cs->beg = beg;
  cs->tid = tid;
  cs->n = len;
  cs->x = calloc(len, sizeof(bam_cs_data_t));
}

void
bl_bamCrossSectionDestruct(bam_cs_t *cs) {
  uint32_t k;
  if(cs == NULL) return;

  for(k=0; k < cs->n; k++) {
    if(cs->x[k].n) FREEMEMORY(NULL, cs->x[k].d);
  }

  FREEMEMORY(NULL, cs->x);
}

void
bl_bamCrossSectionPrint(FILE* dev, bam_cs_t *cs) {
  uint32_t k,i;
  for(k=0; k < cs->n; k++){
    fprintf(dev, "%d\t%c\t", cs->beg+k, seq_nt16_str[cs->x[k].c]);
    for(i=0; i < cs->x[k].n; i++) {
      uint32_t val = cs->x[k].d[i];
      fprintf(dev, "%c", seq_nt16_str[BAM_CS_GET_NT(val)]);
    }
    fprintf(dev, "\n");
  }
}

void
bl_bamCrossSectionAddRef(bam_cs_t *cs, char *refseq, uint64_t reflen) {
  uint32_t k;  
  cs->prev = 15;
  cs->next = 15;

  if(cs->beg > 0) 
    cs->prev = seq_nt16_table[(int)refseq[cs->beg-1]];
  if(cs->beg+cs->n < reflen) 
    cs->next = seq_nt16_table[(int)refseq[cs->beg+cs->n]];

  for(k=0; k < cs->n; k++) {
    cs->x[k].c = seq_nt16_table[(int)refseq[cs->beg+k]];
  }
}

void
bl_bamCrossSectionDumpRef(bam_cs_t *cs) {
  uint32_t k;
  fprintf(stdout, ">%d:%d-%d", cs->tid, cs->beg, cs->beg+cs->n-1);
  for(k=0; k < cs->n; k++) {
    if(k % 80 == 0) fprintf(stdout, "\n");
    fprintf(stdout, "%c", seq_nt16_str[cs->x[k].c]); 
  }
  fprintf(stdout, "\n");
}

void
bl_bamCrossSectionDumpMethylSimple(FILE *dev, bam_cs_t *cs, bam_hdr_t *hdr) {
  
  uint32_t k, start=0, end=cs->n; 
  char *out;

  char *chrom = hdr->target_name[cs->tid];
 
  //adjustments necessary to allow joint CG calling 
  if(cs->beg > 0) start = 1;
  if(cs->n > 0  && cs->beg + cs->n == hdr->target_len[cs->tid]) end = cs->n-1;


  for(k=start; k < end; k++) {
    out = bl_bamCrossSectionMethylSimple(dev, cs, k, cs->n, chrom);
    fprintf(dev, "%s", out);
  }

}

void
bl_bamCrossSectionDumpMethylBuffered(circbuffer_t *cb, bam_cs_t *cs, bam_hdr_t *hdr) {
  uint32_t k, start=0, end=cs->n;
  char *out;

  char *chrom = hdr->target_name[cs->tid];
  char *msg=NULL;

  //adjustments necessary to allow joint CG calling 
  if(cs->beg > 0) start = 1;
  if(cs->n > 0  && cs->beg + cs->n == hdr->target_len[cs->tid]) end = cs->n-1;

  bl_bsprintf(&msg, "[%d,%d]\n", cs->beg, cs->beg+cs->n-1);
  bl_circBufferAddSave(cb, msg, strlen(msg));


  for(k=start; k < end; k++) {
    out = bl_bamCrossSectionMethylSimple(NULL, cs, k, cs->n, chrom);
    if(out) bl_circBufferAddSave(cb, out, strlen(out));
    FREEMEMORY(NULL, out);
  }
  FREEMEMORY(NULL, msg);

}

void 
bl_bamInitIterationBuffer(bl_bam1_buffer_t *b,  uint8_t isthreaded, pthread_mutex_t *mtx) {
  b->n = 0;
  b->i = 0;
  b->maxn = 0;
  b->next = 1;
  b->d = (bam1_t*)calloc(BAM_ITR_BUFSZ, sizeof(bam1_t));
  b->isthreaded = isthreaded;
  b->mtx = mtx;
}

void 
bl_bamDestructIterationBuffer(bl_bam1_buffer_t *b) {
  uint32_t k; 

  for(k=0 ; k < BAM_ITR_BUFSZ; k++) {
    if(b->d[k].l_data) free(b->d[k].data);
  }

  b->n = 0;
  b->mtx = NULL;
  FREEMEMORY(NULL, b->d);
}

char*
bl_bamGetQuerySeqString(bam1_t *b) {

  uint32_t len = b->core.l_qseq; //length of the read.
  uint8_t *q = bam_get_seq(b); 
  char *qseq = (char *)malloc(len+1);

  for(int i=0; i< len ; i++){
    qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
  }
  qseq[len]=0;

  return qseq;
}

char*
bl_bamGetRefSeqString(bam1_t *b, char *ref) {
  uint32_t len = b->core.l_qseq; //length of the read. 
  char *rseq = (char *)malloc(len+1);
  int32_t pos = b->core.pos;
  
  for(int i=0; i< len ; i++){
    rseq[i] = ref[pos+i]; //gets nucleotide id and converts them into IUPAC id.
  }
  rseq[len]=0;

  return rseq;
}

int
bl_sam_itr_next_buffered(samFile *in, hts_itr_t* iter, bl_bam1_buffer_t *buf, bam1_t** b) {
  uint16_t k=0;

  if(buf->i >= buf->n) {

    if(buf->isthreaded) pthread_mutex_lock(buf->mtx);
    buf->i = 0; 
    while (k < BAM_ITR_BUFSZ && (sam_itr_next(in, iter, &buf->d[k]) >= 0)) 
      k++;
    buf->n = k;
    if(buf->isthreaded) pthread_mutex_unlock(buf->mtx); 
  } 


  if(buf->i < buf->n) {
    *b=&buf->d[buf->i];
    buf->i++;
  } else { 
    return -1;
  }

  return 0;
}

void
bl_bamMethylStringWrite(bam1_t *b, char *mseq, bl_bam_methyl_worker_t* w) {

  char *out = NULL; 
  kstring_t str = {0, 0, NULL};

  if(sam_format1(w->hdr, b, &str) == -1) {
    NFO("error writing sam format.\n", NULL);
    exit(EXIT_FAILURE);
  }

  if(!w->cb) {
    //lock out device and report
    if (w->isthreaded) pthread_mutex_lock(w->devmtx);
    fprintf(w->out, "%s\tZM:Z:%s\n", str.s, mseq);
    if (w->isthreaded) pthread_mutex_unlock(w->devmtx);
  } else {
    bl_bsprintf(&out, "%s\tZM:Z:%s\n", str.s, mseq);
    if(out) {
      bl_circBufferAddSave(w->cb, out, strlen(out)); 
    } else {
      NFO("error writing methylstring to buffer.\n", NULL);
      exit(EXIT_FAILURE);
    }
    FREEMEMORY(NULL, out);
  }

  free(str.s);

}

/* bam_cigar_type returns a bit flag with:
 *   bit 1 set if the cigar operation consumes the query
 *   bit 2 set if the cigar operation consumes the reference
 *
 * For reference, the unobfuscated truth table for this function is:
 * BAM_CIGAR_TYPE  QUERY  REFERENCE
 * --------------------------------
 * BAM_CMATCH      1      1
 * BAM_CINS        1      0
 * BAM_CDEL        0      1
 * BAM_CREF_SKIP   0      1
 * BAM_CSOFT_CLIP  1      0
 * BAM_CHARD_CLIP  0      0
 * BAM_CPAD        0      0
 * BAM_CEQUAL      1      1
 * BAM_CDIFF       1      1
 * BAM_CBACK       0      0
 * --------------------------------
 */

char*
bl_bamMethylStringProcess(samFile *in, bam_hdr_t *header,  hts_idx_t *idx, 
    hts_itr_t *iter, faidx_t *fai, char* rseq, bl_bam_methyl_worker_t *w) {
  bl_bam1_buffer_t buf; 
  bam1_t *b= NULL;
  Uint cnt = 0;
  unsigned char nt;


  bl_bamInitIterationBuffer(&buf, w->isthreaded, w->bammtx);

  while (bl_sam_itr_next_buffered(in, iter, &buf, &b) >= 0) {

    if(!(b->core.flag & BAM_FUNMAP)) {
      cnt++;
      uint8_t *XB = bam_aux_get(b, "XB");
      uint64_t xb = 0;


      if(XB) {
        char* tmp = bam_aux2Z(XB);
        if(strcmp(&tmp[3], "CT")==0) xb = 1;
        if(strcmp(&tmp[3], "GA")==0) xb = 2;
      } 
 
      char *mseq = NULL;
      
      uint32_t qpos = 0;
      uint64_t rpos = b->core.pos;
      uint32_t *cigar = bam_get_cigar(b);
      uint8_t *seq = bam_get_seq(b);
      
      if(w->forest) {
        annotationitem_t myannot;
        void **results = NULL;
        uint64_t nelems = 0;
        bl_annotationitemInit(&myannot, 0);
        myannot.chromname = w->hdr->target_name[b->core.tid];
        myannot.start = rpos;
        myannot.end = rpos + bam_cigar2rlen(b->core.n_cigar, cigar)-1;
        
        bl_intervalforestSearch(w->forest, myannot.chromname, &myannot, 
            getlow_annotitem, gethigh_annotitem, &results, &nelems, NULL);
        
        if(nelems == 0) {  
          continue;
        } 
        FREEMEMORY(NULL, results);
      }
     
      int cc = 0, ck;
      
      for(int i=0; i < b->core.n_cigar; i++) {

        int op = cigar[i]&0xf;
        int32_t ln = bam_cigar_oplen(cigar[i]);

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
          //regular operations that consume reference and query here

          for(int32_t j=0; j < ln; j++) { 

            uint8_t trgt = bl_bamIsBisulfiteTarget(seq_nt16_table[(int)rseq[rpos+j]]);

            nt = seq_nt16_str[bam_seqi(seq, qpos+j)];

            if(trgt) {

              cc = bl_bamGetContextCode(rseq, rpos+j, trgt); 
              ck = seq_nt16_C_context_char[cc];

              if(xb == 1) {
                /*
                 *  read maps to fwd
                 */

                if(nt == 'C') {
                  bl_bsprintf(&mseq, "%c", toupper(ck)); 
                } else if (nt == 'T') {
                  bl_bsprintf(&mseq, "%c", ck); 
                } else {
                  bl_bsprintf(&mseq, ".");
                } 

              } else if (xb == 2){
                /* 
                 *  read maps to rev
                 */

                if(nt == 'G') {
                  bl_bsprintf(&mseq, "%c", toupper(ck));
                } else if (nt == 'A') {
                  bl_bsprintf(&mseq, "%c", ck); 
                } else {
                  bl_bsprintf(&mseq, ".");
                }
              }
            } else {
              bl_bsprintf(&mseq, ".");
            }
          }

          rpos+=ln;
          qpos+=ln;

        } else if (op == BAM_CREF_SKIP || op == BAM_CDEL) {
          //irregular skip operation (split read) 
          //and deletions that only consume reference
          rpos += ln;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
          //insertion or SC only consumes read 
          //SC can only occur at the beginning of a read
          bl_strnfill(&mseq, ln, '.');
          qpos += ln;
        }
      }


      bl_bamMethylStringWrite(b, mseq, w);
      FREEMEMORY(NULL, mseq);
    } 
  } 
 
  bl_bamDestructIterationBuffer(&buf);
  return NULL;
}

void
bl_bamCrossSectionProcess(samFile *in, bam_hdr_t *header,  hts_idx_t *idx, 
    hts_itr_t *iter, faidx_t *fai, bam_cs_t *cs, bl_bam_methyl_worker_t *w) {

  uint8_t isthreaded = w->isthreaded;
  uint8_t uniqueonly = w->uniqueonly;
  pthread_mutex_t *mtx = w->bammtx;

  bl_bam1_buffer_t buf;
  bl_bamInitIterationBuffer(&buf, isthreaded, mtx);
  bam1_t *b= NULL;
  /*
#ifndef BAM_CS_BUFFERED_SAM_ITR
b = bam_init1();
while ( sam_itr_next(in, iter, b) >= 0) {
#else 
*/
  while (bl_sam_itr_next_buffered(in, iter, &buf, &b) >= 0) {
    /*#endif */
    if(!(b->core.flag & BAM_FUNMAP)) {
      uint8_t *NH = bam_aux_get(b, "NH");
      uint64_t nh = bam_aux2i(NH);    

      if(!uniqueonly || nh == 1){ 
        //uint8_t *MD = bam_aux_get(b, "MD");
        uint8_t *NM = bam_aux_get(b, "NM");
        uint8_t *XB = bam_aux_get(b, "XB");

        //char *md = bam_aux2Z(MD);
        uint64_t nm = bam_aux2i(NM);

        uint64_t xb = 0;

        if(XB) {
          char* tmp = bam_aux2Z(XB);
          if(strcmp(&tmp[3], "CT")==0) xb = 1;
          if(strcmp(&tmp[3], "GA")==0) xb = 2;
        } 

        uint64_t rpos = b->core.pos;
        uint32_t qpos = 0;

        uint32_t *cigar = bam_get_cigar(b);
        uint8_t *seq = bam_get_seq(b);
        uint8_t *qual = bam_get_qual(b);
        uint8_t rc = (b->core.flag & BAM_FREVERSE) >> 4;
        uint8_t pp = (b->core.flag & BAM_FPAIRED) | (b->core.flag & BAM_FPROPER_PAIR);
        uint8_t mq = b->core.qual;
        uint8_t nq = qual[0];


        for(int i=0; i < b->core.n_cigar; i++) {

          int op = cigar[i]&0xf;
          int32_t ln = bam_cigar_oplen(cigar[i]);


          if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            //regular operations that consume reference and query here

            for(int32_t j=0; j < ln; j++) { 
              //boundary check
              if (rpos+j >= cs->beg+cs->n) break; 
              //proceed
              int32_t qp = qpos+j;
              int8_t nt = bam_seqi(seq, qpos+j);
              nq = qual[qpos+j];
              if(rpos+j >= cs->beg) 
                bl_bamCrossSectionAddCoded(cs, rpos+j, nt, rc, nq, mq, nh, pp, xb, qp, nm);
            }

            rpos+=ln;
            qpos+=ln;

          } else if (op == BAM_CREF_SKIP || op == BAM_CDEL) {
            //irregular skip operation (split read) 
            //and deletions that only consume reference
            if(op == BAM_CDEL) { 
              for(int32_t j=0; j < ln; j++) {
                if (rpos+j >= cs->beg+cs->n) break; 

                if(rpos+j >= cs->beg)
                  bl_bamCrossSectionAddCoded(cs, rpos+j, BAM_CS_NT_D, rc, nq, mq, nh, pp, xb, qpos, nm);
              }
            }

            rpos += ln;
          } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            //insertion or SC only consumes read 
            //SC can only occur at the beginning of a read

            qpos += ln;
          }

        }
      }
    }
  }

  //bam_destroy1(b);
  bl_bamDestructIterationBuffer(&buf);

  return;
}

bl_bam_methyl_worker_t*
bl_bamInitWorker(bl_bam_methyl_worker_t *wk, bl_bam_methyl_master_t *ms,
    FILE *out, bam_info_t *nfo, uint8_t uniqueonly, uint8_t isthreaded,
    pthread_mutex_t* devmtx, pthread_mutex_t *bammtx, circbuffer_t *cb, intervalforest_t *forest){

  wk->out = out; 
  wk->idx = nfo->idx; 
  wk->hdr = nfo->hdr;
  wk->fai = nfo->fai;
  wk->in = nfo->in;
  wk->devmtx = devmtx; 
  //wk->faimtx = faimtx; 
  wk->bammtx = bammtx; 
  wk->ms = ms;
  wk->isthreaded = isthreaded;
  wk->cb = cb;
  wk->forest = forest;
  wk->uniqueonly = uniqueonly;
  return wk;
}

uint32_t 
bl_bamMaster(bl_bam_methyl_master_t *ms, uint32_t *beg, uint32_t *len, uint8_t ovl) {

  uint32_t _beg=0, _len=0, _tid=0; 

  if(ms->isthreaded) pthread_mutex_lock(ms->mtx); 

  _tid = ms->last_tid;    

  if(_tid == -1) {
    if(ms->isthreaded) pthread_mutex_unlock(ms->mtx); 
    return -1;
  }

  _beg = ms->last_beg + ms->last_len;

  if(_beg < ms->hdr->target_len[_tid]) {
    _len = MIN(ms->hdr->target_len[_tid] - _beg, BAM_CS_DEFAULT_REG);
  } else if (_tid+1 < ms->hdr->n_targets) {
    NFO("processing chromosome '%d'\n", _tid);
    _tid++;
    _beg = 0;
    _len = MIN(ms->hdr->target_len[_tid] - _beg, BAM_CS_DEFAULT_REG);
  } else {
    MSG("done.\n");
    _tid=-1;
  }



  ms->last_tid = _tid;
  ms->last_beg = _beg;
  ms->last_len = _len;

  if(ms->isthreaded) pthread_mutex_unlock(ms->mtx);

  //give one nucleotide overlap on each side here in order to
  //call joint CG values in methylcall
  if(_tid != -1 && ovl) {
    if(_beg > 0) { _beg--; _len++;}
    if(_len < ms->hdr->target_len[_tid] - _beg) _len++;
  }
  *beg = _beg;
  *len = _len;

  return _tid;
}



void*
bl_bamCrossSectioMethylWorker(void *args) {
  bl_bam_methyl_worker_t *w = (bl_bam_methyl_worker_t*) args;

  faidx_t *fai = w->fai;
  samFile *in = w->in;
  //uint64_t exitcnt=0;
  uint64_t oldk = -1;
  FILE *out = w->out;
  hts_idx_t* idx = w->idx;
  bam_hdr_t* hdr = w->hdr;


  hts_itr_t *iter=NULL;    
  bam_cs_t cs;
  char *rseq = NULL;
  int rlen = 0; 
  uint32_t k;
  uint32_t beg;
  uint32_t len;


  while((k=bl_bamMaster(w->ms, &beg, &len, 1))!=-1) {

    //loading region assuming that k = tid
    //extended from beg+len-1 to beg+len
    iter = sam_itr_queryi(idx, k, beg, beg+len); 
    if(iter==NULL) return NULL;

    //init cross sections
    bl_bamCrossSectionInit(&cs, beg, len,  k);

    //fetch_seq is not thread safe -> lock mutex
    //if(w->isthreaded) pthread_mutex_lock(w->faimtx);
    if(k != oldk) {
      if(rseq) free(rseq);
      rseq = faidx_fetch_seq (fai, hdr->target_name[k], 0, INT_MAX, &rlen); 
    }

    bl_bamCrossSectionAddRef(&cs, rseq, rlen);   
    bl_bamCrossSectionProcess(in, hdr, idx, iter, fai, &cs, w);

 
    if(!w->cb) {
      //lock out device and report
      if (w->isthreaded) pthread_mutex_lock(w->devmtx);
      bl_bamCrossSectionDumpMethylSimple(out, &cs, hdr);
      fflush(out);
      if (w->isthreaded) pthread_mutex_unlock(w->devmtx);
    } else {
      bl_bamCrossSectionDumpMethylBuffered(w->cb, &cs, hdr);
    }


    //cleanup
    bl_bamCrossSectionDestruct(&cs);
    hts_itr_destroy(iter);  
    oldk = k;

  }

  free(rseq);
  return NULL;
}

void*
bl_bamMethylStringWorker(void *args) {
  bl_bam_methyl_worker_t *w = (bl_bam_methyl_worker_t*) args;

  faidx_t *fai = w->fai;
  samFile *in = w->in;
  //uint64_t exitcnt=0;
  uint64_t oldk = -1;
  hts_idx_t* idx = w->idx;
  bam_hdr_t* hdr = w->hdr;


  hts_itr_t *iter=NULL;    

  char *rseq = NULL;
  int rlen = 0; 
  uint32_t k;
  uint32_t beg;
  uint32_t len;


  while((k=bl_bamMaster(w->ms, &beg, &len, 0))!=-1) {
    iter = sam_itr_queryi(idx, k, beg, beg+len); 
    if(iter==NULL) return NULL;

    if(k != oldk) {
      if(rseq) free(rseq);
      rseq = faidx_fetch_seq (fai, hdr->target_name[k], 0, INT_MAX, &rlen); 
    }

    //fprintf(stderr, "tid:%d, %d-%d\n", k, beg, beg+len-1);
    bl_bamMethylStringProcess(in, hdr, idx, iter, fai, rseq, w);
    
    
    //cleanup 
    hts_itr_destroy(iter);  
    oldk = k;

  }
  free(rseq);
  return NULL;
}


bam_info_t*
bl_bamLoadInfo(bam_info_t *nfo, char *bamfn,  char *fafn) {


  nfo->in = sam_open(bamfn, "r");
  if(nfo->in==NULL) {
    NFO("error opening the bam file '%s'\n", bamfn);
    exit(EXIT_FAILURE);
  }

  if ((nfo->hdr = sam_hdr_read(nfo->in)) == 0) {
    exit(EXIT_FAILURE);
  }

  nfo->idx = sam_index_load(nfo->in, bamfn);
  if(nfo->idx==NULL) {
    NFO("error opening the index file for '%s'\n", bamfn);
    exit(EXIT_FAILURE);
  }

  nfo->fai = fai_load(fafn);
  if(nfo->fai == NULL) {

    NFO("error opening the fasta index file for '%s'\n", fafn);
    exit(EXIT_FAILURE);
  }
  return nfo;
}

void
bl_bamDestructInfo(bam_info_t *nfo) {

  fai_destroy(nfo->fai); 
  bam_hdr_destroy(nfo->hdr);
  hts_idx_destroy(nfo->idx);
  sam_close(nfo->in);

  return;
}


/*------------------------------ bam_initHeader ------------------------------
 *    
 * @brief initalize bam header structure
 * @author Steve Hoffmann 
 *   
 */
bam_hdr_t*
bl_bamGetHeader (samheader_t *head, Uint binno) {
  kstring_t str = { 0, 0, NULL };
  char *text;
  bam_hdr_t *h = NULL;
  
  text = bl_samwriteHeader(head, binno, 0, '\t', '\n'); 

  kputsn(text, strlen(text), &str);
  kputc('\n', &str);
  h = sam_hdr_parse(str.l, str.s);
  
  h->l_text = str.l; 
  h->text = str.s;

  FREEMEMORY(NULL, text);

  return h;
}
/*----------------------------- bl_samSamrec2Bamrec ----------------------------
 *    
 * @brief print one sam record
 * @author Steve Hoffmann 
 *   
 */

bam1_t*
bl_bamSamrec2Bamrec(bam_hdr_t *h, samrec_t *r) {
  Uint i;
  int ret;
  kstring_t *s;
  bam1_t* b;

  b = bam_init1();

  s = (kstring_t*)calloc(1, sizeof(kstring_t));
  ksprintf(s, "%s\t%u\t%s\t%ju\t%u\t%s\t",
      r->qname, r->flag, r->rname, r->pos, r->mapq, r->cigar);
 
  if(r->rnext) {
    ksprintf(s,"%s\t%ju\t%jd\t", r->rnext, r->pnext, r->tlen);
  } else {
    ksprintf(s,"*\t0\t0\t");
  }

  ksprintf(s, "%s\t%s\t", r->seq, r->qual);


  for(i=0; i < r->nooftags; i++) {
    ksprintf(s, "%s", r->tags[i].tag);
    if(i < r->nooftags-1) ksprintf(s,"\t");
  }

  ret = sam_parse1(s, h, b);
  assert(ret >= 0);

  if(ret <0) {
    NFO("parse error in sam2bam: '%s'\n", s->s);
    exit(-1);
  }

  FREEMEMORY(NULL, s->s);
  FREEMEMORY(NULL, s);

  return b;
}

/*----------------------------- bl_bamprintSamrec ------------------------------
 *    
 * @brief print a bam record
 * @author Steve Hoffmann 
 *   
 */

void
bl_bamPrintBamrec (htsFile *fp, samrec_t *rec, bam_hdr_t *hdr, pthread_mutex_t *mtx) {
  bam1_t *b;
  int ret;
  

  b = bl_bamSamrec2Bamrec(hdr, rec);
  pthread_mutex_lock(mtx);
  ret = sam_write1(fp, hdr, b); 
  pthread_mutex_unlock(mtx);

  if(ret < 0) {
    NFO("error writing bam file.\n", NULL);
    exit(EXIT_FAILURE);
  }

  bam_destroy1(b);

  return;
}

/*----------------------------- bam_openFile ------------------------------
 *    
 * @brief print a bam record
 * @author Steve Hoffmann 
 *   
 */
htsFile*
bl_bamOpenFile(char *fn) {
 htsFile *fp;



 if ((fp = sam_open(fn, "wb")) == 0) {
   NFO("error opening bam file '%s'\n", fn);
 }

 hts_set_cache_size(fp, 500000000);

 return fp;
}


