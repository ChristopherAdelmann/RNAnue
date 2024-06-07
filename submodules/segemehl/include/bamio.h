#ifndef BAMITER_H
#define BAMITER_H

/*
 *
 *	bamiter.h
 *  alignment representation
 *  
 *
 */

#include <stdint.h>
#include <inttypes.h>
#include <pthread.h>

#include "filebuffer.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "samio.h"
#include "intervaltree.h"
#define BAM_ITR_BUFSZ 1000


typedef struct{
    uint32_t n;
    uint32_t c; //massive padding for better alignment 
    uint64_t *d;
} bam_cs_data_t;

typedef struct{
    uint8_t prev;
    uint8_t next;
    uint32_t tid;
    uint32_t beg;
    uint32_t n;
    bam_cs_data_t *x;
} bam_cs_t;

/*
 * lower -> higher bits
 *  0:  4-bit : char encoding
 *  4:  1-bit : rc 
 *  5:  8-bit : nucleotide qual (64)
 *  13: 8-bit : mapping qual (64)
 *  21: 8-bit : mm/nh
 *  29: 2-bit : paired
 *  31: 2-bit : conversion protocol
 *  33: 16-bit: query position 
 *  49: 15-bit: number mismatches/nm 
 
 */

#define BAM_CS_NT_MASK ((uint64_t)((1<<4)-1))
#define BAM_CS_RC_MASK ((uint64_t)((1<<1)-1))
#define BAM_CS_NQ_MASK ((uint64_t)((1<<8)-1))
#define BAM_CS_MQ_MASK (((uint64_t)(1<<8)-1))
#define BAM_CS_MM_MASK (((uint64_t)(1<<8)-1))
#define BAM_CS_PP_MASK (((uint64_t)(1<<2)-1))
#define BAM_CS_CP_MASK (((uint64_t)(1<<2)-1))
#define BAM_CS_QP_MASK (((uint64_t)(1<<16)-1)) 
#define BAM_CS_NM_MASK (((uint64_t)(1<<15)-1)) 
//#define BAM_CS_RS_MASK (1<<2)-1 

#define BAM_CS_NT_LBIT 0
#define BAM_CS_RC_LBIT 4
#define BAM_CS_NQ_LBIT 5
#define BAM_CS_MQ_LBIT 13
#define BAM_CS_MM_LBIT 21
#define BAM_CS_PP_LBIT 29
#define BAM_CS_CP_LBIT 31
#define BAM_CS_QP_LBIT 33 
#define BAM_CS_NM_LBIT 49
//#define BAM_CS_RS_LBIT 62

/*
 *  encoding of common nucleotides
 */

#define BAM_CS_NT_A 0x1
#define BAM_CS_NT_C 0x2
#define BAM_CS_NT_G 0x4
#define BAM_CS_NT_T 0x8
#define BAM_CS_NT_N 0xF
#define BAM_CS_NT_D 0x0

#define BAM_CS_DEFAULT_REG 1000000
#define BAM_CS_DEFAULT_COVERAGE 30 
#define BAM_CS_DEFAULT_MAX_COVERAGE 10000
#define BAM_CS_DEFAULT_CIRCBUF_SZ 10000000

#define METHYL_FWD_STRAND 1
#define METHYL_REV_STRAND 2
#define METHYL_BTH_STRAND 0

//char cop = (char) bam_cigar_opchr(cigar[i]);
//fprintf(stderr, "%c",  "=ACMGRSVTWYHKDBN"[bam_seqi(seq, qpos+j)]);

#define BAM_CS_GET_NT(val) (((val) & BAM_CS_NT_MASK) )
#define BAM_CS_GET_RC(val) (((val) >> BAM_CS_RC_LBIT) & BAM_CS_RC_MASK)
#define BAM_CS_GET_NQ(val) (((val) >> BAM_CS_NQ_LBIT) & BAM_CS_NQ_MASK)
#define BAM_CS_GET_MQ(val) (((val) >> BAM_CS_MQ_LBIT) & BAM_CS_MQ_MASK)
#define BAM_CS_GET_MM(val) (((val) >> BAM_CS_MM_LBIT) & BAM_CS_MM_MASK)
#define BAM_CS_GET_PP(val) (((val) >> BAM_CS_PP_LBIT) & BAM_CS_PP_MASK)
#define BAM_CS_GET_CP(val) (((val) >> BAM_CS_CP_LBIT) & BAM_CS_CP_MASK)
#define BAM_CS_GET_QP(val) (((val) >> BAM_CS_QP_LBIT) & BAM_CS_QP_MASK)
#define BAM_CS_GET_NM(val) (((val) >> BAM_CS_NM_LBIT) & BAM_CS_NM_MASK)
#define BAM_CS_GET_RS(val) (((val) >> BAM_CS_RS_LBIT) & BAM_CS_RS_MASK)


typedef struct{

    bam_hdr_t *hdr;
    uint32_t last_tid;
    uint32_t last_beg;
    uint32_t last_len;

    uint8_t isthreaded;
    pthread_mutex_t *mtx;

} bl_bam_methyl_master_t;

typedef struct {
    FILE *out;

    hts_idx_t *idx;
    bam_hdr_t *hdr;
    samFile *in;   
    faidx_t *fai;
    intervalforest_t *forest;
    uint8_t uniqueonly;
    uint8_t isthreaded;
    bl_bam_methyl_master_t *ms;
    pthread_mutex_t *devmtx; 
    pthread_mutex_t *faimtx; 
    pthread_mutex_t *bammtx; 
    circbuffer_t *cb; 

} bl_bam_methyl_worker_t;



typedef struct {
    samFile *in;
    faidx_t *fai;
    hts_idx_t *idx;
    bam_hdr_t *hdr;
} bam_info_t;


typedef struct {
    uint16_t n;
    uint16_t i;
    uint16_t maxn;
    uint8_t next;
    bam1_t *d;
    uint8_t isthreaded;
    pthread_mutex_t *mtx;
} bl_bam1_buffer_t;

bam_info_t* bl_bamLoadInfo(bam_info_t *nfo, char *bamfn,  char *fafn);

bl_bam_methyl_master_t* bl_bamInitMaster(bl_bam_methyl_master_t *ms, 
        bam_hdr_t *hdr, uint8_t isthreaded, pthread_mutex_t *mtx);

bl_bam_methyl_worker_t* bl_bamInitWorker(bl_bam_methyl_worker_t *wk, 
        bl_bam_methyl_master_t *ms, FILE *out, bam_info_t *nfo, uint8_t uniqueonly, 
        uint8_t isthreaded, pthread_mutex_t* devmtx, pthread_mutex_t *bammtx, 
        circbuffer_t *cb, intervalforest_t *forest);
void* bl_bamCrossSectioMethylWorker(void *args);
void* bl_bamMethylStringWorker(void *args);

void bl_bamDestructInfo(bam_info_t *nfo);
void bl_bamPrintBamrec (htsFile *fp, samrec_t *rec, bam_hdr_t *hdr, pthread_mutex_t *mtx) ;
htsFile* bl_bamOpenFile(char *fn); 
bam_hdr_t* bl_bamGetHeader (samheader_t *head, Uint binno);
#endif
