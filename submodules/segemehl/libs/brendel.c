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
 *  brendel.c
 *  Implementation of Volker Brendels EST align algorithm
 *
 *  @author Steve Hoffmann
 *  @email Deletion@bioinf.uni-leipzig.de
 *  @date 01/03/2016 12:42:41 AM CET
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include "basic-types.h"
#include "alignment.h"
#include "memory.h"
#include "mathematics.h"
#include "iupac.h"
#include "biofiles.h"
#include "locus.h"
#include "mapfrag.h"


#define BRENDELCONSTSCORES

/*---------------------------- splicedalignscore -----------------------------
 *    
 * @brief scoring function for splice alignment
 * @author Steve Hoffmann 
 *   
 */
 
double
splicedalignscore(char a, char b)
{

  double sigma = 2.0; // score of match
  double mu = -2.0; // score of a mismatch
  double delta = -4.0; // score of a deletion
	
  //if(a == 'N' && b == '-' ) return 3.0;
  if(a == 'N' || b == 'N' ) return 0.0;

  if(a== '-' || b == '-') return delta;
  if(a == b) return sigma;
  if(a != b) return mu;


  return delta;
}


/*------------------------------ splicedaligndp ------------------------------
 *    
 * @brief solving the spliced alignment problem with full dp
 * @author Steve Hoffmann 
 *   
 */
 
Alignment*
splicedaligndp (char *read, unsigned int m, char *genome, unsigned int n, gene_t **model)
{

  unsigned int u, v, maxu=0, maxv=0;
  char maxstate = 'E', state;
  char acc[] ={'-','-'};
  char explicitmatches = 1;
  Alignment *al;

  gene_t *gene;
  gene = bl_initGene("NULL",'+');

  al = ALLOCMEMORY(NULL, NULL, Alignment, 1);
  initAlignment(al, read, m, 0, genome, n, 0);
  unsigned int exonstart=0, exonlength=0;

  double p_gap = 0.03; // prob of a gap
  double p_don = 0.0005; //previsously 0.01
  double p_acc = 0.0005; //previously 0.01
  double maxprob = 0.0;
  
  double tau_start = log(0.5); // exon state prob at start
  double tau_exex = log((1.0-p_gap)*(1.0-p_don)); //transition exon -> exon
  double tau_inex = log(p_acc*(1.0-p_gap)); // intron -> exon
  double tau_exin = log((1.0-p_gap)*p_don); // exon -> intron
  double tau_inin = log(1.0-p_acc); // intron->intron
  double tau_exde = log(p_gap); // a genomic gap inside an exon
  double tau_inde = log(p_acc*p_gap); // a genomic gap at an acceptor site

  double **E = ALLOCMEMORY(NULL, NULL, double*, n+1); //now n*m 
  double **I = ALLOCMEMORY(NULL, NULL, double*, n+1); //now n*m 

  //initializations
  for(u=0; u < n+1; u++) {
    E[u] = ALLOCMEMORY(NULL, NULL, double, m+1);
    I[u] = ALLOCMEMORY(NULL, NULL, double, m+1);
    I[u][0] = 1.0; // v == 0 -> I[u][v] = 1
    E[u][0] = 1.0;
  }


  for(v=0; v < m+1; v++) {
    I[0][v] = 0.0;
    E[0][v] = 1.0;
  }

  //first recursion step with init transitions
  for(v=1; v < m+1; v++) {
    //recursions
#ifdef BRENDELCONSTSCORES
    double t1 =  MAX(E[0][v] + tau_start, I[0][v] + tau_start) - 4.0;
    double t3 = MAX(E[1][v-1] + tau_exde, I[1][v-1] + tau_inde) -4.0;
    double t2 = MAX(E[0][v-1] + tau_start, I[0][v-1] + tau_start);
    if(genome[0] != 'N' && read[v-1] != 'N') { 
      if(genome[0] == read[v-1]) t2 += 2.0; else t2 -= 2.0;
    }
#else
    double t1 = MAX(E[0][v] + tau_start, I[0][v] + tau_start) + splicedalignscore(genome[0], '-');
    double t3 = MAX(E[1][v-1] + tau_exde, I[1][v-1] + tau_inde) + splicedalignscore('-', read[v-1]);
    double t2 = MAX(E[0][v-1] + tau_start, I[0][v-1] + tau_start) + splicedalignscore(genome[0], read[v-1]);
#endif
    E[1][v] = MAX3(t1, t2, t3);
    I[1][v] = MAX(E[0][v] + (1-tau_start), I[0][v] + (1-tau_start));

    if(E[1][v] >= maxprob || I[1][v] >= maxprob) {
      maxu = 1;
      maxv = v;
      maxprob = MAX(E[1][v], I[1][v]);
    }  
  }

  //recursion 2
  for(u=2; u < n+1; u++) {
    for(v=1; v < m+1; v++) {
      
  
     double epsilon=.0;
     
     if(n > 2 && u < n-2) { 
        acc[0] = genome[u-1];
        acc[1] = genome[u];
      } else {
        acc[0] = '-';
        acc[1] = '-';
      }

      if((acc[0]=='G' && acc[1]=='T') || (acc[0]=='C' && acc[1]=='T')) epsilon = 0.001; 
   
      //recursions
#ifdef BRENDELCONSTSCORES
      double t1 = MAX(E[u-1][v] + tau_exex, I[u-1][v] + tau_inex) - 4.0;
      double t2 = MAX(E[u-1][v-1] + tau_exex, I[u-1][v-1] + tau_inex);
      double t3 = MAX(E[u][v-1] + tau_exde, I[u][v-1] + tau_inde) - 4.0;
      if(genome[u-1] != 'N' && read[v-1] != 'N'){ 
        if(genome[u-1] == read[v-1]) t2 += 2.0; else t2 -= 2.0;
      }
#else
      double t1 = MAX(E[u-1][v] + tau_exex, I[u-1][v] + tau_inex) +  splicedalignscore(genome[u-1], '-');
      double t2 = MAX(E[u-1][v-1] + tau_exex, I[u-1][v-1] + tau_inex) +  splicedalignscore(genome[u-1], read[v-1]) ;
      double t3 = MAX(E[u][v-1] + tau_exde, I[u][v-1] + tau_inde) +  splicedalignscore('-', read[v-1]);
#endif
      E[u][v] = MAX3(t1, t2, t3);
      I[u][v] = MAX(E[u-1][v] + tau_exin + epsilon, I[u-1][v] + tau_inin);

      if(E[u][v] >= maxprob || I[u][v] >= maxprob) {
        maxu = u;
        maxv = v;
        maxprob = MAX(E[u][v], I[u][v]);
        maxstate = (maxprob == E[u][v]) ? 'E' : 'I';
      }
    }
  }

  //backtrace
  u = maxu;
  v = maxv;
  state = maxstate;
  double tau_exex_ = tau_exex;
  double tau_inex_ = tau_inex;
  double tau_exin_ = tau_exin;
  double tau_inin_ = tau_inin;
  double tau_inde_ = tau_inde;

  while( u > 0 && v > 0) {

    if(u==1) {
      tau_exex_ = tau_start;
      tau_inex_ = tau_start;
      tau_exin_ = 1.0-tau_start;
      tau_inin_ = 1.0-tau_start;

    } 
  
     double epsilon=.0;
#ifdef BRENDELCONSTSCORES
     double temp = .0;
     if(genome[u-1] != 'N' && read[v-1] != 'N') { 
       temp = (genome[u-1] == read[v-1]) ? 2.0 : -2.0;
     }
#endif

     if(n > 2 && u < n-2 && u > 1) { 
        acc[0] = genome[u-1];
        acc[1] = genome[u];
      } else {
        acc[0] = '-';
        acc[1] = '-';
      }
      
      //simple score to slightly favor the canonical splice site if
      //available as opposed to non-canonical ones
      if((acc[0]=='G' && acc[1]=='T') || (acc[0]=='C' && acc[1]=='T')) epsilon = 0.001; 

      if(state == 'E') { 
     
#ifdef BRENDELCONSTSCORES
      if(E[u][v] == E[u-1][v-1] + tau_exex_ + temp) {
#else
      if(E[u][v] == E[u-1][v-1] + tau_exex_ + splicedalignscore(genome[u-1], read[v-1])) {
#endif   
        if(explicitmatches)  { 
          if (matchIUPAC(genome[u-1], read[v-1]))
            insertEop(al, Match);
          else
            insertEop(al, Mismatch);
        } else { 
          insertEop(al, Replacement);
        }
        
        state = 'E';
        u = u-1;
        v = v-1;
        exonlength++;
  
      } else

#ifdef BRENDELCONSTSCORES
    if (E[u][v] == E[u-1][v] + tau_exex_ - 4.0) { 
#else
    if (E[u][v] == E[u-1][v] + tau_exex_ + splicedalignscore(genome[u-1], '-')) {
#endif
        state = 'E';
        u = u-1;
        exonlength++;
        insertEop(al, Deletion);
      } else 
              
#ifdef BRENDELCONSTSCORES
      if(E[u][v] == E[u][v-1] + tau_exde - 4.0) {
#else
      if(E[u][v] == E[u][v-1] + tau_exde + splicedalignscore('-', read[v-1])) {
#endif
        state = 'E';
        v = v-1;
        insertEop(al, Insertion);
      } else

#ifdef BRENDELCONSTSCORES
      if (E[u][v] == I[u-1][v] + tau_inex_ - 4.0) { 
#else
      if (E[u][v] == I[u-1][v] + tau_inex_ + splicedalignscore(genome[u-1], '-')) { 
#endif
        state = 'I';
        u = u-1;
        exonstart = u;
        exonlength++;
        bl_addExon(gene, exonstart, exonstart+exonlength-1, '+', "chr");
        exonlength = 0;
        insertEop(al, Deletion);

      } else


#ifdef BRENDELCONSTSCORES
      if(E[u][v] == I[u-1][v-1] + tau_inex_ + temp) {
#else
      if(E[u][v] == I[u-1][v-1] + tau_inex_ + splicedalignscore(genome[u-1], read[v-1])) {
#endif
        state = 'I';

        if(explicitmatches)  { 
          if (matchIUPAC(genome[u-1], read[v-1]))
            insertEop(al, Match);
          else
            insertEop(al, Mismatch);
        } else { 
          insertEop(al, Replacement);
        }
        
        u = u-1;
        v = v-1;
        exonstart = u;
        exonlength++;
        bl_addExon(gene, exonstart, exonstart+exonlength-1, '+', "chr");
        exonlength = 0;

      } else


#ifdef BRENDELCONSTSCORES
      if(E[u][v] == I[u][v-1] + tau_inde_ - 4.0){
#else
      if(E[u][v] == I[u][v-1] + tau_inde_ + splicedalignscore('-', read[v-1])){
#endif
        state = 'I';
        v= v-1;
        exonstart=u;
        exonlength++;
        bl_addExon(gene, exonstart, exonstart+exonlength-1, '+', "chr");
        exonlength=0;
        insertEop(al, Insertion);
      }
    } else {

      if(I[u][v] == E[u-1][v] + tau_exin_ + epsilon) { 
        state = 'E';
        u = u-1;
        insertEop(al, Deletion);

      } else

      if(I[u][v]== I[u-1][v] + tau_inin_) {
        state = 'I';
        u= u-1;
        insertEop(al, Deletion);
      } 
    }


  }

  if(state == 'E' && exonlength > 0) {
      exonstart=u;
      bl_addExon(gene, exonstart, exonstart+exonlength-1, '+', "chr");
      exonlength=0;
  }
  
  al->voff =u;
  al->uoff =v;

  revMeops(al);

  //cleanup
  for(u=0; u < n+1; u++) {
    FREEMEMORY(NULL, E[u]);
    FREEMEMORY(NULL, I[u]);
  }

  FREEMEMORY(NULL, E);
  FREEMEMORY(NULL, I);
  
  *model = gene;
  return al;
}


/*-------------------------- bl_dpsplicedalign2map ---------------------------
 *    
 * @brief return a mapping set for alignments with gene models
 * @author Steve Hoffmann 
 *   
 */
 

mapping_t*
bl_dpsplicealign2map(Alignment *al, gene_t *model, MultiCharSeq *mseq, 
    Uint vpos, Uint vlen, char strand, char *querydesc, char *query, char *qual, 
    Uint ulen, char ismate) {

  Uint n,i;
  Alignment *new;
  mapping_t *mymapping;
  MultiCharSeqAlignment *mcsa=NULL;

  n = bl_getExonNumber(model);
  mymapping = ALLOCMEMORY(NULL, NULL, mapping_t, 1);
  bl_initMapping(mymapping, query, qual, 0, 0); 


  for(i=0; i < n; i++) {

    mcsa = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
    
    initMultiCharSeqAlignment(NULL, mcsa, mseq, 
            vpos, 0, vlen, strand, querydesc, query, qual, ulen);
    
    wrapAlignment(mcsa->al);
    FREEMEMORY(NULL, mcsa->al);

    new = getSubAlignment(al, model->exons[i].start, model->exons[i].end);
    mcsa->al = new;
    bl_addMapFrag(mymapping, mcsa, NULL, ismate, 1);

  }

 
  return mymapping;
}


/*--------------------------- bl_checkSpliceAlign ----------------------------
 *    
 * @brief check a spliced alignment
 * @author Steve Hoffmann 
 *   
 */
 

char
bl_checkSpliceAlign(mapping_t *m) {
  Uint i,n, edist, len, totaledist=0, totallen=0;
  mapfrag_t* frags = bl_getMapFrags (m, &n);

  for(i =0; i < n; i++) {
    edist = bl_getMapFragEdist(&frags[i]);
    len = bl_getMapFragGetUlen(&frags[i]);
    totaledist += edist;
    totallen += len;

    if(len <  8 ) break;
    if(len >= 8 && len <= 10 && edist > 1) break;
    if(len > 10 && len <= 15 && edist > 2) break;
    if(len > 15 && len <= 20 && edist > 2) break;
    if(len > 20 && edist > (Uint)(((float)len)*0.2)) break;
  }

  char val = 1;

  if(i < n) val =0;
  if(totallen < (Uint)(((float)frags[0].mcsa->qrylen)*0.9)) val=0;

  FREEMEMORY(NULL, frags);
  
  return val;
}




#ifdef SPLITALIGNTEST
unsigned char mute=0;
int main(int argc, char**argv) {

  initIUPAC(1, 1);
  Alignment* al;
  gene_t *model;
  
  static char genome[] = "TGTGCTCTCAACCGTACTAGCATTTGCATACGGACAGTACGGCTCCTAATACGTGTAACTATATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAAACCCAAATCGACAATGCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCCATATTATAGGGAGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTTTGATTTTTTGGCAACCCAAAATGGTGGCGGATGAACGAGATGATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCATTGCAACGTTAAATACAGCACAATATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTAATGAGTGCCTCTCGTTCTCTGTCTTATATTACCGCAAACCCAAAAAGACAATACACGACAGAGAGAGAGAGCAGCGGAGATATTTAGATTGCCTATTAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTCTATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAAACCCAAATCGACAATGCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCCATATTATAGGGAGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTTTGATTTTTTGGCAACCCAAAATGGTGGCGGATGAACGAGATGATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCATTGCAACGTTAAATACAGCACAATATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTAATGAGTGCCTCTCGTTCTCTGTCTTATATTACCGCAAACCCAAAAAGACAATACACGACAGAGAGAGAGAGCAGCGGAGATATTTAGATTGCCTATTAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTCTATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAAACCCAAATCGACAATGCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCCATATTATAGGGAGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTTTGATTTTTTGGCAACCCAAAATGGTGGCGGATGAACGAGATGATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCATTGCAACGTTAAATACAGCACAATATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTAATGAGTGCCTCTCGTTCTCTGTCTTATATTACCGCAAACCCAAAAAGACAATACACGACAGAGAGAGAGAGCAGCGGAGATATTTAGATTGCCTATTAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTCTATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAAACCCAAATCGACAATGCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCCATATTATAGGGAGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTTTGATTTTTTGGCAACCCAAAATGGTGGCGGATGAACGAGATGATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCATTGCAACGTTAAATACAGCACAATATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTAATGAGTGCCTCTCGTTCTCTGTCTTATATTACCGCAAACCCAAAAAGACAATACACGACAGAGAGAGAGAGCAGCGGAGATATTTAGATTGCCTATTAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTCTATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAAACCCAAATCGACAATGCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCCATATTATAGGGAGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTTTGATTTTTTGGCAACCCAAAATGGTGGCGGATGAACGAGATGATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCATTGCAACGTTAAATACAGCACAATATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTAATGAGTGCCTCTCGTTCTCTGTCTTATATTACCGCAAACCCAAAAAGACAATACACGACAGAGAGAGAGAGCAGCGGAGATATTTAGATTGCCTATTAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTCTATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAAACCCAAATCGACAATGCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCCATATTATAGGGAGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTTTGATTTTTTGGCAACCCAAAATGGTGGCGGATGAACGAGATGATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCATTGCAACGTTAAATACAGCACAATATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTAATGAGTGCCTCTCGTTCTCTGTCTTATATTACCGCAAACCCAAAAAGACAATACACGACAGAGAGAGAGAGCAGCGGAGATATTTAGATTGCCTATTAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTCTATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAAACCCAAATCGACAATGCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCCATATTATAGGGAGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTTTGATTTTTTGGCAACCCAAAATGGTGGCGGATGAACGAGATGATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCATTGCAACGTTAAATACAGCACAATATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTAATGAGTGCCTCTCGTTCTCTGTCTTATATTACCGCAAACCCAAAAAGACAATACACGACAGAGAGAGAGAGCAGCGGAGATATTTAGATTGCCTATTAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTCTATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAAACCCAAATCGACAATGCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCCATATTATAGGGAGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTTTGATTTTTTGGCAACCCAGTCGACGACCCGTAACCGTACTAGCATTTGCATACGGACAGTACGGCTCCTAATACGTGTAACTAG\0";
  static char read[] = "AACCGTACTAGCATGCATACGGATAGTGGGACCGCTCCTAATACGTGTAACTAGACTGCCTCTCATTCTGTCTTATTTTACCGCAAACCCACGAGATGATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCATTGCAACGTTAAATACAGCACAATATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTAATGAGTGCCTCTCGTTCTCTGTCTTATATTACCGCAAACCCAAAAACTATATAATGACTGCCTCTCATTCTGTCTTATTTTACCGCAAACCCAAAGTCGACGACCCGTAACCGTACTAGCATTTGCATACGGACAGTACGGCTCCTAATACGTGTAACTAGAAAAAAAAAA\0";
  
 
  fprintf(stderr, "starting alignment\n");
  al = splicedaligndp (read, strlen(read), genome, strlen(genome), &model);
  fprintf(stderr, "output\n");
  showAlignModel(al, stderr, model);
  
  //cleanup
  wrapAlignment(al);
  FREEMEMORY(NULL, al);
  bl_wrapGene(model);
  FREEMEMORY(NULL, model);

  return EXIT_SUCCESS;
}

#endif
