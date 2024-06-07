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
 *  mappingqual.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 08/19/2015 10:50:05 AM EDT
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "mapfrag.h"
#include "debug.h"
#include "info.h"
#include "basic-types.h"
#include "memory.h"
#include "mathematics.h"
#include "alignment.h"

double recoverprobs[][4] = {
                        
                        {0.9275  , 0.7365  , 0.131618, 0.0537588}, //evalue 1000 A 80 
                        {0.937173, 0.88199 , 0.595244, 0.301946 }, 
                        {0.927928, 0.899936, 0.737285, 0.474131 }, 
                        {0.947242, 0.915428, 0.720003, 0.482145}, 
                        {0.997669, 0.984942, 0.962289, 0.932011}, 
                        {0.995098, 0.994730, 0.993867, 0.992455}}; 

double betterhits[][4] = {
                            
                        {0.0807  , 0.412474 , 0.430131  , 0.0608373},
                        {0.24324 , 0.268112 , 0.186927  , 0.033561 },
                        {0.146536, 0.28438  , 0.177582  , 0.03283  },
                        {0.185668, 0.334409 , 0.230898  , 0.039658 },
                        {0.314135, 0.078896 , 0.0250372 , 0.00138042},
                        {0.100000, 0.0314964, 0.00305881, 0.000127136}};


extern char decodeEop[];

double decodeQual33(char ascii, char* desc) {
  unsigned int myqual = (unsigned int) ascii;
  double mydoublequal = (double) myqual;

  if(mydoublequal < 33.0) {  
    fprintf(stderr, "\nQUALITY STRING ERROR in %s - %f\n",desc, mydoublequal);
  }
  assert(mydoublequal >= 33.0);

  return mydoublequal-33.0;
}


double 
logreadcoexistence(Uint besterr, Uint maxerr, Uint length) {

  Uint diff; 
  Uint lenidx;

  assert(besterr <= maxerr);

  if(maxerr > besterr+3) { 
    diff = 3;
  } else {
    diff= maxerr-besterr;
  }
    
  if(length <   25) lenidx = 0;
  if(length >=  25) lenidx = 1;
  if(length >=  28) lenidx = 2;
  if(length >=  30) lenidx = 3;
  if(length >=  50) lenidx = 4;
  if(length >= 100) lenidx = 5;


  return log(betterhits[lenidx][diff]);
}

double 
logmissedprobability(Uint maxerr, Uint length) {
  Uint erridx;
  Uint lenidx;

  //length index in matrix
  if(maxerr <=3) { 
    erridx = maxerr;
    if(length <   25) lenidx = 0;
    if(length >=  25) lenidx = 1;
    if(length >=  28) lenidx = 2;
    if(length >=  30) lenidx = 3;
    if(length >=  50) lenidx = 4;
    if(length >= 100) lenidx = 5;
  } else {

    lenidx = 0; 
    erridx = 3; 
    
    if(maxerr == 4) { 
      if(length >=   50) { lenidx = 1; erridx = 2; }
      if(length >=  100) { lenidx = 4; erridx = 2; }
    }
    
    if(maxerr == 5) {
      if(length >=   50) { lenidx = 0; erridx = 2; }
      if(length >=  100) { lenidx = 0; erridx = 0; }
    }

    if(length >= 100 && maxerr == 6) { lenidx = 3; erridx = 3; }
    if(length >= 100 && maxerr == 7) { lenidx = 1; erridx = 3; }
    if(length >= 100 && maxerr == 8) { lenidx = 0; erridx = 3; }
  }

  return log(1.0-recoverprobs[lenidx][erridx]);
}




//for the alignments not seen, i.e. missed, we calculate type2
double type2(Uint l, Uint besterr, Uint maxerr, char *qual, char *desc) {
  Uint i;
  double avg = 0, temp =.0, sum = -INFINITY;

  //calculate an average base-calling error for the read
  for(i=0; i < l; i++) {
    avg += exp(decodeQual33(qual[i],desc)/-4.342945);
  }

  fprintf(stdout, "avg: %f, log(%d): %f\n", avg, l, log(l));
  avg /= (double)l; // average error in log 10 space
  avg = log(avg);


  for(i=besterr; i <= maxerr; i++) { 
    
    //first we calculate the probability of having a true hit with more errors than the one found
    //using the bionomial coefficient for possible substitutions (approx). Theoretically maxerr = l, but
    //we need to find a better upper bound [costs time and doesnt help to calculate much more ]
    temp  = logbinomialapprox(l-besterr, i-besterr) ;
    fprintf(stdout, "\t%d \\in [%d,%d]: logbinom(%d,%d)=%f\n", i, besterr, maxerr, l-besterr, i-besterr, logbinomialapprox(l-besterr, i-besterr));
    //now we "multiply" the average error probability to the power of the additional mismatches (maxerr-besterr)
    temp += avg * ((double)(i-besterr));
    //temp += log(1.0-exp(avg)) * ((double)(i-besterr));
    

    //temp += log(1.0-exp(avg))*(l-i); 
    //temp += avg * ((double)(l-i));

    fprintf(stdout, "\t%d \\in [%d,%d]: %f*(%d-%d)=%f, \n", i, besterr, maxerr, avg, i, besterr, 
        avg*((double)(i-besterr)));

    fprintf(stdout, "\t%d \\in [%d,%d]: %f*(%d-%d)=%f, \n", i, besterr, maxerr, log(1.0-exp(avg)), l, i, 
        log(1.0-exp(avg))*((double)(l-i)));

    //now we "multiply" the probability that such alginments (with maxerr and besterr coexist)
    temp += logreadcoexistence(besterr, i, l);
    fprintf(stdout, "\t%d \\in [%d,%d]: factor=%f\n", i, besterr, maxerr, logreadcoexistence(besterr, i, l));
    //now we "multiply" the probabilty that an alignment is missed (most important because this is the 
    //probability that gives rise to this type of error)
    temp += logmissedprobability(i, l);
    fprintf(stdout, "\t%d \\in [%d,%d]: factor=%f\n", i, besterr, maxerr, logmissedprobability(i, l));

    sum = logadd(sum, temp);
    fprintf(stdout, "\t%d \\in [%d,%d]: temp=%f, sum=%f\n", i, besterr, maxerr, temp, sum);
  }

  return sum;
}


//coverter
double type2mappingset(mappingset_t *set, char *desc, Uint n, karlin_t* stats) {
 
  Uint j, k, len, edist;
  char *qual;
  Uint noofaligns = set->n;
  Uint nooffrags;


  for(j=0; j < noofaligns; j++) {
     
    nooffrags = set->elem[j].n;

    //iterate fragments of each independent alginment
    for(k=0; k < nooffrags; k++) {
    
      qual  = bl_getMapFragQual(&set->elem[j].f[k]);
      edist = bl_getMapFragEdist(&set->elem[j].f[k]);
      len = strlen(qual);
      bl_setMapFragMapQual(&set->elem[j].f[k], type2(len, edist, 8, qual, desc));
      
    }
  }

  return .0;
}


double longestmatchqual(mappingset_t *set, Uint n, karlin_t* stats) {

  Uint j, k, len, lmat;
  //Uint edist;

  Uint noofaligns = set->n;
  Uint nooffrags;
  double lmateval;
  
  for(j=0; j < noofaligns; j++) {
     
    nooffrags = set->elem[j].n;
    //iterate fragments of each independent alginment
    for(k=0; k < nooffrags; k++) {


      len = getUalignlen(set->elem[j].f->mcsa->al);
      //edist = bl_getMapFragEdist(&set->elem[j].f[k]);
      lmat = bl_getMapFragLongestMatch(&set->elem[j].f[k]);

      
      double space = spacemult(len, n, stats->H, stats->K);
      double scale = evaluelog(stats->lambda, stats->K, space, 25); 

      lmateval = evaluelog(stats->lambda, stats->K, space, lmat);
      lmateval /= -1.0*scale; 
       
//      fprintf(stdout, "lmat: %d, E-val:%f/%f=%f, log10:%f\n", lmat, lmateval*scale, scale, lmateval, lmateval/log(10));           
      lmateval = MIN(-0.23026, lmateval);

      bl_setMapFragMapQual(&set->elem[j].f[k], lmateval);
    }
  }

  return .0;
}


//coverter
double type3mappingset(mappingset_t *set, char *desc) {

  Uint j, i, k, u=0, v=0, x=0, uallen, maxj=0;
  Alignment *aln;
  char *qual;
  double **p, *q, sigma, min;
  
  sigma = -INFINITY;
  Uint noofaligns = set->n;
  Uint nooffrags;
  //          unsigned int temp = 0;

  p = ALLOCMEMORY(NULL, NULL, double*, noofaligns);
  q = ALLOCMEMORY(NULL, NULL, double, noofaligns);

  //iterate elements of set (independent alignments)
  for(j=0; j < noofaligns; j++) {

    nooffrags = set->elem[j].n;
    p[j] = ALLOCMEMORY(NULL, NULL, double, nooffrags);
    q[j] = 0;
    
    //iterate fragments of each independent alginment
    for(k=0; k < nooffrags; k++) {

      aln = set->elem[j].f[k].mcsa->al;
      qual = bl_getMapFragQual(&set->elem[j].f[k]);
      uallen = getUalignlen(aln);

      p[j][k] = 0.0;
   //   fprintf(stdout, "alignment [%d,%d]: parsing eopstring of length %d\n", j, k, noofeops);
   //   fprintf(stdout, "%s\n%s\n", eopstring, qual);
      u = 0;
      v = 0;
      for(x=0; x < aln->numofmeops; x++) { 

        switch(decodeEop[aln->meops[x].eop]) {
          case 'R':
            for(i=0; i < aln->meops[x].steps; i++) {
            if(isMatch(aln, u, v)) {  
              p[j][k] += log(1.0 - exp(decodeQual33(qual[u],desc) / -4.342945)); 
              
     //         fprintf(stdout, "Match - qual:%d -> p=%f -> sum=%f\n", qual[u], 
     //           log(1.0 - exp(decodeQual33(qual[u],desc) / -4.342945)), p[j][k]) ;

              // = 1 - log(e) , where e is the error prob
            } else {

              p[j][k] += (decodeQual33(qual[u],desc) / -4.342945) - 1.097; 
                            
       //       fprintf(stdout, "MisMatch - qual:%d -> p=%f -> sum=%f\n", qual[u],
       //           ( decodeQual33(qual[u],desc) / -4.342945), p[j][k]) ;

              // log(e) - log(3) = log(e/3) (the e/3 is doubtful)
            }
            u++;
            v++;
            }
            break;
          case '=':
            for(i=0; i < aln->meops[x].steps; i++) {
              p[j][k] += log(1.0 - exp(decodeQual33(qual[u],desc) / -4.342945)); 
              u++;
              v++;
            }
            break;
          case 'X':
            for(i=0; i < aln->meops[x].steps; i++) {
              p[j][k] += (decodeQual33(qual[u],desc) / -4.342945) - 1.097; 
              u++;
              v++;
            }
            break;
          case 'D':
            for(i=0; i < aln->meops[x].steps; i++) {
            if(u > 0) {
              min = ((decodeQual33(qual[u-1],desc) / -4.342945) - 1.3863);
       //       temp=u-1;
            } else {
              min = ((decodeQual33(qual[0],desc) / -4.342945) - 1.3863); 
       //       temp = 0;
            }

            if(u < uallen && min > ((decodeQual33(qual[u],desc) / -4.342945) - 1.3863)) {
              min = (decodeQual33(qual[u],desc) / -4.342945) - 1.3863;
       //       temp = u;
            }
          
        //    fprintf(stdout, "Deletion - qual:%d -> p=%f -> sum=%f\n", qual[temp],
        //          ( decodeQual33(qual[temp],desc) / -4.342945) -1.3863, p[j][k]) ;

            p[j][k] += min;
            v++;
            }
            break;
          case 'I':
            for(i=0; i < aln->meops[x].steps; i++) {
            p[j][k] += (decodeQual33(qual[u],desc) / -4.342945)-1.3863; 
      
          //  fprintf(stdout, "Insertion - qual:%d -> p=%f -> sum=%f\n", qual[u],
          //        ( decodeQual33(qual[u],desc) / -4.342945) - 1.3863, p[j][k]) ;


            // log(e) see comment on e/3
            u++;
            }
            break;
         default:
            break;

        }
      }
    
      q[j] += p[j][k];
      FREEMEMORY(NULL, qual);
    }


   // get the sigma = sum of p(z | x_v^l) over 
   // all alignments v \in V in log space
   sigma = log(exp(sigma) + exp(q[j]));
   //fprintf(stdout, "sigma=%f\n", sigma);
 }

  maxj = 0;
  for(j=0; j < noofaligns; j++) {

    //iterate fragments of each independent alginment
    //sum up alignment score
    //fprintf(stdout, "mapping qual -------- qual:%f, sigma:%f, post:%f\n", q[j], sigma, q[j]-sigma);
    q[j] -= sigma;
    set->elem[j].mapqual = log(1.0-exp(q[j]));
    set->elem[j].sigma = sigma;
    if (q[j] > q[maxj]) {
      maxj = j;
      set->elem[j].maxqual = q[maxj];
    }
    FREEMEMORY(NULL, p[j]);
  }

  // fprintf(stdout, "maximum mapping qual ------ %f\n", q[maxj]);

  //error prob only of the best hit
  //double ebest = 1-q[maxj];
  //double qbest = -4.342945 * ebest;

  FREEMEMORY(NULL, q);
  FREEMEMORY(NULL, p);


  return .0; //qbest;
}


//for the alignments seen we calucate the type3:wq
double 
type3(Alignment** aln, Uint noofaligns, char *qual) {

  Uint j, i, u=0, v=0, noofeops, uallen, maxj=0;
  double *p, sigma, min;
  char *eopstring;
  sigma = -INFINITY;

  p =ALLOCMEMORY(NULL, NULL, double, noofaligns);
  
  // calculate p(z | x_u^l ) where z is the read, x_u^l 
  // is the reference sequence at pos U=u of length l for 
  // each of the alignments
  // by using the sum of log10 logarithms of qualities over 
  // all l nucleotides (for matches) and the 

  for(j=0; j < noofaligns; j++) { 
    eopstring = getEopString(aln[j]);
    noofeops = strlen(eopstring);
    uallen = getUalignlen(aln[j]);

    for(i=0; i < noofeops; i++) {
      switch(eopstring[i]) {
        case Replacement:
          if(isMatch(aln[i], u, v)) {  
            p[j] += log(1 - exp(qual[u] / -4.342945)); 
            // = log(1) - log(e) , where e is the error prob
          } else {
            p[j] += (qual[u] / -4.342945) - 1.097; 
            // log(e) - log(3) = log(e/3) (the e/3 is doubtful)
          }
          u++;
          v++;
          break;
        
        case Match:
          p[j] += log(1 - exp(qual[u] / -4.342945)); 
          u++;
          v++;
          break;
        case Mismatch:
          p[j] += (qual[u] / -4.342945) - 1.097; 
          u++;
          v++;
          break; 
        case Deletion:
          if(u > 0) {
            min = (qual[u-1] / -4.342945);
          } else {
            min = (qual[0] / -4.342945); 
          }

          if(u < uallen && min > (qual[u] / -4.342945) ) {
            min = (qual[u] / -4.342945);
          }

          p[j] += min;
          v++;
          break;
        case Insertion:
          p[j] += (qual[u] / -4.342945); 
          // log(e) see comment on e/3
          u++;
          break;
        case Skipped:
          break;
          
      }
    }
    
    // get the sigma = sum of p(z | x_v^l) over 
    // all alignments v \in V in log space
    sigma = logadd(sigma, p[j]);
  }

  // get posterior probability p(u | x, z) for 
  // each of the observed alignments by 
  // dividing p(z | x_u^l) / sigma
 
  maxj = 0;
  for(j=0; j < noofaligns; j++) { 
    p[j] -= sigma;
    if (p[j] > p[maxj]) {
      maxj = j;
    }
  }

  //error prob only of the best hit
  double ebest = 1-p[maxj];
  double qbest = -4.342945 * ebest;

  FREEMEMORY(NULL, p);

  return qbest;
}


double
mappingquality(Alignment **aln, Uint noofaligns, char *qual, Uint length, Uint besterr, char *desc) {

    double t3, t2;
    Uint maxerr=3;

    t3 = type3(aln, noofaligns, qual);
    t2 = type2(length, besterr, maxerr, qual, desc);

    return MIN(t2,t3);
}

//coverter
double
playground(mappingset_t *set, Uint n, karlin_t* stats) {
 
  Uint j, k, len, edist, lmat;
  char *qual;
  Uint noofaligns = set->n;
  Uint nooffrags;


  
  fprintf(stdout,"error type 2 --------------- \n");
  for(j=0; j < noofaligns; j++) {
     
    nooffrags = set->elem[j].n;

    //iterate fragments of each independent alginment
    for(k=0; k < nooffrags; k++) {


    
      qual  = bl_getMapFragQual(&set->elem[j].f[k]);
      edist = bl_getMapFragEdist(&set->elem[j].f[k]);
      lmat = bl_getMapFragLongestMatch(&set->elem[j].f[k]);
      len = strlen(qual);
      fprintf(stdout,"[%d,%d]-%d:\n",j, k, edist);
  
      double space = spacemult(len, n, stats->H, stats->K);
      double pval = significance(stats->lambda, stats->K, 1, len-(3*edist));
      double eval = evaluelog(stats->lambda, stats->K, space, len-(3*edist));
      double lmateval = evaluelog(stats->lambda, stats->K, space, lmat);

      space = spacemult(20, n, stats->H, stats->K);
      double base20 = evaluelog(stats->lambda, stats->K, space, 20); 
      
      space = spacemult(25, n, stats->H, stats->K);
      double base25 = evaluelog(stats->lambda, stats->K, space, 25); 
      
      fprintf(stdout, "len:%d, edist:%d, logged evalue:%f, pvalue:%f, lmateval:%f\n", len, edist, eval, pval, lmateval);
      fprintf(stdout, "base20: %f, base25: %f\n", base20, base25);
      fprintf(stdout, "log10(eval-base20): %f, log10(eval-base25): %f\n", (eval-base20)/log(10), (eval-base25)/log(10));
      fprintf(stdout, "log10(eval/base20): %f, log10(eval/base25): %f\n", (eval/base20)/log(10), (eval/base25)/log(10));
      fprintf(stdout, "log10(eval/log(len)): %f\n", (eval/log(len))/log(10));
      fprintf(stdout, "--------\n");

      
      bl_setMapFragMapQual(&set->elem[j].f[k], pval);

   //   bl_setMapFragMapQual(&set->elem[j].f[k], type2(len, edist, 8, qual));
    }
  }

  return .0;
}
