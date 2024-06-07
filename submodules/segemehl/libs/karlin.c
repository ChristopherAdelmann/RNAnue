
/*
 *   segemehl - a read aligner
 *   Copyright (C) 2008-2017  Steve Hoffmann and Stephan Kurtz
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "debug.h"
#include "info.h"
#include "basic-types.h"
#include "memory.h"
#include "mathematics.h"

#define MAXIT 150    /* Maximum number of iterations used in calculating K */

static int karlinpp(void *space, int low, int high, double *pr, 
                     double *lambda, double *K)
{
  int i, j, range, lo, hi, first, last;
  double upval, Sumval, av, sum, *p, *P, *ptrP, *ptr1, *ptr2, *ptr1e, newval,
         Ktmp;

  /* Check that scores and their associated probabilities are valid     */

  if (low >= 0)
  {
    NFO("Lowest score %ld must be negative", low);
    return  -1;
  }
  for (i = range = high - low; i > -low && !pr[i]; --i)
    /* Nothing */ ;
  if (i <= -low)
  {
    MSG("A positive score must be possible");
    return  -2;
  }
  for (sum = 0.0, i = 0; i <= range; sum += pr[i++])
  {
    if (pr[i] < 0.0)
    {
      DBG("Negative probability %.2f not allowed",pr[i]);
      return  -3;
    }
  }
  if (sum < 0.99995 || sum > 1.00005)
  {
    DBGL(3,"Probabilities sum to %.4f. Normalizing.\n", sum);
  }
  
  p = ALLOCMEMORY(space, NULL, double, (Uint)(range+100));
  if(p == NULL)
  {
    return  -4;
  }
  for (Sumval = (double) low, i = 0; i <= range; ++i)
  {
    Sumval += i * (p[i] = pr[i] / sum);
  }
  if(Sumval >= 0.0)
  {
    NFO("Invalid (non-negative) expected score:  %.3f", Sumval);
    return  -5;
  }

  /* Calculate the parameter lambda */

  upval = 0.5;
  do
  {
    upval *= 2;
    ptr1 = p;
    for (sum = 0.0, i = low; i <= high; ++i)
    {
      sum += *ptr1++ * exp (upval * i);
    }
  } while (sum < 1.0);
  for (*lambda = 0.0, j = 0; j <  40; ++j)
  {
    newval = (*lambda + upval) / 2.0;
    ptr1 = p;
    for (sum = 0.0, i = low; i <= high; ++i)
    {
      sum += *ptr1++ * exp (newval * i);
    }
    if (sum > 1.0)
    {
      upval = newval;
    } else
    {
      *lambda = newval;
    }
  }

  /* Calculate the pamameter K */

  ptr1 = p;
  for (av = 0.0, i = low; i <= high; ++i)
  {
    av += *ptr1++ * i * exp (*lambda * i);
  }
  if (low ==  -1 || high ==  1)
  {
    *K = (high ==  1) ? av : Sumval * Sumval / av;
    *K *= 1.0 - exp (-*lambda);
    free (p);
    return 0;  /* Parameters calculated successfully */
  }
  Sumval = 0.0;
  lo = 0;
  hi = 0;
  P = ALLOCMEMORY(space, NULL,double,(Uint) (MAXIT * range + 100));
  if(P == NULL)
  {
    return  -6;
  }
  for (*P = 1.0, sum = 1.0, j =  1; 
       j <=  MAXIT && sum > 0.00001; Sumval += sum /= j++)
  {
    first = last = range;
    for (ptrP = P + (hi += high) - (lo += low); ptrP >= P; *ptrP-- = sum)
    {
      ptr1 = ptrP - first;
      ptr1e = ptrP - last;
      ptr2 = p + first;
      for (sum = 0.0; ptr1 >= ptr1e;)
      {
        sum += *ptr1-- * *ptr2++;
      }
      if (first)
      {
        --first;
      }
      if (( (ptrP - P)) <= range)
      {
        --last;
      }
    }
    for (sum = 0.0, i = lo; i; ++i)
    {
      sum += *++ptrP * exp (*lambda * i);
    }
    for (; i <= hi; ++i)
    {
      sum += *++ptrP;
    }
  }
  if (j >  MAXIT)
  {
    MSG("Value for K may be too large due to insufficient iterations");
    return  -7;
  }
  for (i = low; !p[i - low]; ++i)
    /* Nothing */ ;
  for (j = -i; i < high && j >  1;)
  {
    if (p[++i - low] != 0.0)
    {
      j = gcd (j, i);
    }
  }
  Ktmp = (double) (j * exp (-2 * Sumval));
  *K = Ktmp / (av * (1.0 - exp (-*lambda * j)));

  FREEMEMORY(space, P);
  FREEMEMORY(space, p);
  return 0;  /* Parameters calculated successfully */
}

int karlinunitcostpp(void *space, double *lambda, double *H, double *K)
{
  int ret;
  int match = 1;
  int mismatch = -1;
  double targetid;
  
  double pr[] = {0.75, 0.0, 0.25};
  ret = karlinpp(space, mismatch, match, &pr[0], lambda, K);
  if (ret  != 0) return ret;

  targetid = 0.25 * exp(*lambda*match);
  *H = (*lambda * match * targetid) + (*lambda * mismatch * (1-targetid)); 

  return ret;
}


double significance (double lambda,double K,double multiplier, int score)
{
  double y;

  y = -lambda * score;
  y = K * multiplier * exp (y);
  return exp (-y);
}


double evalue (double lambda, double K, double multiplier, int score)
{
  double y;

  y = -lambda * score;
  y = K * multiplier * exp (y);
  return y;
}

double evaluelog (double lambda, double K, double multiplier, int score)
{
  double y;

  y = -lambda * score;
  y = log(K) + log(multiplier) + y;
  return y;
}


double bitscoreevalue (double lambda,double K,double multiplier, int score)
{
  double y;

  y = -1 * score;
  y = multiplier * pow(2,(y));
  return y;
}

double bitscore(int score, double lambda, double K) {
    return ((lambda*score)-log(K))/log(2);
}

double explength(Uint m, Uint n, double H, double K) {
    return log(m*n*K)/H;
}

double effSubjectLength(Uint m, Uint n, double H, double K) {
  double effl = (double)n - (1*explength(m, n, H, K));
  return effl < 1/K ? 1/K : effl;

}
double effQueryLength(Uint m, Uint n, double H, double K) {
  double effl = (double)m - explength(m, n, H, K);
  return effl < 1/K ? 1/K : effl;
}

double spacemult(Uint m, Uint n, double H, double K) {
  return (double) effSubjectLength(m,n,H,K) * effQueryLength(m,n,H,K);
}

