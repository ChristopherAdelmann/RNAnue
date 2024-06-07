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


#ifndef KARLIN_H
#define KARLIN_H
/*
 *
 *	karlin.h
 *  declaration for karlin.c
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 05/06/2008 08:55:40 PM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: karlin.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/karlin.h $
 */


typedef struct karlin_s {
  double lambda;
  double H;
  double K;
} karlin_t;


int karlinunitcostpp(void *space, double *lambda, double*H, double *K);
double significance (double lambda,double K,double multiplier, int score);
double evalue (double lambda,double K,double multiplier, int score);
double explength(Uint m, Uint n, double H, double K);
double effSubjectLength(Uint m, Uint n, double H, double K);
double effQueryLength(Uint m, Uint n, double H, double K);
double spacemult(Uint m, Uint n, double H, double K);
double evaluelog (double lambda, double K, double multiplier, int score);

#endif
