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
 *  vtprogressbar.c
 *  implementation for a very simple
 *  progress bar
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/05/06 02:11:30 CET
 *  
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: vtprogressbar.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/vtprogressbar.c $
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "basic-types.h"
#include "vtprogressbar.h"


void
cursorInvisible() {
  fprintf(stderr, "%c%c%c%d%c", 27, '[', '?', 25, 'l');
}

void
cursorVisible() {
  fprintf(stderr, "%c%c%c%d%c", 27, '[', '?', 25, 'h');
}



/*---------------------------- initProgressBarVT -----------------------------
 *    
 * initializes a progress bar for VT
 * 
 */
 
void
initProgressBarVT ()
{	
  	fprintf(stderr, "%c%c%c", 27, '[', 's');
    fprintf(stderr, "%c%c%c", 27, '[', 'K');
    return ;
}

/*------------------------------ progressBarVT -------------------------------
 *    
 * a simple progress bar for VT terminals
 * 
 */

void
progressBarVT (char *message, Uint complete, Uint processed, Uint size)
{
    Uint i, bar, percent;
    double p;
	char cur;
    
    if (complete == 0) complete = 1;
    p = ((double) processed)/complete;
    bar = (Uint) (size * p);
    percent = (Uint) (100 * p);
    fprintf(stderr, "[");
    for(i=0; i < size; i++) {
      if(i<=bar) fprintf(stderr, "=");
      else fprintf(stderr," ");
    }
    i = processed % 30;
    if (i<=10) cur = '/'; 
    else if (i<=20)cur='\\';
    else cur='-';
    fprintf(stderr,"]   %d%c(%d)  %s  %c\n", percent, '%',
            processed, message, cur);
    fprintf(stderr,"%c%c%c", 27, '[', 'A');
    return;
}


