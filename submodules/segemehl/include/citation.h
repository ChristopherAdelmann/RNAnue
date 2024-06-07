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
 *
 *	citation.h
 *  Beertime!
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/25/2007 10:50:10 PM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: citation.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/citation.h $
 */

#include <stdlib.h>
#include <time.h>

 char* cite[] = { "\"Beertime!\" (A. Torda)\0", 
                 "\"Ick fahr nur noch die janz jrossen Poette, wa!\" (Apotheker Lenz)\0",
                 "\"Nochn' schoenes Bier verhaften?\" (M. Mosisch)\0",
                 "\"Mahlzeit!\" (Ditsche, Ingo, Schildkroede)\0",
                 "\"Halt die Klappe, ich hab' Feierabend.\" (Schildkroede)\0",
                 "\"Gehen Sie vorsichtig mit dem Begriff der Unendlichkeit um!\" (Shorty)\0",
                 "\"Die Ficker!\" (Thommy)\0",
                 "\"Hab' ich gerade inner Bild gelesen!\" (Bienchen)\0",
                 "\"Tschüss, Herr Kayser!\" (Shorty)\0",
                 "\"Ich bin neu in der Hamburger Schule\" (Tocotronic)\0",
                 "\"Es wäre nicht zum aushalten, wäre er echt\" (Kettcar)\0",
                 "\"Stefan Kurtz uses suffix arrays to fix his bike.\" (A. Torda)\0",
                 "\"Ich hol' jetzt die Hilti!\" (Ein verzweifelter Bauarbeiter)\0",
                 "\"Kaeff'chen?\" (Lars)\0",
                 "\"Wir sind hier nicht in Seattle Dirk!\" (Tocotronic)\0", 
                 "\"Boooooring!\" (David S.)\0"
                 }; 


 unsigned citenumber = 16;

 char* citerand() { 
   Uint r;
   srand(time(NULL)); 
   r = ((unsigned)rand()%((int)(citenumber)));
   //fprintf(stderr, "cite %d\n",r);
   return cite[r];
 }


