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


 #ifndef FALPHABET_H
 #define FALPHABET_H

/*
 * alphabet.h
 * declarations for a flexible alphabet
 *
 *  SVN
 *  Revision of last commit: $Rev: 40 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-09-03 10:42:48 +0200 (Wed, 03 Sep 2008) $
 *
 *  Id: $Id: falphabet.h 40 2008-09-03 08:42:48Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/trunk/libs/sufarray/falphabet.h $
 */

 #include "basic-types.h"

 typedef struct {

	Uint *characters,
		 *mapdomain;
	
	Uint domainsize,
		 mapsize,

		 mappedwildcards,
		 undefsymbol,
		 *symbolmap;

 } FAlphabet;


 
#endif
