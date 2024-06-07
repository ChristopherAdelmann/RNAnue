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


#ifndef MMCHAR_H
#define MMCHAR_H

/*
 *
 *	mmchar.h
 *  declaration for searches manber myers style
 *  on enhanced suffix arrays
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 12/22/06 19:11:16 CET  
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: mmchar.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/sufarray/mmchar.h $
 */
 
 #include "basic-types.h"
 #include "sufarray.h"

 int mmleft(Suffixarray *, char*, Uint, int, int, int);
 int mmright(Suffixarray *, char*, Uint, int, int, int);
 PairSint mmcompare(Suffixarray *, char*, Uint, int, int);
 PairSint mmsearch(Suffixarray *, char*, Uint, Uint, Uint, Uint);

#endif
