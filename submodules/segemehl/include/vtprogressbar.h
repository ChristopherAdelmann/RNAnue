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


#ifndef VTPROGRESSBAR_H
#define VTPROGRESSBAR_H

/*
 * =====================================================================================
 * 
 *       Filename:  vtprogressbar.h
 * 
 *    Description:  header file for a simple vt100 progressbar
 * 
 *        Version:  1.0
 *        Created:  12/07/06 00:11:54 CET
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:  Steve Hoffmann (SH), shoffmann@zbh.uni-hamburg.de
 *        Company:  Center for Bioinformatics, Hamburg
 * 
 * =====================================================================================
 */
#include "basic-types.h"

void progressBarVT (char *message, Uint complete, Uint processed, Uint size);
void cursorInvisible();
void cursorVisible();
void initProgressBarVT();

#endif
