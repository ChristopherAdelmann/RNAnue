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
 *	matealign.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 29.07.2013 18:07:27 CEST  
 *
 */

mappingset_t* bl_greedypairMappingSets(mappingset_t *l, mappingset_t *r);

unsigned int bl_myersMCSA(MultiCharSeqAlignment *mcsa, MultiCharSeq *mseq, 
    unsigned int *enctab, bitvector *peq, bitvector *D, unsigned int maxE);

mappingset_t *
bl_matealign(mappingset_t *l, MultiCharSeq *mseq, char **seqs, char **quals, char *qname, unsigned int m,
    unsigned int *enctab, bitvector *D, unsigned int maxL, char ismate, Suffixarray *s, segemehl_t *nfo);

mappingset_t *
bl_pairMateMapping (mappingset_t *l, mappingset_t *r, unsigned int insertsize);

unsigned int
bl_localMCSA(MultiCharSeqAlignment *mcsa, MultiCharSeq *mseq, int* scores, int indel, int minscore); 
