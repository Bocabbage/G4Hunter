#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Update date: 2020/2/11(unfinished)
########################################################################
"""
    <G4Hunter - a program to search quadruplex-forming regions in DNA.>
    Copyright (C) <2012>  <Bedrat amina  supervised by Dr.Jean Louis Mergny>
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
########################################################################
import argparse
from Bio import SeqIO


class Soft:
    def __init__(self):
        pass

    # Use Bio.SeqIO API to read the fasta file.
    def ReadFile(self,Filein):
        # ListSeq,LHeader=[],[]
        seqs_with_headers = []
        for record in SeqIO.parse(Filein, "fasta") :
            # LHeader.append(record.id)
            # ListSeq.append(record.seq)
            seqs_with_headers.append([record.id, record.seq])
        # return LHeader,ListSeq
        # Modify Note: return List of [Header, Seq]
        return seqs_with_headers
        
    def __BaseScore(self, seq):
        item, bases_score = 0, []
        seq_len = len(seq)
        while (item < seq_len):
            if(seq[item] == 'G' or seq[item] == 'g'):
                # G-base(s) scoring rules:
                # if single G, gives score 1;
                # if GG, gives each G a score 2;
                # the max score for a base is 4
                new_item = item + 1
                new_score = 1
                for increment in range(4):
                    if (
                         new_item < seq_len and 
                         (seq[new_item] == 'G' or seq[new_item] == 'g')
                       ):
                        new_item += 1
                        new_score += 1
                    else:
                        break
                bases_score += [ new_score for x in range(new_item - item)]
                item = new_item

                # if more than 4 continuous G, all of them
                # are given score 4.
                while(item < seq_len and (seq[item] == 'G' or seq[item == 'g'])):
                    bases_score.append(4)
                    item += 1

            elif(seq[item] == 'C' or seq[item] == 'c'):
                # C-base(s) scoring rules:
                # if single C, gives score -1;
                # if CC, gives each G a score -2;
                # the max score for a base is -4
                new_item = item + 1
                new_score = -1
                for increment in range(4):
                    if (
                         new_item < seq_len and 
                         (seq[new_item] == 'C' or seq[new_item] == 'c')
                       ):
                        new_item += 1
                        new_score -= 1
                    else:
                        break
                bases_score += [ new_score for x in range(new_item - item)]
                item = new_item

                # if more than 4 continuous C, all of them
                # are given score -4.
                while(item < seq_len and (seq[item] == 'C' or seq[item == 'c'])):
                    bases_score.append(-4)
                    item += 1
                
            else:
                # A, T, U and N
                bases_score.append(0)
                item += 1

        return seq, bases_score















        


if __name__ == '__main__':
    # Args' Parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--ifile', help='input seq file.')
    parser.add_argument('-o', '--ofile', help='ouput seq file.')
    parser.add_argument('-w', '--wind', default=25, type=int,  
                        help='width of the slide-window.')
    parser.add_argument('-s', '--score', type=float, 
                        help='threshold of G4H score.')
    args = parser.parse_args()

