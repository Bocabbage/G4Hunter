#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Update date: 2020/2/13
# Note: No visualization function which is existed in py2 version.
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
import os
import numpy as np
from Bio import SeqIO


class Soft:
    def __init__(self):
        pass

    def ReadFile(self, Filein):
        '''
            Use Bio.SeqIO API to read the fasta file.
        '''
        # ListSeq,LHeader=[],[]
        seqs_with_headers = []
        for record in SeqIO.parse(Filein, "fasta") :
            # LHeader.append(record.id)
            # ListSeq.append(record.seq)
            seqs_with_headers.append([record.id, record.seq])
        # return LHeader,ListSeq
        # Modify Note: return List of [Header, Seq]
        return seqs_with_headers

    def GFinder(self, infile, k):
        '''
            For every entries in the fasta file, read and call
            __BaseScore and __CalScore to get the mean-scores list
            of every windows in the seq, then packaging them
            as output.
        '''
        seqs_with_headers = self.ReadFile(infile)
        headers_seqs_scores = []
        for seq_with_header in seqs_with_headers:
            bases_score = self.__BaseScore(seq_with_header[1])
            mean_scores_list = self.__CalScore(bases_score, k)
            headers_seqs_scores.append([seq_with_header[0], seq_with_header[1], 
                                             mean_scores_list])
        return headers_seqs_scores


    def GetG4(self, ofile, entry, threshold, k):
        '''
            Compare mean-scores of windows in the input-seq
            with the threshold and return the indexes of 
            G4-windows-long piece as well as writing them 
            to the ofile using __Write api.

            input:
                    ofile     -- the output file discriptor
                    entry     -- [header, seq, seq_mean_scores](output of GFinder)
                    threshold -- cut-off G4Score
                    k         -- window width
        '''
        LG4_idxes = []
        header, seq, seq_mean_scores = entry
        SEQ=">"+header+"\n Start \t End \t Sequence\t Length \t Score\n"
        ofile.write(SEQ)
        for i in range(len(seq_mean_scores)):
            if(seq_mean_scores[i] >= threshold or seq_mean_scores[i] <= -threshold):
                seq_in_win = seq[i:i+k]
                LG4_idxes.append(i)
                self.__Write(ofile, i, k, 0, 0, seq_in_win, k, seq_mean_scores[i])
                ofile.write("\n")
        return LG4_idxes


    def WriteSeq(self, ofile, entry, g4_idx_list, F):
        '''
            Merge the overlapped G4 pieces.
        '''
        header, line, liste = entry
        i,k,I=0,0,0
        a=b=g4_idx_list[i]
        MSCORE=[]
        SEQ=">"+header+"\nStart\tEnd\tSequence\tLength\tScore\tNBR\n"
        ofile.write(SEQ)

        if (len(g4_idx_list) > 1):
            c=g4_idx_list[i+1]
            while (i< len(g4_idx_list) - 2):
                if(c==b+1):
                    k=k+1
                    i=i+1
                else:
                    I=I+1
                    seq=line[a:a+F+k]
                    liste2=self.__BaseScore(seq)
                    self.__Write(ofile, a, k , F, 0, seq ,len(seq) , round(np.mean(liste2),2))
                    MSCORE.append(abs(round(np.mean(liste2),2)))
                    ofile.write("\n")  
                    k=0
                    i=i+1
                    a=g4_idx_list[i]
                b=g4_idx_list[i] 
                c=g4_idx_list[i+1] 
            I=I+1
            seq=line[a:a+F+k+1]
            liste2=self.__BaseScore(seq)
            self.__Write(ofile, a, k , F, 1, seq ,len(seq) , round(np.mean(liste2),2))
            MSCORE.append(abs(round(np.mean(liste2),2)))
            ofile.write("\t")          
            ofile.write(str(I))
            ofile.write("\n")
        else:
            I=I+1
            seq=line[a:a+F]
            self.__Write(ofile, a, 0, F, 0, seq, len(seq) , liste[a])
            MSCORE.append(abs(liste[a]))
            ofile.write("\t")
            ofile.write(str(I))
            ofile.write("\n")   

        return MSCORE 
        

    def __BaseScore(self, seq):
        '''
            take a sequence as input and return the score of every base
        '''
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
                for increment in range(1, 4):
                    if (
                         new_item < seq_len and 
                         (seq[new_item] == 'G' or seq[new_item] == 'g')
                       ):
                        new_item += 1
                        new_score += 1
                    else:
                        break
                bases_score += [ new_score for x in range(new_item - item) ]
                item = new_item

                # if more than 4 continuous G, all of them
                # are given score 4.
                while(item < seq_len and (seq[item] == 'G' or seq[item] == 'g')):
                    bases_score.append(4)
                    item += 1

            elif(seq[item] == 'C' or seq[item] == 'c'):
                # C-base(s) scoring rules:
                # if single C, gives score -1;
                # if CC, gives each G a score -2;
                # the min score for a base is -4
                new_item = item + 1
                new_score = -1
                for increment in range(1, 4):
                    if (
                         new_item < seq_len and 
                         (seq[new_item] == 'C' or seq[new_item] == 'c')
                       ):
                        new_item += 1
                        new_score -= 1
                    else:
                        break
                bases_score += [ new_score for x in range(new_item - item) ]
                item = new_item

                # if more than 4 continuous C, all of them
                # are given score -4.
                while(item < seq_len and (seq[item] == 'C' or seq[item] == 'c')):
                    bases_score.append(-4)
                    item += 1
                
            else:
                # A, T, U and N
                bases_score.append(0)
                item += 1

        return bases_score


    def __CalScore(self, bases_score, k):
        '''
            Calculate the mean bases score of each
            width-k-windows.
            Noted that the sliding step is 1, which
            means overlapping will happen.

            input:
                    bases_score -- list of base score of a seq
                    k           -- window width

            output:
                    mean_scores_list -- list including mean score of windows
        '''
        mean_scores_list = []
        for i in range(len(bases_score) - k + 1):
            k_bases_mean = round(np.mean(bases_score[i:i+k]),2)
            mean_scores_list.append(k_bases_mean)
        return mean_scores_list

    def __Write(self, ofile, i, k, F, X, seq, l, score):
        line = "{} \t {} \t {} \t {} \t {}".format(
                i, i+k+F+X, seq, l, score)
        ofile.write(line)








if __name__ == '__main__':
    # Args' Parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--ifile', help='input seq file.')
    parser.add_argument('-o', '--outdir', help='ouput dir(with / at the end).')
    parser.add_argument('-w', '--wind', default=25, type=int,  
                        help='width of the slide-window.')
    parser.add_argument('-s', '--score', type=float, 
                        help='threshold of G4H score.')
    args = parser.parse_args()

    # Use os-module to replace the origin-version's
    # directly string process
    in_file_name = os.path.basename(args.ifile)
    out_file_dir = args.outdir + "/Results_" + in_file_name.split('.')[0]


    soft = Soft()
    with open("{}W{}_S{}.txt".format(out_file_dir, args.wind, args.score), 'w') as ofile1:
        with open("{}Merge.txt".format(out_file_dir), 'w') as ofile2:
            headers_seqs_scores = soft.GFinder(args.ifile, args.wind)
            for entry in headers_seqs_scores:
                g4_idx_list = soft.GetG4(ofile1, entry, args.score, args.wind)
                if(len(g4_idx_list) > 0):
                    mscore = soft.WriteSeq(
                                            ofile2,
                                            entry,
                                            g4_idx_list, 
                                            args.wind
                                           )




