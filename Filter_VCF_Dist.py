#!/usr/bin/env python
# coding: utf-8

# In[7]:

#Python script usefull to create a bed file with regions corresponding to start and end of scaffolds.
#The dist option can be used to customize the number of nucleotides.
#The resulting bed file can be used for example to remove SVs of SNPs calling around end of scaffolds with SURVIVOR

import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Get blastx of the 20 best insertions of each sequence present in the library.')
parser.add_argument('--genome', required=True, help='input genome')
parser.add_argument('--out', required=True, help='prefix output file')
parser.add_argument('--dist', required=True, help='Minimum allowed distance from end of scaffold')

args = parser.parse_args()

GENOME =  args.genome
OUT = args.out
DIST = args.dist

df = open(OUT,'w')
for record in SeqIO.parse(GENOME, "fasta") :
    end = len(record)
    end_2 = len(record) - int(DIST)
    line_1 = "\t".join([record.id,"0",DIST])
    line_2 = "\t".join([record.id,str(end_2),str(end)])
    df.write(str(line_1))
    df.write('\n')
    df.write(str(line_2))
    df.write('\n')
df.close()
    
    

