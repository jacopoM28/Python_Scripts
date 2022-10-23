#!/usr/bin/env python
# coding: utf-8

# In[23]:



##################################################################################################
##################################################################################################
##                                                                                              ## 
##                                      From Blastx to fasta                                    ##
##                                                                                              ##
##################################################################################################
##################################################################################################
##NB: The input file should be a blastx outfmt 6 with standard format + qseq (Aligned part of query sequence)
#as 13 column
##Example: -outfmt "6 std qseq"

#Python libraries : Biopython, Pandas

#Author : Jacopo Martelossi
#E-mail : jacopo.martelossi2@unibo.it

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='From blastx outfmt 6 to pep fasta')
parser.add_argument('--in_tsv', required=True, help='input blastx outfmt 6 output')
parser.add_argument('--out', required=True, help='output prefix')

args = parser.parse_args()
OUT = args.out
IN = args.in_tsv

Sequences = []
df = pd.read_csv(IN,sep="\t",header=None)
for index, row in df.iterrows():
    Pep = SeqRecord(Seq(row[12].replace("-","")), id = row[0], 
                    description = "Start=" + str(row[6]) + " End=" + str(row[7]) )
    Sequences.append(Pep)
SeqIO.write(Sequences,"%s_Blastx.fa" %OUT, "fasta")


# In[15]:


Pep

