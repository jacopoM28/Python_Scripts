#!/usr/bin/env python
# coding: utf-8

# In[30]:


import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse

#Parisng arguments
parser = argparse.ArgumentParser(description='From TIDEHunter tsv to fasta')
parser.add_argument('--in_tsv', required=True, help='Input TIDEHunter output file')
parser.add_argument('--out', required=True, help='prefix output fasta file')

args = parser.parse_args()
OUT = args.out
IN = args.in_tsv

#------------------------------------------MAIN--------------------------------------------------#
df = pd.read_csv(IN,sep="\t",header=None)
df.set_axis(['readName', 'repN', 'copyNum', 'readLen', 'start',
             'end','consLen','aveMatch','fullLen','subPos','consSeq'], axis=1, inplace=True)

#Empty list to store tandem repeats sequences
Tandems = []
for ind in df.index:
    #Keep only consensus longer than 50bp
    if df['consLen'][ind] >= 50 :
        #Create repeat name
        name = df['readName'][ind] + "_" + df['repN'][ind]
        #Store tandem repeat sequences
        seq = df['consSeq'][ind].strip()
        #Create description
        descript = "copyNum=" + str(df['copyNum'][ind]) + " " + "readLen=" + str(df['readLen'][ind])
        #Create SeqRecord object
        record = SeqRecord(Seq(seq),id=name, name="Tandem repeat",
                   description=descript)
        Tandems.append(record)
#------------------------------------------FINISHED--------------------------------------------------#
#Write everything to a fasta file
SeqIO.write(Tandems, "%s_Tandem.Repeats.fasta" %OUT, "fasta")
#Remove redundancy in the output
subprocess.run('cd-hit-est -T 0 -i %s_Tandem.Repeats.fasta -o %s_Tandem.Repeats_nr.fasta -d 0' % (OUT, OUT), shell=True)

