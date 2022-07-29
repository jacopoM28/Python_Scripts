#!/usr/bin/env python
# coding: utf-8

# In[3]:


##################################################################################################
##################################################################################################
##                                                                                              ## 
##                                  MERGE CONSENSUS LIBRARIES                                   ##
##                                                                                              ##
##################################################################################################
##################################################################################################

#Usefull script to merge multiple, possible redundant (i.e. with duplicated entries) consensus libraries.
#The script will read each file in the input directory (assumed to be a fasta file) and will search for
#sequences with the same name before the classification of the element (format similar to "Lapu_rnd-1_64#DNA").
#For each redundant entry it will keep the classified one or the longest one if both are unclassified.
#It is possible to provide a manually curated library assumed to be correct that will recive maximum priority
#(i.e. its consensus sequences will always be selected). Sequences coming from this library will have the MANUAL
#suffix

##IMPORTANT!!!!!!
#This script assume that all redundant sequences have exatly the same name before the classification,
#also in the manually curated library (i.e >sequence_name#Classification). Moreover the input directory
#must contain only libraries that you want to analyze.

#Python libraries : Biopython, Pandas

#Author : Jacopo Martelossi
#E-mail : jacopo.martelossi2@unibo.it

#------------------------------------------Libraries--------------------------------------------#

import os
import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo

#-----------------------------------Parsing arguments-------------------------------------------#

parser = argparse.ArgumentParser(description='Merge consensus TE libraries')
parser.add_argument('--path', required=True, help='path to consensus libraries, with final slash ')
parser.add_argument('--manual', help='name of the manually curated library that will recive maximum priority. Do not specify the path, it is assumed to be the same of others libraries')
parser.add_argument('--out', required=True, help='output prefix')

args = parser.parse_args()
PATH = args.path
OUT = args.out
if args.manual is None :
    CONDITION = 0
else :
    CONDITION = 1
    MANUAL = args.manual
                    
#Empty lists to store consensus sequences and consensus names without classification (to find index only)
Seqs = []
Check_List = []

#---------------------------------------------MAIN---------------------------------------------#
                    
#Read all fasta file in a directory and selection of "best" consensus sequences
file = os.listdir(PATH)

#Do not read manually curated library if present
if CONDITION == 1 :
    file.remove(MANUAL) 
    
for filename in file :
        in_file = os.path.join(PATH, filename)
        with open(in_file) as handle :
            for record in SeqIO.parse(handle, "fasta") :
                name = record.id.split("#")[0]
                #Check if the consensus is already in the list, if not add it
                if name not in Check_List :
                    Seqs.append(record)
                    Check_List.append(name)
                else :
                    my_index = Check_List.index(name)
                    #if an already present consensus is found in another library, keep the classified one
                    if Seqs[my_index].id.split("#")[1] == "Unknown" and record.id.split("#")[1] != "Unknown" :
                        Seqs[my_index] = record
                    #If both are "Unknown" keep the longest one
                    elif Seqs[my_index].id.split("#")[1] == "Unknown" and record.id.split("#")[1] == "Unknown" :
                        if len(record.seq) >= len(Seqs[my_index].seq) :
                            Seqs[my_index] = record
                        
#Read manually curated consensus and keep them
if CONDITION == 1 :
    with open(os.path.join(PATH,MANUAL)) as handle :
        for record in SeqIO.parse(handle, "fasta") :
            name = record.id.split("#")[0] 
            if name not in Check_List :
                record.description = "MANUAL"
                Seqs.append(record)
                Check_List.append(name)
            else : 
                my_index = Check_List.index(name)
                Seqs[my_index] = record
                Seqs[my_index].description = "MANUAL"

#---------------------------------------------OUT---------------------------------------------#

#Write output fasta file
SeqIO.write(Seqs, "%s_Merged.Consensus.fa" %OUT , "fasta")

