#!/usr/bin/env python

#################################################################################
#################################################################################
##                       TE Annotatio summary script                           ##
#################################################################################
#################################################################################

#Script to summarize the results of a TE annotation and "quality" of a consensus library
#Six metrics are currently implemented :
#           1. Mean length of consensus sequences
#           2. Number of annotated insertions
#           3. Number of short annotated insertions (as reported by RepeatCraft (https://academic.oup.com/bioinformatics/article/35/6/1051/5079332)
#           4. Number of "Known" insertions
#           5. Mean length of insertions
#           6. Number of significant Blastx results of the 10% longest insertions against a reference TE - derived protein db

import os
import sys
import csv
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MafftCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

parser = argparse.ArgumentParser(description='Get summary statistics of a repeat library and annotation in fasta format')
parser.add_argument('--lib', required=True, help='Input TE library')
parser.add_argument('--genome', required=True, help='input genome')
parser.add_argument('--out', required=True, help='prefix output file')
parser.add_argument('--gff',required=True, help="gff3 file of annotated repeats. PS. First line of the file will be SKIPPED")
parser.add_argument('--TE_prot',required=True, help="TE derived protein library")
parser.add_argument('--num_threads',required=True, help="Number of cores for blastx")
parser.add_argument('--Lib_Name',required=True, help="Library name to use in the ouput table")

args = parser.parse_args()

#Parsing required argument
GENOME =  args.genome
LIBRARY = args.lib
OUT = args.out
GFF = args.gff
TE_ProtLib = args.TE_prot
N_threads = int(args.num_threads)
LibName = args.Lib_Name

mydict = {'Name' : [], "Cons.Len" : [],
       'N.ins' : [], "Known.Ins" : [],
       'N.Short.ins' : [], "Mean.Ins.Len" : [],
       'N.Blastx' : []}

################################################################################
#------------------------Mean library consensus length-------------------------#
################################################################################

i = 0
Length = 0
for seq_record in SeqIO.parse(LIBRARY, "fasta") :
    Length = Length + len(seq_record)
    i = i + 1
Mean_Len = Length/i                                                             #Calculating mean consensus length

################################################################################
#-----------------------GFF working (point 2,3,4,5)----------------------------#
################################################################################

Known_ins = 0
Nshort = 0

Scaffold = []
Start = []
End = []
Name = []
InsLen = []

with open(GFF) as f:
    first_line = f.readline()
    for line in f:                                                          
        InsLen.append(int(line.split("\t")[4]) - int(line.split("\t")[3]))        #Insertion length
        InsClass = line.split("\t")[2]                                            #Classification
        ID = line.split(";")[2].split("=")[1]                                     #Family name
        Scaffold.append(line.split("\t")[0])                                      #Scaffold name
        Start.append(line.split("\t")[3])                                         #Start coordinates
        End.append(line.split("\t")[4])                                           #End coordinates
        Name.append("#".join([ID,InsClass]))                                      #Concatenate Fam and Class
        InType = line.split(";")[3].split("=")[1]                                 #Insertion type (Fragmented or not)                                                   #Ins length
        if InsClass != "Unknown" :
            Known_ins += 1                                                        #Known insertions counter
        if InType == "T" :
            Nshort += 1                                                           #Fragmented insertions counter

Mean_InsLen = sum(InsLen)/len(InsLen)
Ins = len(InsLen)

##############################################################################
#-----------------Extract longest insertions and run Blastx------------------#
##############################################################################

#Dataframe with coordinates of insertions
Ins_df = pd.DataFrame({"scaffold" : Scaffold, "start" : Start, "end" : End, 
                          "TE" : Name, "Len" : InsLen})
Nrow = round(len(Ins_df)*0.10)                                                  #Calculate the 10% of the total insertions                                                   
Ins_df = Ins_df.sort_values(by=['Len'],ascending=False)                         #Sort by Insertion Length
Ins_df = Ins_df.head(Nrow)                                                      #keep only the 10% opf longest insertions
Ins_df= Ins_df.sort_index(ascending=True)                                       #Reorder by index
Ins_df.to_csv("%s_Insertions.bed" %OUT, sep="\t", index=False, header = False)          #Write bed file


with open("%s_Insertions.fa" %OUT, 'w') as extracted_faa :                           #Extraction of insertion copies 
    subprocess.run(["bedtools", "getfasta", "-name","-fi", 
                    GENOME , "-bed" , "%s_Insertions.bed" %OUT], 
                   stdout=extracted_faa, stderr=subprocess.DEVNULL)
    
db = NcbimakeblastdbCommandline(dbtype="prot",                                  #Create blastdb of TE proteins
                                   input_file= TE_ProtLib)
stdout, stderr = db()

blastn_cline = NcbiblastxCommandline(query = "%s_Insertions.fa" %OUT, num_threads = N_threads, #Blastx of insertions again TE proteins
                                     db=TE_ProtLib, evalue=1e-5, outfmt=6, 
                                     out="%s_Blast.out" %OUT,max_target_seqs = 1, max_hsps = 1)
stdout, stderr = blastn_cline()

nBlastxHit = 0                                                                    #Count number of hits
with open("%s_Blast.out" %OUT) as f:
    for line in f:
        nBlastxHit += 1
        
##############################################################################
#--------------------Append results and create ouput df----------------------#
##############################################################################

mydict['Name'].append(LibName)
mydict['Cons.Len'].append(Mean_Len)
mydict['N.Blastx'].append(nBlastxHit)
mydict['Known.Ins'].append(Known_ins)
mydict['N.ins'].append(Ins)
mydict['N.Short.ins'].append(Nshort)
mydict['Mean.Ins.Len'].append(Mean_InsLen)

df = pd.DataFrame(mydict)
df.to_csv("%s_Summary.tsv" %OUT, sep="\t", index=False, header = True)
