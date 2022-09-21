#!/usr/bin/env python

#################################################################################
#################################################################################
##                                                                             ##   
##                         TE Annotation summary script                        ##
##                                                                             ##
#################################################################################
#################################################################################

#Script to summarize the results of a TE annotation and "quality" of a consensus library
#Seven metrics are currently implemented :
#           1. Mean length of consensus sequences
#           2. Number of annotated insertions
#           4. Number of "Known" insertions
#           5. Mean length of insertions
#           6. Number of significant Blastx results of the 10% longest insertions against a reference TE - derived protein db
#           7. Number of significant Blastx results on consensuxs sequences
##IMPORTANT: No header must be present in the GFF file or the script will rise an error "list indes out of range" in line 104

#Dependecies necessary : Blast, Bedtools Trimal
#Python libraries : Biopython, Pandas

#Author : Jacopo Martelossi
#E-mail : jacopo.martelossi2@unibo.it

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

print("\n")
print("#################################################")
print("#EVALUATING THE LIBRARY AND THE ANNOTATION...   #")
print("#################################################")
print("\n")

mydict = {'Name' : [], "Cons.Len" : [],
       'N.ins' : [], "Known.Ins" : [], "Mean.Ins.Len" : [],
       'N.ins_Blastx' : [], "N.cons_Blastx" : []} 

################################################################################
#------------------------Mean library consensus length-------------------------#
################################################################################

print("\n")
print("Calculating mean consensus length...")

#Calculating mean consensus length
i = 0
Length = 0
for seq_record in SeqIO.parse(LIBRARY, "fasta") :
    Length = Length + len(seq_record)
    i = i + 1
Mean_Len = Length/i                                                             

print("...Done")

################################################################################
#-----------------------GFF working (point 2,3,4,5)----------------------------#
################################################################################

print("\n")
print("Creating summary statistics of the annotation...")
Known_ins = 0

Scaffold = []
Start = []
End = []
Name = []
InsLen = []

with open(GFF) as f:
    first_line = f.readline()
    for line in f:                                                         
        #Insertion length
        InsLen.append(int(line.split("\t")[4]) - int(line.split("\t")[3]))
        #Classification
        InsClass = line.split("\t")[2]
        #Family name
        ID = line.split(";")[2].split("=")[1]
        #Scaffold name
        Scaffold.append(line.split("\t")[0])
        #Start coordinates
        Start.append(line.split("\t")[3])
        #End coordinates
        End.append(line.split("\t")[4])
        #Concatenate Fam and Class
        Name.append("#".join([ID,InsClass]))
        #Known insertions counter
        if InsClass != "Unknown" :
            Known_ins += 1                                                                                                    
#Mean ins len
Mean_InsLen = sum(InsLen)/len(InsLen)
#N. ins
Ins = len(InsLen)

print("...Done")

##############################################################################
#-----------------Extract longest insertions and run Blastx------------------#
##############################################################################

#----------------------------Blastx of insertions----------------------------#

print("\n")
print("Searching for TE protein similarities in the insertions...")

#Dataframe with coordinates of insertions
Ins_df = pd.DataFrame({"scaffold" : Scaffold, "start" : Start, "end" : End, 
                          "TE" : Name, "Len" : InsLen})
#Calculate the 10% of the total insertions
Nrow = round(len(Ins_df)*0.10)                                                       
#Sort by Insertion Length
Ins_df = Ins_df.sort_values(by=['Len'],ascending=False)                         
#keep only the 10% opf longest insertions
Ins_df = Ins_df.head(Nrow)
#Reorder by index
Ins_df= Ins_df.sort_index(ascending=True)
#Write bed file
Ins_df.to_csv("%s_Insertions.bed" %OUT, sep="\t", index=False, header = False)

#Extraction of insertion copies
with open("%s_Insertions.fa" %OUT, 'w') as extracted_faa :                            
    subprocess.run(["bedtools", "getfasta", "-name","-fi", 
                    GENOME , "-bed" , "%s_Insertions.bed" %OUT], 
                   stdout=extracted_faa, stderr=subprocess.DEVNULL)

#Create blastdb of TE proteins
db = NcbimakeblastdbCommandline(dbtype="prot",                           
                                   input_file= TE_ProtLib)
stdout, stderr = db()

#Blastx of insertions again TE proteins
blastn_cline = NcbiblastxCommandline(query = "%s_Insertions.fa" %OUT, num_threads = N_threads,
                                     db=TE_ProtLib, evalue=1e-5, outfmt=6, 
                                     out="%s_Ins_Blast.out" %OUT,max_target_seqs = 1, max_hsps = 1)
stdout, stderr = blastn_cline()

#Count number of hits
nBlastxHit = 0
with open("%s_Ins_Blast.out" %OUT) as f:
    for line in f:
        nBlastxHit += 1

print("...Done")
        
#----------------------Blastx of consensus sequences-------------------------#
print("\n")
print("Searching for TE protein similarities in the consensus...")

blastn_cline = NcbiblastxCommandline(query = LIBRARY, num_threads = N_threads,
                                     db=TE_ProtLib, evalue=1e-5, outfmt=6, 
                                     out="%s_Cons_Blast.out" %OUT,max_target_seqs = 1, max_hsps = 1)
stdout, stderr = blastn_cline()

#Count number of hits
nBlastxHit_cons = 0
with open("%s_Cons_Blast.out" %OUT) as f:
    for line in f:
        nBlastxHit_cons += 1

print("...Done")
        
##############################################################################
#--------------------Append results and create ouput df----------------------#
##############################################################################

print("\n")
print("Collecting all results...")

mydict['Name'].append(LibName)
mydict['Cons.Len'].append(Mean_Len)
mydict['N.ins_Blastx'].append(nBlastxHit)
mydict['N.cons_Blastx'].append(nBlastxHit_cons)
mydict['Known.Ins'].append(Known_ins)
mydict['N.ins'].append(Ins)
mydict['Mean.Ins.Len'].append(Mean_InsLen)

df = pd.DataFrame(mydict)
df.to_csv("%s_Summary.tsv" %OUT, sep="\t", index=False, header = True)

print("...Done")
