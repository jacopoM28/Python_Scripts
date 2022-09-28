#!/usr/bin/env python
# coding: utf-8

# In[ ]:


##################################################################################################
##################################################################################################
##                                                                                              ## 
##                                        BLAST-EXTEND-EXTRACT                                  ##
##                                                                                              ##
##################################################################################################
##################################################################################################
#Author : Jacopo Martelossi
#E-mail : jacopo.martelossi2@unibo.it


#Automatic Blast-Extend-Extract process. Raw consensus libraries are blasted against the genome
#a user-defined number of hits are extracted and aligned together with the raw consensus. Finally,
#all positions in which more than the 99% of sequences have gaps are removed. All optional arguments 
#have default values. Flags allow to personalize the process, if you want to change more parameters
#for blast and/or mafft you should change the corresponding lines. By default a maximum of 50 sequences 
#are extracted to build up consensus, if you want to change this behaviour look at line 140

#V2: Removed max_hsps and max_target_seqs; added an awk command to kee only top 50 blast hits of raw cons
#V3: Added merging of overlapping or close (<100bp) blast hits of raw consensus

##IMPORTANT!!!!!!
#Consensus sequences are created with emboss cons, gaps are coded as <n>, all sequences MUST be
#manually curated!

#Dependecies necessary : Blast, Bedtools, Mafft, Trimal, emboss.
#Python libraries : Biopython, Pandas

#Everything can be installed through Conda!

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
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

parser = argparse.ArgumentParser(description='Blast-Extend-Extract script')
parser.add_argument('--genome', required=True, help='Input genome fasta file')
parser.add_argument('--lib', required=True, help='raw consensus library')
parser.add_argument('--out', required=True, help='prefix output file')
parser.add_argument('--mafft_strategy', choices=['einsi', 'auto'],help='Mafft_Strategy', required=True, type=str)
parser.add_argument('--cutoff', help='Minimum length cutoff; default = 50')
parser.add_argument('--blast_evalue', help='Blast evalue treshold; default = 10e-10')
parser.add_argument('--blast_query_cov', help='query coverage cutoff value; default = 80')
parser.add_argument('--min_Blast_Hits', help='Minum number of blast hits; default = 10')
parser.add_argument('--num_threads', help='Number of cpus for Blast and mafft; default = 1')
parser.add_argument('--extension', help='number of bp to extend blast hits; default = 2000pb')
parser.add_argument('--blast_identity', help='Minimum identity treshold between consensus and each copy on the genome to be considered a significant hit; default 80')

args = parser.parse_args()

#Parsing required argument
GENOME =  args.genome
LIBRARY = args.lib
OUT = args.out
MAFFT_OPT = args.mafft_strategy

#Parsing arguments with default values

if args.blast_identity is None :
    IDENTITY = 80
else :
    IDENTITY = args.blast_identity
if args.cutoff is None:
    CUTOFF = 50
else :
    CUTOFF = args.cutoff
if args.blast_evalue is None :
    EVALUE = 10e-10
else:
    EVALUE = args.blast_evalue
if args.blast_query_cov is None :
    QUERY_COV = 80
else :
    QUERY_COV = args.blast_query_cov
if args.min_Blast_Hits is None :
    MIN_HITS = 10
else :
    MIN_HITS = args.min_Blast_Hits
if args.num_threads is None :
    NTHREADS = 1
else :
    NTHREADS = args.num_threads
if args.extension is None :
    EXTENSION = 2000
else :
    EXTENSION = args.extension

print("\n")
print("#######################################")
print("#STARTING BLAST-EXTEND-EXTRACT PROCESS#")
print("#######################################")
print("\n")

os.mkdir("intermediate_files")

#######################################################################################################    
##-----------------------------Cleaning raw library based on cut-off value----------------------------#
#######################################################################################################

print("Filtering raw consensus library based on a minimum length of :", CUTOFF)

output_name = OUT + "_" + str(CUTOFF) + "pb"
New_Sequences = []

with open(LIBRARY) as handle :
    for record in SeqIO.parse(handle, "fasta") :
        if len(record) >= int(CUTOFF) :
            New_Sequences.append(record)
SeqIO.write(New_Sequences, "%s.fasta" %output_name ,"fasta")

print("----> Your filtered raw consensus library has :", len(New_Sequences), "sequences")
print("########DONE")

##########################################################################################################
##------------------------------------------Blast--------------------------------------------------------#
##########################################################################################################

print("Blasting filtered consensus library with a minimum query coverage of :", QUERY_COV)

db = NcbimakeblastdbCommandline(dbtype="nucl",
                                   input_file= GENOME)
stdout, stderr = db()
blastn_cline = NcbiblastnCommandline(query = "%s.fasta" %output_name, num_threads = NTHREADS,
                                     db=GENOME ,qcov_hsp_perc = QUERY_COV,
                                     evalue=EVALUE, outfmt=6, out="./intermediate_files/Blast.out", perc_identity = IDENTITY)
print("---> Blast command line :", blastn_cline,".")
print("########DONE")

stdout, stderr = blastn_cline()

bashCommand = "awk '{ if (++count[$1] <= 50) print $0 }' ./intermediate_files/Blast.out> ./intermediate_files/Blast.out.filtered" #Change the default 50 to keep more or less sequences
process = subprocess.Popen(bashCommand, stdout=subprocess.PIPE,shell=True)
output, error = process.communicate()

#########################################################################################################
##------------------------------------FROM BLAST TO BED-------------------------------------------------#
#########################################################################################################
query_count = {}
query = []
scaffold = []
start = []
end = []
strand = []
evalue = []

with open("./intermediate_files/Blast.out.filtered") as tsv:
    for line in csv.reader(tsv, delimiter="\t") :
        query.append(line[0])
        scaffold.append(line[1])
        start_coord = int(line[8])
        end_coord = int(line[9])
        evalue.append(line[10])
        #Switch coordinates for hits on the negative strand
        if start_coord > end_coord :
            start_coord = int(line[9])
            end_coord = int(line[8])
            strand.append("-")
        else :
            strand.append("+")
        start.append(start_coord)
        end.append(end_coord)
        if not line[0] in query_count :
            query_count[line[0]] = 1
        else :
            query_count[line[0]] += 1

#Creation of a dataframe with coordinates
Blast_bed = pd.DataFrame({"scaffold" : scaffold, "start" : start, "end" : end, "query" : query,
                          "evalue" : evalue ,"strand" : strand})

#Filter transposon based on number of hits
for key,value in query_count.items() :
    if value < int(MIN_HITS) :
        indexNames = Blast_bed[ Blast_bed['query'] == key ].index
        Blast_bed.drop(indexNames , inplace=True)    

#Saving filtered dataframe in a bed file
Blast_bed.to_csv("./intermediate_files/Blast.bed", sep="\t", index=False, header = False)

print("----> After filtering you are analyzing :", Blast_bed['query'].nunique(), "consensus with more than :", 
      int(MIN_HITS) -1, "hits on the genome."  )
print("########DONE")
#########################################################################################################
##---------------------------------------------EXTEND and EXTRACT---------------------------------------#
#########################################################################################################
print("Extending the blast hits to", int(EXTENSION), "pb and filtering based on number of blast hits")

subprocess.run(["samtools", "faidx", GENOME],stderr=subprocess.DEVNULL)

#Extend blast hits
with open("./intermediate_files/Blast_Extended.bed", 'w') as Extended_Hits :
    subprocess.run(["bedtools", "slop", "-i", "./intermediate_files/Blast.bed",
                    "-g", "%s.fai" %GENOME, "-b", str(EXTENSION)], 
                   stdout=Extended_Hits, stderr=subprocess.DEVNULL)
    
#Sort the bed file
with open("./intermediate_files/Blast_Extended.sorted.bed", 'w') as Extended_sorted_Hits :
    subprocess.run(["bedtools", "sort", "-i", "./intermediate_files/Blast_Extended.bed"], 
                   stdout=Extended_sorted_Hits, stderr=subprocess.DEVNULL)
    
#---------------------Merge overlapping hits coming from the same consensus-----------------------------# 
df = pd.read_csv("intermediate_files/Blast_Extended.sorted.bed",sep="\t",header = None)
df.columns = ['scaffold', 'start', 'end', 'name','evalue','strand']

#Merge scaffold and name columns. This allow us to have unique identifiers
df["scaffold"] = df["name"] + "_MERGE_" + df["scaffold"]
df.sort_values(['scaffold', 'start'])
df.to_csv("./intermediate_files/Blast_Extended.tmp",sep="\t",index=False, header = False)

#Sort again
with open("./intermediate_files/Blast_Extended.sorted.tmp", 'w') as Extended_sorted_Hits :
    subprocess.run(["bedtools", "sort", "-i", "./intermediate_files/Blast_Extended.tmp"], 
                   stdout=Extended_sorted_Hits, stderr=subprocess.DEVNULL)

#Merge overlapping or closer than 100bp intervals
with open("./intermediate_files/Blast_Extended.sorted.merged.bed", 'w') as Extended_merged_Hits :
    subprocess.run(["bedtools", "merge", '-d', '100' , '-s', "-c","4,5,6", "-o", "distinct,mean,distinct", "-i",
                    "./intermediate_files/Blast_Extended.sorted.tmp"], 
                   stdout=Extended_merged_Hits, stderr=subprocess.DEVNULL)

#Remove the name ID from the first column
bashCommand = "sed -i 's/^.*_MERGE_//' ./intermediate_files/Blast_Extended.sorted.merged.bed"
process = subprocess.Popen(bashCommand, stdout=subprocess.PIPE,shell=True)
output, error = process.communicate()

bashCommand = '''awk -v OFS="\t" -F"\t" '{print$1,$2,$3,$5,$6,$7}' ./intermediate_files/Blast_Extended.sorted.merged.bed > ./intermediate_files/Blast_Extended.sorted.merged_FINAL.bed'''
process = subprocess.Popen(bashCommand, stdout=subprocess.PIPE,shell=True)
output, error = process.communicate()

with open("./intermediate_files/BlastExtendend.fasta", 'w') as extracted_faa :
    subprocess.run(["bedtools", "getfasta", "-name", "-s","-fi", 
                    GENOME , "-bed" , "./intermediate_files/Blast_Extended.sorted.merged_FINAL.bed"], 
                   stdout=extracted_faa, stderr=subprocess.DEVNULL)

#Parsing each transposon in a different multifasta and add original consensus    
os.mkdir("Extracted_Sequences")

TransposonsID = [] 

for record in SeqIO.parse("./intermediate_files/BlastExtendend.fasta", "fasta"):
    ID = record.id.split("#")
    TransposonsID.append(ID[0])
TransposonsID = set(TransposonsID)

for element in TransposonsID :
    Sequences = []
    for record in SeqIO.parse("./intermediate_files/BlastExtendend.fasta", "fasta"):
        if record.id.split("#")[0] == element :
            #Susbtituting : with _ in output fasta files for downstream analyses
            record.id = record.id.replace(":", "_")
            record.description = record.description.replace(":", "_")
            Sequences.append(record)
    SeqIO.write(Sequences,"Extracted_Sequences/%s.fasta" %element, "fasta")
    
print("########DONE")

#########################################################################################################
#---------------------------------------------ALIGN-----------------------------------------------------#
#########################################################################################################

print("Aligning with mafft and", MAFFT_OPT, "strategy.")

os.mkdir("Alignments")

directory = './Extracted_Sequences'

i = 1
for filename in os.listdir(directory):
    #Check that we are using only fasta files and setting mafft parameters
    if filename.endswith(".fasta") :
        print(filename, "(", i, "/", Blast_bed['query'].nunique(), ")")
        in_file = os.path.join(directory, filename)
        mafft_cline = MafftCommandline(input=in_file)
        mafft_cline.thread = int(NTHREADS)
        mafft_cline.adjustdirection = True
        if MAFFT_OPT == 'einsi' :
            mafft_cline.genafpair = True
            mafft_cline.maxiterate = 1000
        stdout, stderr = mafft_cline()
        i += 1
        
        with open("./Alignments/%s.mafft" %filename, "w") as output:
            output.write(stdout)
            align = AlignIO.read("./Alignments/%s.mafft" %filename, "fasta")
            
print("########DONE")
            
############################################################################################################
#----------------------------------------------CLEAN-------------------------------------------------------#
############################################################################################################

print("Cleaning alignments...")
print("\n")

os.mkdir("Cleaned_Alignments")
directory = './Alignments/'

for filename in os.listdir(directory):
    if filename.endswith(".mafft") :
        in_file = os.path.join(directory, filename)
        with open(in_file, "r+") as handle :
            subprocess.run(["trimal", "-gt", "0.05", "-in", in_file, "-out", "./Cleaned_Alignments/%s.trimmed" %filename ],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
############################################################################################################
#-------------------------------------------Build up consensus---------------------------------------------#
############################################################################################################

print("Building up consensus...")
print("\n")

os.mkdir("Cons_Alignments")
directory = './Cleaned_Alignments/'

for filename in os.listdir(directory):
    if filename.endswith(".trimmed") :
        in_file = os.path.join(directory, filename)
        seq_name = filename.split(".")[0]
        print(seq_name)
        with open(in_file, "r+") as handle :
            subprocess.run(["cons", "-name", "%s_cons" %seq_name, "-sequence" , in_file, "-outseq",  "./Cons_Alignments/%s.cons" %filename, "-plurality", "2" ],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

print("\n")

print("###############################################")
print("#EVERYTHING DONE, HAVE A NICE MANUAL CURATION!#")
print("###############################################")
