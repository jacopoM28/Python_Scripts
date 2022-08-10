#!/usr/bin/env python
# coding: utf-8

# In[ ]:


##################################################################################################
##################################################################################################
##                                                                                              ## 
##                                  FULL-AUTO BLAST-EXTEND-EXTRACT                              ##
##                                                                                              ##
##################################################################################################
##################################################################################################

#Automatic Blast-Extend-Extract process. Raw consensus libraries are blasted against the genome
#a user-defined number of hits are extracted and aligned. Poorly aligned regions are removed with 
#trimal and a consensus sequences is generated with CIAlign. All optional arguments 
#have default values. Flags allow to personalize the process, if you want to change more parameters
#for blast, mafft and/or number of hits to align, you should change the corresponding lines. 
#By default a maximum of 30 sequences are extracted to build up consensus, if you want to change this behaviour 
#look at line 151

##IMPORTANT!!!!!!
#This version of the script implements CIAalign for building up consensus, for hardly trimmed sequences it works very well.
#However, don't use it if you want to manually curate consensus sequences, for that is better :
#https://github.com/jacopoM28/Python_Scripts/blob/main/Blast_Extend_Extract_v2.0.py).
#This script can be run mutiple times to improve quality of consensus sequences.

#Dependecies necessary : Blast, Bedtools, Mafft, Trimal, CIAlign.
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
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

parser = argparse.ArgumentParser(description='Automatic Blast-Extend-Extract script')
parser.add_argument('--genome', required=True, help='Input genome fasta file')
parser.add_argument('--lib', required=True, help='raw consensus library')
parser.add_argument('--out', required=True, help='prefix output file')
parser.add_argument('--cutoff', help='Minimum length cutoff; default = 50')
parser.add_argument('--blast_evalue', help='Blast evalue treshold; default = 10e-10')
parser.add_argument('--blast_query_cov', help='query coverage cutoff value; default = 70')
parser.add_argument('--min_Blast_Hits', help='Minum number of blast hits; default = 5')
parser.add_argument('--num_threads', help='Number of cpus for Blast and mafft; default = 1')
parser.add_argument('--extension', help='number of bp to extend blast hits; default = 1000pb')
parser.add_argument('--blast_identity', help='Minimum identity treshold between consensus and each copy on the genome to be considered a significant hit; default 70')
parser.add_argument('--gt', help='1 - (fraction of sequences with a gap allowed) default = 0.6')
parser.add_argument('--st', help='Minimum average similarity allowed default = 0.6')

args = parser.parse_args()

#Parsing required argument
GENOME =  args.genome
LIBRARY = args.lib
OUT = args.out

#Parsing arguments with default values

if args.blast_identity is None :
    IDENTITY = 70
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
    QUERY_COV = 70
else :
    QUERY_COV = args.blast_query_cov
if args.min_Blast_Hits is None :
    MIN_HITS = 5
else :
    MIN_HITS = args.min_Blast_Hits
if args.num_threads is None :
    NTHREADS = 1
else :
    NTHREADS = args.num_threads
if args.extension is None :
    EXTENSION = 1000
else :
    EXTENSION = args.extension
if args.gt is None :
    GT = 0.6
else :
    GT = args.gt
if args.st is None :
    ST = 0.6
else :
    ST = args.gt
    
print("\n")
print("#################################################")
print("#STARTING AUTOMATIC BLAST-EXTEND-EXTRACT PROCESS#")
print("#################################################")
print("\n")
    
os.mkdir("intermediate_files")

#######################################################################################################    
##-----------------------------Cleaning raw library based on cut-off value----------------------------#
#######################################################################################################

print("Filtering raw consensus library based on a minimum length of :", CUTOFF)

output_name = OUT + "_" + str(CUTOFF) + "pb"
New_Sequences = []

#Remove sequences shorter than cutoff value
with open(LIBRARY) as handle :
    for record in SeqIO.parse(handle, "fasta") :
        if len(record) >= int(CUTOFF) :
            New_Sequences.append(record)
SeqIO.write(New_Sequences, "./intermediate_files/%s.fasta" %output_name ,"fasta")

print("----> Your filtered raw consensus library has :", len(New_Sequences), "sequences")
print("########DONE")


##########################################################################################################
##------------------------------------------Blast--------------------------------------------------------#
##########################################################################################################

print("Blasting filtered consensus library with a minimum query coverage of :", QUERY_COV)

#Create Blastdb
db = NcbimakeblastdbCommandline(dbtype="nucl",
                                   input_file= "Lapu.genomic.fa")
stdout, stderr = db()
#Blastn command
blastn_cline = NcbiblastnCommandline(query = "./intermediate_files/%s.fasta" %output_name, num_threads = NTHREADS,
                                     db="Lapu.genomic.fa" ,qcov_hsp_perc = QUERY_COV,
                                     evalue=EVALUE, outfmt=6, out="./intermediate_files/Blast.out", 
                                     perc_identity = IDENTITY)

stdout, stderr = blastn_cline()

#Keep only 30 hits (based on evalue)
bashCommand = "awk '{ if (++count[$1] <= 30) print $0 }' ./intermediate_files/Blast.out> ./intermediate_files/Blast.out.filtered"
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

#Extract them
with open("./intermediate_files/BlastExtendend.fasta", 'w') as extracted_faa :
    subprocess.run(["bedtools", "getfasta", "-name+", "-s","-fi", 
                    GENOME , "-bed" , "./intermediate_files/Blast_Extended.bed"], 
                   stdout=extracted_faa, stderr=subprocess.DEVNULL)
    
#Parsing each transposon in a different multifasta 
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

    
#########################################################################################################
#---------------------------------------------ALIGN-----------------------------------------------------#
#########################################################################################################

os.mkdir("Alignments")

directory = './Extracted_Sequences'

i = 1
for filename in os.listdir(directory):
    #Check that we are using only fasta files and setting up mafft parameters
    if filename.endswith(".fasta") :
        print(filename, "(", i, "/", Blast_bed['query'].nunique(), ")")
        in_file = os.path.join(directory, filename)
        mafft_cline = MafftCommandline(input=in_file)
        mafft_cline.thread = int(NTHREADS)
        mafft_cline.adjustdirection = True
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

os.mkdir("Cleaned_Alignments")
directory = './Alignments/'

for filename in os.listdir(directory):
    if filename.endswith(".mafft") :
        in_file = os.path.join(directory, filename)
        upper = []
        #Make all sequences to uppercases for compatibility with Trimal
        for record in SeqIO.parse(in_file, "fasta"):
            record.seq = record.seq.upper()
            upper.append(record)
        SeqIO.write(upper, "%s" %in_file , "fasta")
        #Remove poorly aligned positions
        subprocess.run(["trimal", "-gt", str(GT), "-st", str(ST), "-in", in_file, "-out", "./Cleaned_Alignments/%s.trimmed" %filename ],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
print("Building up consensus...")
print("\n")

############################################################################################################
#------------------------------------------------Consensus-------------------------------------------------#
############################################################################################################

os.mkdir("Cons_Alignments")
directory = './Cleaned_Alignments/'

for filename in os.listdir(directory):
    if filename.endswith(".trimmed") :
        in_file = os.path.join(directory, filename)
        seq_name = filename.split(".")[0]
        print(seq_name)
        #Extracting annotation of the elements
        with open(in_file) as f:
            first_line = f.readline()
            clas = first_line.split("__")[0].replace(">","").split("#")[1]
            print(clas)
        #Build up consensus with CIAlign and ignorig gaps (They have been already delated by Trimal)
        subprocess.run(["CIAlign", "--infile", in_file, "--outfile_stem", "./Cons_Alignments/%s" %seq_name,
                        "--make_consensus", "--consensus_name", "%s#%s" % (seq_name,clas),
                        "--consensus_type", "majority_nongap"],
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        

#Putting all consensus together
directory = "./Cons_Alignments"

consensus = []
for filename in os.listdir(directory):
    if filename.split("_")[-1] == "consensus.fasta" and filename.split("_")[-2] != "with" :
        in_file = os.path.join(directory, filename)
        for record in SeqIO.parse(in_file, "fasta") :
            consensus.append(record)
        SeqIO.write(consensus, "%s_Consensus.fa" %OUT , "fasta")        

print("\n")
print("All extended consensus sequences can be found in : %s_Consensus.fa" %OUT)
        
print("\n")
print("#################################################")
print("#################EVERYTHING DONE#################")
print("#################################################")
print("\n")
