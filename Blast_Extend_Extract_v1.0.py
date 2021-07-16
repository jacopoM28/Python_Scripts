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
#all positions in which more than the 95% of sequences have gaps are removed. All optional arguments 
#have default values. Flags allow to personalize the process, if you want to change more parameters
#for blast and/or mafft you should change the corresponding lines.

#Dependecies necessary : Blast, Bedtools, Mafft, Trimal.
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
parser.add_argument('--hits_to_keep', help='N. of hits to keep for blast; default = 25')
parser.add_argument('--blast_query_cov', help='query coverage cutoff value; default = 80')
parser.add_argument('--min_Blast_Hits', help='Minum number of blast hits; default = 10')
parser.add_argument('--num_threads', help='Number of cpus for Blast and mafft; default = 1')
parser.add_argument('--extension', help='number of bp to extend blast hits; default = 2000pb')

args = parser.parse_args()

#Parsing required argument
GENOME =  args.genome
LIBRARY = args.lib
OUT = args.out
MAFFT_OPT = args.mafft_strategy

#Parsing arguments with default values

if args.cutoff is None:
    CUTOFF = 50
else :
    CUTOFF = args.cutoff
if args.blast_evalue is None :
    EVALUE = 10e-10
else:
    EVALUE = args.blast_evalue
if args.hits_to_keep is None :
    NHITS = 25
else :
    NHITS = args.hits_to_keep
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

#Print error and exit if number of blast hits to keep is lesser than minimum number of required hits
if int(NHITS) < int(MIN_HITS) :
    print("\n")
    print("ERROR : Number of blast hits to keep cannot be lower than minium number of required hits!", file = sys.stderr)
    sys.exit()

print("\n")
print("#######################################")
print("#STARTING BLAST-EXTEND-EXTRACT PROCESS#")
print("#######################################")
print("\n")

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

print("Blasting filtered consensus library keeping", NHITS, "hits with a minimum query coverage of :", QUERY_COV)

db = NcbimakeblastdbCommandline(dbtype="nucl",
                                   input_file= GENOME)
stdout, stderr = db()
blastn_cline = NcbiblastnCommandline(query = "%s.fasta" %output_name, max_hsps=1, num_threads = NTHREADS,
                                     max_target_seqs = NHITS, db=GENOME ,qcov_hsp_perc = QUERY_COV,
                                     evalue=EVALUE, outfmt=6, out="%s_Blast.out" %OUT)
print("---> Blast command line :", blastn_cline,".")
print("########DONE")

stdout, stderr = blastn_cline()

#########################################################################################################
##----------------------------------------------EXTEND--------------------------------------------------#
#########################################################################################################

print("Extending the blast hits to", int(EXTENSION), "pb and filtering based on number of blast hits")

scaffold_len = {}

#Storing informations about length of scaffolds...
for record in SeqIO.parse(GENOME, "fasta"):
    scaffold_len[record.id] = int(len(record.seq))


query_count = {}
query = []
scaffold = []
start = []
end = []
strand = []
random_line = []

#Extend
with open("%s_Blast.out" %OUT) as tsv:
    for line in csv.reader(tsv, delimiter="\t") :
        random_line.append("1")
        query.append(line[0])
        scaffold.append(line[1])
        start_coord = int(line[8])
        end_coord = int(line[9])
        #if starting coordinates > ending means that the transposon is present in the opposite strand. 
        #We can just switch starting and ending coordinates and add negative strand feature
        if start_coord > end_coord :
            start_coord = int(line[9])
            end_coord = int(line[8])
            strand.append("-")
        else :
            strand.append("+")
        #If the extensin goes out the beggining of the scaffold, stop at the 1st base
        if start_coord - int(EXTENSION) < 1 :
            start.append(0)
        else :
            start.append(start_coord - int(EXTENSION))
        #If the extension goes out the end of the scaffold just use the end of the scaffold
        if end_coord + int(EXTENSION) > scaffold_len[line[1]] :
            end.append(scaffold_len[line[1]])
        else :
            end.append(end_coord + int(EXTENSION))
        if not line[0] in query_count :
            query_count[line[0]] = 1
        else :
            query_count[line[0]] += 1

#Creation of a dataframe with extended coordinates
Blast_bed = pd.DataFrame({"scaffold" : scaffold, "start" : start, "end" : end, "query" : query, 
                          "strand" : strand})

#Filter transposon based on number of hits
for key,value in query_count.items() :
    if value < int(MIN_HITS) :
        indexNames = Blast_bed[ Blast_bed['query'] == key ].index
        Blast_bed.drop(indexNames , inplace=True)    

#Saving filtered dataframe in a bed file
Blast_bed.to_csv("%s_BlastExtendend.bed" %OUT, sep="\t", index=False, header = False)

print("----> After filtering you are analyzing :", Blast_bed['query'].nunique(), "consensus with more than :", 
      int(MIN_HITS) -1, "hits on the genome."  )
print("########DONE")

####################################################################################################
#-------------------------------------------------EXTRACT------------------------------------------#
####################################################################################################

print("Extracting filtered extended blast hits...")

#Extraction of extended hits with bedtools getfasta. Hits on the opposite strand are reversed complemented
with open("%s_BlastExtendend.fasta" %OUT, 'w') as extracted_faa :
    subprocess.run(["bedtools", "getfasta", "-name", "-s","-fi", 
                    GENOME , "-bed" , "%s_BlastExtendend.bed" %OUT], 
                   stdout=extracted_faa, stderr=subprocess.DEVNULL)

#Parsing each transposon in a different multifasta and add original consensus    
os.mkdir("Extracted_Sequences")

TransposonsID = [] 

for record in SeqIO.parse("%s_BlastExtendend.fasta" %OUT, "fasta"):
    ID = record.id.split("#")
    TransposonsID.append(ID[0])
TransposonsID = set(TransposonsID)

for element in TransposonsID :
    Sequences = []
    for record in SeqIO.parse("%s_BlastExtendend.fasta" %OUT, "fasta"):
        if record.id.split("#")[0] == element :
            #Susbtituting : with _ in output fasta files for downstream analyses
            record.id = record.id.replace(":", "_")
            record.description = record.description.replace(":", "_")
            Sequences.append(record)
    for record2 in SeqIO.parse("%s.fasta" %output_name, "fasta") :
        if record2.id.split("#")[0] == element :
               Sequences.append(record2)
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
        if MAFFT_OPT == 'einsi' :
            mafft_cline.genafpair = True
            mafft_cline.maxiterate = 1000
            mafft_cline.adjustdirection = True
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
            

print("###############################################")
print("#EVERYTHING DONE, HAVE A NICE MANUAL CURATION!#")
print("###############################################")

