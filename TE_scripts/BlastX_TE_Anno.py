#!/usr/bin/env python
# coding: utf-8

# In[8]:


#!/usr/bin/env python
# coding: utf-8

# In[ ]:

#################################################################################
#################################################################################
##                                                                             ##   
##                       Blastx based TE classification                        ##
##                                                                             ##
#################################################################################
#################################################################################

#This script will blast all sequences of a TE library against a genome and will search
#for copies with similarities to TE-Related proteins. By default an e-value of 1e-10
#is used for blastn and an e-value of 1e-05 for Blastx. Moreover only 20 insertions longer than 
#100 bp will be screened. If all significant Blastx are against
#elements with the same annotation, the annotation is simply transfered to the consensus.
#If not, the script will cechk if they belong at least at the same TE class.


import os
import sys
import csv
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastxCommandline

parser = argparse.ArgumentParser(description='Get blastx of the 20 best insertions of each sequence present in the library.')
parser.add_argument('--lib', required=True, help='Input TE library')
parser.add_argument('--genome', required=True, help='input genome')
parser.add_argument('--num_threads',required=True, help="Number of cores for blastx")
parser.add_argument('--out', required=True, help='prefix output file')
parser.add_argument('--anno', action='store_true', default=False)
parser.add_argument('--TE_prot',required=True, help="TE derived protein library")

args = parser.parse_args()

GENOME =  args.genome
LIBRARY = args.lib
N_threads = int(args.num_threads)
OUT = args.out
MODE = args.anno
TE_ProtLib = args.TE_prot

os.mkdir("intermediate_files")

##################################################################################################
#---------------------------Running Blastn and reformat output-----------------------------------#
##################################################################################################

db = NcbimakeblastdbCommandline(dbtype="nucl",
                                   input_file= GENOME)
stdout, stderr = db()

blastn_cline = NcbiblastnCommandline(query=LIBRARY, num_threads=N_threads,
                                     db=GENOME,
                                     evalue=1e-10, outfmt=6, out="./intermediate_files/Blast.out")
stdout, stderr = blastn_cline()

scaffold = []
start = []
end = []
name = []

df = pd.read_csv("./intermediate_files/Blast.out", sep="\t",header = None)
#Remove insertions shorter than 100bp
df.drop(df[df[3] < 100].index, inplace=True)
#Sort by name, alignment length and bitscore
df = df.sort_values([0, 3,11], ascending=[True,False,False ])
#Keep only first 20 hits
df = df.groupby(df[0]).head(20)
df.reset_index(drop=True, inplace=True)

for index, row in df.iterrows():
    name.append(row[0])
    scaffold.append(row[1])
    if row[8] < row [9] :
        start.append(row[8])
        end.append(row[9])
    else :
        start.append(row[9])
        end.append(row[8])

Blast_bed = pd.DataFrame({"scaffold" : scaffold, "start" : start, "end" : end, "name" : name})
Blast_bed.to_csv("./intermediate_files/Blast.bed", sep="\t", index=False, header = False)

###################################################################################################
#--------------------------------Extract hits and run Blastx--------------------------------------#
###################################################################################################

subprocess.run(["samtools", "faidx", GENOME],stderr=subprocess.DEVNULL)

with open("./intermediate_files/Blast.fasta", 'w') as extracted_faa :
    subprocess.run(["bedtools", "getfasta", "-nameOnly","-fi", 
                    GENOME , "-bed" , "./intermediate_files/Blast.bed"], 
                   stdout=extracted_faa, stderr=subprocess.DEVNULL)
    
db = NcbimakeblastdbCommandline(dbtype="prot",                           
                                   input_file= TE_ProtLib)

stdout, stderr = db()

blastn_cline = NcbiblastxCommandline(query = "./intermediate_files/Blast.fasta", num_threads = 2,
                                     db=TE_ProtLib, evalue=1e-5, outfmt=6, 
                                     out="./intermediate_files/Ins_Blastx.out",max_target_seqs = 1, max_hsps = 1)
stdout, stderr = blastn_cline()

##################################################################################################
#---------------------Automatic protein level classification if requested------------------------#
##################################################################################################

#------------------------------------------MAIN FUNCTION-----------------------------------------#

#Check annotation of a list of elements. If all annotations are the same the function will directly return 
#the annotation. If they are different it will check for concardant class-level classification

def ckeckList(lst):
  
    ele = lst[0]
    chk = True
    ClassList = []
    #If the element has only one hit just take that one as annotation
    if len(lst) == 1 :
        anno = ele
    else :
        # Comparing each element with first item...
        for item in lst:
            ClassList.append(item.split("/")[0])
            if ele != item:
                chk = False
                ClassList.append(item.split("/")[0])
        #If the second item is already different from the first one, try to check the annotation at the 
        #class-level
        if chk == False :
            ele = ClassList[0]
            for item in ClassList :
                #If also class-level classification is incosistent annotate it as an Unknown
                if ele != item :
                    anno = "Unkown"
                #If class-level classification is cosistent annotate it th the class level
                else : 
                    anno = ele
        #If all element of the list have the same annotation use it!
        else :
            anno = ele
    return anno

#---------------------------------------------END-----------------------------------------------#

#Old and new classification that will be used to classify consensus sequences
Old_Anno = []
New_Anno = []

df = pd.read_csv("./intermediate_files/Ins_Blastx.out", sep="\t",header = None)
#Add a last empty line to the dataframe. Necessary to run the function also in the last line of the original 
#file in same special cases.
df = df.append(pd.Series(), ignore_index = True)

element = df.iloc[0].tolist()[0]
#List on which the ckeckList function will be run.
Element_Anno = []

for index, row in df.iterrows() :
        if row[0] != element :
            anno = ckeckList(Element_Anno)
            Old_Anno.append(element)
            New_Anno.append(element.split("#")[0] + "#" + anno)
            #Check if it is the last, fake line of the file. If yes, stop
            if pd.isna(row[0]) == True :
                break
            else :
                element = row[0]
                Element_Anno = [row[1].split("#")[1]]
        else :
            Element_Anno.append(row[1].split("#")[1])

Classification = pd.DataFrame({"Old name": Old_Anno, "New name": New_Anno })  
Classification.to_csv("Consensus_Classification.tsv", sep="\t", index=False)

if MODE == True :
#Reclassify consensus sequences
    mapper = Classification.set_index("Old name").to_dict()["New name"]
    records = list(SeqIO.parse(LIBRARY, "fasta"))
    for i in records:
        i.id = mapper.get(i.id,i.id)
        i.description = mapper.get(i.description,i.description)
    with open("%s-classified.fa" %OUT, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

