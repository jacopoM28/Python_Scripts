#!/usr/bin/env python
# coding: utf-8

# In[162]:


##################################################################################################
##################################################################################################
##                                                                                              ## 
##                       Extract ubiqutous single copy busco genes                              ##
##                                                                                              ##
##################################################################################################
##################################################################################################
#The script must be run in a directory in which are stored the output dir of Busco v.5
#IMPORTANT: The directories bust be named ${specie_name}_*_Busco. In lines 20 and 30  you can change the
#name of the reference dataset (depending on the name the Busco output directory looks different)

from Bio import SeqIO
import pandas as pd
import os, glob

os.mkdir("Ubiqutous_SingleCopy_Genes")
Sequences = []
for directory in os.listdir() :
    if directory.endswith("_Busco") :
        tsv_path = directory + '/run_metazoa_odb10/full_table.tsv'
        df = pd.read_csv(tsv_path,comment='#',header=None,sep="\t")
        Single_Copy = list(df.loc[df[1] == 'Complete'][0])
        Sequences.append(Single_Copy)
Ubiqoutous_Genes_List = list(set.intersection(*map(set, Sequences)))

Ubiqoutous_Genes_Fasta = {}
for directory in os.listdir() :
    if directory.endswith("_Busco") :
        specie_name = directory.split("_")[0]   
        fasta_fold = directory + "/run_metazoa_odb10/busco_sequences/single_copy_busco_sequences"
        fasta_paths = glob.glob(os.path.join(fasta_fold, '*.faa'))
        for fasta_path in fasta_paths :
            gene = fasta_path.split("/")[-1].split(".")[0]
            for seq_record in SeqIO.parse(fasta_path, "fasta"):
                if gene in Ubiqoutous_Genes :
                    seq_record.id = specie_name
                    seq_record.description = ""
                    seq_record.name = ""
                    if gene in Ubiqoutous_Genes_Fasta :
                        Ubiqoutous_Genes_Fasta[gene].append(seq_record)
                    else :
                        Ubiqoutous_Genes_Fasta[gene] = [seq_record]

for key, value in Ubiqoutous_Genes_Fasta.items() :
    SeqIO.write(value,"Ubiqutous_SingleCopy_Genes/%s.fasta" %key, "fasta")
    

