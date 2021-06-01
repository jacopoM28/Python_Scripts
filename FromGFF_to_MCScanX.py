#!/usr/bin/env python

######################################################################
#           GFF and Gene Family parser for MCScanX                   #
######################################################################
#Usefull script to prepare the GFF file necessary for MCScanX from a #
#standard NCBI GFF3 file.                                            #
######################################################################
 
import argparse
import os
import pandas as pd

#-------------------------Importing arguments------------------------#

parser = argparse.ArgumentParser(description='Python script usefull to reformat an NCBI gff3 into a gff file suitable for MCScanX')
#GFF3 input file
parser.add_argument('--in_gff', required=True, help='The gff file.')
#Prefix for output file
parser.add_argument('--in_chr2acc', required=True, help='Mapping file to rename chromosomes')
#Output prefix
parser.add_argument('--out', required=True, help='prefix for output files.')

args = parser.parse_args()
IN_GFF = args.in_gff
IN_MAPPING = args.in_chr2acc
OUT = args.out

#------------------------------MAIN---------------------------------#

#Open mapping file and store mapping information in a dictionary
mapping_dict = {}
with open(IN_MAPPING) as mapping_file:
    for line in mapping_file:
        if not line.startswith('#') :
            value = line.split("\t")[0]
            key = line.split("\t")[1].rstrip()
            mapping_dict[key] = value

Chromosomes = []
gene_ids = []
start_coords = []
ending_coords = []


#Opening original gff file...
with open(IN_GFF) as gcf:
    first_line=gcf.readline()
    #Checking the first line of the file...
    if not first_line == "##gff-version 3\n" :
        print("ERROR : The input GFF is not a standard NCBI GFF3")
        print("EXITING...")
    #Parsing the gff file
    else :
        for line in gcf:
            if not line.startswith('#') :
                if line.split("\t")[2] == "gene" :
                    chromosome = line.split("\t")[0]
                    start_gene = line.split("\t")[3]
                    end_gene = line.split("\t")[4]
                    id_gene = line.split("\t")[8].split(";")[0].replace("ID","gene").replace("gene-","")
                    gene_type = line.split("\t")[8].split(";")[5].split("=")[1].rstrip()
                    #Checking that the gene is a protein coding gene...
                    if gene_type == "protein_coding" :
                    #Checking that the the annotation is related to a chromosome and not to an unplaced scaffold...
                    	if chromosome in mapping_dict :
                        	Chromosomes.append(chromosome)
                        	gene_ids.append(id_gene)
                        	start_coords.append(start_gene)
                        	ending_coords.append(end_gene)
                
#Check that everything has worked...
length = len(Chromosomes)
if not all(len(lst) == length for lst in [Chromosomes, gene_ids, start_coords, ending_coords]) :
    print("Something has gone wrong....Errors occured while formatting the GFF.")
    print("Ensure that is a standard NCBI GFF3")
else :
    #Creating the reformatted gff...
    reformatted_GFF = pd.DataFrame(
    {"#chr":Chromosomes,
     "#gene_id":gene_ids,
     "#start_coord":start_coords,
     "#ending_coord":ending_coords
    })

#Changing the name of the chromosomes based on the mapping file...
reformatted_GFF = reformatted_GFF.replace({"#chr": mapping_dict})

#---------------------Printing out the resulting file--------------------#

reformatted_GFF.to_csv(path_or_buf=open('%s.gff' % OUT, 'w'),sep="\t", header = False,index=False)

