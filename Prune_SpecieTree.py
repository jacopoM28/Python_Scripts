#!/usr/bin/env python

##################################################################################################################
#                             Prune Specie tree based on gene family occurency                                   #
#The script take as input a directory supposing to contain a gene family with one fasta file for each            #
#specie. As for Orthofinder, the name of the file must be: 'SpecieName.{pep,fa,fasta}' (e.g R.philippinarum.pep) #
#Species names must be in completely agreement with specie tree tip labels. For each gene family the script      #
#will prune the specie tree based on fasta file names.                                                           #
#Usefull for automatize within gene families Orthofinder runs with user supplied specie tree                     #                                               
##################################################################################################################

#Libraries needed : ETE3

import shutil
from ete3 import Tree
import pandas as pd
import os
from pathlib import Path
import sys
import argparse

parser = argparse.ArgumentParser(description='Usefull script to automatically prune a species tree based on a set of fasta files, where each one represent one OTU.')
parser.add_argument('--in_dir', required=True, help='Input directory where to screen for fasta files. Files must ended with .pep, .fa or .fasta.')
parser.add_argument('--in_tree', required=True, help='Input specie tree.')
parser.add_argument('--out', required=True, help='output prefix.')
args = parser.parse_args()

DIR =  args.in_dir
TREE = args.in_tree
OUT = args.out

Species_Proteins = []

#Reading protein files and storing informationa bout species...

for root, dirs, files in os.walk(DIR):
    for file in files :
        filename, extension = os.path.splitext(file)
        if extension == '.pep' or extension == '.fa' or extension == ".fasta" :
            Species_Proteins.append(filename)
         
OTUs = []

#Reading tree file and storing leaf informations...

t = Tree(TREE)
for node in t.traverse() :
    if node.is_leaf() :
        OTUs.append(node.name)

Data_not_Tree = list(set(Species_Proteins).difference(OTUs))
Tree_not_Data = list(set(OTUs).difference(Species_Proteins))
To_Keep = set(OTUs).intersection(Species_Proteins)

#Checking that al species for which we have proteins are present in the specie tree...

if len(Data_not_Tree) > 0 :
    print("\n")
    print("ERROR : Some species for which you have proteins are missing from the specie tree", file = sys.stderr)
    print("Please check it!", file = sys.stderr)
    sys.exit()
if len(Tree_not_Data) == 0 :
    print("You have proteins for all yours species...exiting")
    shutil.copy(TREE, DIR)
else :
#Pruning the tree...
    print("in ", "your tree", "is/are missing the following specie/s :", Tree_not_Data)
    print("They will be pruned from the specie tree...")

    t.prune(To_Keep)
    t.write(format=1, outfile="%s.nwk" %OUT)

