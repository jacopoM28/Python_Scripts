#!/usr/bin/env python

############################################################
#   Parser for create MCScanX gene family file             #
############################################################
 
# Simple script usefull to create a gene family file for MCScanX. Multiple gene family files can be appended in order to analyze all families together. 
#NB: The input gene family file must have one gene per line and the gene name must be at the end of the line separated by a pipe character ("|")!

#Input family file example :
#C.gigas|CDS_000111222.3_XP_9992221|gene=LOC1111678
#C.gigas|CDS_050112123.3_XP_9234781|gene=LOC1111679

#Output file :
#C1qDC	gene=LOC1111678	gene=LOC1111679

import argparse
import os

parser = argparse.ArgumentParser(description='Python script usefull to reformat a gene family file - one line per gene - to a format suitable for MCScanX - all genes in one line.')
#Gene family input file
parser.add_argument('--in_file', required=True, help='The input gene family file.')
#Gene Family name
parser.add_argument('--name', required=True, help='Name of the gene family. It will be the firs element of the line')
#Output prefix
parser.add_argument('--out', required=True, help='prefix for output files.')

args = parser.parse_args()
IN_FILE = args.in_file
NAME = args.name
OUT = args.out

genes = [NAME]
with open(IN_FILE) as file:
    for line in file :
        gene = line.split("|")[2].rstrip()
        genes.append(gene)

with open('%s.MCScanX_GeneFamily' % OUT, 'w') as file:
    file.write('\t'.join(genes))

