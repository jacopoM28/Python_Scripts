#######################################################################################################
#                             From Blast to GFF3                                                      #
#The script will convert a Blast output to a gff3, by default it will use blast 6  out format,        #
#but each GFF3 field can be filled based on custom blast output.                                      #
#User can use both query and subject coordinates, just pay attention to select the right columns      #                                                           
####################################################################################################### 
#!/usr/bin/env python

import pandas as pd
import argparse
import os

##############################################################################
#------------------Parsing arguments and setting variables-------------------#
##############################################################################

parser = argparse.ArgumentParser(description='Convert blast output to gff3.')
parser.add_argument('--in_file', required=True, help='Input blast output.')
parser.add_argument('--query_field', nargs="?", help='Query or subject field number',const="1")
parser.add_argument('--source', required=True, help='Hitted database')
parser.add_argument('--start_coord', nargs="?", help='Start of alignment in subject or query field number',const="7")
parser.add_argument('--end_coord', nargs="?", help='End of alignment in query or subject field number',const="8")
parser.add_argument('--score', nargs="?", help='Score field',const="11")
parser.add_argument('--attribute', nargs="?", help='Attributes field numbers. It can accept multiple numbers separetd by a comma',const="2")
parser.add_argument('--out', required=True, help='Output file name.')

args = parser.parse_args()

IN_FILE = args.in_file
QUERY = int(args.query_field) - 1
SOURCE = args.source
START_COORD = int(args.start_coord) - 1
END_COORD = int(args.end_coord) - 1
SCORE = int(args.score) - 1
ATTRIBUTE = args.attribute
#Checking number of field for the Attribute field
if "," in ATTRIBUTE :
    ATTRIBUTE = ATTRIBUTE.split(',')
    ATTRIBUTE = [int(i) for i in ATTRIBUTE]
    ATTRIBUTE = [ i - 1 for i in ATTRIBUTE]
else :
    ATTRIBUTE = int(ATTRIBUTE) - 1
OUT = args.out

#########################################################################
#-------------------------------MAIN------------------------------------#
#########################################################################

queries = []
source = []
feature = []
feature_start = []
feature_end = []
score = []
strand = []
phase = []
attribute = []
with open(IN_FILE, 'r') as file:
    for l in file :
        sl = l.split('\t')
        queries.append(sl[QUERY])
        source.append(SOURCE)
        feature.append("Blast_hit")
        #Check the strand of the hit, if starting coordinates are smaller than ending use + strand
        if int(sl[START_COORD]) < int(sl[END_COORD]) :
            feature_start.append(sl[START_COORD])
            feature_end.append(sl[END_COORD])
            strand.append("+")
        #If ending coordinates are smaller than starting ones, reverse them in the output and specify - strand
        else :
            feature_start.append(sl[END_COORD])
            feature_end.append(sl[START_COORD])
            strand.append("-")
        score.append(sl[SCORE])
        phase.append(".")
        #joining all attributes values
        if ATTRIBUTE > 1 :
            tmp = []
            for i in ATTRIBUTE :
                tmp.append(sl[i])
            tmp = ';'.join(tmp)
            attribute.append(tmp)
        else :
            attribute.append(sl[ATTRIBUTE])
     
gff3_out = pd.DataFrame(
    {'#query': queries,
     '#source': source,
     '#feature': feature,
     '#feature_start':feature_start,
     '#feature_end':feature_end,
     '#score':score,
     '#strand':strand,
     '#phase':phase,
     '#attribute':attribute
    })
with open(OUT, 'a') as f:
    f.write('#gff3\n')
    gff3_out.to_csv(f,index=False,sep="\t")