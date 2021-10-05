#!/usr/bin/env python
# coding: utf-8

# In[175]:


####################################################################################
####################################################################################
#                                                                                  #
#                       REFORMAT REPEATMASKER DIVERGENCE                           #
#                                                                                  # 
####################################################################################
####################################################################################

#Author : Jacopo Martelossi
#E-mail : jacopo.martelossi2@unibo.it


#This script collapse the RepeatMasker divergence table at Class level. Genome size of the organism must
#be provided by the user. You can directly parse the resulting table to R for nice plotting

#RUNNING EXAMPLE

#1. perl calcDivergenceFromAlign.pl -s <output.div> <input.alig>
#2. tail -n 72 <div_file> > <Output_Distance_Table>
#3. python RepLandscape_Class.py  --div_file <Distance_Table> --genome_size <Genome_size> --out <Output_Name>

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Reformat repeat landascape to class level')
parser.add_argument('--div_file', required=True, help='Input file obtained from the last 72 lines of the div file obtained with the CalcDiv_fromAlign.pl script')
parser.add_argument('--genome_size', required=True, help='genome size in bp')
parser.add_argument('--out', required=True, help='Output file name.')

args = parser.parse_args()

IN_FILE = args.div_file
GENOME_SIZE = int(args.genome_size)
OUT = args.out

df_New = pd.DataFrame(columns = ["div", "DNA","LINE", "SINE", "RC", "MITE", "LTR","Unknown", "Satellite"])
df = pd.read_csv(IN_FILE,sep=" ", header=0)

for index, row in df.iterrows() :
    df_New.loc[index, "div"] = row[0]
    df_New = df_New.fillna(0)
    for cell in df.columns :
        if cell.split("/")[0] == "DNA" :
            df_New.loc[index,"DNA" ] = df_New["DNA"][index] + float(df[cell][index])
        elif cell.split("/")[0] == "LINE" :
            df_New.loc[index,"LINE" ] = df_New["LINE"][index] + float(df[cell][index])
        elif cell.split("/")[0] == "SINE" :
            df_New.loc[index,"SINE" ] = df_New["SINE"][index] + float(df[cell][index])
        elif cell.split("/")[0] == "RC" :
            df_New.loc[index,"RC" ] = df_New["RC"][index] + float(df[cell][index])
        elif cell.split("/")[0] == "MITE" :
            df_New.loc[index,"MITE" ] = df_New["MITE"][index] + float(df[cell][index])
        elif cell.split("/")[0] == "LTR" :
            df_New.loc[index,"LTR" ] = df_New["LTR"][index] + float(df[cell][index])
        elif cell.split("/")[0] == "Unknown" :
            df_New.loc[index,"Unknown" ] = df_New["Unknown"][index] + float(df[cell][index])
        elif cell.split("/")[0] == "Satellite" :
            df_New.loc[index,"Satellite" ] = df_New["Satellite"][index] + float(df[cell][index])

df_New = df_New/GENOME_SIZE * 100
df_New["div"] = range(1,len(df_New)+1)
df_New = df_New[0:51]
df_New.to_csv("%s_RepLandscape_Class.tsv" %OUT, sep ="\t", index=False)

