#!/usr/bin/env python
# coding: utf-8

# In[59]:


import pandas as pd

inFile = "PTHR10071_PfamProsite_Hits.txt"
Family = inFile.split("_")[0]
Specie_count={Family:{}} 

with open("PTHR10071_PfamProsite_Hits.txt") as file:
    for line in file:
        Specie = line.split("|")[0]
        if not Specie in Specie_count[Family] :
            Specie_count[Family][Specie] = int(1)
        else :
            Specie_count[Family][Specie] = Specie_count[Family][Specie] + 1
pd.DataFrame.from_dict(Specie_count,orient='index')


# In[17]:


########################################################################################################
#  Convert a set of txt files in which each line correspond to a sequence in a table                   #
#  in which each line correspond to an original txt file and each column to the number of              #
#  sequences present for a specif specie in the correspondig file. The table is the equivalent of      #
#  the Orthofinder Gene.count.tsv output table.                                                        #
#                                                                                                      #
# NB: It's a very specifi scripts, not generalized. Each sequence must have the specie name at the     #
#   beggining of the line, separated by the rest from a pipe. E.g :                                    #
#                                                                                                      #
#                   H.robusta|KB097620.1_cds_ESN93303.1_19906_1|locus_tag=HELRODRAFT_89090             #
#                                                                                                      #
#  Moreover shouldn't be any other .txt files in any  downstream folders (I know it's really terrible!)#
########################################################################################################

import pandas as pd
import os
from pathlib import Path

#Creation of an empty dictionary
Specie_count={}

#Starting iterating through the files...
for root, dirs, files in os.walk(os.getcwd()):
    for file in files :
        filename, extension = os.path.splitext(file)
        
        #if files extension is ".txt" open the file...
        if extension == '.txt' :
            with open(file) as file_2 :
                Family = filename.split("_")[0] 
                i = 0
                for line in file_2:
                    Specie = line.split("|")[0]
                    #If it's the first line create the first key value for the first appearing specie...
                    if i == 0 :
                        Specie_count[Family] = {Specie:int(1)}
                    
                    #else check if Specie key already exist...
                    else :
                        if not Specie in Specie_count[Family] :
                            Specie_count[Family][Specie] = int(1)
                        else :
                            Specie_count[Family][Specie] = Specie_count[Family][Specie] + 1
                    i += 1

                    #Create a table from the ensted dictionary
df = pd.DataFrame.from_dict(Specie_count,orient='index')
#Substitute NaN values with 0s
df = df.fillna(0)
#Remove decimal 0s
df = df.astype(int)
print(df)
df.to_csv("Gene.count.tsv",sep =" ")

