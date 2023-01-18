#!/usr/bin/env python
# coding: utf-8

# In[93]:


from pyfaidx import Fasta
import argparse

parser = argparse.ArgumentParser(description='Fix vcf entries files adding padding bases. Usefull for SNIFFLES insertions.')
parser.add_argument('--vcf', required=True, help='input vcf')
parser.add_argument('--genome', required=True, help='input genome')

args = parser.parse_args()
GENOME =  args.genome
VCF = args.vcf

genome = Fasta(GENOME)
with open(VCF) as f:
    for item in f:
        if item[0] == '#':
            line = item.rstrip()
            print(line)
        else :
            line = item.rstrip().split('\t')
            svtype = line[7].split(";")[1].split("=")[1]
            chrom = line[0]
            pos = int(line[1])
            ref = line[3]
            alt = line[4]
            if ref[0] == alt[0] :
                print('\t'.join(line))
            else :
                if svtype == "INS" :
                    nuc = genome[str(chrom)][pos - 1]
                    line[3] = nuc.seq
                    line[4] = nuc.seq + line[4]
                    print('\t'.join(line))


# In[92]:


print("#SNIFFLES_DEL",genome["NC_047559.1"][28813:28818]) #SNIFFLES DEL
print("#BSV_INS", genome["NC_047559.1"][34141:34148]) #PBSV INS
print("#SNIFFLES_INS", genome["NC_047559.1"][577839:577845]) #SNIFFLES INS
print("#PBSV_DEL", genome["NC_047559.1"][709547:709560]) #PBSV DEL

#SNIFFLES DEL == Correct

