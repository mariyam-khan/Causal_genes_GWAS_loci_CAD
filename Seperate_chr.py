from numpy import *
import numpy as np
import pandas as pd
import sys
from collections import Counter
import os

"""


You can give as input 
 
1. 1st argument, exposure data which you used to extract the GWAS summary data from 
the MR-Base package 
2. 2nd argument GWAS summary data 
3. 3rd argument integer value (usually 2)


This function will take the exposure and output files, 
segregate both files per chromosome and additionally per position.


You will have multiple files saved in your initial directory, 
for each chromosome and for different positions of the variants on the chromosome. 

As an example, with argument 2, you will have all SNPs on a particular chromosome  
on positions 13________, in one file and on the same chromosome, position 91________ in 
a different file (into batches of SNPs sharing only the first two digits in their base pair position). 
The input and outcome would both have comma as a separator for the .csv files. Using this function 
you can approximately segregate data and then manually check for exceptions. 


To use this function, make sure your output and exposure data are in the format you need for the MRBase package.


you can give the exposure data as the 1st argument, the outcome data as the 2nd argument and 
the position (integer value) as the 3rd argument. This function will take the exposure and output file, 
then separate the file per chromosome and additionally per position. As an example
if you give the position argument to be 2, the code will segregate each chromosome 
into batches of SNPs sharing only the first two values in their base pair position.


"""

file_EX = sys.argv[1]
file_EY = sys.argv[2]
position = sys.argv[3]
path = os.path.dirname(file_EX)
name_EX = os.path.splitext(os.path.basename(file_EX))[0]
name_EY = os.path.splitext(os.path.basename(file_EY))[0]
Data_out_all = pd.read_csv(file_EY, sep=',')
Data_exp = pd.read_csv(file_EX, sep=',')

Data_out_all = Data_out_all[["SNP", "effect_allele.outcome", "other_allele.outcome", "beta.outcome", "pval.outcome", "chr", "pos"]]
unique_chr = Counter(Data_out_all['chr']).keys()
no_chr = len(unique_chr)
get_chr = Counter(Data_out_all['chr']).most_common(no_chr)

for i in range(no_chr):
    Dataframe = (Data_out_all.loc[Data_out_all['chr'] == get_chr[i][0]])
    Dataframe.reset_index(drop=True, inplace=True)
    Dataframe = Dataframe[Dataframe['pval.outcome'] <= 5E-8]
    Dataframe.reset_index(drop=True, inplace=True)

    pd.options.mode.chained_assignment = None
    Dataframe['pos'] = Dataframe['pos'].astype('string')
    pd.options.mode.chained_assignment = None
    Dataframe['pos1'] = Dataframe['pos'].str[0:int(position)]
    unique_pos = Counter(Dataframe['pos1']).keys()
    no_pos = len(unique_pos)
    get_pos = Counter(Dataframe['pos1']).most_common(no_pos)
    for k in range(no_pos):
        Data_out = (Dataframe.loc[Dataframe['pos1'] == get_pos[k][0]])
        Data_out.reset_index(drop=True, inplace=True)
        Data_exp = Data_exp.loc[Data_exp['SNP'].isin(np.intersect1d(Data_exp.SNP, Data_out.SNP))]
        Data_out = Data_out.loc[Data_out['SNP'].isin(np.intersect1d(Data_out.SNP, Data_exp.SNP))]
        Data_exp.reset_index(drop=True, inplace=True)
        Data_out.reset_index(drop=True, inplace=True)
        if isinstance(Data_out, pd.DataFrame):
            if not Data_out.empty:
                Data_exp = Data_exp.loc[Data_exp['SNP'].isin(np.intersect1d(Data_exp.SNP, Data_out.SNP))]
                Data_out = Data_out.loc[Data_out['SNP'].isin(np.intersect1d(Data_out.SNP, Data_exp.SNP))]
                Data_exp.reset_index(drop=True, inplace=True)
                Data_out.reset_index(drop=True, inplace=True)
                Data_out = Data_out[
                    ["SNP", "effect_allele.outcome", "other_allele.outcome", "beta.outcome",
                     "pval.outcome", "chr", "pos"]]
                Data_exp.to_csv(path + "/" + name_EX + "_" + str(get_chr[i][0]) + "_" + str(get_pos[k][0]) + ".csv", sep=",",
                                index=False)
                Data_out.to_csv(path + "/" + name_EY + "_" + str(get_chr[i][0]) + "_" + str(get_pos[k][0]) + ".csv", sep=",",
                                index=False)
            else:
                sys.exit('Error Message : Dataset empty.')
        else:
            sys.exit('Error Message : Dataset empty.')

