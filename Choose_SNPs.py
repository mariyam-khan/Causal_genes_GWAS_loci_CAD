from numpy import *
import numpy as np
import pandas as pd
import sys
from collections import Counter
import os

"""
you can give the exposure data as the 1st argument, the outcome data as the 2nd argument and 
the position  of the SNP around which you wish to get SNPs within 1Mb distance (integer value) as the 3rd argument. 
This function will take the exposure and output file,  take the SNP with it's position and return SNPs within +- 100000
Kb distance then give .
"""

file_EX = sys.argv[1]
file_EY = sys.argv[2]
chromosome = int(sys.argv[3])
position = int(sys.argv[4])
path = os.path.dirname(file_EX)
name_EX = os.path.splitext(os.path.basename(file_EX))[0]
name_EY = os.path.splitext(os.path.basename(file_EY))[0]
Data_out_all = pd.read_csv(file_EY, sep=',')
Data_exp = pd.read_csv(file_EX, sep=',')

Data_out_all = Data_out_all[
    ["SNP", "effect_allele.outcome", "other_allele.outcome", "beta.outcome", "pval.outcome", "chr", "pos"]]
Data_out = (Data_out_all.loc[Data_out_all['chr'] == chromosome])
Data_out.reset_index(drop=True, inplace=True)
Data_out = Data_out[Data_out['pval.outcome'] <= 5E-8]
Data_out.reset_index(drop=True, inplace=True)
if isinstance(Data_out, pd.DataFrame):
    if not Data_out.empty:
        lower = int(position) - 50000
        upper = int(position) + 50000

        Data_out = Data_out[(Data_out['pos'] >= lower) & (Data_out['pos'] <= upper)]
        Data_out.reset_index(drop=True, inplace=True)
        Data_exp = Data_exp.loc[Data_exp['SNP'].isin(np.intersect1d(Data_exp.SNP, Data_out.SNP))]
        Data_out = Data_out.loc[Data_out['SNP'].isin(np.intersect1d(Data_out.SNP, Data_exp.SNP))]
        Data_exp.reset_index(drop=True, inplace=True)
        Data_out.reset_index(drop=True, inplace=True)
        if isinstance(Data_out, pd.DataFrame):
            if not Data_out.empty:
                Data_exp.to_csv(path + "/" + name_EX + "_" + str(chromosome) + "_" + str(position) + ".csv", sep=",",
                                index=False)
                Data_out.to_csv(path + "/" + name_EY + "_" + str(chromosome) + "_" + str(position) + ".csv", sep=",",
                                index=False)
            else:
                sys.exit('Error Message : Dataset is empty, no SNPs shared with the exposure data.')
        else:
            sys.exit('Error Message : Dataset is empty, no SNPs shared with the exposure data.')
else:
    sys.exit('Error Message : Dataset is empty.There were no SNPs that passed the p-value threshold.')