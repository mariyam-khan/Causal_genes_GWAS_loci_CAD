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
These two files should be the output of either Seperate_chr.py or Choose_SNPs.py, 
with comma as a separator and the output would be a file with the format required for the
causal analysis saved in the same input directory with suffix prepared.csv

"""

file_EX = sys.argv[1]
file_EY = sys.argv[2]
path = os.path.dirname(file_EX)

name_EX = os.path.splitext(os.path.basename(file_EX))[0]
name_EY = os.path.splitext(os.path.basename(file_EY))[0]

Data_out = pd.read_csv(file_EY, sep=',')
Data_exp = pd.read_csv(file_EX, sep=',')

Data_exp = Data_exp.loc[Data_exp['SNP'].isin(np.intersect1d(Data_exp.SNP, Data_out.SNP))]
Data_out = Data_out.loc[Data_out['SNP'].isin(np.intersect1d(Data_out.SNP, Data_exp.SNP))]
Data_exp.reset_index(drop=True, inplace=True)
Data_out.reset_index(drop=True, inplace=True)

unique_genes = Counter(Data_exp['exposure']).keys()
no_genes = len(unique_genes)
most_common = Counter(Data_exp['exposure']).most_common(no_genes)
selected_columns = Data_exp['effect_allele.exposure']
new_df = selected_columns.copy()
o = Data_out['SNP'].values
for m in range(len(Data_out['SNP'].values)):
    value = o[m]
    idx = Data_exp[Data_exp['SNP'] == value].index.values
    string_allele = new_df[idx[0]]
    string_allele = string_allele.replace("'", "")
    a = Data_out.loc[m, 'effect_allele.outcome']
    b = Data_out.loc[m, 'other_allele.outcome']
    if string_allele != Data_out.loc[m, 'effect_allele.outcome']:
        Data_out.loc[m, 'beta.outcome'] = - Data_out.loc[m, 'beta.outcome']
        Data_out.loc[m, 'effect_allele.outcome'] = b
        Data_out.loc[m, 'other_allele.outcome'] = a

        Data_out.loc[m, 'SNP'] = value + "_" + b + "_" + a
    else:
        Data_out.loc[m, 'SNP'] = value + "_" + a + "_" + b
genes = []
for j in range(no_genes):
    genes.append(most_common[j][0])
selected_columns2 = Data_out['beta.outcome']
new_df2 = selected_columns2.copy()
column_names = genes
Data_out = Data_out.sort_values(by='pval.outcome', ascending=True)
Data_out.reset_index(drop=True, inplace=True)
row_names = list(Data_out.SNP)
matrix = np.zeros((len(row_names), len(column_names)))
R_data = pd.DataFrame(matrix, columns=column_names, index=row_names, dtype='object')
s_data = pd.DataFrame(np.zeros(len(row_names)), columns=['outcome'], index=row_names, dtype='object')

for n in range(no_genes):
    name = genes[n]
    exp = (most_common[n][0])
    Data = Data_exp[Data_exp['exposure'] == exp]
    Data.reset_index(drop=True, inplace=True)
    o = Data['SNP'].values
    for j1 in range(int(most_common[n][1])):
        value = o[j1] + "_" + Data.loc[j1, 'effect_allele.exposure'] + "_" + Data.loc[j1, 'other_allele.exposure']
        R_data.at[value, name] = Data.loc[j1, 'beta.exposure']
        idx = Data_out[Data_out['SNP'] == value].index.values
        s_data.loc[value, 'outcome'] = new_df2[idx[0]]

R_data.replace(np.nan, 0)
s_data.replace(np.nan, 0)
df_R = pd.DataFrame(data=R_data)
df_s = pd.DataFrame(data=s_data)
RandS = pd.concat([df_R, df_s], axis=1)
RandS.index.name = 'SNPs'
RandS.to_csv(path + "/" + name_EX + "_prepared.csv", sep=",")
