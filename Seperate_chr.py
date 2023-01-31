from numpy import *
import numpy as np
import pandas as pd
import sys
from collections import Counter

"""

you can give the exposure data as the 1st argument, the outcome data as the 2nd argument and 
the position (integer value) as the 3rd argument. This function will take the exposure and output file, 
then separate the file per chromosome and additionally per position. As an example
if you give the position argument to be 2, the code will segregate each chromosome 
into batches of SNPs sharing only the first two values in their base pair position.

"""

file_EX = sys.argv[1]
file_EY = sys.argv[2]
position = sys.argv[3]
Data_out_all = pd.read_csv(file_EY, sep=',')
Data_exp = pd.read_csv(file_EX, sep=',')

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
    Dataframe['pos1'] = Dataframe['pos'].str[0:position]
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
                Data_exp.to_csv(file_EX + "_" + str(get_chr[i][0]) + "_" + str(get_pos[k][0]) + ".csv", sep=",",
                                index=False)
                Data_out.to_csv(file_EY + "_" + str(get_chr[i][0]) + "_" + str(get_pos[k][0]) + ".csv", sep=",",
                                index=False)
