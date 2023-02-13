from scipy.linalg import lstsq
from numpy import *
import numpy as np
import pandas as pd
import sys


"""
getting the causal genes for cases the user can supply the LD-matrix.


To run the fil, type in the command line python3 run_MVMR.py "/home/user/exposure_outcome.csv" "/home/user/ld.csv"


where: 

"/home/user/exposure_outcome.csv"

sys.argv[1] is the first argument when you run the file and should be the exposure_outcome.csv file 
containing the SNPs to exposure effects and SNPs to outcome effects.


"/home/user/ld.csv"

sys.argv[2] is the second argument i.e. ld.csv file for the LD matrix.
Please make sure the ordering of the SNPs is same as in the exposure_outcome file


The output after running this file would be a .csv with results from the methods (depending
 on the dimensions, least-squares, generalized method of moments and ratio method)
saved as .csv file in the same directory as the one given for the exposure_outcome.csv file. 

i.e. for this example "/home/user/exposure_outcome_results.csv.csv"

"""

file_EXEY = sys.argv[1]
file_EE = sys.argv[2]
Data = pd.read_csv(file_EXEY, sep=',')
cov_EE = np.genfromtxt(file_EE, delimiter=' ')
for (columnName, columnData) in Data.iteritems():
    if (Data[columnName] == 0).all():
        Data.drop(columnName, axis=1, inplace=True)
snps = Data.loc[:, 'SNPs'].values
if len(snps) != 1:
    cov_data = pd.DataFrame(cov_EE, columns=snps, index=snps, dtype='float')
    upper_tri = cov_data.where(np.triu(np.ones(cov_data.shape), k=1).astype(bool))
    to_drop = [column for column in upper_tri.columns if any(upper_tri[column] == 1.0)]
    cov_data = cov_data.drop(labels=to_drop, axis=1)
    cov_data = cov_data.drop(labels=to_drop, axis=0)
    Data = Data[Data.SNPs.isin(cov_data.index.values)]
    Data.reset_index(drop=True, inplace=True)
    cov_EE = cov_data.values
    snps_new = Data.loc[:, 'SNPs'].values
    no_snps = len(snps_new)
    no_genes = len(Data.columns) - 2
    df = Data.loc[:, ~Data.columns.isin(['SNPs', 'outcome'])]
    gene_names = column_names = list(df.columns.values)
    covEY = Data.loc[:, 'outcome'].values
    covEX = df.values
    if no_snps == no_genes:
        if no_snps == 1:
            b_est = covEY / covEX
            d1 = {'gene': gene_names,
                  'Causal Estimate Ratio method': b_est}
            df1 = pd.DataFrame(data=d1)
            df1 = df1.set_index('gene')
        else:
            b_est = np.linalg.solve(covEX, covEY)
            d1 = {'gene': gene_names,
                  'Causal Estimate': b_est}
            df1 = pd.DataFrame(data=d1)
            df1 = df1.set_index('gene')
    else:
        if no_snps > no_genes:
            b_est2, res, rnk, sy = lstsq(covEX, covEY)
            b_est1 = np.linalg.inv(covEX.T @ np.linalg.inv(cov_EE) @ covEX) @ (covEX.T @ np.linalg.inv(cov_EE) @ covEY)
            d1 = {'gene': gene_names,
                  'Causal Estimate Least Squares': b_est2, 'Causal Estimate GMM': b_est1}
            df1 = pd.DataFrame(data=d1)
            df1 = df1.set_index('gene')
        else:
            sys.exit('Error Message : You require at least as many instruments as exposures to run this analysis.')

else:
    no_genes = len(Data.columns) - 2
    if len(snps) >= no_genes:
        df = Data.loc[:, ~Data.columns.isin(['SNPs', 'outcome'])]
        gene_names = column_names = list(df.columns.values)
        covEY = Data.loc[:, 'outcome'].values
        covEX = df.values
        b_est = covEY / covEX
        d1 = {'gene': gene_names,
              'Causal Estimate Ratio method': b_est[0]}
        df1 = pd.DataFrame(data=d1)
        df1 = df1.set_index('gene')
    else:
        sys.exit('Error Message : You require at least as many instruments as exposures to run this analysis.')
df1.to_csv(file_EXEY + "_results.csv", sep=",", float_format='%g')

