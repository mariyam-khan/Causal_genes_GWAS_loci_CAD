from scipy.linalg import lstsq
from numpy import *
import numpy as np
import pandas as pd
import sys
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
pandas2ri.activate()

file_EXEY = sys.argv[1]
Data = pd.read_csv(file_EXEY, sep=',')
for (columnName, columnData) in Data.iteritems():
    if (Data[columnName] == 0).all():
        Data.drop(columnName, axis=1, inplace=True)
snps = Data.loc[:, 'SNPs'].values
if len(snps) != 1:
    cov = r("TwoSampleMR::ld_matrix")(snps)
    cov_data = pd.DataFrame(cov, columns=snps, index=snps, dtype='float')
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
    df = Data.loc[:, ~Data.columns.isin(['SNPs', 'CAD'])]
    gene_names = column_names = list(df.columns.values)
    covEY = Data.loc[:, 'CAD'].values
    covEX = df.values
    if no_snps == no_genes:
        if no_snps == 1:
            b_est = covEY / covEX
            d1 = {'gene': gene_names,
                  'Causal Estimate Ratio method': b_est}
            df1 = pd.DataFrame(data=d1)
        else:
            b_est = np.linalg.solve(covEX, covEY)
            d1 = {'gene': gene_names,
                  'Causal Estimate': b_est}
            df1 = pd.DataFrame(data=d1)
    else:
        b_est2, res, rnk, sy = lstsq(covEX, covEY)
        b_est1 = np.linalg.inv(covEX.T @ np.linalg.inv(cov_EE) @ covEX) @ (covEX.T @ np.linalg.inv(cov_EE) @ covEY)
        d1 = {'gene': gene_names,
              'Causal Estimate Least Squares': b_est2, 'Causal Estimate GMM': b_est1}
        df1 = pd.DataFrame(data=d1)
    df1 = df1.set_index('gene')
else:
    df = Data.loc[:, ~Data.columns.isin(['SNPs', 'CAD'])]

    gene_names = column_names = list(df.columns.values)
    covEY = Data.loc[:, 'CAD'].values
    covEX = df.values
    b_est = covEY / covEX
    d1 = {'gene': gene_names,
          'Causal Estimate Ratio method': b_est[0]}
    df1 = pd.DataFrame(data=d1)

df1.to_csv(file_EXEY + "_results.csv", sep=",", float_format='%g')