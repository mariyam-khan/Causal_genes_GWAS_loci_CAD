from numpy import *
import numpy as np
import pandas as pd
import sys
import os
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

ro.r['options'](warn=-1)

base = importr('base')
base.warnings()

if len(sys.argv) < 5:
    print("You failed to either provide the minimum or maximum number of parameters required to run this code. Please "
          " read the documentation for more details")
    sys.exit(1)
else:
    file_EX1 = sys.argv[1]
    file_EY1 = sys.argv[2]
    path_original = os.path.dirname(file_EX1)
    name_EX = os.path.splitext(os.path.basename(file_EX1))[0]
    name_EY = os.path.splitext(os.path.basename(file_EY1))[0]
    if len(sys.argv) == 4:
        LD_threshold_lower1 = float(sys.argv[3])
        LD_threshold_upper1 = 1.0
        if not isinstance(LD_threshold_lower1, float):
            print(
                "Provide a valid value for the lower LD threshold. Please "
                " read the documentation for more details.")
            sys.exit(1)

        if not LD_threshold_lower1 <= 0.5 and LD_threshold_lower1 >= 0.0:
            print(
                "Provide threshold in the range, 0.0 >= LD_threshold_lower <= 0.5. Please "
                " read the documentation for more details.")
            sys.exit(1)
    elif len(sys.argv) == 5:
        LD_threshold_lower1 = float(sys.argv[3])
        LD_threshold_upper1 = float(sys.argv[4])
        if not (isinstance(LD_threshold_lower1, float) and isinstance(LD_threshold_lower1, float)):
            print(
                "Provide a valid value for the lower/upper LD threshold. Please "
                "read the documentation for more details.")
            sys.exit(1)

        if not (0.5 >= LD_threshold_lower1 >= 0.0) and (LD_threshold_upper1 <= 1.0 and
                                                        LD_threshold_lower1 >= 0.0):
            print(
                "Provide threshold in the range, 0.0 >= LD_threshold_lower <= 0.5 and 0.0 >= LD_threshold_lower <= 1.0."
                "Please"
                " read the documentation for more details.")
            sys.exit(1)

    else:
        LD_threshold_lower1 = 0.0
        LD_threshold_upper1 = 1.0


def SNPs_LDrange(file_EX, file_EY, LD_threshold_lower, LD_threshold_upper):
    Data_out = pd.read_csv(file_EY, sep=',')
    Data_exp = pd.read_csv(file_EX, sep=',')
    pvalue = 5E-8
    Data_out = Data_out[Data_out['pval.outcome'] <= pvalue]
    Data_out.reset_index(drop=True, inplace=True)
    if not Data_out.empty:
        snps = list(Data_out.loc[:, 'SNP'].values)
        if len(snps) != 1:
            cov_data = r("TwoSampleMR::ld_matrix")(snps, with_alleles=False, pop="EUR")

            pandas2ri.activate()
            pd_from_r_df = ro.conversion.rpy2py(cov_data)
            cov_data = pd.DataFrame(pd_from_r_df, columns=r('colnames')(cov_data),
                                    index=r('colnames')(cov_data), dtype='float')

            cov_data = pd.DataFrame.abs(cov_data)
            upper_tri = cov_data.where(np.triu(np.ones(cov_data.shape), k=1).astype(bool))
            to_drop = [column for column in upper_tri.columns if
                       any(upper_tri[column] >= int(LD_threshold_upper)) or any(upper_tri[column]
                                                                                <= int(
                           LD_threshold_lower))]
            cov_data = cov_data.drop(labels=to_drop, axis=1)
            cov_data = cov_data.drop(labels=to_drop, axis=0)

            Data_out = Data_out[Data_out['SNP'].isin(cov_data.index.values)]

            Data_exp1 = Data_exp.loc[
                Data_exp['SNP'].isin(np.intersect1d(Data_exp.SNP, Data_out.SNP))]
            Data_out1 = Data_out.loc[
                Data_out['SNP'].isin(np.intersect1d(Data_out.SNP, Data_exp.SNP))]
            Data_exp1.reset_index(drop=True, inplace=True)
            Data_out1.reset_index(drop=True, inplace=True)
            pandas2ri.deactivate()
        else:
            Data_exp1 = Data_exp.loc[
                Data_exp['SNP'].isin(np.intersect1d(Data_exp.SNP, Data_out.SNP))]
            Data_out1 = Data_out.loc[
                Data_out['SNP'].isin(np.intersect1d(Data_out.SNP, Data_exp.SNP))]
            Data_exp1.reset_index(drop=True, inplace=True)
            Data_out1.reset_index(drop=True, inplace=True)
        Data_exp1.to_csv(
            path_original + "/" + name_EX + "_" + "_pruned.csv",
            sep=",",
            index=False)
        Data_out1.to_csv(
            path_original + "/" + name_EY + "_pruned.csv",
            sep=",",
            index=False)
        print("New Dataset has been saved in the directory " + path_original)
    else:
        sys.exit('Error Message : Dataset is empty. There were no SNPs that passed the p-value threshold.')


SNPs_LDrange(file_EX1, file_EY1, LD_threshold_lower1, LD_threshold_upper1)
