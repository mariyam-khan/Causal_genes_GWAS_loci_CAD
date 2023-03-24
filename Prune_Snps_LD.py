from numpy import *
import numpy as np
import pandas as pd
import sys
import os
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from pathlib import Path

ro.r['options'](warn=-1)
base = importr('base')
base.warnings()

"""This function takes as input the 1.exposure data file and 2.outcome data file. 

where exposure/outcome data are given as: 

"/home/user/file.csv"

3. 3rd argument chromosome (chromosome of the SNP around which you want to get the data for analysis)

Additionally you can give as input:

4. 4th argument position. (position of the SNP around which you want to get the data for analysis)

5. 5th argument LD_threshold_lower , where we will remove SNPs which are in LD lower than LD_threshold_lower with the 
   lead SNP

6. 6th argument pvalue threshold, where we will remove SNPs whose pvalue is higher in the GWAS outcome data
   than the mentioned pvalue threshold.

If you do not give the  4th argument, the lead SNP (SNP with the lowest p-value in the GWAS data is chosen)

If you do not give the 5th and 6th arguments then default values are:

LD_threshold_lower = 0.0
pvalue = 5E-8



"""


def check_valid(user_input, check):
    out = None
    allowed_extensions = [".csv"]
    if check == "path":
        if not (Path(user_input).exists() and Path(user_input).is_file and Path(
                user_input).suffix in allowed_extensions):
            print(
                "Provide a valid path with .csv extension")
            sys.exit(1)
        else:
            out = user_input
    elif check == "LD":
        if not (0.5 >= user_input >= 0.0):
            print(
                "Provide threshold in the range, 0.0 >= LD_threshold_lower <= 0.5. Please "
                " read the documentation for more details.")
            sys.exit(1)
        else:
            out = user_input

    elif check == "pvalue":
        out = user_input

    elif check == "chromosome":
        if not user_input <= 1 and user_input >= 21:
            print(
                "Please provide a valid chromosome.")
        else:
            out = user_input
    elif check == "position":
        out = user_input
    else:
        print("No parameter given")
    return out


try:
    file_EX1 = sys.argv[1]
    file_EY1 = sys.argv[2]
    file_EX1 = check_valid(file_EX1, "path")
    file_EY1 = check_valid(file_EY1, "path")
    chromosome1 = int(sys.argv[3])
    chromosome1 = check_valid(chromosome1, "chromosome")
    try:
        position1 = float(sys.argv[4])
        position1 = check_valid(position1, "position")
    except ValueError:
        print(
            "Provide an integer value for position. Please "
            " read the documentation for more details.")
        sys.exit(1)
    except IndexError:
        position1 = None

    try:
        LD_threshold_lower1 = float(sys.argv[5])
        LD_threshold_lower1 = check_valid(LD_threshold_lower1, "LD")
    except ValueError:
        print(
            "Provide a valid value for the lower LD threshold. Please "
            " read the documentation for more details.")
        sys.exit(1)
    except IndexError:
        LD_threshold_lower1 = 0.0
    try:
        pvalue1 = float(sys.argv[6])
        pvalue1 = check_valid(pvalue1, "pvalue")
    except ValueError:
        print(
            "Provide a valid value for pvalue threshold. Please "
            " read the documentation for more details.")
        sys.exit(1)
    except IndexError:
        pvalue1 = 5E-8

except ValueError:
    print(
        "Provide a valid value for chromosome. Please "
        " read the documentation for more details.")
    sys.exit(1)

except IndexError:
    print("You failed to either provide the minimum or maximum number of parameters required to run this code. Please "
          " read the documentation for more details")
    sys.exit(1)

path_original = os.path.dirname(file_EX1)
name_EX = os.path.splitext(os.path.basename(file_EX1))[0]
name_EY = os.path.splitext(os.path.basename(file_EY1))[0]
LD_threshold_upper = 1.0


def SNPs_LDrange(file_EX, file_EY, chromosome, position, LD_threshold_lower, pvalue):
    Data_out = pd.read_csv(file_EY, sep=',')
    Data_exp = pd.read_csv(file_EX, sep=',')
    Data_out = (Data_out.loc[Data_out['chr'] == chromosome])
    Data_out.reset_index(drop=True, inplace=True)
    Data_out = Data_out[Data_out['pval.outcome'] <= pvalue]
    Data_out.reset_index(drop=True, inplace=True)

    Data_out = Data_out.sort_values(by='pval.outcome', ascending=True)
    Data_out.reset_index(drop=True, inplace=True)

    if position is None:
        Data_out_new = Data_out
    else:
        idx = Data_out[Data_out['pos'] == position].index.values
        df1 = Data_out.iloc[idx].to_frame('row')
        df2 = Data_out.drop([idx])
        Data_out_new = pd.concat([df1, df2], axis=0)
    if not Data_out_new.empty:
        snps = list(Data_out_new.loc[:, 'SNP'].values)
        if len(snps) != 1:
            cov_data = r("TwoSampleMR::ld_matrix")(snps, with_alleles=False, pop="EUR")

            pandas2ri.activate()
            pd_from_r_df = ro.conversion.rpy2py(cov_data)
            cov_data = pd.DataFrame(pd_from_r_df, columns=r('colnames')(cov_data),
                                    index=r('colnames')(cov_data), dtype='float')

            cov_data = pd.DataFrame.abs(cov_data)
            upper_tri = cov_data.where(np.triu(np.ones(cov_data.shape), k=1).astype(bool))
            to_drop = [column for column in upper_tri.columns if
                       any(upper_tri[column] >= float(LD_threshold_upper))]
            cov_data = cov_data.drop(labels=to_drop, axis=1)
            cov_data = cov_data.drop(labels=to_drop, axis=0)

            upper_tri = cov_data.where(np.triu(np.ones(cov_data.shape), k=1).astype(bool))
            to_drop = [column for column in upper_tri.columns if
                       upper_tri[column][0] <= float(LD_threshold_lower)]
            cov_data = cov_data.drop(labels=to_drop, axis=1)
            cov_data = cov_data.drop(labels=to_drop, axis=0)

            Data_out_new = Data_out_new[Data_out_new['SNP'].isin(cov_data.index.values)]

            Data_exp1 = Data_exp.loc[
                Data_exp['SNP'].isin(np.intersect1d(Data_exp.SNP, Data_out.SNP))]
            Data_out1 = Data_out.loc[
                Data_out_new['SNP'].isin(np.intersect1d(Data_out_new.SNP, Data_exp.SNP))]
            Data_exp1.reset_index(drop=True, inplace=True)
            Data_out1.reset_index(drop=True, inplace=True)
            pandas2ri.deactivate()
        else:
            Data_exp1 = Data_exp.loc[
                Data_exp['SNP'].isin(np.intersect1d(Data_exp.SNP, Data_out_new.SNP))]
            Data_out1 = Data_out.loc[
                Data_out_new['SNP'].isin(np.intersect1d(Data_out_new.SNP, Data_exp.SNP))]
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


SNPs_LDrange(file_EX1, file_EY1, chromosome1, position1, LD_threshold_lower1, pvalue1)
