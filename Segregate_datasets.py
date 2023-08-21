from numpy import *
import numpy as np
import pandas as pd
import sys
from collections import Counter
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

"/home/user/exposure_file_name.csv"
"/home/user/output_file_name.csv"



(This function takes files which have exposure and outcome data across Chromosomes and divides the dataset into sub-datasets
with each sub-dataset having locus-specific data)


Additionally you can give as input:

3. 3rd argument position, where we will remove SNPs which are at distance(radius) greater than position from the 
   lead SNP.
4. 4th argument LD_threshold_lower, where we will remove SNPs which are in LD lower than LD_threshold_lower with the 
   lead SNP
 
5. 5th argument pvalue threshold, where we will remove SNPs whose pvalue is higher in the GWAS outcome data
   than the mentioned pvalue threshold.
   
If you do not give the 3rd, 4th and 5th arguments then default values are:

position =  500000
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

    elif check == "distance":
        if not (1000000 > user_input > 25000):
            print(
                "Provide distance in the range, 25000 > distance < 1000000. Please "
                " read the documentation for more details.")
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
    else:
        print("No parameter given")
    return out


try:
    file_EX1 = sys.argv[1]
    file_EY1 = sys.argv[2]
    file_EX1 = check_valid(file_EX1, "path")
    file_EY1 = check_valid(file_EY1, "path")
    try:
        distance1 = int(sys.argv[3])
        distance1 = check_valid(distance1, "distance")
    except ValueError:
        print(
            "Provide a non-zero integer value for distance. Please "
            " read the documentation for more details.")
        sys.exit(1)
    except IndexError:
        distance1 = 500000
    try:
        LD_threshold_lower1 = float(sys.argv[4])
        LD_threshold_lower1 = check_valid(LD_threshold_lower1, "LD")
    except ValueError:
        print(
            "Provide a valid value for the lower LD threshold. Please "
            " read the documentation for more details.")
        sys.exit(1)
    except IndexError:
        LD_threshold_lower1 = 0.0
    try:
        pvalue1 = float(sys.argv[5])
        pvalue1 = check_valid(pvalue1, "pvalue")
    except ValueError:
        print(
            "Provide a valid value for pvalue threshold. Please "
            " read the documentation for more details.")
        sys.exit(1)
    except IndexError:
        pvalue1 = 5E-8
except IndexError:
    print("You failed to either provide the minimum or maximum number of parameters required to run this code. Please "
          " read the documentation for more details")
    sys.exit(1)

path_original = os.path.dirname(file_EX1)
name_EX = os.path.splitext(os.path.basename(file_EX1))[0]
name_EY = os.path.splitext(os.path.basename(file_EY1))[0]
LD_threshold_upper = 1.0


def get_datasets(file_EX, file_EY, distance, LD_threshold_lower, pvalue):
    prepare = True
    Data_outcome = pd.read_csv(file_EY, sep=',')
    Data_exposure = pd.read_csv(file_EX, sep=',')

    Data_outcome = Data_outcome[
        ["SNP", "effect_allele.outcome", "other_allele.outcome", "beta.outcome", "pval.outcome", "chr", "pos"]]
    Data_outcome = Data_outcome[Data_outcome['pval.outcome'] <= pvalue]
    Data_outcome.reset_index(drop=True, inplace=True)

    if not Data_outcome.empty:

        Data_exposure = Data_exposure.loc[
            Data_exposure['SNP'].isin(np.intersect1d(Data_exposure.SNP, Data_outcome.SNP))]
        Data_outcome = Data_outcome.loc[Data_outcome['SNP'].isin(np.intersect1d(Data_outcome.SNP, Data_exposure.SNP))]
        Data_exposure.reset_index(drop=True, inplace=True)
        Data_outcome.reset_index(drop=True, inplace=True)

        if not Data_outcome.empty:

            unique_chr = Counter(Data_outcome['chr']).keys()
            no_chr = len(unique_chr)
            get_chr = Counter(Data_outcome['chr']).most_common(no_chr)

            for i in range(no_chr):
                Dataframe = (Data_outcome.loc[Data_outcome['chr'] == get_chr[i][0]])
                Dataframe.reset_index(drop=True, inplace=True)
                Dataframe = Dataframe.sort_values(by='pval.outcome', ascending=True)
                Dataframe.reset_index(drop=True, inplace=True)
                chromosome = get_chr[i][0]
                if not Data_outcome.empty:
                    while not Dataframe.empty:

                        # choose SNPs within 'distance' around the lead SNP ######################
                        position = Dataframe.loc[0, 'pos']
                        lower = int(position) - int(distance)
                        upper = int(position) + int(distance)
                        Data_out_pos = Dataframe[(Dataframe['pos'] >= lower) & (Dataframe['pos'] <= upper)]
                        Data_out_pos.reset_index(drop=True, inplace=True)
                        Dataframe = Dataframe.loc[~Dataframe['SNP'].isin(Data_out_pos['SNP'])].dropna()
                        Dataframe.reset_index(drop=True, inplace=True)

                        Data_out_pos = Data_out_pos.loc[
                            Data_out_pos['SNP'].isin(np.intersect1d(Data_out_pos.SNP, Data_exposure.SNP))]
                        Data_out_pos.reset_index(drop=True, inplace=True)
                        Data_exp_pos = Data_exposure.loc[
                            Data_exposure['SNP'].isin(np.intersect1d(Data_exposure.SNP, Data_out_pos.SNP))]
                        Data_exp_pos.reset_index(drop=True, inplace=True)
                        ##########################################################################

                        # prune SNPs according to LD ############################################
                        Data_exp_ld, Data_out_ld = SNPs_LDrange(Data_exp_pos, Data_out_pos, LD_threshold_lower)

                        ##########################################################################

                        if isinstance(Data_out_ld, pd.DataFrame):
                            if not Data_out_ld.empty:
                                path_new = os.path.join(path_original, "Datasets")
                                create_directory_if_not_exists(path_new)
                                Data_exp_ld.to_csv(
                                    path_new + "/" + name_EX + "_" + str(chromosome) + "_" + str(position) + ".csv",
                                    sep=",",
                                    index=False)
                                Data_out_ld.to_csv(
                                    path_new + "/" + name_EY + "_" + str(chromosome) + "_" + str(position) + ".csv",
                                    sep=",",
                                    index=False)
                                print("New Datasets have been saved in the directory " + path_new)
                                path_RandS = path_new
                                if prepare:
                                    Data_preparation(Data_exp_ld, Data_out_ld, path_RandS, chromosome, position)
                                else:
                                    continue
                            else:

                                continue
                        else:
                            continue
                else:
                    continue

            else:
                sys.exit('Error Message : Dataset is empty. There were no SNPs that passed the p-value threshold.')

        else:
            sys.exit('Error Message : Dataset is empty, no SNPs shared with the exposure data.')

    else:
        sys.exit('Error Message : Dataset is empty. There were no SNPs that passed the p-value threshold.')


def create_directory_if_not_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f" ")


def SNPs_LDrange(Data_exp, Data_out_pos, LD_threshold_lower):
    print("LD_threshold_lower", LD_threshold_lower)
    print("LD_threshold_upper", LD_threshold_upper)
    snps = list(Data_out_pos.loc[:, 'SNP'].values)
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
                   upper_tri[column] <= float(LD_threshold_lower)]
        cov_data = cov_data.drop(labels=to_drop, axis=1)
        cov_data = cov_data.drop(labels=to_drop, axis=0)
        Data_out_pos = Data_out_pos[Data_out_pos['SNP'].isin(cov_data.index.values)]

        Data_exp1 = Data_exp.loc[
            Data_exp['SNP'].isin(np.intersect1d(Data_exp.SNP, Data_out_pos.SNP))]
        Data_out1 = Data_out_pos.loc[
            Data_out_pos['SNP'].isin(np.intersect1d(Data_out_pos.SNP, Data_exp.SNP))]
        Data_exp1.reset_index(drop=True, inplace=True)
        Data_out1.reset_index(drop=True, inplace=True)
        pandas2ri.deactivate()
    else:
        Data_exp1 = Data_exp.loc[
            Data_exp['SNP'].isin(np.intersect1d(Data_exp.SNP, Data_out_pos.SNP))]
        Data_out1 = Data_out_pos.loc[
            Data_out_pos['SNP'].isin(np.intersect1d(Data_out_pos.SNP, Data_exp.SNP))]
        Data_exp1.reset_index(drop=True, inplace=True)
        Data_out1.reset_index(drop=True, inplace=True)
    return Data_exp1, Data_out1


def Data_preparation(Data_exp, Data_out, pathRS, chromosome, position):
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
    pathRS_new = os.path.join(pathRS, "Data_prepared")
    create_directory_if_not_exists(pathRS_new)
    RandS.to_csv(pathRS_new + "/" + name_EX + "_" + str(chromosome) + "_" + str(position) + "_prepared.csv",
                 sep=",")
    print("Prepared Datasets have been saved in the directory " + pathRS_new)


get_datasets(file_EX1, file_EY1, distance1, LD_threshold_lower1, pvalue1)
