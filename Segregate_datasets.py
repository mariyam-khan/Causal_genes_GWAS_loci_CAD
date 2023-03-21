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

ro.r['options'](warn=-1)

base = importr('base')
base.warnings()
"""This function takes as input 


"""


if len(sys.argv) < 3 or len(sys.argv) > 6:
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
        distance1 = int(sys.argv[3])
        if not isinstance(distance1, int):
            print(
                "Provide a non-zero integer value for distance. Please "
                " read the documentation for more details.")
            sys.exit(1)
        if not distance1 < 1000000 and distance1 > 25000:
            print(
                "Provide distance in the range, 25000 > distance < 1000000. Please "
                " read the documentation for more details.")
            sys.exit(1)

    else:
        distance1 = 500000

    if len(sys.argv) == 5:
        LD_threshold_lower1 = float(sys.argv[4])
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
    elif len(sys.argv) == 6:
        LD_threshold_lower1 = float(sys.argv[4])
        LD_threshold_upper1 = float(sys.argv[5])
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
        LD_threshold_lower1 = None
        LD_threshold_upper1 = None


def get_datasets(file_EX, file_EY, distance, LD_threshold_lower, LD_threshold_upper):
    prepare = True
    Data_outcome = pd.read_csv(file_EY, sep=',')
    Data_exposure = pd.read_csv(file_EX, sep=',')

    Data_outcome = Data_outcome[
        ["SNP", "effect_allele.outcome", "other_allele.outcome", "beta.outcome", "pval.outcome", "chr", "pos"]]
    Data_outcome = Data_outcome[Data_outcome['pval.outcome'] <= 5E-8]
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
                        print(distance)
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
                        if LD_threshold_upper:
                            Data_exp_ld, Data_out_ld = SNPs_LDrange(Data_exp_pos, Data_out_pos, LD_threshold_lower,
                                                                    LD_threshold_upper)
                        else:
                            Data_exp_ld = Data_exp_pos.loc[
                                Data_exp_pos['SNP'].isin(np.intersect1d(Data_exp_pos.SNP, Data_out_pos.SNP))]
                            Data_out_ld = Data_out_pos.loc[
                                Data_out_pos['SNP'].isin(np.intersect1d(Data_out_pos.SNP, Data_exp_pos.SNP))]
                            Data_exp_pos.reset_index(drop=True, inplace=True)
                            Data_out_pos.reset_index(drop=True, inplace=True)
                        ##########################################################################
                        if isinstance(Data_out_ld, pd.DataFrame):
                            if not Data_out_ld.empty:
                                # Path("/Datasets").mkdir(parents=True, exist_ok=True)
                                path_new = os.path.join(path_original, "Datasets")

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


def SNPs_LDrange(Data_exp, Data_out_pos, LD_threshold_lower, LD_threshold_upper):
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
                   any(upper_tri[column] >= int(LD_threshold_upper)) or any(upper_tri[column]
                                                                            <= int(
                       LD_threshold_lower))]
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
    RandS.to_csv(pathRS_new + "/" + name_EX + "_" + str(chromosome) + "_" + str(position) + "_prepared.csv",
                 sep=",")
    print("Prepared Datasets have been saved in the directory " + pathRS_new)


get_datasets(file_EX1, file_EY1, distance1, LD_threshold_lower1, LD_threshold_upper1)
