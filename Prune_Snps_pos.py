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

5. 5th argument distance_threshold, where we will remove SNPs which are in LD lower than LD_threshold_lower with the 
   lead SNP

6. 6th argument pvalue threshold, where we will remove SNPs whose pvalue is higher in the GWAS outcome data
   than the mentioned pvalue threshold.

If you do not give the  4th argument, the lead SNP (SNP with the lowest p-value in the GWAS data is chosen)

If you do not give the 5th and 6th arguments then default values are:

distance_threshold = 500000
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
        distance1 = int(sys.argv[5])
        distance1 = check_valid(distance1, "distance")
    except ValueError:
        print(
            "Provide a non-zero integer value for distance. Please "
            " read the documentation for more details.")
        sys.exit(1)
    except IndexError:
        distance1 = 500000

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


def SNPs_posrange(file_EX, file_EY, chromosome, position, distance, pvalue):
    Data_out = pd.read_csv(file_EY, sep=',', low_memory=False)
    Data_exp = pd.read_csv(file_EX, sep=',', low_memory=False)

    Data_out = (Data_out.loc[Data_out['chr'] == chromosome])
    Data_out.reset_index(drop=True, inplace=True)
    Data_out = Data_out[Data_out['pval.outcome'] <= pvalue]
    Data_out.reset_index(drop=True, inplace=True)

    Data_out = Data_out.sort_values(by='pval.outcome', ascending=True)
    Data_out.reset_index(drop=True, inplace=True)
    if position is None:
        position = Data_out.at[0, 'pos']
    index = Data_out[Data_out['pos'] == position].index
    if index.empty:
        error_message = (f"Value {position} not found in the 'pos' column. "
                         "Please provide a valid position.")
        sys.exit(error_message)
    if not Data_out.empty:

        # choose SNPs within 'distance' around the lead SNP ######################

        lower = int(position) - int(distance)
        upper = int(position) + int(distance)
        Data_out_pos = Data_out[(Data_out['pos'] >= lower) & (Data_out['pos'] <= upper)]
        Data_out_pos.reset_index(drop=True, inplace=True)
        Data_out_pos = Data_out_pos.loc[
            Data_out_pos['SNP'].isin(np.intersect1d(Data_out_pos.SNP, Data_exp.SNP))]
        Data_out_pos.reset_index(drop=True, inplace=True)
        Data_exp_pos = Data_exp.loc[
            Data_exp['SNP'].isin(np.intersect1d(Data_exp.SNP, Data_out_pos.SNP))]
        Data_exp_pos.reset_index(drop=True, inplace=True)
        ##########################################################################

        Data_exp_pos.to_csv(
            path_original + "/" + name_EX + "_" + "pruned.csv",
            sep=",",
            index=False)
        Data_out_pos.to_csv(
            path_original + "/" + name_EY + "_pruned.csv",
            sep=",",
            index=False)
        print("New Dataset has been saved in the directory " + path_original)
    else:
        sys.exit('Error Message : Dataset is empty. There were no SNPs that passed the p-value threshold.')


SNPs_posrange(file_EX1, file_EY1, chromosome1, position1, distance1, pvalue1)
