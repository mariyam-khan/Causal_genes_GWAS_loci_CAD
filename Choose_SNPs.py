from numpy import *
import numpy as np
import pandas as pd
import sys
import os
from pathlib import Path

"""This function takes as input 



1. 1st argument, exposure data which you used to extract the GWAS summary data from 
the MR-Base package 

2. 2nd argument GWAS summary data 

3. 3rd argument chromosome 


Additionally can be given as input:

4. 4th argument position. 

5. distance, where we will remove SNPs which are at distance(radius) greater than position from the 
   lead SNP.

6. 6th argument pvalue threshold

If you do not give the 4th argument, the lead SNP (SNP with the lowest p-value in the GWAS data is chosen)

If you do not give the 5th  and 6th arguments then default are:

position =  5000000
pvalue = 5E-8



Once you  run this file, you will get SNPs on the chromosome (integer given as argument for chromosome number) 1Mb 
around the position of the lead SNP you gave as argument for position, saved in exposure and outcome data .csv files. 
These files  will be saved in the same directory as the original files with the suffix of the chromosome and position 
appended to them.

To use this function, make sure your output and exposure data are in the format you need for the MRBase package.

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


def SNPs_pos_range(file_EX, file_EY, chromosome, position, distance, pvalue ):
    Data_out_all = pd.read_csv(file_EY, sep=',')
    Data_exp = pd.read_csv(file_EX, sep=',')

    Data_out_all = Data_out_all[
        ["SNP", "effect_allele.outcome", "other_allele.outcome", "beta.outcome", "pval.outcome", "chr", "pos"]]
    Data_out = (Data_out_all.loc[Data_out_all['chr'] == chromosome])
    Data_out.reset_index(drop=True, inplace=True)
    Data_out = Data_out[Data_out['pval.outcome'] <= pvalue]
    Data_out.reset_index(drop=True, inplace=True)
    Data_out = Data_out.sort_values(by='pval.outcome', ascending=True)
    Data_out.reset_index(drop=True, inplace=True)
    if position is None:
        position = Data_out.iloc[0, 'pos']
    if isinstance(Data_out, pd.DataFrame):
        if not Data_out.empty:
            lower = int(position) - 500000
            upper = int(position) + 500000

            Data_out = Data_out[(Data_out['pos'] >= lower) & (Data_out['pos'] <= upper)]
            Data_out.reset_index(drop=True, inplace=True)
            Data_exp = Data_exp.loc[Data_exp['SNP'].isin(np.intersect1d(Data_exp.SNP, Data_out.SNP))]
            Data_out = Data_out.loc[Data_out['SNP'].isin(np.intersect1d(Data_out.SNP, Data_exp.SNP))]
            Data_exp.reset_index(drop=True, inplace=True)
            Data_out.reset_index(drop=True, inplace=True)
            if isinstance(Data_out, pd.DataFrame):
                if not Data_out.empty:
                    Data_exp.to_csv(path + "/" + name_EX + "_" + str(chromosome) + "_" + str(position) + ".csv",
                                    sep=",",
                                    index=False)
                    Data_out.to_csv(path + "/" + name_EY + "_" + str(chromosome) + "_" + str(position) + ".csv",
                                    sep=",",
                                    index=False)
                else:
                    sys.exit('Error Message : Dataset is empty, no SNPs shared with the exposure data.')
            else:
                sys.exit('Error Message : Dataset is empty, no SNPs shared with the exposure data.')
    else:
        sys.exit('Error Message : Dataset is empty.There were no SNPs that passed the p-value threshold.')
