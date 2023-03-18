#
#  This script selects all significant SNPs for the GTEx files.
#
from datetime import datetime
from os import path, makedirs
from glob import glob
from pandas import read_csv

def main(input_data_path="./", version="test_01", threshold=0.05):
    """ Docstring for main.

    :input_data_path: string, relative path to the input data.
    :version: string, version of this analysis.
    :threshold: float, significance threshold for p-values.

    """
    print("#    Select SNPs from {input_data_path} at p-values below {threshold}")
    # Start time:
    t0 = datetime.now()
    print("Starting at:", t0)
    out_path=input_data_path+"output/"+version+"/"
    if not path.isdir(out_path):
        makedirs(out_path)
    for a_file in glob(input_data_path + "out*.csv"):
        print(a_file)
        df=read_csv(a_file)
        print(df.shape)
        dft=df[df["pval.outcome"]<threshold]
        print(dft.shape)
        dft=dft.dropna(subset=["effect_allele.outcome"])
        df_selection=dft.dropna(subset=["other_allele.outcome"])
        print(df_selection.shape)
        print(df_selection["effect_allele.outcome"])
        print(df_selection["other_allele.outcome"])
        a_file_name=a_file.split("/")[-1]
        df_selection.to_csv(out_path+"significant_snips_from"+a_file_name)
    print(" ### TODO")

    print("Done.")
    # End time:
    t1 = datetime.now()
    print(t1)
    print("It took:", t1 - t0)
    print("End.")

if __name__ == "__main__":
    main(input_data_path="../data/GTEx_data/", version="snps_e_minus8_outcome_exists", threshold=1.0*10**(-8))
