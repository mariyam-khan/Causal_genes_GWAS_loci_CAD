# Causal_genes_GWAS_loci_CAD
  Identification of causal genes at GWAS loci with pleiotropic gene regulatory effects using instrumental variable sets
  

1. run_MVMR.py
  This is the file for getting the causal genes for cases the user can supply the LD-matrix. This file can be downloaded and run on the terminal as :
  
          python3 run_MVMR.py "/home/user/exposure_outcome.csv" "/home/user/ld.csv"
  
  
After python3, the first argument should be the exposure_outcome.csv file containing the SNPs to exposure effects and SNPs to outcome effects. The first column is always the SNP ID's (the ID is irrelevant if you are providing the LD-matrix but should still be filled with default values), the last column is always the SNPs to outcome effect. Every other column in between is treated as an exposure variable. The separator to be used is comma. As an example:



SNPs,AS3MT,SFXN2,CAD <br />
rs11191416,-0.5,0.37,0.079 <br />
rs7098825,-0.34,0.0,0.078 <br />
rs17115100,-0.4,0.54,0.05 




The second argument is ld.csv file for the LD matrix. Please make sure the ordering of the SNPs is same as in the exposure_outcome file. This has the format:


1.0 0.9 0.8 <br />
0.9 1.0 0.7 <br />
0.8 0.7 1.0 <br />



The delimiter is space. 
The results are saved as .csv file in the same directory as the one given for the exposure_outcome.csv file. 
The code automatically prunes for SNPs in perfect-LD and keeps only the first occurring SNPs. If you would like to keep the most significant SNPs amongst perfect LD SNPs then please order the SNPs in decreasing order of significance in both .csv files.


2. MVMR_withoutLD.py
This file has the same purpose as run_MVMR.py except here you have to run the file  on the terminal as : 


          python3 "/home/user/exposure_outcome.csv" 


i.e. without the file with the LD-matrix. The LD-matrix is generated using the TwoSampleMR function ld_matrix. Please make sure all SNPs belong to the LD panel.


# Guide to software requirements

To use this code, download the code files and ensure that you have the dependencies explained below. To run the code files, you have to type:

          python3 run_MVMR.py "/home/user/exposure_outcome.csv" "/home/user/ld.csv"

Here, python3 should refer to at least Python 3.5 and depends on your specific installation of Python.

Before running the code files, check whether you have the following requirements and install them if necessary:

          Python 3.5 or later

You can check this by typing 'python' (or a more specific command as explained above) in the command line. For further support, in particular how to install python please visit https://www.python.org/.


This Python version has the packages numpy (version 1.11.0 or later), scipy, pandas, sys and rpy2 (version 2.9.4 or later).

You can check this by starting this python version (check it especially carefully if you have multiple Python versions on your system) and typing
    
          import numpy
          import pandas 
          import rpy2


You would also need statistical programming language R 3.2.0 or later

You can check this by typing 'R' in the command line. To install R, please visit https://www.r-project.org/.
If you cannot provide the LD-matrix file, you need to have the package TwoSampleMR installed.



For a detailed explainations, check out the [code vigenette](Code_Vigenette.pdf)
