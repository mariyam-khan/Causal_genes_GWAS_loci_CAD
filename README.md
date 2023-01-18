# Causal_genes_GWAS_loci_CAD
  Identification of causal genes at GWAS loci with pleiotropic gene regulatory effects using instrumental variable sets
  
  
1. Results.xlsx has all results chromosome-wise for our analysis using the three GWAS studies with the Starnet and Gtex exposure data. This also includes
   analysis for the causal tissues.
2. Starnet.py is the code for the starnet exposure data, you need to give the exposure data and outcome data, it will divide the files by chromosome
   and position to give you the causal genes at each loci and position.
3. Gtex.py is the same file for the gtex exposure data.
4. Images has 3 components, images of all simulations, images of gene expression levels tissue wise comparision for intresting genes and finally gene          expression plots against genotype data. 

1. run_MVMR.py
  This is the file for getting the causal genes for cases the user can supply the LD-matrix. This file can be downloaded and run on the terminal as :
  
          python3 "/home/user/exposure_outcome.csv" "/home/user/ld.csv"
  
  
After python3, the first argument should be the exposure_outcome.csv file containing the SNPs to exposure effects and SNPs to outcome effects. The first column is always the SNP id's (the ID is irrelevant if you are providing the LD-matrix but should still be filled with default values), the last column is always the SNPs to outcome effect. Everything column in between is treated as an exposure. The separator to be used is comma. As an example:



SNPs,AS3MT,SFXN2,CAD

rs11191416,-0.5,0.37,0.079 

rs7098825,-0.34,0.0,0.078 

rs17115100,-0.4,0.54,0.05 




The second argument is ld.csv file for the LD matrix. Please make sure the ordering of the SNPs is same as the exposure_outcome file. This has the format:


1 0.98852 0.96565 0.98

0.98852 1 0.965852 0.994254 

0.96565 0.965852 1 



The delimiter is space. 
The results are saved as a .csv file in the same directory as the one given for the exposure_outcome.csv file. 
The code automatically prunes for SNPs in perfect-LD and keeps only the first occurring SNPs. If you would like to keep the most significant SNPs amongst perfect LD SNPs then please order the SNPs in decreasing order of significance in both .csv files.


2. MVMR_withoutLD.py
This file has the same purpose as run_MVMR.py except here you have to run the file  on the terminal as : 


          python3 "/home/user/exposure_outcome.csv" 


i.e. without the file with the LD-matrix. The LD-matrix is generated using the TwoSapleMR function ld_matrix. Please make sure all SNPs belong to the LD panel.
