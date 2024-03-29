# Causal genes at GWAS loci for Coronary Artery Disease
  Identification of causal genes at GWAS loci with pleiotropic gene regulatory effects using instrumental variable sets.

  
# 1. Guide to software requirements


To use this code, download the code files and ensure that you have the dependencies explained below. To run the code files, you have to type:

        python3 run_MVMR.py "/home/user/file.csv" "/home/user/ld.csv"

Here, python3 should refer to at least **Python 3.5** and depends on your specific installation of Python.

Before running the code files, check whether you have the following requirements and install them if necessary:


        Python 3.5 or later
        
You can check this by typing ’python’ (or a more specific command as explained above) in the command line. For further support, in particular how to
install python please visit https://www.python.org/.

This Python version has the packages numpy (version 1.11.0 or later), scipy, pandas, sys and rpy2 (version 2.9.4 or later).You can check this by starting this python version (check it especially carefully if you have multiple Python versions on your system) and typing


        import numpy
        import pandas
        import rpy2
        import sys
        import scipy


You would also need statistical programming language **R 3.2.0** or later. You can check this by typing ’R’ in the command line. To install R, please visit
https://www.r-project.org/. If you need to use GWAS summary data for your exposure data (gene expression data) or if you need to get the LD-matix of *Snps* in your data, you can use the package TwoSampleMR (https://mrcieu.github.io/TwoSampleMR/). It uses the IEU GWAS database to obtain data automatically, and you can install and call them in R using the following commands


        install.packages("devtools")
        remotes::install_github("MRCIEU/TwoSampleMR")
        remotes::install_github("MRCIEU/MRInstruments")
and

        library(devtools)
        library(TwoSampleMR)
        library(MRInstruments)



In a situation where some packages are not compatible, you can make sure to have exactly these versions installed:

        pandas 1.5.3
        rpy2 3.5.13
        numpy 1.21
        R  version: 4.3.1

        
# 2. Guide to data requirements

To estimate the causal effect, the minimum information required is as follows:


- A. *Snps* to exposure effect (exposure can be *expression of a gene*).

- B. *Snps* to outcome effect (outcome can be a disease like *Coronary Artery Disease*).


For the causal analysis, A. and B.  need to be in one file as described in **Section 5.1**


- C. Lastly, optional is LD-matrix of the *Snps* (*Snps* in the *Snps*-exposure/outcome data). This is optional because, *run\_MVMR.py* has in-built functionality to run the analysis *without the user providing this LD-matrix* **(Section 5.1)**. This is not optional, in the case you want to use *run\_MVMR\_LD.py* as this function allows the user to *specify their own LD-matrix* for their data **(Section 5.2)**.


# 3. Steps to running the code

**Case 1:** If you already have the data in the format containing the *Snps* to exposure effects and *Snps* to outcome effects **(Section 5)**. 

   - **Case 1.1 :** If you do not provide an LD-matrix.
   
     - **Step 1 (Section 5.1) :** `run_MVMR.py` 
   
   
   - **Case 1.2 :** If you wish to provide your own LD-matrix.
   
     - **Step 1 (Section 5.2)** : `run_MVMR_LD.py`, LD-matrix can be generated in R as described in **Section 5.2.1**
     
   
 
*In Case 1, the analysis is finished after these steps and you will be returned the *.csv* files of estimates of the causal parameters by different methods and their standard errors.*



     
**Case 2:** If you have the dataset from GTEx and wish to prepare the datasets to run the analysis. 

  - **Step 1 (Section 4.1) :** Download exposure data from GTEx https://gtexportal.org/home/datasets
  
  - **Step 2 (Section 4.2) :** Align the *Snp* id's in the GTEx data with corresponding rs id's `Match_rs_id.py`.

  - **Step 3 (Section 4.3.1) :** Extract outcome data using TwoSampleMR package. 
  
  - **Step 4 (Section 4.3.2) :** Harmonize the outcome and exposure data using TwoSampleMR package. 

  - **Step 4 (Section 4.4) :** From the harmonized dataset, get data on a specific chromosome and within *chosen* distance and LD-range. You can                                        either 

    - **Step 4.1 (Section 4.4.1) :**   Choose your *Snp* and get data within *chosen* distance using `Prune_Snps_pos.py`.
    - **Step 4.2 (Section 4.4.2) :**   Choose your *Snp* and get data within *chosen* LD-range using `Prune_Snps_LD.py`.
    
    You can run both Step 4.1 and Step 4.2 in any order to get *Snps* within *chosen* distance and *chosen* LD-range.

    Alternate
    
- **Step 4  (Section 4.4.3) :**   For entire harmonized data, get datasets of each chromosome and *Snps* within *chosen* distance and *Snps*                                              within        *chosen* LD-range using `Segregate_datasets.py`.

- **Step 5 (Section 5) :** Run the causal analysis, as in Case 1.
  
  Alternate
   
- **Step 5  (Section 5.2) :** In case you want to run the causal analysis for all datasets in a directory (as generated by Step 4, Section 4.5.4) simultaneously,                                      you can use `MVMR_directory.py`



# 4. Data and Preparation


We have used this method to estimate the causal effect of *genes* which are shared on a locus on outcome *Coronary Artery Disease* using summary statistics from *genome wide association studies (GWAS)*. We used for GWAS summary data, *ebi-a-GCST003116* with trait as coronary artery disease, from the year 2015. The summary statistics were downloaded from the **TwoSampleMR Package**  for a European population.

## 4.1 Download exposure data from GTEx

For the summary data on eQTL analysis, we have used STARNET for association analysis from instruments (*Snps*) to exposures (*genes*). But since this
data is not publicly available, exposure data fom GTEx can be download from https://gtexportal.org/home/datasets. Here in the Single-Tissue cis-QTL
data, you can download the full summary statistics of the cis-eQTLs mapped in European-American subjects. You can check the alignment of the effect allele
from https://www.gtexportal.org/home/faq#interpretEffectSize.



When you download exposure data from Gtex, you would have data in the format
such as:

 ![Format of exposure data from GTEx](gtex_exp_data.png) 

## 4.2 Align the *Snp*-ID in the GTEx data with corresponding rs-ID

Since the data you downloaded is not in format for extraction of GWAS summary data using the **TwoSampleMR Package**, you will need to get the GTEx exposure data in **TwoSampleMR Package** format. You would need to align the variant_id with rs-ID of the *Snps* using the GRCh37 build. For this, you would need to download the *Snps* annotation file from GTEx. This file can be found in the *Reference* section of the GTEx portal:

 ![Location of the GTEx *Snp* annotation file](Gtex_annotation_des.png) 

 The following image gives the format of this annotation file as downloaded from GTEx:

 ![Format of annotation data from GTEx](gtex_annotation_file.png) 
 
 You can then run the function below


        python3 Match_rs_id.py "/home/user/Exposure.csv" "/home/user/annotation.csv" 

Here, *Exposure.csv* is the gene expression data for a specific tissue you downloaded from GTEx and *annotation.csv* is the annotation file used to align the *Snp*-ID in the GTEx data with corresponding rs-ID. You will be returned *.csv* file with the annotated rs-ID’s in the same directory as your *Exposure.csv* file. The following image shows you how the annotated output of  `Match_rs_id.py` looks like:

 ![Format of annotated file returned from `Match_rs_id.py`](gtex_annotated.png) 

In summary you need to ensure that the GTEx data is in the format required by the **TwoSampleMR Package**.

To extract outcome data for the study of interest, we used the **TwoSampleMR Package** (package for performing Mendelian randomization using GWAS
summary data, https://mrcieu.github.io/TwoSampleMR/). Note that if you are using GTEx for exposure datasets, you have to make sure that this minimum information is provided for the extraction of the GWAS summary data from the MRBase package:


        SNP - rs ID.
        beta - The effect size (binary traits log(OR)).
        se - The standard error of the effect size.
        effect_allele - allele of SNP which has the effect marked in beta.
	
If you have used the function `Match_rs_id.py`, you would automatically have the data is format for **TwoSampleMR Package**, if not you need to ensure that the following mapping is done from GTEx:


 	    'phenotype_id': 'exposure',
	    'variant_id': 'variant_id',
	    'rs_id_dbSNP151_GRCh38p7': 'SNP',
	    'maf': 'eaf.exposure',
	    'ref': 'other_allele.exposure',
	    'alt': 'effect_allele.exposure',
	    'slope': 'beta.exposure',
	    'slope_se': 'se.exposure',
	    'pval_beta': 'pval.exposure',
	    'chr': 'chr',
	    'variant_pos': 'pos'

Here the value 'SNP' needs to be the rs-ID of the *Snps* using the GRCh37 build as required in the  **TwoSampleMR Package**.

## 4.3.1 Extract outcome data using TwoSampleMR package 

Now if you have the exposure data in the correct format, you can get the GWAS summary data as follows:


The `available_outcomes` function returns a table of all the available studies in the database. Each study has a unique ID.

        library(TwoSampleMR)
        available_outcomes <- available_outcomes()
        View(available_outcomes)
        
        
After this you can extract the outcome data for any study using this code that specifies these steps in R:

        id_outcome <- ebi-a-GCST005195
        outcome_data <- extract_outcome_data(exposure_data$SNP, id_outcome)
        write.csv(outcome_data, "Outcome.csv" ,row.names = FALSE)
        write.csv(exposure_data, "Exposure.csv" ,row.names = FALSE)
        

## 4.3.2 Harmonize the outcome and exposure data using  **TwoSampleMR Package** package

You can now harmonize the exposure and outcome datasets using the following command in the  **TwoSampleMR Package**. Make sure that you use the option to remove palindromes. The harmonize function will not remove the paliindromic *Snps* but it will flag them as TRUE which we will then use to remove them in our MVMR analysis.

	dat <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data, action = 3)

 The following image shows you what the dataset recieved after using the harmonise_data outcome for the exposure data of the Aorta tissue from GTEx.

 ![Format of harmonised file returned from `harmonise_data`](gtex_harmonised.png) 
		
		

## 4.4 From the GTEx data, get the data on a specific chromosome and within *chosen* distance and *chosen* LD-range



Now that you have the harmonised data, you can choose a chromosome and genetic variants within *chosen* distance of each other and within *chosen* Ld-range, to run the causal analysis. 


For all the following functions, we use *harmonised dataset* to mean the output of the harmonisation function in **TwoSampleMR** package. Please make sure you have these columns in your harmonised dataset:


		     gene name or ID : 'exposure',
		     rs_id of SNP GRCh37 : 'SNP',
		     Marker if the SNP is plaindromic or not : 'palindromic',
		     other allele in exposure data : 'other_allele.exposure',
		     other allele in outcome data : 'other_allele.outcome',
		     effect allele in exposure data: 'effect_allele.exposure',
		     effect allele in outcome data: 'effect_allele.outcome',
		     effect size for exposure : 'beta.exposure',
		     effect size for outcome : 'beta.outcome',
		     pvalue of the effect in exposure : 'pval.exposure',
		     pvalue of the effect in outcome : 'pval.outcome',
		     chromosome of the SNP : 'chr' or 'chr.x',
		     position of the SNP : 'pos' or 'pos.x'

       

The following image shows a short pipeline of functions in Section 4.4.1 (`Prune_Snps_pos.py`) , 4.4.2 (`Prune_Snps_LD.py`) and their usage for MVMR analysis in Section 5.1 (`run_MVMR.py`) , when you do not provide the LD-matrix for the analysis and Section 5.2 (`run_MVMR_LD.py`) , when you provide an LD-matrix for the analysis.

  ![short pipeline](Pipeline2.png) 
  
### 4.4.1  Choose your *Snp* and get exposure and outcome data within *chosen* distance

To make sure that you have data from the same locus and of *Snps* within *chosen* distance (ex. 1Mb) of each other for the MVMR analysis, you can choose a *Snp* and run the function `Prune_Snps_pos.py`. Once you run this, you will get *Snps* (that pass a certain p-value threshold, ex. 5E-8) on the chromosome (integer given as argument for chromosome number) within *chosen* distance around the position of the *Snp* (specified by the argument for position of this *Snp*), saved in *.csv* file.  In the following image, you see a harmonized file *Aor_no_palindromes.csv* as returned using the harmonisation fuction in the  **TwoSampleMR Package** and the new dataset for the causal analysis is saved in the same directory in the new folder *Datasets*.

![Directory craeted for new dataset in Prune_Snps_pos.py](git_pic3.png)

To achieve this dataset you need:

Compulsory:

- Harmonized Dataset/Path to harmonized dataset.

- Chromosome on which your chosen *Snp* lies.

Additionally (optional),

- Position (int, ex 109817192) on which your chosen *Snp* lies. (Default: Position of lead *Snp*)

- Distance (int, example 250000 for 0.25 Mb radius), *Snps* within this distance around the chosen/lead *Snp* will be included for a given dataset. 
  default: 500000

- P-value threshold (float, example 5E-8), *Snps* with p-value lower than this (in the GWAS data), will not be included.
  default: 5E-8
 


As an example, if you run


        python3 Prune_Snps_pos.py "/home/user/data.csv"  1 109817192 250000 5E-8
        
        
Here the  *Snp* has position 109817192 on Chromosome 1 and you wish to have *Snps* around this *Snp* within 0.25Mb of distance. You will then get
the same data\_1\_109817192.csv files with SNPs which are significant (p − value ≤ 5E-08) and 0.25Mb around 109817192 on
chromosome 1. Please make sure you give the p-value in the format ()E-(), distance as int (250000, not 250000.0) and position as int (109817192, not 109817192.0).

The following image gives as example of the output file of the `Prune_Snps_pos.py` function.

![Output exposure dataset of Prune_Snps_pos.py](pos_prep_tog_outcome.png)



### 4.4.2  Choose your *Snp* and get exposure and outcome data within *chosen* LD-range

To make sure that you have data from the same locus and have *Snps* above *chosen* LD-threshold (ex. 0.01) excluded from the MVMR analysis, you can choose a *Snp* and run the function `Prune_Snps_LD.py`. Once you run this, you will get *Snps* (that pass a certain p-value threshold, ex. 5E-8) on the chromosome (integer given as argument for chromosome number) within *chosen* distance around the position of the *Snp* (specified by the argument for position of this *snp*), saved in  *.csv* files. These files will be saved in the new folder (datasets) within the same directory as the original files with the suffix of the chromosome and position appended to them.


Please give as input:

Compulsory:

- Harmonized Dataset/Path to harmonized dataset.

- Chromosome on which your chosen *Snp* lies.

Additionally (optional),

- Position (int, ex 109817192) on which your chosen *Snp* lies. (Default: Position of lead *Snp*)

- Lower LD threshold (float, example 0.01), *Snps* with LD less than or equal to 0.01 with lead *Snp* will not be included. 
  default: 0.01

- P-value threshold (float, example 5E-8), *Snps* with p-value lower than this (in the GWAS data), will not be included.
  default: 5E-8


As an example, if you run


        python3 Prune_Snps_LD.py "/home/user/data.csv" 1 109817192 0.01 5E-8
        
        
Here the  *Snp* has position 109817192 on Chromosome 1 and you wish to have *Snps* around this *Snp* within LD greather than or equal 0.01. You will then get the same data\_1\_109817192.csv file with SNPs which are significant (p − value ≤ 5E-08) and within LD >= 0.01 around 109817192 on chromosome 1. Please make sure you give the p-value in the format ()E-(), LD as float (1.0, not 1) and position as int (109817192, not 109817192.0).

The following image gives as example of the output file of the `Prune_Snps_LD.py` function.


![Output exposure dataset of Prune_Snps_LD.py](LD_prep_tog_output.png)



### 4.4.3  For entire exposure/outcome data, get datasets of each chromosome and *Snps* within *chosen* distance and *chosen* LD-range

The following image shows a short pipeline of functions in 4.4.3 (`Segregate_datasets.py`) and their usage for MVMR analysis in Section 5.3 (`MVMR_directory.py`) for causal analysis.

  ![short pipeline 2](Pipeline3.png)

If there is no specific locus where you wish to perform the analysis but rather the entire exposure and outcome data, we have a code `Segregate_datasets.py` which can output the data needed for MVMR analysis, segregated per chromosome and "hotspots" per chromosome.

A *hotspot* is characterized by the lead *Snp* (the *Snp* which is the most GWAS significant) and *Snps* around it which are either a. certain *distance* around the lead *Snp* (example 0.5 Mb) or b. are in a minimum amount of LD with the lead *Snp* *(lower LD threshold)*.


Please give as input:

Compulsory:

- Harmonized Dataset/Path to harmonized dataset.

  Additionally (optional),

- Distance (int, example 250000 for 0.25 Mb radius), *Snps* within this distance around the lead *Snp* will be chosen for a given dataset. 
  default: 500000

- Lower LD threshold (float, example 0.01), *Snps* with LD less than or equal to 0.01 with lead *Snp* will not be included. 
  default: 0.01

- P-value threshold (float, example 5E-8), *Snps* with p-value lower than this (in the GWAS data), will not be included.
  default: 5E-8



For the following code


        python3 Segregate_datsets.py "/home/user/Exposure.csv" 500000 0.01 1.0


You will have multiple files saved in your initial directory, for each chromosome and for different *hotspots*. To use this function, make sure your harmonized data are in the format you got from the **TwoSampleMR Package**  If you do not provide any distance, the default distance is taken to be 500000, while the default lower LD threshold is 0.0.


In the following images, you see the data saved in directories *Datasets*,  and in the second image, an example output dataset on on Chromosome 1 and and *hotspot* around the lead GWAS *Snp* at position 109817192.

![Seperate Exposure and outcome datasets saved for each chromosome and hotpot, in the directory Datasets.](git_pic2.png)
![Data prepared in format required for causal analysis, saved for each chromosome and hotpot, in the directory Data prepared.](out_segdata.png)


# 5. Scripts for causal analysis


There are two major files for running the causal analysis.

## 5.1  If you do not provide an LD-matrix.

          python3 run_MVMR.py "/home/user/file.csv"
          
After python3, the argument should be the file.csv file containing the *Snps* to exposure effects and SNPs to outcome effects. The first column is
always the SNP ID’s, the last column is always the SNPs to outcome effect. Every other column in between is treated as an exposure variable. The separator to be used is comma. As an example:


        SNPs,gene1,gene2,outcome
        rs11191416,0.5,0.37,0.079
        rs7098825,0.34,0.0,0.078
        rs17115100,0.4,0.54,0.05

## 5.2  If you provide an LD-matrix

        python3 run_MVMR_LD.py "/home/user/file.csv" "/home/user/ld.csv"


After python3, the first argument should be the file.csv file containing the *Snps* to exposure effects and *Snps* to outcome effects. The first column is
always the *Snp* ID’s (the ID is irrelevant if you are providing the LD-matrix but should still be filled with default values), the last column is always the *Snps* to outcome effect. Every other column in between is treated as an exposure variable. The separator to be used is comma. As an example:


        SNPs,gene1,gene2,outcome
        rs11191416,0.5,0.37,0.079
        rs7098825,0.34,0.0,0.078
        rs17115100,0.4,0.54,0.05



The second argument is ld.csv file for the LD-matrix. Please make sure the ordering of the SNPs is same as in the exposure file.csv file. This has the format:


        1.0,0.9,0.8
        0.9,1.0,0.7
        0.8,0.7,1.0


The delimiter is comma. The results are saved as .csv file in the same directory as the one given for the exposure file.csv file. The code automatically
prunes for *Snps* in perfect-LD and keeps only the first occurring *Snps*. If you would like to keep the most significant *Snps* amongst perfect LD *Snps* then
please order the *Snps* in decreasing order of significance in both .csv files.


### 5.2.1  If you provide an LD-matrix

The LD-matrix can be generated using the **TwoSampleMR function** `ld_matrix`. Please make sure all *Snps* belong
to the LD panel. To get the LD-matrix, run the following commands in R



      snps <- list(’rs7776079’,’rs36049381’,’rs9349379’)
      matrix <- ld_matrix(snps, with_alleles = FALSE, pop = "EUR")


In snps, add the *SNps* you wish to get the LD-matrix for.


The following code will save the LD-matrix in the format required for our analysis


        write.table(matrix,file= "ld.csv",sep=" ",quote=F,col.names=F,row.names=F)


## 5.3  If you want to run the causal analysis for all datasets in a given directory.

In this case you do not need to provide the LD-matrix and run the causal analysis by the foloowing code:

        python3 MVMR_directory.py "/home/Data_prepared/"
    
Here Data\_prepared is the directory containing multiple datasets in the format as given in **Section 5.1**
        
        
        
# 6. Error messages


You can get the following errors while using the code:



        Error Message : You require at least as many instruments as exposures to run this analysis.


In this case you cannot run the analysis. The code prunes for SNPs in perfect LD so it may happen that you have more *Snps* in the dataset
you provided but they end up getting removed in the pruning and you get this error.


        R[write to console]: The following variants are not present in the LD reference panel rs28789513. 


In this case you should remove the mentioned *Snps* (rs28789513) from the dataset. You can only get this error if you do not provide
an LD-matrix and it is generated using the R-package TwoSampleMR.


        rpy2.rinterface_lib.embedded.RRuntimeError: Error in matrix(as.numeric(res), nrow(res), ncol(res)) : 
        non-numeric matrix extent

This is a random error with the ld matrix function and re-running ususally resolves it.


Packages errors: In a situation where some packages are not compatible, you can make sure to have exactly these versions installed:

        pandas 1.5.3
        rpy2 3.5.13
        numpy 1.21
        R  version: 4.3.1


# 7. Comparison to other methods


You can compare your results to other methods in the MVMR community like TWMR https://www.nature.com/articles/s41467-019-10936-0, MVMR
https://onlinelibrary.wiley.com/doi/10.1002/gepi.21758, and functions in the MR-Base package https://mrcieu.github.io/TwoSampleMR/.
For our simulations, we used TWMR as from their github page (code found under https://github.com/eleporcu/TWMR), MRBase with the exposure/
outcome data and analysis after harmonization as metioned in https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html#multivariable-mr.
As for MVMR, we supplied data as mentioned in their vignette and then ran the analysis https://cran.r-project.org/web/packages/MendelianRandomization/vignettes/Vignette_MR.pdf.



# 8. FAQ's

1. Why are the estimates for the causal parameter different using Least-squares and GMM?
 - This can be the case if the determinant of the covaraince matrix of the instruments is very small. In these cases,, from our simulations, Least-squares 
   estimates tend to be more relaible. If you want, you can also prune for *Snps*  which are in very high LD. You can set a threshold for the maximum LD      theshold such that, out of the *Snps* in LD higher than a certain value (example 0.95), only the *Snp* with the lowest GWAS significance will be kept      and others will be removed from the analysis. By doing this, the determinant of LD-matrix will not be as small (ideally not smaller than 0.01) and the      GMM estimator will inch closer to the Least-squares estimate.
