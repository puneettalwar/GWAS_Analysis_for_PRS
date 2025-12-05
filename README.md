# GWAS Data Processing for Polygenic score (PRS) computation
GWAS Pre-processing and QC Pipeline for downstream PRS computation

*** Warning! This Shiny App has not been validated enough! There is no warranty for the app! ***

This repository contains R scripts and documentation for the pre-processing and quality control (QC) pipeline used in Genome-Wide Association Studies (GWAS). 

The aim is to provide quick and simple pipeline in R for GWAS data processing suitable for downstream PRS computation and analysis.

The pipeline follows community standards and best practices, including sample and SNP-level filtering, relatedness checks, population structure inference, and data preparation for PRS computation.



**GWAS_Data_Harmonize script is a R shiny app** 

- It performs merging of two datasets along with sex check and filtering (removal of sex chromosomes, non-rsID SNPs, and ambiguous SNPs.

- It also removes duplicates based on genotype call rates.

- If a dataset has participants from multiple studies, then upload it as Data2 and specify prefix identifier from excel file in the tab

- The app requires the datasets to be named 'Data1' and 'Data2

- The datasets must be accompanied by an Excel sheets name Data1 and Data2 with sample information containing the following columns: 

     - Ref_RegNumber, Index, Sample_ID, Call_Rate, Gender, Gender_code
 
     - Gender_code refers to the actual gender of the participants
 
- You must specify the project directory path (e.g., /media/home/media/Data)

- You must also provide a tools directory path (e.g., /media/icaro/T7/CRC_Genetic_Data/Tools)

- Upload the required input datasets in the 'Upload' tab.

- Run harmonization by clicking the 'Harmonize Data' button.

- Wait for confirmation in the 'Status' tab.

**Note: GWAS_Data_merging script is same as GWAS_Data_Harmonize script but without GUI (not a shiny code)**

**The Final_GWAS_all_QC_preprocessing_script script calls other scripts QC1-QC4**

- GWAS QC1 : Basic QC after data merging and filtering
- GWAS QC2 : Using 1000 genome data to match allele frequency & check population stratification
- GWAS QC3 : Preimputation QC and VCF Check
- GWAS QC4 : Post imputation QC
