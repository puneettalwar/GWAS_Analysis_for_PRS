#*********************************************************************************
# GWAS QC1-4 : 
# R 4.3.1
# by Puneet 27.05.2025
#*********************************************************************************

################################################################################################
# Before using this script, comment the configuration part in individual scripts (ex. GWAS QC1)
################################################################################################

# Required installations (system): Samtools, BCFtools, HTSlib
# URLs:
# https://www.htslib.org/download/
# https://gist.github.com/adefelicibus/f6fd06df1b4bb104ceeaccdd7325b856

## Required R packages
# install.packages(c("tidyverse","rstudioapi","data.table","plinkQC","getopt","ggplot2","GGally","dplyr","purrr","ggpubr"))

# Load necessary libraries
library(readxl)
library(rstudioapi) 
library(tidyverse)
library(data.table)
library(plinkQC)
library(ggplot2)
library(ggpubr)
library(GGally)
library(dplyr)
library(purrr)
library(glue)
library(stringr)

#*********************************************************************************
# GWAS QC1 : Basic QC after data merging and filtering [GWAS_QC1.R]
#*********************************************************************************

#========================= CONFIGURATION =============================#

# Tool and data paths
TOOLS_DIR <- "/media/icaro/T7/CRC_Genetic_Data/Tools"
PROJECT_DIR <- "/media/icaro/T7"
QC_DIR <- "PD_PRS_GWAS_QC/COF/Data_QC"
DATA_DIR <- "PD_PRS_GWAS_QC/COF"

# Sample name (prefix of binary PLINK files)
DATA_NAME <- "COF_Data"

#=====================================================================#

source(file.path(PROJECT_DIR, "CRC_Genetic_Data/Codes/Puneet/Final_GWAS_QC_scripts/GWAS_QC1.R"))

#***********************************************************************************************
# GWAS QC2 : Using 1000 genome data to match allele frequency & check population stratification [GWAS_QC2_1000G_generic.R]
#***********************************************************************************************

#========================= CONFIGURATION =============================#

# Tool and data paths
TOOLS_DIR <- "/media/icaro/T7/CRC_Genetic_Data/Tools"
PROJECT_DIR <- "/media/icaro/T7"
QC_DIR <- "PD_PRS_GWAS_QC/COF/Data_QC"
DATA_DIR <- "PD_PRS_GWAS_QC/COF"
PATH_TO_GENOME1000 <- "/media/icaro/T7/CRC_Genetic_Data/common_files/1000genome"

# Sample name (prefix of binary PLINK files)
DATA_NAME <- "COF_Data"

#=====================================================================#

source(file.path(PROJECT_DIR, "CRC_Genetic_Data/Codes/Puneet/Final_GWAS_QC_scripts/GWAS_QC2_1000G.R"))

#*********************************************************************************
# GWAS QC3 : Preimputation QC and VCF Check [GWAS_QC3_preimputation.R]
#*********************************************************************************

#========================= CONFIGURATION =============================#

TOOLS_DIR <- "/media/icaro/T7/CRC_Genetic_Data/Tools"
PROJECT_DIR <- "/media/icaro/T7"
DATA_DIR <- "/media/icaro/T7/PD_PRS_GWAS_QC/COF"
COMMON_FILES <- "/media/icaro/T7/CRC_Genetic_Data/common_files"
INPUT_PREFIX <- "final_QC2_common_flipped_ref_outlier"

REF_LEGEND_FILE <- file.path(COMMON_FILES, "1000genome/1000GP_Phase3_combined.legend")
REF_FASTA <- file.path(COMMON_FILES, "1000genome/human_g1k_v37.fasta")
CHECKVCF_SCRIPT <- file.path(TOOLS_DIR, "checkVCF.py")

# Derived paths
VCF_DIR <- file.path(DATA_DIR, "Data_vcf_check")
PREIMP_DIR <- file.path(VCF_DIR, "preimputation_qc")
IMPUTE_INPUT_DIR <- file.path(PREIMP_DIR, "Input_imputation")
PLINK_EXEC <- file.path(TOOLS_DIR, "plink1.9")

#=====================================================================#

source(file.path(PROJECT_DIR, "CRC_Genetic_Data/Codes/Puneet/Final_GWAS_QC_scripts/GWAS_QC3_preimputation.R"))

#*********************************************************************************
# GWAS : Run Imputation on the Michigan Server
#*********************************************************************************

# Upload your chr wise vcf.gz files to the Michigan imputation server https://imputationserver.sph.umich.edu/. 
# Further guidance at https://imputationserver.readthedocs.io/en/latest/.

# These are the settings used:

# Ref panel: HRC r1.1 2016
# rsq filter: OFF (instead 0.3 can be used but then remove "--exclude-if-info 'R2 < 0.3' " in the QC4 script)
# Phasing: Eagle v2.4
# Population: EUR
# Mode: Quality control and imputation

# Download the data into Data_Imputed folder

#*********************************************************************************
# GWAS QC4 : Post imputation QC
#*********************************************************************************

# Use bash terminal to unzip Michigan output files 
#***************************************************************
# Use 7zip (as recommended) [sudo apt install p7zip-full p7zip-rar]
# cd to the folder "/media/icaro/T7/PD_PRS_GWAS_QC/COF/Data_Imputed"
# unzip Michigan output files using the password received from Michigan server in Email


# Password : XChk02xP}?bYRi
# 
# for chr in $(seq 1 22)
# do
# 7z e chr_$chr.zip
# done

#========================= CONFIGURATION =============================#

TOOLS_DIR <- "/media/icaro/T7/CRC_Genetic_Data/Tools"
PROJECT_DIR <- "/media/icaro/T7"
DATA_DIR <- "/media/icaro/T7/PD_PRS_GWAS_QC/COF"
COMMON_FILES <- "/media/icaro/T7/CRC_Genetic_Data/common_files"
DATASET_NAME <- "COF"  # Change for each dataset

VCF_DIR <- file.path(DATA_DIR, "Data_Imputed")
OUTPUT_DIR <- file.path(DATA_DIR, "Post_imputation_QC")
HRC_VCF <- file.path(COMMON_FILES, "HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz")

#=====================================================================#

source(file.path(PROJECT_DIR, "CRC_Genetic_Data/Codes/Puneet/Final_GWAS_QC_scripts/GWAS_QC4_Post_imputation.R"))




