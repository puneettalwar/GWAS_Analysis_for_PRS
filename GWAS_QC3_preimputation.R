#*********************************************************************************
# GWAS QC3 : Preimputation QC and VCF Check
# R 4.3.1
#*********************************************************************************

# Load necessary libraries
library(rstudioapi)
library(tidyverse)
library(data.table)
library(glue)

#***************************************************************
# Set user variables in configuration sections
#***************************************************************

#========================= CONFIGURATION =============================#

# TOOLS_DIR <- "/media/icaro/T7/CRC_Genetic_Data/Tools"
# PROJECT_DIR <- "/media/icaro/T7"
# DATA_DIR <- "/media/icaro/T7/PD_PRS_GWAS_QC/COF"
# COMMON_FILES <- "/media/icaro/T7/CRC_Genetic_Data/common_files"
# INPUT_PREFIX <- "final_QC2_common_flipped_ref_outlier"
# 
# REF_LEGEND_FILE <- file.path(COMMON_FILES, "1000genome/1000GP_Phase3_combined.legend")
# REF_FASTA <- file.path(COMMON_FILES, "1000genome/human_g1k_v37.fasta")
# CHECKVCF_SCRIPT <- file.path(TOOLS_DIR, "checkVCF.py")
# 
# # Derived paths
# VCF_DIR <- file.path(DATA_DIR, "Data_vcf_check")
# PREIMP_DIR <- file.path(VCF_DIR, "preimputation_qc")
# IMPUTE_INPUT_DIR <- file.path(PREIMP_DIR, "Input_imputation")
# PLINK_EXEC <- file.path(TOOLS_DIR, "plink1.9")

#=====================================================================#

# Set environment
Sys.setenv(PATH = paste(TOOLS_DIR, Sys.getenv("PATH"), sep = ":"))
setwd(DATA_DIR)

# Create required directories
dir.create(VCF_DIR, showWarnings = FALSE)
dir.create(PREIMP_DIR, showWarnings = FALSE)

# Set working directory to preimputation
setwd(PREIMP_DIR)

# Copy PLINK input files
file.copy(file.path(DATA_DIR,"Data_merge_kgb", paste0(INPUT_PREFIX, ".bed")), ".")
file.copy(file.path(DATA_DIR,"Data_merge_kgb", paste0(INPUT_PREFIX, ".bim")), ".")
file.copy(file.path(DATA_DIR,"Data_merge_kgb", paste0(INPUT_PREFIX, ".fam")), ".")
file.copy(PLINK_EXEC, "plink1.9")

# Frequency calculation
system(glue("./plink1.9 --bfile {INPUT_PREFIX} --out {INPUT_PREFIX}_freq --freq --keep-allele-order"))


#***************************************************************
# Run Will Rayner's QC tool before CheckVCF tool
#***************************************************************

# https://www.chg.ox.ac.uk/~wrayner/tools/
# https://github.com/huw-morris-lab/imputation
# https://imputationserver.readthedocs.io/en/latest/prepare-your-data/
# wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
# https://www.chg.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz

## perl script to check plink .bim files against HRC/1000G 
## strand, id names, positions, alleles, ref/alt assignment

# Run HRC-1000G QC tool

# Takes 15 minutes
# 1000G reference with EUR population [output will be in the working directory] 

system(glue(
  "perl {TOOLS_DIR}/HRC-1000G-check-bim.pl ",
  "-b {PREIMP_DIR}/{INPUT_PREFIX}.bim ",
  "-f {PREIMP_DIR}/{INPUT_PREFIX}_freq.frq ",
  "-r {REF_LEGEND_FILE} -g -p EUR -v -t 0.2 -n"
))

# Use tool generated bash script for QC
# Modify 'Run-plink.sh' to use local plink1.9

SCRIPT_PATH <- "Run-plink.sh"
if (file.exists(SCRIPT_PATH)) {
  SCRIPT_LINES <- readLines(SCRIPT_PATH)
  SCRIPT_LINES <- gsub("\\bplink\\b", "./plink1.9", SCRIPT_LINES)
  writeLines(SCRIPT_LINES, SCRIPT_PATH)
  system("sh Run-plink.sh")
}

# As there is no chr 23 script will give error [Error: All variants excluded.]
# This Will Rayne script splits your data into separate chromosomes but use the combined file to create vcf

# Correct REF/ALT swaps using reference
system(glue(
  "plink2 --bfile {INPUT_PREFIX}-updated ",
  "--ref-from-fa --fa {REF_FASTA} ",
  "--make-bed --out {INPUT_PREFIX}_updated_ref"
))

# Convert to VCF
system(glue(
  "plink2 --bfile {INPUT_PREFIX}_updated_ref ",
  "--recode vcf-iid bgz --out Final_Data_QC3_Preimputation"
))

#*****************************
# CHECK VCF BEFORE IMPUTATION
#*****************************

setwd(VCF_DIR)

# https://github.com/zhanxw/checkVCF?tab=readme-ov-file
# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/

# Check VCF using checkVCF tool [Works with python2 only]
# python checkVCF.py -r hs37d5.fa.gz -o test $your_VCF  - USE human_g1k_v37.fasta [https://imputationserver.readthedocs.io/en/latest/prepare-your-data/]
# output will be saved in log file
# It's done if message appears - "No error found by checkVCF.py, thank you for cleanning VCF file."
# In case of issue check the step that requires ACTION.

system(glue(
  "python2.7 {CHECKVCF_SCRIPT} -r {REF_FASTA} -o Final_Data_Checkvcf_Preimputation ",
  "preimputation_qc/Final_Data_QC3_Preimputation.vcf.gz"
))

# Split VCF by chromosome
setwd(PREIMP_DIR)
dir.create(IMPUTE_INPUT_DIR, showWarnings = FALSE)
system("tabix -p vcf Final_Data_QC3_Preimputation.vcf.gz")

for (CHR in 1:22) {
  CMD <- sprintf(
    "bcftools view -r %d Final_Data_QC3_Preimputation.vcf.gz -Oz -o %s/chr%d.vcf.gz",
    CHR, IMPUTE_INPUT_DIR, CHR
  )
  system(CMD)
}

cat("\n✔ FINAL QC3 COMPLETE!\nYour files are ready for imputation in the Input_imputation folder\n")

#*********************************************************************************
# QC3 COMPLETE !!! 
# RUN Imputation
#*********************************************************************************