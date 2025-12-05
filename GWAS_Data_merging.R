#*********************************************************************************
# GWAS: Data filtering, sorting & merging using Data1 and Data2
# R 4.3.1
#*********************************************************************************

# Data1 refers to original Dataset 1 (copy raw data and rename COF_Data_87 to Data1)
# Data2 refers to original Dataset 2 (copy raw data and rename COF+_48 to Data2)
# Add sheets Data1 and Data2 in the CRC_DNA_consolidated_data.xlsx (copy COF_Data_87 sheet and rename it as Data1)
# Copy Datasets to Data_Merge folder
# Update line number 84 in the script
# Check for duplicates and remove them (modify the code)

#========================= CONFIGURATION =============================#

TOOLS_DIR        <- "/media/icaro/T7/CRC_Genetic_Data/Tools"
PROJECT_DIR      <- "/media/icaro/T7"
QC_DIR           <- file.path(PROJECT_DIR, "PD_PRS_GWAS_QC")
DATA_DIR         <- "/media/icaro/T7/PD_PRS_GWAS_QC/COF"
MERGE_DIR        <- file.path(DATA_DIR, "Data_Merge")
EXCEL_FILE       <- file.path(PROJECT_DIR, "CRC_Genetic_Data/CRC_DNA_consolidated_data.xlsx")
DATA_NAME        <- "COF_Data"

#=====================================================================#

# Load libraries

library(tidyverse)
library(data.table)
library(readxl)
library(glue)

# Set environment & working directory
Sys.setenv(PATH = paste(TOOLS_DIR, Sys.getenv("PATH"), sep = ":"))
setwd(QC_DIR)

# Create required directory
if (!dir.exists(MERGE_DIR)) dir.create(MERGE_DIR, recursive = TRUE)

#============================================================== #
# Run sex check  
# If Data1 has data from only one study with no duplicates
#============================================================== #

# Convert .PED TO .BED
system("plink1.9 --file Data_Merge/Data1 --make-bed --out Data_Merge/Data1")

# Update .fam file with the sex information

fam <- read.table(file.path(MERGE_DIR, "Data1.fam"), header = FALSE)
colnames(fam) <- c("FID", "IID", "FatherID", "MotherID", "Sex", "Phenotype")

excel_data <- read_excel(EXCEL_FILE, sheet = "Data1") %>%
  rename(IID = Sample_ID, NewSex = Gender_code)

fam_updated <- fam %>%
  left_join(excel_data %>% select(IID, NewSex), by = "IID") %>%
  mutate(Sex = NewSex) %>%
  select(-NewSex)

write.table(fam_updated, file.path(MERGE_DIR, "Data1.fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# Check sex 

system("plink1.9 --impute-sex --bfile Data_Merge/Data1 --make-bed --out Data_Merge/Data1_sex_check")
sexcheck <- read.table(file.path(MERGE_DIR, "Data1_sex_check.sexcheck"), header = TRUE)
sex_mismatch <- subset(sexcheck, STATUS == "PROBLEM")

write.table(sex_mismatch[, c("FID", "IID")], file.path(MERGE_DIR, "sex_mismatch.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

system(glue("plink1.9 --bfile {MERGE_DIR}/Data1 --remove {MERGE_DIR}/sex_mismatch.txt --make-bed --out {MERGE_DIR}/Data1_filtered"))

#===============================================================#
# FILTER Data2 and run sex check
# If Data2 has data from more than one study with no duplicates 
#==============================================================#

DATA2_DIR <- file.path(DATA_DIR, "Data2")
if (!dir.exists(DATA2_DIR)) dir.create(DATA2_DIR)

data2_info <- read_excel(EXCEL_FILE, sheet = "Data2")

filtered_data2 <- data2_info %>%
  filter(!is.na(Ref_RegNumber)) %>%
  filter(!grepl("^COF", Ref_RegNumber, ignore.case = TRUE) | Call_Rate < 0.9)

write.table(filtered_data2[, c("Index", "Sample_ID")],
            file.path(DATA2_DIR, "filtered_non_data2.txt"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

system(glue("plink1.9 --remove {DATA_DIR}/Data2/filtered_non_data2.txt --file {DATA_DIR}/Data2/Data2 --make-bed --out {MERGE_DIR}/Data2_filtered"))


# Update .fam file with the sex information

fam2 <- read.table(file.path(MERGE_DIR, "Data2_filtered.fam"), header = FALSE)
colnames(fam2) <- c("FID", "IID", "FatherID", "MotherID", "Sex", "Phenotype")

excel_data2 <- read_excel(EXCEL_FILE, sheet = "Data2") %>%
  rename(IID = Sample_ID, NewSex = Gender_code)

fam2_updated <- fam2 %>%
  left_join(excel_data2 %>% select(IID, NewSex), by = "IID") %>%
  mutate(Sex = NewSex) %>%
  select(-NewSex)

write.table(fam2_updated, file.path(MERGE_DIR, "Data2_filtered.fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)

system(glue("plink1.9 --impute-sex --bfile {MERGE_DIR}/Data2_filtered --make-bed --out {MERGE_DIR}/Data2_sex_check"))

#=========================#
#   SNP Harmonization     
#=========================#

data1_bim <- fread(file.path(MERGE_DIR, "Data1_filtered.bim"))
data2_bim <- fread(file.path(MERGE_DIR, "Data2_filtered.bim"))

filter_snps <- function(df) {
  df %>%
    filter(V1 > 0 & V1 < 23) %>%
    filter(str_detect(V2, "^r")) %>%
    filter(!(V5 == "G" & V6 == "C")) %>%
    filter(!(V5 == "C" & V6 == "G")) %>%
    filter(!(V5 == "T" & V6 == "A")) %>%
    filter(!(V5 == "A" & V6 == "T"))
}

data1_bim <- filter_snps(data1_bim) %>% arrange(V2)
data2_bim <- filter_snps(data2_bim) %>% arrange(V2)

common_snps <- intersect(data1_bim$V2, data2_bim$V2)

data1_bim <- data1_bim %>% filter(V2 %in% common_snps) %>% arrange(V2)
data2_bim <- data2_bim %>% filter(V2 %in% common_snps) %>% arrange(V2)

exact_match <- data1_bim$V1 == data2_bim$V1 &
  data1_bim$V2 == data2_bim$V2 &
  data1_bim$V4 == data2_bim$V4 &
  ((data1_bim$V5 == data2_bim$V5 & data1_bim$V6 == data2_bim$V6) |
     (data1_bim$V5 == data2_bim$V6 & data1_bim$V6 == data2_bim$V5))

final_snps <- data1_bim[exact_match, ]
unmatched_snps <- data1_bim[!exact_match, ]

write.table(unmatched_snps, file.path(MERGE_DIR, "unmatched_snps.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(final_snps$V2, file.path(MERGE_DIR, "snps_qc.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(final_snps[, .(V2, V5)], file.path(MERGE_DIR, "new_ref_alleles.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

system(glue("plink1.9 --bfile {MERGE_DIR}/Data1_filtered --extract {MERGE_DIR}/snps_qc.txt --make-bed --out {MERGE_DIR}/Data1_final"))
system(glue("plink1.9 --keep-allele-order --bfile {MERGE_DIR}/Data2_filtered --extract {MERGE_DIR}/snps_qc.txt --make-bed --out {MERGE_DIR}/Data2_final --reference-allele {MERGE_DIR}/new_ref_alleles.txt"))

#=====================================================#
# FINAL MERGE if there are no duplicate samples 
#====================================================#

system("plink1.9 --bfile Data_Merge/Data1_final --bmerge Data_Merge/Data2_final.bed Data_Merge/Data2_final.bim Data_Merge/Data2_final.fam --make-bed --out Data_Merge/Merged_Data1_Data2 --reference-allele Data_Merge/new_ref_alleles.txt")

#=========================#
#  REMOVE DUPLICATE IID   
#=========================#

# Load sample IDs, sample names, and call rates from both Data1 and Data2
data1_samples <- read_excel(EXCEL_FILE, sheet = "Data1") %>%
  select(IID = Sample_ID, Sample_Name, CallRate1 = Call_Rate)

data2_samples <- read_excel(EXCEL_FILE, sheet = "Data2") %>%
  select(IID = Sample_ID, Sample_Name, CallRate2 = Call_Rate)

# Find duplicates based on Sample_ID
dup_by_id <- inner_join(data1_samples, data2_samples, by = "IID")

# Find duplicates based on Sample_Name (not already found by IID)
dup_by_name <- inner_join(data1_samples, data2_samples, by = "Sample_Name") %>%
  filter(!(IID.x %in% dup_by_id$IID & IID.y %in% dup_by_id$IID))

# Combine duplicates
duplicates_combined <- bind_rows(
  dup_by_id %>%
    transmute(IID1 = IID, IID2 = IID, CR1 = CallRate1, CR2 = CallRate2),
  dup_by_name %>%
    transmute(IID1 = IID.x, IID2 = IID.y, CR1 = CallRate1, CR2 = CallRate2)
)

# Decide which to remove based on lower call rate
remove_from_data1 <- duplicates_combined %>%
  filter(CR1 < CR2) %>%
  pull(IID1) %>%
  unique()

remove_from_data2 <- duplicates_combined %>%
  filter(CR1 >= CR2) %>%  # Includes equal call rates
  pull(IID2) %>%
  unique()

# Write to files for plink removal
write.table(data.frame(FID = remove_from_data1, IID = remove_from_data1),
            file.path(MERGE_DIR, "dups_data1.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(data.frame(FID = remove_from_data2, IID = remove_from_data2),
            file.path(MERGE_DIR, "dups_data2.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


# Load original Data2 sheet
data2_info <- read_excel(EXCEL_FILE, sheet = "Data2") %>%
  select(Index, Sample_ID)

# Filter rows where Sample_ID is in the list of removed duplicates
remove_dups_data2_info <- data2_info %>%
  filter(Sample_ID %in% remove_from_data2)

# Write to text file: index and sample ID of removed duplicates
write.table(remove_dups_data2_info,
            file.path(MERGE_DIR, "remove_dups_data2.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Filter out duplicates with lower call rates
system(glue("plink1.9 --bfile {MERGE_DIR}/Data1_final --remove {MERGE_DIR}/remove_dups_data1.txt --make-bed --out {MERGE_DIR}/Data1_final_nodups"))
system(glue("plink1.9 --bfile {MERGE_DIR}/Data2_final --remove {MERGE_DIR}/remove_dups_data2.txt --make-bed --out {MERGE_DIR}/Data2_final_nodups"))

#==================================================#
# FINAL MERGE after duplicate removal
#=================================================#

system(glue("plink1.9 --bfile {MERGE_DIR}/Data1_final_nodups --bmerge {MERGE_DIR}/Data2_final_nodups.bed {MERGE_DIR}/Data2_final_nodups.bim {MERGE_DIR}/Data2_final_nodups.fam --make-bed --out {MERGE_DIR}/{DATA_NAME} --reference-allele {MERGE_DIR}/new_ref_alleles.txt"))

