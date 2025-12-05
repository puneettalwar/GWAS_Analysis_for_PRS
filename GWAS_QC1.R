#*********************************************************************************
# GWAS QC1 : Basic QC after data merging and filtering
# R 4.3.1
#*********************************************************************************

# Required installations (system): Samtools, BCFtools, HTSlib
# URLs:
# https://www.htslib.org/download/
# https://gist.github.com/adefelicibus/f6fd06df1b4bb104ceeaccdd7325b856

# install.packages(c("tidyverse","rstudioapi","data.table","plinkQC","getopt","ggplot2","GGally","dplyr","purrr","ggpubr"))

# Load libraries
library(rstudioapi) 
library(tidyverse)
library(data.table)
library(plinkQC)
library(ggplot2)
library(ggpubr)
library(readxl)

#************************* Modify only configuration part & Run****************#

#========================= CONFIGURATION =============================#

# # Tool and data paths
# TOOLS_DIR <- "/media/icaro/T7/CRC_Genetic_Data/Tools"
# PROJECT_DIR <- "/media/icaro/T7"
# QC_DIR <- "PD_PRS_GWAS_QC/COF/Data_QC"
# DATA_DIR <- "PD_PRS_GWAS_QC/COF"
# 
# # Sample name (prefix of binary PLINK files)
# DATA_NAME <- "COF_Data"

#=====================================================================#

# Set PATH for PLINK
Sys.setenv(PATH = paste(TOOLS_DIR, Sys.getenv("PATH"), sep = ":"))

# Set working directory
setwd(PROJECT_DIR)

# Create QC directories if missing
dir.create(file.path(DATA_DIR, "Data_QC"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(PROJECT_DIR, DATA_DIR))


############  QC steps [Missingness, Genotyping rate, MAF & HWE]

# Perform QC filtering steps
system(sprintf("plink1.9 --bfile Data_merge/%s --mind 0.2 --make-bed --out Data_QC/%s_qc_mind", DATA_NAME, DATA_NAME))
system(sprintf("plink1.9 --bfile Data_QC/%s_qc_mind --geno 0.05 --maf 0.01 --hwe 1E-4 --make-bed --out Data_QC/%s_filtered", DATA_NAME, DATA_NAME))

################### Check cryptic Relatedness ######################
# LD pruning
system(sprintf("plink1.9 --bfile Data_QC/%s_filtered --indep-pairphase 1000 10 0.2 --out Data_QC/%s_filtered_pruned", 
               DATA_NAME, DATA_NAME))

# Relatedness check
system(sprintf("plink1.9 --bfile Data_QC/%s_filtered --extract Data_QC/%s_filtered_pruned.prune.in --genome --min 0.2 --make-bed --out Data_QC/%s_filtered_pruned_related", 
               DATA_NAME, DATA_NAME, DATA_NAME))

genome <- read.table(sprintf("Data_QC/%s_filtered_pruned_related.genome", DATA_NAME), header = TRUE)
relatives <- genome %>% filter(PI_HAT > 0.5) # Katya used 0.5 and Pouya used 0.25 [It should be 0.25]
#head(relatives)
write.table(relatives, sprintf("Data_QC/relatives_PIHAT_0.5.txt"), row.names = FALSE, quote = FALSE)

# Read .fam
sample_info <- read.table(sprintf("Data_QC/%s_filtered_pruned_related.fam", DATA_NAME), header = FALSE)
colnames(sample_info) <- c("FID", "IID", "PID", "MID", "SEX", "PHENOTYPE")

# Join sex info for relatives
relatives_with_sex <- relatives %>%
  left_join(sample_info %>% select(FID, IID, SEX), by = c("FID1" = "FID", "IID1" = "IID")) %>%
  rename(SEX1 = SEX) %>%
  left_join(sample_info %>% select(FID, IID, SEX), by = c("FID2" = "FID", "IID2" = "IID")) %>%
  rename(SEX2 = SEX)

# Save related individuals (remove one of each pair)
pihat_related <- data.frame(
  FID = relatives$FID1,
  IID = relatives$IID1
)

write.table(pihat_related, "Data_QC/pihat_pruned_related.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

########### Plot cryptic relatedness ######################

# PCA for visualization
# Note: genome output is not used for pca but only for relatedness check
system(sprintf("plink1.9 --bfile Data_QC/%s_filtered --pca --out Data_QC/%s_filtered", DATA_NAME, DATA_NAME))

eigenvec <- read.table(sprintf("Data_QC/%s_filtered.eigenvec", DATA_NAME), header = FALSE)
colnames(eigenvec) <- c("FID", "IID", paste0("PC", 1:(ncol(eigenvec) - 2)))

# Flag related individuals
rel_IIDs <- unique(c(relatives$IID1, relatives$IID2))
eigenvec$related <- ifelse(eigenvec$IID %in% rel_IIDs, "Related", "Unrelated")

# Plot PCA
plot_related <- ggplot(eigenvec, aes(x = PC1, y = PC2, color = related)) +
  geom_point() +
  theme_classic2() +
  scale_color_manual(values = c("Unrelated" = "grey", "Related" = "red")) +
  labs(title = "PCA Plot with Relatives Highlighted", x = "PC1", y = "PC2")

ggsave("Data_QC/plot_related.pdf", plot = plot_related)

#############  Check Duplicates  (by high genetic relatedness) #############

# Check .genome output file for pihat>0.9 [Don't use output file in the next step]

system(sprintf("plink1.9 --bfile Data_QC/%s_filtered --extract Data_QC/%s_filtered_pruned.prune.in --genome --min 0.9 --make-bed --out Data_QC/%s_filtered_pruned_duplicates", 
               DATA_NAME, DATA_NAME, DATA_NAME))

genome <- read.table(sprintf("Data_QC/%s_filtered_pruned_duplicates.genome", DATA_NAME), header = TRUE)
duplicates <- genome %>% filter(PI_HAT > 0.9)
write.table(duplicates, "Data_QC/relatives_PIHAT_0.9.txt", row.names = FALSE, quote = FALSE)

################## Final file after QC ######################

system(sprintf("plink1.9 --bfile Data_QC/%s_filtered --remove Data_QC/pihat_pruned_related.txt --make-bed --out Data_QC/%s_final_qc1_complete", 
               DATA_NAME, DATA_NAME))
# 596198 variants and 92 people 

#========================= plinkQC Section ==============================#

dir.create("Data_QC/plinkQC", showWarnings = FALSE)

package.dir <- find.package('plinkQC')
path2plink <- file.path(TOOLS_DIR, "plink1.9")
indir <- file.path(PROJECT_DIR, DATA_DIR, "Data_QC")
qcdir <- file.path(indir,"plinkQC")
name <- paste0(DATA_NAME, "_final_qc1_complete")

source(file.path(TOOLS_DIR, "../Codes/Puneet/plinkQC.R"))

# Optionally run these lines interactively:
fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name, path2plink=path2plink,
                            verbose=TRUE, interactive=TRUE, showPlinkOutput=FALSE)
overview_marker <- overviewPerMarkerQC(fail_markers, interactive = TRUE)

dev.cur()
while (!is.null(dev.list())) dev.off()

cat("\n FINAL QC2 COMPLETE!\nYour file is ready for QC2: Run GWAS_QC2_1000G \n")

#*********************************************************************************
# QC1 COMPLETE !!! 
# NEXT: RUN GWAS_QC2_1000G SCRIPT
#*********************************************************************************
