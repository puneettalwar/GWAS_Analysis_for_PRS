#*********************************************************************************
# GWAS QC2 : Using 1000 genome data to match allele frequency & check population stratification
# R 4.3.1
# by Puneet 27.05.2025
#*********************************************************************************

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

#========================= CONFIGURATION =============================#

# # Tool and data paths
# TOOLS_DIR <- "/media/icaro/T7/CRC_Genetic_Data/Tools"
# PROJECT_DIR <- "/media/icaro/T7"
# QC_DIR <- "PD_PRS_GWAS_QC/COF/Data_QC"
# DATA_DIR <- "PD_PRS_GWAS_QC/COF"
# PATH_TO_GENOME1000 <- "/media/icaro/T7/CRC_Genetic_Data/common_files/1000genome"
# 
# # Sample name (prefix of binary PLINK files)
# DATA_NAME <- "COF_Data_94subs"

#=====================================================================#

# Set PATH for PLINK
Sys.setenv(PATH = paste(TOOLS_DIR, Sys.getenv("PATH"), sep = ":"))
# system("plink1.9 --help")
# system("plink2 --help")

# Set working directory
setwd(PROJECT_DIR)

# Create QC directories if missing
dir.create(file.path(DATA_DIR, "Data_merge_kgb"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(PROJECT_DIR, DATA_DIR))

#*********************************************************************************
# PCA Plot with 1000 genome data
#*********************************************************************************

# This is part of the second script of Pouya/Katya ["0_GENOME1000_QC"]
# All the required 1000 genome files are already present in "/media/icaro/T7/CRC_Genetic_Data/common_files/1000genome"

# library(ggplot2)
# library(GGally)

# Create QC directory if it doesn't exist
# if (!dir.exists("Data_merge_kgb")) dir.create("Data_merge_kgb")

# Load data
pcs <- read.table(file.path(PATH_TO_GENOME1000, "plink2.eigenvec"), header = FALSE)
pop <- read.table(file.path(PATH_TO_GENOME1000, "integrated_call_samples_v3.20130502.ALL.panel"), header = TRUE)

# Sorting and resetting index
pop <- pop %>% arrange(sample)
pcs <- pcs %>% arrange(V2)  # Assuming V2 is the column with sample ID

# Plot with base ggplot2: PC1 vs PC2 colored by population
df <- data.frame(pc1 = pcs[,3], pc2 = pcs[,4], pop = pop$pop,super_pop = pop$super_pop)

plot_kgb_pop <- ggplot(df, aes(x = pc1, y = pc2, color = pop)) +
  geom_point(size = 1.5) +
  theme_minimal()

ggsave("Data_merge_kgb/plot_kgb_pop.pdf", plot = plot_kgb_pop)

# Pair plot equivalent: PC1 vs PC2
ggpairs(df, columns = 1:2, mapping = aes(color = pop))

# 2nd against 3rd PC
df1 <- data.frame(pc2 = pcs[,3], pc3 = pcs[,4], pop = pop$pop)
ggpairs(df1, columns = 1:2, mapping = aes(color = pop))


# Plot with base ggplot2: PC1 vs PC2 colored by super population

plot_kgb_super_pop <- ggplot(df, aes(x = pc1, y = pc2, color = super_pop)) +
  geom_point(size = 1.5) +
  theme_classic2()

ggsave("Data_merge_kgb/plot_kgb_super_pop.pdf", plot = plot_kgb_super_pop)

#*********************************************************************************
# Merging 1000 genome data with processed data after QC 
#*********************************************************************************
# This is the third script of Pouya/Katya ["1_KGB+GWA"]

# DATA_DIR  <- "/media/icaro/T7/PD_PRS_GWAS_QC"
# setwd(DATA_DIR )

# library(data.table)
# library(dplyr)
# library(purrr)

# Create QC directory if it doesn't exist
# if (!dir.exists("Data_merge_kgb")) dir.create("Data_merge_kgb")
# Read the .bim file (no header)
qc1_final_bim <- read.table(sprintf("Data_QC/%s_final_qc1_complete.bim", DATA_NAME),header = FALSE)
dim (qc1_final_bim)
# 508811 
head(qc1_final_bim)

# Save the second column (SNP IDs) to a text file
write.table(qc1_final_bim$V2,"Data_merge_kgb/snp_qc1_final.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

system(sprintf("plink1.9 --vcf %s/kgb.vcf.gz --extract Data_merge_kgb/snp_qc1_final.txt --make-bed --out Data_merge_kgb/kgb_filtered", PATH_TO_GENOME1000))

# Load .bim file
kgb_bim <- fread("Data_merge_kgb/kgb_filtered.bim", header = FALSE)
head(kgb_bim)

# Remove duplicated SNPs
# kgb_bim <- kgb_bim[!duplicated(kgb_bim[[2]])]
kgb_bim <- kgb_bim[!duplicated(kgb_bim$V2)]

# Save SNP IDs to file
write.table(kgb_bim$V2, "Data_merge_kgb/snps_common.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Run PLINK to extract common SNPs
system(sprintf("plink1.9 --bfile Data_QC/%s_final_qc1_complete --extract Data_merge_kgb/snps_common.txt --make-bed --out Data_merge_kgb/qc1_final_QC_common --noweb",DATA_NAME))

# Load filtered BIM files
# kgb_bim <- fread(file.path(DATA_DIR , "kgb_filtered.bim"), header = FALSE)
# gwa_bim <- fread(file.path(DATA_DIR , "qc1_final_QC_common.bim"), header = FALSE)

kgb_bim <- fread("Data_merge_kgb/kgb_filtered.bim", header = FALSE)
qc1_final_bim <- fread("Data_merge_kgb/qc1_final_QC_common.bim", header = FALSE)

# Remove unmatched SNPs
common_snps <- qc1_final_bim$V2
kgb_bim <- kgb_bim[kgb_bim$V2 %in% common_snps]

# Sort and reset
kgb_bim <- kgb_bim[order(kgb_bim$V2), ]
qc1_final_bim <- qc1_final_bim[order(qc1_final_bim$V2), ]

# Check allele alignment
allele_match <- qc1_final_bim$V5 == kgb_bim$V5 & qc1_final_bim$V6 == kgb_bim$V6
allele_swap  <- qc1_final_bim$V5 == kgb_bim$V6 & qc1_final_bim$V6 == kgb_bim$V5

flip_function <- function(allele) {
  allele <- toupper(allele)
  return(switch(allele,
                "A" = "T",
                "T" = "A",
                "C" = "G",
                "G" = "C",
                "N"))  # default return for unrecognized input
}

allele_flip <- qc1_final_bim$V5 == map_chr(kgb_bim$V5, flip_function) &
  qc1_final_bim$V6 == map_chr(kgb_bim$V6, flip_function)

allele_flip_swap <- qc1_final_bim$V5 == map_chr(kgb_bim$V6, flip_function) &
  qc1_final_bim$V6 == map_chr(kgb_bim$V5, flip_function)

# Count
sum(allele_match)
sum(allele_swap)
sum(allele_flip)
sum(allele_flip_swap)

# SNPs to exclude
valid <- allele_match | allele_swap | allele_flip | allele_flip_swap
excludedSNP <- qc1_final_bim[!valid, ]

write.table(excludedSNP$V2, "Data_merge_kgb/snps_exclude.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# SNPs to flip
flipped_snps <- qc1_final_bim[allele_flip | allele_flip_swap, ][[2]]

write.table(flipped_snps, "Data_merge_kgb/snps_flipped.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Apply exclusions and flipping via PLINK
system("plink1.9 --bfile Data_merge_kgb/kgb_filtered --exclude Data_merge_kgb/snps_exclude.txt --make-bed --out Data_merge_kgb/kgb_filtered_common")
system("plink1.9 --bfile Data_merge_kgb/qc1_final_QC_common --exclude Data_merge_kgb/snps_exclude.txt --flip Data_merge_kgb/snps_flipped.txt --make-bed --out Data_merge_kgb/qc1_final_QC_common_flipped")

# Reference allele file
kgb_bim <- fread("Data_merge_kgb/kgb_filtered_common.bim", header = FALSE)
ref_allele <- data.frame(SNP = kgb_bim$V2, REF = kgb_bim$V5)

write.table(ref_allele, "Data_merge_kgb/snps_kgb_reference.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# write.table(ref_allele, file = file.path(PROJECT_DIR, DATA_DIR , "Data_merge_kgb/snps_kgb_reference.txt"),
#             quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Final PLINK call
system("plink1.9 --bfile Data_merge_kgb/qc1_final_QC_common_flipped --reference-allele Data_merge_kgb/snps_kgb_reference.txt --make-bed --out Data_merge_kgb/qc1_final_QC_common_flipped_ref")
# 592636

##################### Read and Merge Population and .fam Data #####################

# DATA_DIR  <- "/media/icaro/T7/PD_PRS_GWAS_QC"
#setwd(PROJECT_DIR,DATA_DIR )

# Read population panel
pop_merge <- read.table(sprintf("%s/integrated_call_samples_v3.20130502.ALL.panel", PATH_TO_GENOME1000), header = FALSE)
#head(pop_merge)

# Rename pop_merge columns to match content
setnames(pop_merge, c("sample", "pop", "super_pop", "gender"))
#head(pop_merge)
#dim(pop_merge)

# Read fam file
qc1_fam <- fread("Data_merge_kgb/qc1_final_QC_common_flipped_ref.fam", header = FALSE)
#qc1_fam <- fread(file.path(PROJECT_DIR,DATA_DIR , "Data_merge_kgb/qc1_final_QC_common_flipped_ref.fam"), header = FALSE)
#head(gwa_fam)

# Determine number of samples in .fam file
num_samples <- nrow(qc1_fam)

# Create data frame for new samples
pop_data_fam <- data.frame(
  sample = qc1_fam$V2,
  pop = rep("diff", num_samples),
  super_pop = rep("Eur", num_samples)
)
#head(pop_data_fam)

# Add missing columns to pop_data_fam
pop_data_fam$gender <- NA

# Ensure columns match and are in the same order
pop_data_fam <- pop_data_fam[, names(pop_merge)]
#head(pop_data_fam)
#dim(pop_data_fam)

# Merge with 1000 Genome population info
pop_merge <- rbind(pop_merge, pop_data_fam)
#head(pop_merge)
#dim(pop_merge)

##################### Subset for Europe and Write "europe.txt" #####################

eur <- pop_merge %>% filter(super_pop == "EUR")
# dim(eur)
europe_txt <- eur %>% select(sample) %>% mutate(sample2 = sample)
fwrite(europe_txt, "Data_merge_kgb/europe.txt", sep = "\t", col.names = FALSE)
#fwrite(europe_txt, file.path(PROJECT_DIR,DATA_DIR , "Data_merge_kgb/europe.txt"), sep = "\t", col.names = FALSE)

system("plink1.9 --bfile Data_merge_kgb/kgb_filtered_common --keep Data_merge_kgb/europe.txt --make-bed --out Data_merge_kgb/europe_subset --keep-allele-order")
system("plink1.9 --bfile Data_merge_kgb/europe_subset --freq --out Data_merge_kgb/plink_eur --keep-allele-order")
system("plink1.9 --bfile Data_merge_kgb/kgb_filtered_common --freq --out Data_merge_kgb/plink_kgb --keep-allele-order")
system("plink1.9 --bfile Data_merge_kgb/qc1_final_QC_common_flipped_ref --out Data_merge_kgb/plink_qc1_New --freq --keep-allele-order --noweb")

##################### Load Allele Frequencies and Sort  #####################

qc1_new_freq <- fread(file.path(PROJECT_DIR, DATA_DIR, "Data_merge_kgb/plink_qc1_New.frq"))
kgb_freq <- fread(file.path(PROJECT_DIR, DATA_DIR, "Data_merge_kgb/plink_kgb.frq"))
eur_freq <- fread(file.path(PROJECT_DIR, DATA_DIR, "Data_merge_kgb/plink_eur.frq"))

qc1_new_freq <- qc1_new_freq %>% arrange(desc(SNP))
eur_freq <- eur_freq %>% arrange(desc(SNP))
kgb_freq <- kgb_freq %>% arrange(desc(SNP))

##################### Compare Alleles ##################### 

sum(qc1_new_freq$A1 == eur_freq$A1)
sum(qc1_new_freq$A2 == eur_freq$A2)

sum(qc1_new_freq$A1 == kgb_freq$A1)
sum(qc1_new_freq$A2 == kgb_freq$A2)

##################### Plot MAF Distributions ##################### 

maf_plot <- ggplot() +
  geom_density(aes(qc1_new_freq$MAF, color = "final cohort")) +
  geom_density(aes(kgb_freq$MAF, color = "Genome 1000")) +
  geom_density(aes(eur_freq$MAF, color = "Europe subset")) +
  labs(x = "MAF", y = "Density") +
  theme_classic2() +
  scale_color_manual(values = c("blue", "green", "red")) +
  theme(legend.title = element_blank())

ggsave("Data_merge_kgb/plot_maf_kgb+data.pdf", plot = maf_plot, width = 8, height = 6)

##################### Joint Plots (Correlation Scatterplots) #####################

#library(GGally)

allele_freq <- data.frame(Data_final = qc1_new_freq$MAF, Europe = eur_freq$MAF)
plot_freq_eur <- GGally::ggpairs(allele_freq)
ggsave("Data_merge_kgb/plot_frq_eur+data.pdf", plot = plot_freq_eur)

allele_freq2 <- data.frame(Data_final = qc1_new_freq$MAF, Genome1000 = kgb_freq$MAF)
plot_freq_kgb <- GGally::ggpairs(allele_freq2)
ggsave("Data_merge_kgb/plot_frq_kgb+data.pdf", plot = plot_freq_kgb)

##################### Outlier Removal #####################

outliers <- eur_freq[(qc1_new_freq$MAF > eur_freq$MAF + 0.2) | 
                       (qc1_new_freq$MAF < eur_freq$MAF - 0.2), ]

eur_freq_new <- eur_freq[!(qc1_new_freq$MAF > eur_freq$MAF + 0.2 |
                             qc1_new_freq$MAF < eur_freq$MAF - 0.2), ]
qc1_new_freq_new <- qc1_new_freq[!(qc1_new_freq$MAF > eur_freq$MAF + 0.2 |
                                                   qc1_new_freq$MAF < eur_freq$MAF - 0.2), ]

# Save outliers

write.table(outliers$SNP, 
            file = file.path(PROJECT_DIR,DATA_DIR, "Data_merge_kgb/snps_OutliersRemoval.txt"), 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)


system("plink1.9 --bfile Data_merge_kgb/kgb_filtered_common --exclude Data_merge_kgb/snps_OutliersRemoval.txt --make-bed --out Data_merge_kgb/kgb_filtered_common_outlier")
system("plink1.9 --bfile Data_merge_kgb/europe_subset --exclude Data_merge_kgb/snps_OutliersRemoval.txt --reference-allele Data_merge_kgb/snps_kgb_reference.txt --make-bed --out Data_merge_kgb/europe_subset_outlier")
system("plink1.9 --bfile Data_merge_kgb/qc1_final_QC_common_flipped_ref --exclude Data_merge_kgb/snps_OutliersRemoval.txt --reference-allele Data_merge_kgb/snps_kgb_reference.txt --make-bed --out Data_merge_kgb/final_QC2_common_flipped_ref_outlier")

# system("plink2 --bfile Data_merge_kgb/final_QC2_common_flipped_ref_outlier --ref-from-fa --fa /media/icaro/T7/CRC_Genetic_Data/common_files/1000genome/human_g1k_v37.fasta --make-bed --out Data_merge_kgb/final_QC2_common_flipped_ref_outlier_ref")
# correcting the REF/ALT swaps using bcftools [same as using Plink2]
## system("bcftools +fixref Data_vcf_check/final_QC2_common_flipped_ref_outlier.vcf.gz -Ov -o Data_vcf_check/test_swapped2.vcf.gz -- -d -f /media/icaro/T7/CRC_Genetic_Data/common_files/1000genome/human_g1k_v37.fasta -m flip")

##################### Merge Bfiles #####################

system("plink1.9 --bfile Data_merge_kgb/kgb_filtered_common_outlier --bmerge Data_merge_kgb/final_QC2_common_flipped_ref_outlier.bed Data_merge_kgb/final_QC2_common_flipped_ref_outlier.bim Data_merge_kgb/final_QC2_common_flipped_ref_outlier.fam --make-bed --out Data_merge_kgb/qc1_kgb_merged_final --reference-allele Data_merge_kgb/snps_kgb_reference.txt")
system("plink2 --bfile Data_merge_kgb/qc1_kgb_merged_final --pca --out Data_merge_kgb/qc1_kgb_merged_final --keep-allele-order")

#*****************************************************************************************
# Freq check after outlier removal
#*****************************************************************************************
system("plink1.9 --bfile Data_merge_kgb/final_QC2_common_flipped_ref_outlier --out Data_merge_kgb/final_QC2_common_flipped_ref_outlier_freq --freq --keep-allele-order")
system("plink1.9 --bfile Data_merge_kgb/europe_subset_outlier --out Data_merge_kgb/europe_subset_outlier_freq --freq --keep-allele-order")

qc1_freq_nooutlier <- fread(file.path(PROJECT_DIR,DATA_DIR , "Data_merge_kgb/final_QC2_common_flipped_ref_outlier_freq.frq"))
eur_freq_outlier <- fread(file.path(PROJECT_DIR,DATA_DIR , "Data_merge_kgb/europe_subset_outlier_freq.frq"))

qc1_freq_nooutlier <- qc1_freq_nooutlier %>% arrange(desc(SNP))
eur_freq_outlier <- eur_freq_outlier %>% arrange(desc(SNP))

sum(qc1_freq_nooutlier$A1 == eur_freq_outlier$A1)
sum(qc1_freq_nooutlier$A2 == eur_freq_outlier$A2)

##################### Plot MAF Distributions ##################### 

maf_plot <- ggplot() +
  geom_density(aes(qc1_freq_nooutlier$MAF, color = "final cohort")) +
  geom_density(aes(kgb_freq$MAF, color = "Genome 1000")) +
  geom_density(aes(eur_freq_outlier$MAF, color = "Europe subset")) +
  labs(x = "MAF", y = "Density") +
  theme_classic2() +
  scale_color_manual(values = c("blue", "green", "red")) +
  theme(legend.title = element_blank())

ggsave("Data_merge_kgb/plot_maf_kgb+data_outlier.pdf", plot = maf_plot, width = 8, height = 6)

#dev.off()

##################### Joint Plots (Correlation Scatterplots) #####################

#library(ggpubr)

allele_freq <- data.frame(Data_final = qc1_freq_nooutlier$MAF, Europe = eur_freq_outlier$MAF)

# Create scatter plot with correlation coefficient
plot_freq_eur <- ggscatter(allele_freq, x = "Europe", y = "Data_final",
                           add = "reg.line",         # Add regression line
                           conf.int = TRUE,          # Add confidence interval
                           cor.coef = TRUE,          # Add correlation coefficient
                           cor.method = "pearson",   # Can be "pearson", "spearman", or "kendall"
                           xlab = "MAF in Europe", 
                           ylab = "MAF in Data_final",
                           color = "blue") +
  theme_pubr()

ggsave("Data_merge_kgb/plot_frq_eur+data_outlier.pdf", plot = plot_freq_eur)
#dev.off()

# allele_freq2 <- data.frame(Data_final = qc1_freq_nooutlier$MAF, Genome1000 = kgb_freq$MAF)
# plot_freq_kgb <- GGally::ggpairs(allele_freq2)
# ggsave("Data_merge_kgb/plot_frq_kgb+data_outlier.pdf", plot = plot_freq_kgb)

##################### Final Population Merge and PCA Plot ##################### 

# Read population panel
pop_merge <- read.table(sprintf("%s/integrated_call_samples_v3.20130502.ALL.panel", PATH_TO_GENOME1000), header = TRUE)
#head(pop_merge)

# Rename pop_merge columns to match content
setnames(pop_merge, c("sample", "pop", "super_pop", "gender"))
head(pop_merge)
dim(pop_merge)

qc1_fam <- fread(file.path(PROJECT_DIR,DATA_DIR , "Data_merge_kgb/final_QC2_common_flipped_ref_outlier.fam"), header = FALSE)
head(qc1_fam)
dim(qc1_fam)

# Add gender column (NA for now), and build pop_fam data.table
pop_fam <- data.table(sample = qc1_fam$V2, 
                      pop = "OLD", 
                      super_pop = "OLD", 
                      gender = NA)


# Combine using rbindlist with fill=TRUE
# pop_merge <- rbind(pop_merge, pop_fam, fill=TRUE)
pop_merge <- rbindlist(list(pop_merge, pop_fam), fill = TRUE)
dim(pop_merge)
head(pop_merge)

pop_merge$pop[!(pop_merge$super_pop %in% c("EUR", "OLD"))] <- "OTHERS"

pop_merge <- pop_merge %>% arrange(desc(sample))
head(pop_merge)
dim(pop_merge)
#2504

pcs_merged <- fread(sprintf("%s/plink2.eigenvec", PATH_TO_GENOME1000), header = TRUE)
head(pcs_merged)
setnames(pcs_merged, "#FID", "sample")

pcs_merged <- pcs_merged %>% arrange(desc(sample))
dim(pcs_merged)
#2504

names(pcs_merged)
names(pop_merge)

# Merge PC and population
df_merge <- merge(pcs_merged, pop_merge, by = "sample")
head(df_merge)
dim(df_merge)
#2504

plot_kgb_data_superpop <- ggplot(df_merge, aes(x = PC1, y = PC2, color = super_pop)) +
  geom_point(alpha = 0.7) +
  theme_classic2() +
  labs(title = "PCA kgb+data by Super Population")

ggsave("Data_merge_kgb/plot_kgb_data_superpop.pdf", plot = plot_kgb_data_superpop)

plot_kgb_data_pop <- ggplot(df_merge, aes(x = PC1, y = PC2, color = pop)) +
  geom_point() +
  theme_classic2() +
  labs(title = "PCA kgb+data by Population")

ggsave("Data_merge_kgb/plot_kgb_data_pop.pdf", plot = plot_kgb_data_pop)
ggsave("Data_merge_kgb/plot_kgb_data_pop.png", plot = plot_kgb_data_pop, width = 8, height = 6)

cat("\n FINAL QC2 COMPLETE!\nYour file is ready for preimputation QC: Run GWAS_QC3_Preimputation\n")

#*********************************************************************************
# QC2 COMPLETE !!! 
# RUN GWAS_QC3_Preimputation
#*********************************************************************************