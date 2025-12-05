#*********************************************************************************
# GWAS QC4 : Post imputation QC
# R 4.3.1
# Updated: 2025-07-15
#*********************************************************************************

# Load R packages
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(glue)

#========================= CONFIGURATION =============================#

# TOOLS_DIR <- "/media/icaro/T7/CRC_Genetic_Data/Tools"
# PROJECT_DIR <- "/media/icaro/T7"
# DATA_DIR <- "/media/icaro/T7/PD_PRS_GWAS_QC/COF"
# COMMON_FILES <- "/media/icaro/T7/CRC_Genetic_Data/common_files"
# DATASET_NAME <- "COF"  # Change for each dataset
# 
# VCF_DIR <- file.path(DATA_DIR, "Data_Imputed")
# OUTPUT_DIR <- file.path(DATA_DIR, "Post_imputation_QC")
# HRC_VCF <- file.path(COMMON_FILES, "HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz")

#=====================================================================#

# Set environment
Sys.setenv(PATH = paste(TOOLS_DIR, Sys.getenv("PATH"), sep = ":"))

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUTPUT_DIR)

#***************************************************************
# unzip Michigan output files
#***************************************************************
# Use bash terminal to unzip Michigan output files 
#***************************************************************
# cd to the folder "/media/icaro/T7/PD_PRS_GWAS_QC/COF/Data_Imputed"
# unzip Michigan output files using the password received from Michigan server in Email
# Use 7zip (as recommended) [sudo apt install p7zip-full p7zip-rar]

# Password : XChk02xP}?bYRi
# 
# for chr in $(seq 1 22)
# do
# 7z e chr_$chr.zip
# done

#***************************************************************
# Concatenate and QC VCFs
#***************************************************************
# Concatenate VCFs [Takes ~15 mins for N=90]
system(glue("vcf-concat {VCF_DIR}/*dose.vcf.gz | bgzip -c > {OUTPUT_DIR}/{DATASET_NAME}_Imputed.vcf.gz"))

# Preprocessing using PLINK2
# imputation quality is stored in the R2 field, not in a generic INFO field
# system("plink2 --vcf Post_imputation_QC/COF_Imputed.vcf.gz --geno 0.05 --maf 0.01 --hwe 1E-4 --exclude-if-info 'INFO < 0.3' --recode vcf bgz --out Post_imputation_QC/COF_Imputed_QC")

system(glue("
plink2 --vcf {OUTPUT_DIR}/{DATASET_NAME}_Imputed.vcf.gz \\
       --geno 0.05 --maf 0.01 --hwe 1E-4 \\
       --exclude-if-info 'R2 < 0.3' \\
       --recode vcf bgz \\
       --out {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC"))

# Index the file with tabix
system(glue("tabix -p vcf {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC.vcf.gz"))

# Annotate with CHROM:POS 
system(glue("
bcftools annotate -I '%CHROM:%POS' \\
    {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC.vcf.gz \\
    -o {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr.vcf.gz -Oz"))

#***************************************************************
# Allele Frequency Comparison with HRC
#***************************************************************

# Allele frequency calculation
system(glue("
plink2 --vcf {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr.vcf.gz \\
       --freq --keep-allele-order \\
       --out {OUTPUT_DIR}/plink2_{DATASET_NAME}"))

# Remove duplicates
freq <- fread(glue("{OUTPUT_DIR}/plink2_{DATASET_NAME}.afreq"))
colnames(freq) <- c("#CHROM", "SNP", "REF", "ALT", "Data_ALT_FREQS", "OBS_CT")
freq <- freq[!duplicated(SNP)]

# Read HRC summary, create SNP column and remove duplicates
hrc_summary <- fread(cmd = paste("zcat", HRC_VCF), skip = 49)
hrc_summary[, SNP := paste0(`#CHROM`, ":", POS)]
hrc_summary <- hrc_summary[!duplicated(SNP)]

# Filter and parse HRC AF
hrc_filtered <- hrc_summary %>%
  filter(SNP %in% freq$SNP) %>%
  filter(!is.na(INFO)) %>%
  mutate(AF = str_match(INFO, "AF=([0-9eE\\.\\-]+)")[,2] %>% as.numeric())

# Merge and QC
mergedFreq <- inner_join(freq, hrc_filtered, by = c("SNP", "REF", "ALT"))
outliers <- mergedFreq %>% filter(abs(Data_ALT_FREQS - AF) > 0.2)
mergedFreq_clean <- mergedFreq[!(Data_ALT_FREQS > AF + 0.2 | Data_ALT_FREQS < AF - 0.2),]

#*******************************************************************
# Plotting before and after cleaning based on allele freq difference
#*******************************************************************

plot_AF_density <- ggplot(mergedFreq, aes(x = Data_ALT_FREQS)) +
  geom_density(aes(color = "Dataset")) +
  geom_density(aes(x = AF, color = "HRC")) +
  labs(title = "Allele Frequency Distribution", x = "Allele Frequency") +
  theme_classic2()
ggsave(glue("{OUTPUT_DIR}/Allele_Freq_Dist.png"), plot_AF_density)

plot_AF_hist <- ggplot(mergedFreq, aes(x = Data_ALT_FREQS)) +
  geom_histogram(binwidth = 0.01, fill = "blue", alpha = 0.5, aes(y = after_stat(count))) +
  geom_histogram(aes(x = AF), binwidth = 0.01, fill = "red", alpha = 0.5) +
  labs(title = "Allele Frequency Histogram", x = "Allele Frequency") +
  theme_classic2()
ggsave(glue("{OUTPUT_DIR}/plot_AF_hist.png"), plot_AF_hist)

plot_AF_scatter <- ggplot(mergedFreq, aes(x = Data_ALT_FREQS, y = AF)) +
  geom_point(color = "#1f77b4") +
  labs(title = "AF Comparison", x = "Dataset", y = "HRC") +
  theme_classic2()
ggsave(glue("{OUTPUT_DIR}/plot_AF_scatter.png"), plot_AF_scatter)

# Cleaned plots
AF_cleaned_density <- ggplot(mergedFreq_clean, aes(x = Data_ALT_FREQS)) +
  geom_density(aes(color = "Data")) +
  geom_density(aes(x = AF, color = "HRC")) +
  labs(title = "Cleaned Allele Frequency Distribution") +
  theme_classic2()
ggsave(glue("{OUTPUT_DIR}/AF_cleaned_density.png"), AF_cleaned_density)

AF_cleaned_scatter <- ggplot(mergedFreq_clean, aes(x = Data_ALT_FREQS, y = AF)) +
  geom_point(color = "#1f77b4") +
  labs(x = "Dataset ALT Freqs", y = "HRC ALT Freqs") +
  theme_classic2()
ggsave(glue("{OUTPUT_DIR}/AF_cleaned_scatter.png"), AF_cleaned_scatter)

#***************************************************************
# Final Dataset Extraction
#***************************************************************
write.table(mergedFreq_clean$SNP, file = glue("{OUTPUT_DIR}/snps_HRC_extract.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Final extraction [Number is higher due to multi-allelic SNPs]
system(glue("
plink2 --vcf {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr.vcf.gz \\
       --extract {OUTPUT_DIR}/snps_HRC_extract.txt \\
       --recode vcf bgz \\
       --out {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr_filter_HRC"))

system(glue("bcftools index {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr_filter_HRC.vcf.gz"))

# Add rsid column
system(glue("
bcftools annotate -a {HRC_VCF} -c CHROM,POS,REF,ALT,ID -Oz \\
  -o {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr_filter_HRC_rsid.vcf.gz \\
  {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr_filter_HRC.vcf.gz"))

system(glue("bcftools index {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr_filter_HRC_rsid.vcf.gz"))

# Create Plink files
system(glue("
plink1.9 --vcf {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr_filter_HRC_rsid.vcf.gz \\
         --make-bed --out {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr_filter_HRC_rsid"))

# Remove multi-allelic SNPs
# Check for multi-allelic SNPs
system(glue("
bcftools query -f '%ID\\n' {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr_filter_HRC_rsid.vcf.gz | \\
sort | uniq -d > {OUTPUT_DIR}/multi-allelic_rsid.txt"))

system(glue("
plink1.9 --bfile {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr_filter_HRC_rsid \\
         --exclude {OUTPUT_DIR}/multi-allelic_rsid.txt \\
         --make-bed --out {OUTPUT_DIR}/{DATASET_NAME}_Imputed_QC_chr_clean_HRC_rsid"))

cat("\n✔ FINAL QC4 COMPLETE!\nYour files are ready for PRS computation in the Post_imputation_QC folder\n")

#*********************************************************************************
# QC4 COMPLETE !!! 
# RUN PRS computation
#*********************************************************************************