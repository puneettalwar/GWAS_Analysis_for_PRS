library(shiny)
library(readxl)
library(dplyr)
library(stringr)
library(data.table)
library(glue)

options(shiny.maxRequestSize = 1000 * 1024^2)  # Allow up to 1000 MB upload
options(shiny.save.state = FALSE)

ui <- fluidPage(
  titlePanel("GWAS Data Harmonization Pipeline: Merging and Filtering"),
  sidebarLayout(
    sidebarPanel(
      textInput("project_dir", "Choose Project directory", value = NULL),
      textInput("tools_dir", "Path to Tools Folder (with plink executable):", value = NULL),      
      h4("Data1 Input Type"),
      radioButtons("data1_type", "Choose Data1 format:", choices = c("PED/MAP", "BED/BIM/FAM"), selected = "PED/MAP"),
      conditionalPanel(
        condition = "input.data1_type == 'PED/MAP'",
        fileInput("data1_ped", "Upload Data1 .ped"),
        fileInput("data1_map", "Upload Data1 .map")
      ),
      conditionalPanel(
        condition = "input.data1_type == 'BED/BIM/FAM'",
        fileInput("data1_bed", "Upload Data1 .bed"),
        fileInput("data1_bim", "Upload Data1 .bim"),
        fileInput("data1_fam", "Upload Data1 .fam")
      ),
      fileInput("sample1", "Data1 Sample Metadata (.xlsx)"),
      checkboxInput("sex_check_data1", "Run sex check for Data1", value = TRUE),
      
      h4("Data2 Input Type"),
      radioButtons("data2_type", "Choose Data2 format:", choices = c("PED/MAP", "BED/BIM/FAM"), selected = "PED/MAP"),
      conditionalPanel(
        condition = "input.data2_type == 'PED/MAP'",
        fileInput("data2_ped", "Upload Data2 .ped"),
        fileInput("data2_map", "Upload Data2 .map")
      ),
      conditionalPanel(
        condition = "input.data2_type == 'BED/BIM/FAM'",
        fileInput("data2_bed", "Upload Data2 .bed"),
        fileInput("data2_bim", "Upload Data2 .bim"),
        fileInput("data2_fam", "Upload Data2 .fam")
      ),
      fileInput("sample2", "Data2 Sample Metadata (.xlsx)"),
      textInput("prefix_pattern", "Prefix to Include in Ref_RegNumber from Data2:", value = "^COF"),
      checkboxInput("sex_check_data2", "Run sex check for Data2", value = TRUE),
      actionButton("run_qc", "Run Data Harmonization Pipeline")
    )
    ,
    mainPanel(
      tabsetPanel(
        tabPanel("Status",
                 textOutput("status")), # Show status message
                 
        tabPanel("Instructions",
                          div(style = "font-size: 18px;",  # Adjust as needed
                          h3("Usage Instructions", style = "font-size: 24px;"),
                          tags$ol(
                            tags$li("This app works on Ubuntu (Linux) and R>4.3.1. "),
                            tags$li("It performs merging of two datasets along with sex check and filtering (removal of sex chromosomes, non-rsID SNPs, and ambiguous SNPs)."),
                            tags$li("It also removes duplicates based on genotype call rates."),
                            tags$li("If a dataset has participants from multiple studies, then upload it as Data2 and specify prefix identifier from excel file in the tab"),
                            tags$li("The app requires the datasets to be named 'Data1' and 'Data2'."),
                            tags$li(HTML("The datasets must be accompanied by an Excel sheets name Data1 and Data2 with sample information containing the following columns: <br>
                                          <b>Ref_RegNumber, Index, Sample_ID, Call_Rate, Gender, Gender_code</b>. <br>
                                          <b>Gender_code</b> refers to the actual gender of the participants.")),
                            tags$li("You must specify the project directory path (e.g., /media/home/media/Data)."),
                            tags$li("You must also provide a tools directory path (e.g., /media/icaro/T7/CRC_Genetic_Data/Tools)."),
                            tags$li("Upload the required input datasets in the 'Upload' tab."),
                            tags$li("Run harmonization by clicking the 'Harmonize Data' button."),
                            tags$li("Wait for confirmation in the 'Status' tab.")
                          )
                        )
                      ),
        tabPanel("Sex Mismatch", tableOutput("sex_mismatch")),
        tabPanel("Filtered Data2", dataTableOutput("filtered_data2")),
        tabPanel("Duplicate Samples", tableOutput("duplicate_samples"))
        )
      )
    )
  )

server <- function(input, output, session) {

    # Define reactive values at the top of the server function
  paths <- reactiveValues(
    project_dir = NULL,
    output_dir = NULL,
    data1_base = NULL,
    data2_base = NULL,
    merged_prefix = NULL
  )

  # Update them when input$project_dir changes
  
  # observeEvent(input$project_dir, {
  #   req(input$project_dir)
  #   
  #   # Reinitialize paths every time project_dir changes
  #   paths$project_dir <- normalizePath(input$project_dir, mustWork = FALSE)
  #   paths$output_dir <- file.path(paths$project_dir, "Data_Merge")
  #   dir.create(paths$output_dir, showWarnings = FALSE, recursive = TRUE)
  # 
  #   paths$data1_base <- file.path(paths$project_dir, "Data1")
  #   paths$data2_base <- file.path(paths$project_dir, "Data2")
  #   paths$merged_prefix <- file.path(paths$project_dir, "Merged_Final")
  #   
  #   # Debug message to confirm correct creation
  #   message("Creating Data_Merge in: ", paths$output_dir)
  # })
  
  observeEvent(input$project_dir, {
    req(input$project_dir)
    
    message("Triggered input$project_dir: ", input$project_dir)
    
    # Normalize and validate the input path
    new_project_dir <- normalizePath(input$project_dir, mustWork = FALSE)
    
    # Only proceed if different from current
    if (is.null(paths$project_dir) || paths$project_dir != new_project_dir) {
      paths$project_dir <- new_project_dir
      paths$output_dir <- file.path(paths$project_dir, "Data_Merge")
      dir.create(paths$output_dir, showWarnings = FALSE, recursive = TRUE)
      
      paths$data1_base <- file.path(paths$project_dir, "Data1")
      paths$data2_base <- file.path(paths$project_dir, "Data2")
      paths$merged_prefix <- file.path(paths$project_dir, "Merged_Final")
      
      message("Creating Data_Merge in: ", paths$output_dir)
    }
  })
  
  
  observeEvent(input$run_qc, {
    output$status <- renderText({"Running Data Harmonization Pipeline..."})
    
    isolate({
      req(input$sample1, input$sample2)
      req(input$tools_dir)
      
      # Update PATH with plink executable
      tools_path <- normalizePath(input$tools_dir)
      Sys.setenv(PATH = paste(tools_path, Sys.getenv("PATH"), sep = ":"))
      
      # Clear old data
      unlink(paths$data1_base, recursive = TRUE, force = TRUE)
      unlink(paths$data2_base, recursive = TRUE, force = TRUE)
      dir.create(paths$data1_base, showWarnings = FALSE, recursive = TRUE)
      dir.create(paths$data2_base, showWarnings = FALSE, recursive = TRUE)
      
      unlink(paths$merged_prefix, recursive = TRUE, force = TRUE)
      
      # Optionally validate that plink is accessible
      plink_test <- system2("plink1.9", "--version", stdout = TRUE, stderr = TRUE)
      if (any(grepl("not found", plink_test))) {
        showModal(modalDialog(
          title = "Error",
          "PLINK executable not found in the specified tools directory.",
          easyClose = TRUE
        ))
        return()
      }
      
      dir.create(paths$data1_base, showWarnings = FALSE)
      dir.create(paths$data2_base, showWarnings = FALSE)
      
      #============================#
      #  Save Data1 Files
      #============================#
      if (input$data1_type == "PED/MAP") {
        req(input$data1_ped, input$data1_map)
        file.copy(input$data1_ped$datapath, file.path(paths$data1_base, "Data1.ped"))
        file.copy(input$data1_map$datapath, file.path(paths$data1_base, "Data1.map"))
        system(glue("plink1.9 --file {paths$data1_base}/Data1 --make-bed --out {paths$data1_base}/Data1"))
      } else {
        req(input$data1_bed, input$data1_bim, input$data1_fam)
        file.copy(input$data1_bed$datapath, file.path(paths$data1_base, "Data1.bed"))
        file.copy(input$data1_bim$datapath, file.path(paths$data1_base, "Data1.bim"))
        file.copy(input$data1_fam$datapath, file.path(paths$data1_base, "Data1.fam"))
      }
      
      #============================#
      #  Save Data2 Files
      #============================#
      if (input$data2_type == "PED/MAP") {
        req(input$data2_ped, input$data2_map)
        file.copy(input$data2_ped$datapath, file.path(paths$data2_base, "Data2.ped"))
        file.copy(input$data2_map$datapath, file.path(paths$data2_base, "Data2.map"))
        system(glue("plink1.9 --file {paths$data2_base}/Data2 --make-bed --out {paths$data2_base}/Data2"))
      } else {
        req(input$data2_bed, input$data2_bim, input$data2_fam)
        file.copy(input$data2_bed$datapath, file.path(paths$data2_base, "Data2.bed"))
        file.copy(input$data2_bim$datapath, file.path(paths$data2_base, "Data2.bim"))
        file.copy(input$data2_fam$datapath, file.path(paths$data2_base, "Data2.fam"))
      }
      
      #============================#
      #  Load Sample Metadata
      #============================#
      file.copy(input$sample1$datapath, file.path(paths$output_dir, "data1_sample.xlsx"))
      file.copy(input$sample2$datapath, file.path(paths$output_dir, "data2_sample.xlsx"))
      sample1 <- read_excel(file.path(paths$output_dir, "data1_sample.xlsx"))
      sample2 <- read_excel(file.path(paths$output_dir, "data2_sample.xlsx"))
      
      #============================#
      #  Update Data1 FAM
      #============================#
      fam1 <- fread(file.path(paths$data1_base, "Data1.fam"), sep = " ", header = FALSE)
      dim(fam1)
      colnames(fam1) <- c("FID", "IID", "FatherID", "MotherID", "Sex", "Phenotype")
      fam1 <- fam1 %>%
        left_join(sample1 %>% rename(IID = Sample_ID, NewSex = Gender_code), by = "IID") %>%
        mutate(Sex = NewSex) %>% select(-NewSex)
      fwrite(fam1, file.path(paths$data1_base, "Data1.fam"), sep = " ", col.names = FALSE)

      #============================#
      #  Sex Check Data1
      #============================#

      if (input$sex_check_data1) {
        system(glue("plink1.9 --impute-sex --bfile {paths$data1_base}/Data1 --make-bed --out {paths$data1_base}/Data1_sex_check"))
        sexcheck1 <- fread(file.path(paths$data1_base, "Data1_sex_check.sexcheck"))
        mismatches <- sexcheck1[STATUS == "PROBLEM"]
        output$sex_mismatch <- renderTable({ mismatches })
        fwrite(mismatches[, .(FID, IID)], file.path(paths$data1_base, "sex_mismatch.txt"), col.names = FALSE, sep = "\t")
        system(glue("plink1.9 --bfile {paths$data1_base}/Data1 --remove {paths$data1_base}/sex_mismatch.txt --make-bed --out {paths$data1_base}/Data1_filtered"))
      } else {
        file.copy(file.path(paths$data1_base, "Data1.bed"), file.path(paths$data1_base, "Data1_filtered.bed"))
        file.copy(file.path(paths$data1_base, "Data1.bim"), file.path(paths$data1_base, "Data1_filtered.bim"))
        file.copy(file.path(paths$data1_base, "Data1.fam"), file.path(paths$data1_base, "Data1_filtered.fam"))
        output$sex_mismatch <- renderTable(NULL)
      }
      
      #============================#
      #  Filter & Update Data2
      #============================#
      filtered2 <- sample2 %>% 
        filter(!is.na(Ref_RegNumber)) %>%
        filter(!grepl(input$prefix_pattern, Ref_RegNumber, ignore.case = TRUE) | Call_Rate < 0.9)
      
      # Rows NOT in filtered2 (i.e., rows that passed the filtering)
      filtered2_out <- anti_join(sample2, filtered2)
      
      # Show in DataTable
      output$filtered_data2 <- renderDataTable({
        filtered2_out
      })
      
      fwrite(filtered2[, c("Index", "Sample_ID")],
             file.path(paths$output_dir, "filtered_non_data2.txt"), col.names = FALSE, sep = "\t")
      
      system(glue("plink1.9 --remove {paths$output_dir}/filtered_non_data2.txt --bfile {paths$data2_base}/Data2 --make-bed --out {paths$data2_base}/Data2_filtered"))
      
      fam2 <- fread(file.path(paths$data2_base, "Data2_filtered.fam"), sep = " ", header = FALSE)
      
      colnames(fam2) <- c("FID", "IID", "FatherID", "MotherID", "Sex", "Phenotype")
      fam2 <- fam2 %>%
        left_join(sample2 %>% rename(IID = Sample_ID, NewSex = Gender_code), by = "IID") %>%
        mutate(Sex = NewSex) %>% select(-NewSex)
      fwrite(fam2, file.path(paths$data2_base, "Data2_filtered.fam"), sep = " ", col.names = FALSE)
      
      #============================#
      #  Sex Check Data2
      #============================#
      # system(glue("plink1.9 --impute-sex --bfile {paths$data2_base}/Data2_filtered --make-bed --out {paths$data2_base}/Data2_sex_check"))
      
      if (input$sex_check_data2) {
        system(glue("plink1.9 --impute-sex --bfile {paths$data2_base}/Data2_filtered --make-bed --out {paths$data2_base}/Data2_sex_check"))
        sexcheck2 <- fread(file.path(paths$data2_base, "Data2_sex_check.sexcheck"))
        mismatches <- sexcheck2[STATUS == "PROBLEM"]
        output$sex_mismatch <- renderTable({ mismatches })
        fwrite(mismatches[, .(FID, IID)], file.path(paths$data2_base, "sex_mismatch.txt"), col.names = FALSE, sep = "\t")
        system(glue("plink1.9 --bfile {paths$data1_base}/Data1 --remove {paths$data2_base}/sex_mismatch.txt --make-bed --out {paths$data1_base}/Data1_filtered"))
      } else {
        file.copy(file.path(paths$data2_base, "Data2.bed"), file.path(paths$data2_base, "Data2_filtered.bed"))
        file.copy(file.path(paths$data2_base, "Data2.bim"), file.path(paths$data2_base, "Data2_filtered.bim"))
        file.copy(file.path(paths$data2_base, "Data2.fam"), file.path(paths$data2_base, "Data2_filtered.fam"))
        output$sex_mismatch <- renderTable(NULL)
        }
      
      #============================#
      #  Harmonize SNPs
      #============================#
      bim1 <- fread(file.path(paths$data1_base, "Data1_filtered.bim"))
      bim2 <- fread(file.path(paths$data2_base, "Data2_filtered.bim"))
      
      filter_snps <- function(df) {
        df %>%
          filter(V1 > 0 & V1 < 23) %>%
          filter(str_detect(V2, "^r")) %>%
          filter(!(V5 == "G" & V6 == "C")) %>%
          filter(!(V5 == "C" & V6 == "G")) %>%
          filter(!(V5 == "T" & V6 == "A")) %>%
          filter(!(V5 == "A" & V6 == "T"))
      }
      
      bim1 <- filter_snps(bim1)
      bim2 <- filter_snps(bim2)
      
      common_snps <- intersect(bim1$V2, bim2$V2)
      bim1 <- bim1[V2 %in% common_snps]
      bim2 <- bim2[V2 %in% common_snps]
      
      exact_match <- bim1$V1 == bim2$V1 &
        bim1$V2 == bim2$V2 &
        bim1$V4 == bim2$V4 &
        ((bim1$V5 == bim2$V5 & bim1$V6 == bim2$V6) |
           (bim1$V5 == bim2$V6 & bim1$V6 == bim2$V5))
      
      final_snps <- bim1[exact_match]
      unmatched_snps <- bim1[!exact_match, ]
      
      write.table(unmatched_snps, file.path(paths$output_dir, "unmatched_snps.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      fwrite(as.data.frame(final_snps$V2), file.path(paths$output_dir, "snps_qc.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
      write.table(final_snps[, .(V2, V5)], file.path(paths$output_dir, "new_ref_alleles.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      system(glue("plink1.9 --bfile {paths$data1_base}/Data1_filtered --extract {paths$output_dir}/snps_qc.txt --make-bed --out {paths$data1_base}/Data1_final"))
      system(glue("plink1.9 --keep-allele-order --bfile {paths$data2_base}/Data2_filtered --extract {paths$output_dir}/snps_qc.txt --make-bed --out {paths$data2_base}/Data2_final --reference-allele {paths$output_dir}/new_ref_alleles.txt"))
      
      #============================#
      #  Remove Duplicates
      #============================#
      dup1 <- sample1 %>% select(IID = Sample_ID, Sample_Name, CallRate1 = Call_Rate)
      dup2 <- sample2 %>% select(IID = Sample_ID, Sample_Name, CallRate2 = Call_Rate)
      
      dup_id <- inner_join(dup1, dup2, by = "IID")
      dup_name <- inner_join(dup1, dup2, by = "Sample_Name") %>%
        filter(!(IID.x %in% dup_id$IID & IID.y %in% dup_id$IID))
      
      dups <- bind_rows(
        dup_id %>% transmute(IID1 = IID, IID2 = IID, CR1 = CallRate1, CR2 = CallRate2),
        dup_name %>% transmute(IID1 = IID.x, IID2 = IID.y, CR1 = CallRate1, CR2 = CallRate2)
      )
      output$duplicate_samples <- renderTable({ dups })
      
      to_remove_1 <- dups %>% filter(CR1 < CR2) %>% pull(IID1) %>% unique()
      to_remove_2 <- dups %>% filter(CR1 >= CR2) %>% pull(IID2) %>% unique()
      
      fwrite(data.frame(FID = to_remove_1, IID = to_remove_1), file.path(paths$output_dir, "dups_data1.txt"), col.names = FALSE, sep = "\t")
      fwrite(data.frame(FID = to_remove_2, IID = to_remove_2), file.path(paths$output_dir, "dups_data2.txt"), col.names = FALSE, sep = "\t")
      
      
      # Load original Data1 sheet
      data1_info <- sample1 %>% select(Index, Sample_ID)
      
      # Filter rows where Sample_ID is in the list of removed duplicates
      remove_dups_data1_info <- data1_info %>% filter(Sample_ID %in%  to_remove_1)
      
      # Write to text file: index and sample ID of removed duplicates
      write.table(remove_dups_data1_info,
                  file.path(paths$output_dir, "remove_dups_data1.txt"),
                  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      # Load original Data2 sheet
      data2_info <- sample2 %>% select(Index, Sample_ID)
      
      # Filter rows where Sample_ID is in the list of removed duplicates
      remove_dups_data2_info <- data2_info %>% filter(Sample_ID %in%  to_remove_2)
      
      # Write to text file: index and sample ID of removed duplicates
      write.table(remove_dups_data2_info,
                  file.path(paths$output_dir, "remove_dups_data2.txt"),
                  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      system(glue("plink1.9 --bfile {paths$data1_base}/Data1_final --remove {paths$output_dir}/remove_dups_data1.txt --make-bed --out {paths$data1_base}/Data1_nodups"))
      system(glue("plink1.9 --bfile {paths$data2_base}/Data2_final --remove {paths$output_dir}/remove_dups_data2.txt --make-bed --out {paths$data2_base}/Data2_nodups"))
      
      #============================#
      #  Final Merge
      #============================#
      system(glue("plink1.9 --bfile {paths$data1_base}/Data1_nodups --bmerge {paths$data2_base}/Data2_nodups.bed {paths$data2_base}/Data2_nodups.bim {paths$data2_base}/Data2_nodups.fam --make-bed --out {paths$merged_prefix} --reference-allele {paths$output_dir}/new_ref_alleles.txt"))
      
      output$status <- renderText({
        "✅ Check Console Output.Data Harmonization Complete. Proceed to QC1-4 pipeline."
      })
    })
  })
}

shinyApp(ui = ui, server = server)
