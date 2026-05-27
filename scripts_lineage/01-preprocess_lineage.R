#' Title: Preprocessing of DARLIN data
#' Author: Jarning Gau
#' Date: 2025-03-31
#' Input:
#'   - raw *.mtx data from spatio_DARLIN pipeline
#' (https://github.com/JarningGau/spatio_DARLIN)
#'   - cell barcodes passing mRNA QC
#' Output:
#'   - <clone table by cell>
#'   - <clone table by clone>
#'     Only cells passing mRNA QC and having dominant BC were included.
#'     We do not perform editing and homoplasy QC here.
#'   - <spatial distribution of allele>
#' Filtering strategy
#'   1. mRNA QC
#'   2. dominant QC (cells with dominant lineage BCs) top 1 > rest
#'   3. editing QC
#'   4. homoplasy QC (normalized frequency in allele bank < 1e-5)
setwd("scripts_lineage")

library(tidyverse)
library(glue)

source("../R/darlin_helpers.R")
source("../R/helpers.R")

allele_freq_cutoff <- 1e-5
out_dir <- "../data/lineage/"
data_dir <- "/mnt/e/Spatio_DARLIN_data/S3000"
allele_bank_dir <- "../data/allele_bank/"
mrna_qc_out_dir <- "../data/mRNA/"

main <- function(sample_name, locus_name) {
  print_with_time(glue("Processing {sample_name}_{locus_name} ..."))
  ## Step 1: Load raw data
  mtx_dir <- file.path(
    data_dir,
    "lineage/processed/outs",
    glue("{sample_name}_{locus_name}"),
    "cellbin"
  )
  png_fn <- file.path(
    data_dir,
    "mRNA/processed",
    sample_name,
    "he_roi_small.png"
  )
  seu <- load_darlin_data(matrix_path = mtx_dir, png_path = png_fn)
  fn <- file.path(
    out_dir,
    sample_name,
    glue("Seurat_{locus_name}.qs")
  )
  save_qs(seu, fn)

  ## Step 2: Generate clone tables
  fn <- file.path(
    data_dir,
    "lineage/processed/outs",
    glue("{sample_name}_{locus_name}"),
    "cellbin",
    "features_allele.tsv"
  )
  X_df <- generate_clone_table(seu, fn, sample_name, locus_name)
  X_df <- X_df %>% filter(!is.na(clone_allele))

  ## Step 3: mRNA QC
  fn <- file.path(
    mrna_qc_out_dir,
    sample_name,
    "barcodes_HQ.txt"
  )
  cells_mrna <- readLines(fn)
  X_df <- X_df %>%
    mutate(is_pass_mRNA_QC = cell_bc %in% cells_mrna, .before = 1)

  ## Step 4: Allele Bank QC
  fn <- file.path(
    allele_bank_dir,
    glue("allele_bank_Gr_{locus_name}.csv.gz")
  )
  allele_bank <- read_csv(fn, show_col_types = F)
  allele_bank <- allele_bank[c("allele", "sample_count", "normalized_count")]
  colnames(allele_bank)[1] <- "clone_allele"

  X_df <- X_df %>%
    left_join(allele_bank, by = "clone_allele") %>%
    mutate(
      sample_count = replace_na(sample_count, 0),
      normalized_count = replace_na(normalized_count, 0)
    )

  ## Step5: Dominant QC
  second_max <- function(x) ifelse(length(x) > 1, sort(x, decreasing = TRUE)[2], 0)

  X_df <- X_df %>%
    group_by(cell_id) %>%
    mutate(
      nUMI_per_cell = sum(nUMI),
      prop_per_cell = nUMI / nUMI_per_cell, .before = 2
    ) %>%
    ungroup() %>%
    arrange(cell_id, desc(nUMI)) %>%
    group_by(cell_id) %>%
    mutate(is_dominant = prop_per_cell > second_max(prop_per_cell), .before = 1) %>%
    ungroup() %>%
    mutate(
      is_rare = normalized_count < allele_freq_cutoff,
      is_edited = !is.na(clone_allele) & clone_allele != "[]",
      .before = 1
    )

  ## save tmp data
  fn <- file.path(
    out_dir,
    sample_name,
    glue("clone_table_{locus_name}_by_cell.tsv")
  )
  write_tsv(X_df, fn)

  ## Step 6: QC summary
  log_fn <- file.path(
    out_dir,
    sample_name,
    glue("clone_table_{locus_name}_stat.tsv")
  )
  report_QC(
    cells.mrna = cells_mrna,
    clone.table.file = fn,
    sample.ID = sample_name,
    locus.use = locus_name,
    log.file = log_fn
  )

  ## Step 7: Filtering & Reshape
  # for each clone_bc, wrap the cell_bc into the same column (named CB), seperated by ','
  # only for dominant signal
  clone_df <- X_df %>%
    filter(!is.na(clone_allele) & is_pass_mRNA_QC & is_dominant) %>%
    group_by(clone_allele) %>%
    reframe(
      allele = clone_allele,
      UMI_count = sum(nUMI),
      CB = paste(cell_bc, collapse = ","),
      CB_N = n(),
      normalized_count = ifelse(!is.na(normalized_count), normalized_count, 0),
      sample_count = ifelse(!is.na(sample_count), sample_count, 0)
    ) %>%
    select(-clone_allele) %>%
    mutate(
      sample = sample_name,
      mouse = unlist(strsplit(sample_name, split = "_"))[1],
      .before = CB
    ) %>%
    arrange(desc(allele)) %>%
    distinct()

  ## save data
  fn <- file.path(
    out_dir,
    sample_name,
    glue("clone_table_{locus_name}_by_clone.tsv")
  )
  write_tsv(clone_df, fn)

  print_with_time(glue("Done!"))
}

#### Main ####
samples <- c(
  "L0927_Brain",
  "L126_Brain_s1",
  "L126_Brain_s2",
  "L126_Brain_s3"
)
arrays <- c("CA", "TA", "RA")

for (sample in samples) {
  for (array in arrays) {
    main(sample, array)
  }
}
