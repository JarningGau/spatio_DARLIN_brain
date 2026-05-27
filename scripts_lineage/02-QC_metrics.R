#' Title: QC metrics of DARLIN data
#' Author: Jarning Gau
#' Date: 2025-03-31
#'
setwd("scripts_lineage")

library(tidyverse)
library(glue)

source("../R/darlin_helpers.R")

out_dir <- "../data/lineage/"
mrna_qc_out_dir <- "../data/mRNA/"
freq_cutoff <- 1e-5
max_sample_n <- 1

#### QC helpers ####

qc_summary_1 <- function(sample_name, locus_name) {
  clone_file <- file.path(
    out_dir,
    sample_name,
    glue("clone_table_{locus_name}_by_cell.tsv")
  )
  X_df <- read_tsv(clone_file, show_col_types = FALSE)

  num_alleles <- X_df %>%
    pull(clone_allele) %>%
    unique() %>%
    length()

  num_cells <- X_df %>%
    pull(cell_bc) %>%
    unique() %>%
    length()

  data.frame(
    sample = sample_name,
    locus = locus_name,
    num_alleles = num_alleles,
    num_cells = num_cells
  )
}

qc_summary_2 <- function(sample_name, locus_name, freq_cutoff = 1e-5, max_sample_n = 3) {
  cells_mrna_file <- file.path(
    mrna_qc_out_dir,
    sample_name,
    "barcodes_HQ.txt"
  )
  cells_mrna <- readLines(cells_mrna_file)

  if (locus_name != "Merged") {
    clone_file <- file.path(
      out_dir,
      sample_name,
      glue("clone_table_{locus_name}_by_cell.tsv")
    )
    X_df <- read_tsv(clone_file, show_col_types = FALSE)
    data_stat <- AlleleStepwiseQCByCell(
      X_df,
      cells_mrna,
      freq_cutoff,
      max_sample_n
    )
  } else {
    X_df_CA <- read_tsv(
      file.path(
        out_dir,
        sample_name,
        "clone_table_CA_by_cell.tsv"
      ),
      show_col_types = FALSE
    )
    X_df_TA <- read_tsv(
      file.path(
        out_dir,
        sample_name,
        "clone_table_TA_by_cell.tsv"
      ),
      show_col_types = FALSE
    )
    X_df_RA <- read_tsv(
      file.path(
        out_dir,
        sample_name,
        "clone_table_RA_by_cell.tsv"
      ),
      show_col_types = FALSE
    )
    data_stat <- MergedStepwiseQCByCell(
      X_df_CA,
      X_df_TA,
      X_df_RA,
      cells_mrna,
      freq_cutoff,
      max_sample_n
    )
  }

  data_stat$sample <- sample_name
  data_stat$locus <- locus_name

  data_stat
}

qc_summary_3 <- function(sample_name, locus_name, freq_cutoff = 1e-5, max_sample_n = 3) {
  cells_mrna_file <- file.path(
    mrna_qc_out_dir,
    sample_name,
    "barcodes_HQ.txt"
  )
  cells_mrna <- readLines(cells_mrna_file)

  clone_file <- file.path(
    out_dir,
    sample_name,
    glue("clone_table_{locus_name}_by_cell.tsv")
  )
  X_df <- read_tsv(clone_file, show_col_types = FALSE)

  data_stat <- AlleleStepwiseQCByAllele(
    X_df,
    cells_mrna,
    freq_cutoff,
    max_sample_n
  )
  data_stat$sample <- sample_name
  data_stat$locus <- locus_name

  data_stat
}

#### Main ####

samples <- c(
  "L0927_Brain",
  "L126_Brain_s1",
  "L126_Brain_s2",
  "L126_Brain_s3"
)
arrays <- c("CA", "TA", "RA")

## QC summary 1: detected cells and alleles
for (sample_name in samples) {
  data_stat <- lapply(
    arrays,
    function(locus_name) {
      qc_summary_1(sample_name, locus_name)
    }
  ) %>%
    do.call(rbind, .)

  fn <- file.path(
    out_dir,
    sample_name,
    "QC_summary_1.tsv"
  )
  write_tsv(data_stat, fn)
}

## QC summary 2: fraction of cells left after stepwise QC
for (sample_name in samples) {
  data_stat <- lapply(
    c(arrays, "Merged"),
    function(locus_name) {
      qc_summary_2(sample_name, locus_name)
    }
  ) %>%
    do.call(rbind, .)

  fn <- file.path(
    out_dir,
    sample_name,
    "QC_summary_2.tsv"
  )
  write_tsv(data_stat, fn)
}

## QC summary 3: fraction of alleles left after stepwise QC
for (sample_name in samples) {
  data_stat <- lapply(
    arrays,
    function(locus_name) {
      qc_summary_3(sample_name, locus_name)
    }
  ) %>%
    do.call(rbind, .)

  fn <- file.path(
    out_dir,
    sample_name,
    "QC_summary_3.tsv"
  )
  write_tsv(data_stat, fn)
}
