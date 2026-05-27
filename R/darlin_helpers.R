library(Seurat)
library(tidyverse)
library(glue)

HOME <- here::here()
source(glue("{HOME}/R/bmk_io.R"))
source(glue("{HOME}/R/helpers.R"))

load_darlin_data <- function(matrix_path, png_path) {
  seu <- CreateS1000Object(
    matrix_path = glue(matrix_path),
    png_path = glue(png_path),
    min.cells = 1, min.features = 1, spot_radius = 0.001
  )
  seu <- ZoomImageCoords(seu)
  return(seu)
}

generate_clone_table <- function(seu, features_file, sample.ID, locus.use) {
  allele_df <- read_tsv(features_file, show_col_types = F)
  allele_df <- subset(allele_df, !is.na(allele)) %>%
    select(-type) %>%
    setNames(c("clone_bc", "clone_allele"))
  allele_df <- allele_df[!duplicated(allele_df$clone_bc), ]
  allele_df <- allele_df %>%
    mutate(clone_id = paste0("clone_", row_number()))

  X <- LayerData(seu, layer = "counts")
  X_df <- summary(X) %>%
    data.frame() %>%
    setNames(c("i", "j", "nUMI")) %>%
    as_tibble()

  X_df <- X_df %>%
    mutate(
      cell_bc = colnames(X)[j],
      clone_bc = rownames(X)[i]
    ) %>%
    select(-i, -j) %>%
    left_join(allele_df, by = "clone_bc")

  X_df <- X_df %>%
    mutate(
      library = sample.ID,
      cell_id = paste0(library, "_", cell_bc),
      templete = locus.use
    ) %>%
    select(nUMI, library, cell_id, cell_bc, clone_id, clone_bc, clone_allele, templete)

  return(X_df)
}


report_QC <- function(cells.mrna, clone.table.file, sample.ID, locus.use, log.file = NULL) {
  # Read data
  X_df <- read_tsv(clone.table.file, show_col_types = F)

  # Calculate QC metrics
  cell_QC <- c(
    X_df %>% filter(is_pass_mRNA_QC) %>% select(cell_id) %>% distinct() %>% nrow(),
    X_df %>% filter(is_pass_mRNA_QC & is_dominant) %>% select(cell_id) %>% distinct() %>% nrow(),
    X_df %>% filter(is_pass_mRNA_QC & is_dominant & is_edited) %>% select(cell_id) %>% distinct() %>% nrow(),
    X_df %>% filter(is_pass_mRNA_QC & is_dominant & is_edited & is_rare) %>% select(cell_id) %>% distinct() %>% nrow()
  )
  cell_QC_pct <- round(cell_QC / length(cells.mrna) * 100, 2)

  clone_QC <- c(
    X_df %>% filter(is_pass_mRNA_QC) %>% select(clone_allele) %>% distinct() %>% nrow(),
    X_df %>% filter(is_pass_mRNA_QC & is_dominant) %>% select(clone_allele) %>% distinct() %>% nrow(),
    X_df %>% filter(is_pass_mRNA_QC & is_dominant & is_edited) %>% select(clone_allele) %>% distinct() %>% nrow(),
    X_df %>% filter(is_pass_mRNA_QC & is_dominant & is_edited & is_rare) %>% select(clone_allele) %>% distinct() %>% nrow()
  )

  # Format report messages
  report <- c(
    "=========================================================================",
    glue("Sample: {sample.ID}    Locus: {locus.use}"),
    "Cell-level QC metrics:",
    glue("Cells passing mRNA QC:                   {length(cells.mrna)}"),
    glue("Cells with detected allele:              {cell_QC[1]} ({cell_QC_pct[1]}%)"),
    glue("Cells with dominant allele:              {cell_QC[2]} ({cell_QC_pct[2]}%)"),
    glue("Cells with dominant edited allele:       {cell_QC[3]} ({cell_QC_pct[3]}%)"),
    glue("Cells with dominant edited rare allele:  {cell_QC[4]} ({cell_QC_pct[4]}%)"),
    "-------------------------------------------------------------------------",
    "Clone-level QC metrics:",
    glue("Total clones:                             {clone_QC[1]}"),
    glue("Clones with dominant allele:              {clone_QC[2]}"),
    glue("Clones with dominant edited allele:       {clone_QC[3]}"),
    glue("Clones with dominant edited rare allele:  {clone_QC[4]}"),
    "========================================================================="
  )

  # Print report

  if (is.null(log.file)) {
    for (msg in report) {
      message(msg)
    }
  } else {
    writeLines(report, log.file)
  }
}


AlleleStepwiseQCByCell <- function(X_df, cells.mrna, freq.cutoff = 1e-4, max.sample.N = 3) {
  # cells after stepwise QC
  cells.bc <- X_df %>%
    pull(cell_bc) %>%
    unique() %>%
    intersect(., cells.mrna)
  cells.dominant <- X_df %>%
    filter(is_dominant) %>%
    pull(cell_bc) %>%
    unique() %>%
    intersect(., cells.mrna)
  cells.edited <- X_df %>%
    filter(is_dominant & is_edited) %>%
    pull(cell_bc) %>%
    unique() %>%
    intersect(., cells.mrna)
  cells.rare <- X_df %>%
    mutate(is_rare = normalized_count < freq.cutoff & sample_count <= max.sample.N) %>%
    filter(is_dominant & is_edited & is_rare) %>%
    pull(cell_bc) %>%
    unique() %>%
    intersect(., cells.mrna)
  # store results
  data.tmp <- data.frame(
    cells = sapply(list(cells.mrna, cells.bc, cells.dominant, cells.edited, cells.rare), length),
    step = c("mRNA", "barcoding", "dominant", "edited", "rare")
  )
  data.tmp <- data.tmp %>%
    mutate(frac = cells / length(cells.mrna))
  return(data.tmp)
}


MergedStepwiseQCByCell <- function(X_df.CA, X_df.TA, X_df.RA, cells.mrna, freq.cutoff = 1e-4, max.sample.N = 3) {
  cells.bc <- c(
    X_df.CA %>% pull(cell_bc) %>% unique(),
    X_df.TA %>% pull(cell_bc) %>% unique(),
    X_df.RA %>% pull(cell_bc) %>% unique()
  ) %>%
    unique() %>%
    intersect(., cells.mrna)

  cells.dominant <- c(
    X_df.CA %>% filter(is_dominant) %>% pull(cell_bc) %>% unique(),
    X_df.TA %>% filter(is_dominant) %>% pull(cell_bc) %>% unique(),
    X_df.RA %>% filter(is_dominant) %>% pull(cell_bc) %>% unique()
  ) %>%
    unique() %>%
    intersect(., cells.mrna)

  cells.edited <- c(
    X_df.CA %>% filter(is_dominant & is_edited) %>% pull(cell_bc) %>% unique(),
    X_df.TA %>% filter(is_dominant & is_edited) %>% pull(cell_bc) %>% unique(),
    X_df.RA %>% filter(is_dominant & is_edited) %>% pull(cell_bc) %>% unique()
  ) %>%
    unique() %>%
    intersect(., cells.mrna)

  cells.rare <- c(
    X_df.CA %>%
      mutate(is_rare = normalized_count < freq.cutoff & sample_count <= max.sample.N) %>%
      filter(is_dominant & is_edited & is_rare) %>% pull(cell_bc) %>% unique(),
    X_df.TA %>%
      mutate(is_rare = normalized_count < freq.cutoff & sample_count <= max.sample.N) %>%
      filter(is_dominant & is_edited & is_rare) %>% pull(cell_bc) %>% unique(),
    X_df.RA %>%
      mutate(is_rare = normalized_count < freq.cutoff & sample_count <= max.sample.N) %>%
      filter(is_dominant & is_edited & is_rare) %>% pull(cell_bc) %>% unique()
  ) %>%
    unique() %>%
    intersect(., cells.mrna)

  # store results
  data.tmp <- data.frame(
    cells = sapply(list(cells.mrna, cells.bc, cells.dominant, cells.edited, cells.rare), length),
    step = c("mRNA", "barcoding", "dominant", "edited", "rare")
  )
  data.tmp <- data.tmp %>%
    mutate(frac = cells / length(cells.mrna))
  return(data.tmp)
}

AlleleStepwiseQCByAllele <- function(X_df, cells.mrna, freq.cutoff = 1e-4, max.sample.N = 3) {
  # cells after stepwise QC
  clones.mrna <- X_df %>%
    filter(cell_bc %in% cells.mrna) %>%
    pull(clone_bc) %>%
    unique()
  clones.dominant <- X_df %>%
    filter(cell_bc %in% cells.mrna) %>%
    filter(is_dominant) %>%
    pull(clone_bc) %>%
    unique()
  cells.edited <- X_df %>%
    filter(cell_bc %in% cells.mrna) %>%
    filter(is_dominant & is_edited) %>%
    pull(clone_bc) %>%
    unique()
  cells.rare <- X_df %>%
    filter(cell_bc %in% cells.mrna) %>%
    mutate(is_rare = normalized_count < freq.cutoff & sample_count <= max.sample.N) %>%
    filter(is_dominant & is_edited & is_rare) %>%
    pull(clone_bc) %>%
    unique()
  # store results
  data.tmp <- data.frame(
    clones = sapply(list(clones.mrna, clones.dominant, cells.edited, cells.rare), length),
    step = c("mRNA", "dominant", "edited", "rare")
  )
  data.tmp <- data.tmp %>%
    mutate(frac = clones / length(clones.mrna))
  return(data.tmp)
}
