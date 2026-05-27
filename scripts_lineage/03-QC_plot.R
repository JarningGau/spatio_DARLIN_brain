#' Title: QC plots of DARLIN data
#' Author: Jarning Gau
#' Date: 2025-03-31
#'
setwd("scripts_lineage")

library(tidyverse)
library(glue)
library(qs)
library(Seurat)
library(cowplot)

source("../R/helpers.R")
source("../R/spatial_plot_utils.R")

out_dir <- "../data/lineage"
mrna_qc_out_dir <- "../data/mRNA"

do_fragment_plot <- TRUE
do_edit_plot <- TRUE
do_spatial_plot <- TRUE

scale_bar_um <- 500
um_per_pixel <- 6.8
scale_bar_pixels <- scale_bar_um / um_per_pixel

#### Plot helpers ####

make_plot_data <- function(darlin_data) {
  data_stat <- lapply(darlin_data$clone_allele, function(allele) {
    allele <- unlist(strsplit(allele, split = ","))
    data.frame(
      allele = allele,
      start = sub("([0-9]+)_.*", "\\1", allele) %>% as.integer(),
      end = sub(".*_([0-9]+).+", "\\1", allele) %>% as.integer(),
      type = ifelse(grepl("del", allele), "deletion", "insertion")
    )
  }) %>%
    do.call(rbind, .)

  data_stat
}

editing_site_plot <- function(data_stat, locus_name, display_name) {
  cut_sites <- seq(19, 270, 27)

  data_stat %>%
    pivot_longer(cols = 2:3, names_to = "position", values_to = "value") %>%
    filter(!(type == "insertion" & position == "end")) %>%
    filter(!is.na(value)) %>%
    ggplot(aes(value, fill = type)) +
    annotate(
      "rect",
      xmin = cut_sites,
      xmax = cut_sites + 7,
      ymin = 0,
      ymax = Inf,
      alpha = 0.5,
      fill = "grey"
    ) +
    geom_histogram(binwidth = 2) +
    labs(
      x = "Editing site (bp)",
      y = glue("{display_name}\n\nEvents"),
      title = locus_name
    ) +
    theme(
      legend.position = c(0.9, 1),
      legend.justification = c(1, 1),
      legend.background = element_blank(),
      legend.title = element_blank()
    )
}

fragment_len_plot <- function(darlin_data, locus_name, display_name) {
  orig_len <- c(
    "CA" = 276,
    "RA" = 275,
    "TA" = 275
  )
  darlin_data$seq_len <- sapply(darlin_data$clone_bc, str_length)

  ggplot(darlin_data, aes(seq_len)) +
    geom_histogram(binwidth = 1, fill = "orange") +
    geom_vline(
      xintercept = orig_len[locus_name],
      linewidth = 1.5,
      color = "black",
      alpha = 0.2
    ) +
    labs(
      x = "Fragment length (bp)",
      y = glue("{display_name}\n\nFragments"),
      title = locus_name
    )
}

spatial_qc_plot <- function(seu, locus_name, pt_alpha = 0.5) {
  SpatialDimPlot(
    seu,
    group.by = "orig.ident",
    pt.size.factor = 5,
    image.alpha = 0,
    alpha = pt_alpha
  ) +
    AddScaleBar(seu, scale_bar_pixels, position = "bottom-right") +
    scale_fill_manual(values = "blue") +
    ggtitle(locus_name) +
    NoLegend() +
    theme(
      plot.title = element_text(size = 25, hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
    )
}

#### Main ####

samples <- c(
  "L0927_Brain",
  "L126_Brain_s1",
  "L126_Brain_s2",
  "L126_Brain_s3"
)
arrays <- c("CA", "TA", "RA")
locus_list <- c(arrays, "mRNA")

print_with_time("Starting QC plot script")

## Fragment length distribution
if (do_fragment_plot) {
  print_with_time("Fragment length distribution plots - start")
  for (sample_name in samples) {
    print_with_time(glue("Fragment length plots: processing {sample_name}"))
    p_list <- lapply(arrays, function(locus_name) {
      clone_file <- file.path(
        out_dir,
        sample_name,
        glue("clone_table_{locus_name}_by_cell.tsv")
      )
      darlin_data <- read_tsv(clone_file, show_col_types = FALSE)
      p <- fragment_len_plot(darlin_data, locus_name, sample_name)
      if (locus_name != "CA") {
        p <- p + ylab("")
      }
      p
    })
    p <- plot_grid(plotlist = p_list, nrow = 1)

    fn <- file.path(
      out_dir,
      sample_name,
      "QC_frag_len_dist.png"
    )
    ggsave(fn, p, width = 10, height = 2.5, units = "in", dpi = 200)
    print_with_time(glue("Fragment length plots: wrote {fn}"))
  }
  print_with_time("Fragment length distribution plots - done")
}

## Editing sites
if (do_edit_plot) {
  print_with_time("Editing site plots - start")
  for (sample_name in samples) {
    print_with_time(glue("Editing site plots: processing {sample_name}"))
    p_list <- lapply(arrays, function(locus_name) {
      clone_file <- file.path(
        out_dir,
        sample_name,
        glue("clone_table_{locus_name}_by_cell.tsv")
      )
      darlin_data <- read_tsv(clone_file, show_col_types = FALSE)
      darlin_data <- darlin_data %>% filter(clone_allele != "[]")
      data_stat <- make_plot_data(darlin_data)
      p <- editing_site_plot(data_stat, locus_name, sample_name)
      if (locus_name != "CA") {
        p <- p + ylab("")
      }
      p
    })
    p <- plot_grid(plotlist = p_list, nrow = 1)

    fn <- file.path(
      out_dir,
      sample_name,
      "QC_editing_sites.png"
    )
    ggsave(fn, p, width = 10, height = 2.5, units = "in", dpi = 200)
    print_with_time(glue("Editing site plots: wrote {fn}"))
  }
  print_with_time("Editing site plots - done")
}

## Spatial distribution
if (do_spatial_plot) {
  print_with_time("Spatial distribution plots - start")
  for (sample_name in samples) {
    print_with_time(glue("Spatial distribution plots: processing {sample_name}"))
    p_list <- lapply(locus_list, function(locus_name) {
      if (locus_name == "mRNA") {
        fn <- file.path(
          mrna_qc_out_dir,
          sample_name,
          glue("{sample_name}_Seurat_raw.qs")
        )
        seu <- qs::qread(fn)
        pt_alpha <- 0.1
      } else {
        fn <- file.path(
          out_dir,
          sample_name,
          glue("Seurat_{locus_name}.qs")
        )
        seu <- qs::qread(fn)
        pt_alpha <- 0.5
      }
      spatial_qc_plot(seu, locus_name, pt_alpha)
    })
    p <- plot_grid(plotlist = p_list, ncol = 4)

    fn <- file.path(
      out_dir,
      sample_name,
      "QC_spatial_dist.png"
    )
    ggsave(fn, p, width = 14, height = 4, units = "in", dpi = 200)
    print_with_time(glue("Spatial distribution plots: wrote {fn}"))
  }
  print_with_time("Spatial distribution plots - done")
}

print_with_time("QC plot script completed")
