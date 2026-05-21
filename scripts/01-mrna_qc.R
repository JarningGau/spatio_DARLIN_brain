setwd("scripts")

library(tidyverse)
library(glue)
# reticulate::use_python("../.pixi/envs/default/bin/python")
# reticulate::py_config()

source("../R/bmk_io.R")
source("../R/helpers.R")
source("../R/QC_plot_utils.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: pixi run Rscript scripts/01-mrna_qc.R <sample_name> <do_rotate>\n",
    "  Example: pixi run Rscript scripts/01-mrna_qc.R L0927_Brain TRUE",
    call. = FALSE
  )
}
sample_name <- args[1L]
do_rotate <- as.logical(toupper(args[2L]))
if (is.na(do_rotate)) {
  stop("do_rotate must be TRUE or FALSE, got: ", args[2L], call. = FALSE)
}

out_dir <- "../data/mRNA/"

data_dir <- "/mnt/e/Spatio_DARLIN_data/S3000/mRNA/processed"
sample_dir <- file.path(data_dir, sample_name)
mtx_dir <- file.path(sample_dir, "cell_split/mtx")
image_file <- file.path(sample_dir, "he_roi_small.png")


seu <- CreateS1000Object(
  matrix_path = mtx_dir,
  png_path = image_file,
  min.cells = 1,
  min.features = 1,
  spot_radius = 0.001
)
seu$sample_name <- sample_name
seu <- ZoomImageCoords(seu)
if (do_rotate) {
  seu <- RotateSeurat(seu)
}
fn <- file.path(out_dir, sample_name, glue("{sample_name}_Seurat_raw.qs"))
save_qs(seu, fn)

## QC plots
QCVlnPlot(
  seu,
  title = sample_name,
  save.file = file.path(out_dir, sample_name, "QC", "QC_Vln_UMI.png"),
  metric = "UMI"
)
QCVlnPlot(
  seu,
  title = sample_name,
  save.file = file.path(out_dir, sample_name, "QC", "QC_Vln_Gene.png"),
  metric = "Gene"
)
QCSpatialPlot(seu,
  title = sample_name,
  save.file = file.path(out_dir, sample_name, "QC", "QC_spatial_UMI.png"),
  metric = "UMI"
)
QCSpatialPlot(seu,
  title = sample_name,
  save.file = file.path(out_dir, sample_name, "QC", "QC_spatial_Gene.png"),
  metric = "Gene"
)

seu <- subset(seu, nFeature_Spatial >= 100)
fn <- file.path(out_dir, sample_name, glue("{sample_name}_Seurat_qc.qs"))
save_qs(seu, fn)
fn <- file.path(out_dir, sample_name, "barcodes_HQ.txt")
writeLines(Cells(seu), fn)
