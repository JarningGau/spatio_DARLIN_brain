setwd("scripts")

library(tidyverse)
library(Seurat)
library(glue)

source("../R/helpers.R")
source("../R/themes.R")
source("../R/cluster.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop(
    "Usage: pixi run Rscript scripts/02-cluster.R <sample_name>\n",
    "  Example: pixi run Rscript scripts/02-cluster.R L0927_Brain",
    call. = FALSE
  )
}
sample_name <- args[1L]

## paths
out_dir <- "../data/mRNA/"
## parameters
npcs <- 50
n_downsample <- 200

fn <- file.path(out_dir, sample_name, glue("{sample_name}_Seurat_qc.qs"))
seu <- qs::qread(fn)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu[["Spatial"]]@layers$scale.data <- NULL
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:npcs)
seu <- FindClusters(seu, resolution = 3)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:npcs)

p <- DimPlot(seu,
  reduction = "umap",
  group.by = "Spatial_snn_res.3",
  label = TRUE
) + NoLegend() + NoAxes()
fn <- file.path(out_dir, sample_name, "Louvain_DimPlot.png")
ggsave(fn, p, width = 10, height = 10, units = "in", dpi = 300)

Idents(seu) <- "Spatial_snn_res.3"
seu_ds <- subset(seu, downsample = n_downsample)
seu_ds <- CreateSeuratObject(
  counts = LayerData(seu_ds, layer = "counts"),
  meta.data = seu_ds@meta.data
)
seu_ds <- NormalizeData(seu_ds)
Idents(seu_ds) <- "Spatial_snn_res.3"
all_markers <- FindAllMarkers(seu_ds, only.pos = TRUE, logfc.threshold = 0.2)
all_markers <- subset(all_markers, p_val_adj < 1e-6)
fn <- file.path(out_dir, sample_name, "Louvain_all_markers.tsv")
write_tsv(all_markers, fn)

top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()
fn <- file.path(out_dir, sample_name, "Louvain_top_markers.tsv")
write_tsv(top_markers, fn)

seu <- iter_cluster(
  seu,
  graph.name = "Spatial_snn",
  cluster.name = "Spatial_snn_res.3",
  max.round = 1,
  num.degs.cutoff = 1
)

seu[["Spatial"]]@layers$data <- NULL
fn <- file.path(out_dir, sample_name, glue("{sample_name}_Seurat_cluster.qs"))
save_qs(seu, fn)

Idents(seu) <- "Spatial_snn_res.3"
cells <- CellsByIdentities(seu)
clusters <- levels(Idents(seu))
for (cluster in clusters) {
  p1 <- DimPlot(seu,
    cells.highlight = cells[[cluster]],
    cols.highlight = c("red", "grey"),
    pt.size = 0.1,
    sizes.highlight = 0.1
  ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle(glue("Cluster {cluster}")) +
    NoLegend() +
    NoAxes() +
    theme(plot.title = element_text(hjust = .5, face = "bold"))

  p2 <- SpatialDimPlot(seu,
    cells.highlight = cells[[cluster]], pt.size.factor = 5, image.alpha = 0,
    cols.highlight = c("red", "grey")
  ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle(glue("Cluster {cluster}")) +
    NoLegend() +
    theme(plot.title = element_text(hjust = .5, face = "bold"))

  fn <- file.path(
    out_dir, sample_name,
    "Louvain_spatial", glue("{cluster}.png")
  )
  base_dir <- dirname(fn)
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
  }
  p <- p1 + p2
  ggsave(fn, p, width = 8, height = 4, units = "in", dpi = 300)
}
