## Generate PCA embeddings for clone2vec analysis
setwd("scripts")

library(tidyverse)
library(Seurat)
library(glue)

source("../R/helpers.R")

### Parameters
out_dir <- "../data/mRNA/"

### No batch correction
seu_list <- list(
  L126_Brain_s1 = qs::qread("../data/mRNA/consensus/L126_Brain_s1_seurat_object.qs"),
  L126_Brain_s2 = qs::qread("../data/mRNA/consensus/L126_Brain_s2_seurat_object.qs"),
  L126_Brain_s3 = qs::qread("../data/mRNA/consensus/L126_Brain_s3_seurat_object.qs")
)

for (i in seq_along(seu_list)) {
  sn <- names(seu_list)[i]
  seu_list[[i]] <- RenameCells(seu_list[[i]], add.cell.id = sn)
}
head(Cells(seu_list$L126_Brain_s1))

seu <- merge(seu_list[[1]], seu_list[-1])
seu <- NormalizeData(seu)
VariableFeatures(seu) <- readLines("../data/mRNA/consensus/variable_features.txt")
seu <- ScaleData(seu)
seu <- RunPCA(seu)
LayerData(seu, layer = "scale.data") <- NULL
seu <- RunUMAP(seu, reduction = "pca", dims = 1:50)

p <- DimPlot(seu, reduction = "umap", group.by = "consensus_cluster", label = TRUE) +
  NoLegend() + NoAxes()
fn <- file.path(out_dir, "latent_embeddings", "umap_consensus_clusters.png")
base_dir <- dirname(fn)
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}
ggsave(fn, p, width = 5, height = 5, units = "in", dpi = 300)

p <- DimPlot(seu,
  reduction = "umap", group.by = "sample_name",
  raster = FALSE, pt.size = 0.1, alpha = 0.1
) +
  NoLegend() +
  NoAxes() +
  ggsci::scale_color_d3()
fn <- file.path(out_dir, "latent_embeddings", "umap_sample_name.png")
ggsave(fn, p, width = 5, height = 5, units = "in", dpi = 300)

emb <- Embeddings(seu, reduction = "pca")
emb <- round(emb, 3)
head(emb)
emb <- as.data.frame(emb) %>%
  mutate(cell_id = rownames(.), .before = 1)
fn <- file.path(out_dir, "latent_embeddings", "L126_Brain_all_pca_embeddings.tsv")
write_tsv(emb, fn)
R.utils::gzip(fn, overwrite = TRUE)

##
seu <- qs::qread("../data/mRNA/consensus/L0927_Brain_seurat_object.qs")
emb <- Embeddings(seu, reduction = "pca")
emb <- round(emb, 3)
head(emb)
emb <- as.data.frame(emb) %>%
  mutate(cell_id = rownames(.), .before = 1)
fn <- file.path(out_dir, "latent_embeddings", "L0927_Brain_pca_embeddings.tsv")
write_tsv(emb, fn)
R.utils::gzip(fn, overwrite = TRUE)
