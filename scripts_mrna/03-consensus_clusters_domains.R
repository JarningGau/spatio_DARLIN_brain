setwd("scripts_mrna")

library(tidyverse)
library(Seurat)
library(glue)
library(ggVennDiagram)

source("../R/niche.R")
source("../R/helpers.R")

### Parameters
out_dir <- "../data/mRNA/"

# Step 1: uneven genes
group_by <- "iter_clusters.r1"
gini_threshold <- 0.1

# Step 2: expression matrices
expr_cutoff <- 0
cells_cutoff <- 0

# Step 3: NMF
nmf_k <- 50
nmf_cores <- 5

# Step 4: optimal clusters
k_range <- 10:100
cluster_weight <- 0.7

# Step 5–6: consensus clustering
consensus_k <- 75

# Step 7–8: spatial domains
softmax_k <- 75
softmax_tau <- 0.3
affiliation_threshold <- 0.05
niche_radius <- 12 # pixels; 1 pixel = 6.8 um
domain_K <- 200
domain_Kh <- 40

### Step 1: Find Uneven Genes

seu_list <- list(
  L0927_Brain = qs::qread("../data/mRNA/L0927_Brain/L0927_Brain_Seurat_cluster.qs"),
  L126_Brain_s1 = qs::qread("../data/mRNA/L126_Brain_s1/L126_Brain_s1_Seurat_cluster.qs"),
  L126_Brain_s2 = qs::qread("../data/mRNA/L126_Brain_s2/L126_Brain_s2_Seurat_cluster.qs"),
  L126_Brain_s3 = qs::qread("../data/mRNA/L126_Brain_s3/L126_Brain_s3_Seurat_cluster.qs")
)

seu_list <- lapply(seu_list, NormalizeData)
sapply(seu_list, dim)

uneven_genes_list <- lapply(seu_list, function(xx) {
  FindUnevenGenes(xx, group.by = group_by)
})
head(uneven_genes_list[[1]])

features <- lapply(uneven_genes_list, function(xx) {
  xx %>%
    filter(gini_residual >= gini_threshold) %>%
    pull(gene)
})

p <- ggVennDiagram(features,
  label_alpha = 0,
  edge_alpha = 0.5,
  set_name_size = 5
) +
  scale_fill_gradient(low = "lightblue", high = "steelblue") +
  theme_void() +
  ggtitle("Overlap of Uneven Genes Across Samples")
fn <- file.path(out_dir, "consensus", "uneven_genes_venn.png")
base_dir <- dirname(fn)
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}
ggsave(fn, p, width = 8, height = 6, units = "in", dpi = 300, bg = "white")

features_stat <- UpSetR::fromList(features)
rownames(features_stat) <- unique(unlist(features))
table(rowSums(features_stat))
features_use <- rownames(features_stat)[rowSums(features_stat) >= length(seu_list)]
cat("Number of common uneven genes:", length(features_use), "\n")

fn <- file.path(out_dir, "consensus", "variable_features.txt")
writeLines(features_use, fn)

### Step 2: Calculate Expression Matrices
pct_expr_mats <- lapply(seq_along(seu_list), function(i) {
  xx <- seu_list[[i]]
  Idents(xx) <- group_by
  mat <- CalPctExprMatrix(xx, features_use, expr_cutoff, cells_cutoff)
  sn <- names(seu_list)[i]
  colnames(mat) <- paste0(sn, "-", colnames(mat))
  return(mat)
})
head(colnames(pct_expr_mats[[1]]))
pct_expr_mat <- do.call(cbind, pct_expr_mats)
dim(pct_expr_mat)

### Step 3: Non-negative Matrix Factorization
nmf_result <- RunNMF(pct_expr_mat, k = nmf_k, cores = nmf_cores)
W <- nmf_result$W
H <- nmf_result$H
cat("NMF completed.\nW dimensions:", dim(W), "\nH dimensions:", dim(H), "\n")
## Remove unexpected program
# H <- H[-which.max(W["Sltm", ]), ]

### Step 4: Consensus Clustering Analysis
# Calculate correlation matrix
cor_mat <- cor(H)

# Find optimal number of clusters
optimal_result <- FindOptimalClusters(cor_mat, k.range = k_range, w = cluster_weight)

# Plot clustering statistics
p1 <- ggplot(optimal_result$stats, aes(factor(k), score)) +
  geom_boxplot() +
  labs(
    x = "Number of Clusters", y = "Score",
    title = "Clustering Score (weighted sum of purity and separation)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(optimal_result$stats, aes(factor(k), purity)) +
  geom_boxplot() +
  labs(x = "Number of Clusters", y = "Purity", title = "Clustering Purity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(optimal_result$stats, aes(factor(k), separation)) +
  geom_boxplot() +
  labs(x = "Number of Clusters", y = "Separation", title = "Clustering Separation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p <- p1 / p2 / p3
fn <- file.path(out_dir, "consensus", "consensus_clustering_statistics.png")
ggsave(fn, p, width = 12, height = 8, units = "in", dpi = 300)

### Step 5: Consensus Heatmap Visualization
K <- consensus_k
sections <- sapply(strsplit(rownames(cor_mat), "-"), function(xx) xx[1])
clusters <- cutree(optimal_result$dendrogram, k = K)
clusters <- paste0("C", clusters)

# Create heatmap
fn <- file.path(out_dir, "consensus", glue("consensus_clustering_heatmap_k{K}.png"))
png(fn, width = 12, height = 10, units = "in", res = 300)
heatmap_result <- CreateConsensusHeatmap(cor_mat, sections, clusters, K)
dev.off()

### Step 6: Spatial Pattern Visualization

# Map clusters back to Seurat objects
clusters_df <- data.frame(
  sub.id = names(cutree(optimal_result$dendrogram, k = K)),
  cluster = paste0("C", cutree(optimal_result$dendrogram, k = K)),
  stringsAsFactors = FALSE
)
head(clusters_df)

for (i in seq_along(seu_list)) {
  sn <- names(seu_list)[i]
  seu_list[[i]]$consensus_cluster <- plyr::mapvalues(
    x = paste0(sn, "-", seu_list[[i]][[group_by, drop = TRUE]]),
    from = clusters_df$sub.id,
    to = clusters_df$cluster,
    warn_missing = FALSE
  )
  seu_list[[i]]$consensus_cluster <- factor(
    seu_list[[i]]$consensus_cluster,
    levels = paste0("C", 1:K)
  )
  Idents(seu_list[[i]]) <- "consensus_cluster"
}

for (cluster in paste0("C", 1:K)) {
  message(glue("Plotting {cluster} ..."))
  p <- PlotSpatialCluster(seu_list, cluster, "consensus_cluster",
    section.names = names(seu_list)
  )
  fn <- file.path(out_dir, "consensus", paste0("spatial_clusters_k", K), paste0(cluster, ".png"))
  base_dir <- dirname(fn)
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
  }
  ggsave(fn, p, width = 16, height = 5, units = "in", dpi = 100)
}

## Step7: Spatial Domain Analysis
softmax_result <- compress_to_A_softmax(cor_mat, k = softmax_k, tau = softmax_tau)
A <- softmax_result$A
A[A < affiliation_threshold] <- 0
A <- t(apply(A, 1, function(x) x / sum(x)))
dim(A)
dim(cor_mat)

## r is pixels radius, 1 pixel = 6.8 um
niche_radius <- 20
niche_mat_list <- lapply(seq_along(seu_list), function(i) {
  sn <- names(seu_list)[i]
  niche_mat <- CNM_wrapper(
    seu = seu_list[[i]],
    A = A,
    section.prefix = paste0(sn, "-"),
    r = niche_radius,
    anno_col = group_by
  )
  rownames(niche_mat) <- paste0(sn, "-", rownames(niche_mat))
  niche_mat
})

niche_mat <- do.call(rbind, niche_mat_list)
dim(niche_mat)
head(rownames(niche_mat))
head(colnames(niche_mat))

niche_mat_norm <- t(apply(niche_mat, 1, function(z) z / sum(z)))
dim(niche_mat_norm)

# Identify spatial domains
domain_K <- 400
domain_Kh <- 30
domain_result <- IdentifySpatialDomains(niche_mat_norm, K = domain_K, Kh = domain_Kh)

# Visualize domain correlation heatmap
fn <- file.path(out_dir, "consensus", glue("consensus_domain_heatmap_k{domain_result$Kh}.png"))
png(fn, width = 12, height = 10, units = "in", res = 300)
heatmap_domains <- CreateConsensusHeatmap(
  domain_result$cor.mat,
  sections = rep("Domain", ncol(domain_result$cor.mat)),
  clusters = paste0("D", domain_result$hierarchical_clusters),
  K = domain_result$Kh
)
dev.off()

### Step 8: Spatial Domain Visualization
# Add domain annotations to Seurat objects
for (i in seq_along(seu_list)) {
  sn <- names(seu_list)[i]
  cells_ordered <- rownames(seu_list[[i]]@meta.data)
  cells_ordered_add_prefix <- paste0(sn, "-", cells_ordered)
  z <- domain_result$domain_assignments[cells_ordered_add_prefix]
  names(z) <- cells_ordered
  seu_list[[i]]$spatial_domain <- z
}

# Plot spatial domains
domain_list <- paste0("D", 1:domain_result$Kh)

for (domain in domain_list) {
  message(glue("Plotting {domain} ..."))
  p <- PlotSpatialCluster(seu_list, domain, "spatial_domain",
    section.names = names(seu_list)
  )
  fn <- file.path(out_dir, "consensus", paste0("spatial_domains_k", domain_result$Kh), paste0(domain, ".png"))
  base_dir <- dirname(fn)
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
  }
  ggsave(fn, p, width = 16, height = 5, units = "in", dpi = 100)
}

## save seurat objects
for (i in seq_along(seu_list)) {
  sn <- names(seu_list)[i]
  seu <- seu_list[[i]]
  seu[["Spatial"]]@layers$data <- NULL
  fn <- file.path(out_dir, "consensus", paste0(sn, "_seurat_object.qs"))
  save_qs(seu, fn)
}
