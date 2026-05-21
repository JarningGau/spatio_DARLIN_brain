#' Spatial Niche Analysis Functions
#'
#' This script contains refactored functions for spatial niche analysis
#' including consensus clustering, spatial neighborhood analysis, and domain identification.
#'
#' @author Jarning
#' @date 2025-12-15

# Required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(RSpectra)
  library(ComplexHeatmap)
  library(ggsci)
  library(grid)
  library(circlize)
  library(dendextend)
  library(RcppML)
  library(dbscan)
  library(pbapply)
  library(ineq)
  library(ggplot2)
  library(dplyr)
  library(plyr)
  library(UpSetR)
})

# Global variable bindings to avoid R CMD check warnings
utils::globalVariables(c(
  "mean_expr", "gini", "gini_residual", "gene", "k", "score", "purity", "separation", "."
))

#' Find genes with uneven expression patterns across clusters
#'
#' @param seu Seurat object
#' @param group.by Column name for grouping cells (default: "seurat_clusters")
#' @param assay Assay to use (default: "Spatial")
#' @param min.mean.expr Minimum mean expression threshold (default: 1e-6)
#' @param gini.threshold Threshold for Gini residual (default: 0.1)
#' @return Data frame with gene statistics including Gini coefficients
#' @export
FindUnevenGenes <- function(seu,
                            group.by = "seurat_clusters",
                            assay = "Spatial",
                            min.mean.expr = 1e-6,
                            gini.threshold = 0.1) {
  # Calculate average expression per group
  avg.expr <- AverageExpression(seu, group.by = group.by, assays = assay)[[1]]
  mean.expr <- rowMeans(avg.expr)
  gini.scores <- apply(avg.expr, 1, ineq::Gini)

  result <- data.frame(
    gene = rownames(avg.expr),
    gini = gini.scores,
    mean_expr = mean.expr,
    stringsAsFactors = FALSE
  )

  # Add log mean expression and filter
  result <- result %>%
    mutate(log_mean = log10(mean_expr + 1e-8)) %>%
    filter(!is.na(gini))

  # Fit Gini curve for genes above threshold
  df <- result %>% filter(mean_expr > min.mean.expr)

  if (nrow(df) < 10) {
    warning("Too few genes above expression threshold for curve fitting")
    return(result)
  }

  # Fit nonlinear model
  tryCatch(
    {
      fit <- nls(
        gini ~ a / (1 + (mean_expr / b)^c) + d,
        data = df,
        start = list(a = 0.7, b = 0.01, c = 1, d = 0.2),
        control = list(maxiter = 500, warnOnly = TRUE)
      )

      result$gini_fitted <- predict(fit, newdata = result)
      result$gini_residual <- result$gini - result$gini_fitted
      result$gini_zscore <- scale(result$gini_residual) %>% as.numeric()
    },
    error = function(e) {
      warning("Failed to fit Gini curve: ", e$message)
      result$gini_fitted <- NA
      result$gini_residual <- NA
      result$gini_zscore <- NA
    }
  )

  return(result)
}

#' Calculate percentage expression matrix for given genes and clusters
#'
#' @param seu Seurat object
#' @param genes.use Genes to calculate percentage for
#' @param expr.cutoff Expression cutoff (default: 0)
#' @param cells.cutoff Minimum cells per cluster (default: 0)
#' @param assay Assay to use (default: "Spatial")
#' @param layer Layer to use (default: "data")
#' @return Matrix of percentage expression
#' @import Seurat
#' @export
CalPctExprMatrix <- function(seu,
                             genes.use,
                             expr.cutoff = 0,
                             cells.cutoff = 0,
                             assay = "Spatial",
                             layer = "data") {
  genes.present <- genes.use[genes.use %in% Features(seu, assay = assay)]

  clusters <- Idents(seu)
  cluster.factor <- factor(clusters)
  cluster.levels <- levels(cluster.factor)

  if (length(genes.present) > 0) {
    expr.mat <- LayerData(seu, assay = assay, layer = layer)[
      genes.present, ,
      drop = FALSE
    ]
    binary.expr <- expr.mat > expr.cutoff

    pct.present <- matrix(
      NA,
      nrow = length(genes.present),
      ncol = length(cluster.levels),
      dimnames = list(genes.present, cluster.levels)
    )
    for (cl in cluster.levels) {
      cells.in.cluster <- which(cluster.factor == cl)
      if (length(cells.in.cluster) < cells.cutoff) next
      pct.present[, cl] <- rowMeans(
        binary.expr[, cells.in.cluster, drop = FALSE]
      )
    }

    na.cols <- is.na(colSums(pct.present))
    pct.present <- pct.present[, !na.cols, drop = FALSE]
    keep.cols <- colnames(pct.present)
  } else {
    keep.cols <- cluster.levels
    pct.present <- NULL
  }

  pct.expr.mat <- matrix(
    0,
    nrow = length(genes.use),
    ncol = length(keep.cols),
    dimnames = list(genes.use, keep.cols)
  )
  if (length(genes.present) > 0 && ncol(pct.present) > 0) {
    idx <- match(genes.present, genes.use)
    pct.expr.mat[idx, ] <- as.matrix(pct.present)
  }

  return(pct.expr.mat)
}

#' Run Non-negative Matrix Factorization
#'
#' @param matrix Input matrix
#' @param k Number of factors
#' @param cores Number of cores to use (default: 5)
#' @param seed Random seed (default: 1024)
#' @param verbose Verbose output (default: FALSE)
#' @return List with W and H matrices
#' @export
RunNMF <- function(matrix, k = 5, cores = 5, seed = 1024, verbose = FALSE) {
  # Set number of threads
  RcppML::setRcppMLthreads(cores)

  # Run NMF
  model <- RcppML::nmf(matrix, k = k, verbose = verbose, seed = seed)

  # Format results
  H <- model$h
  rownames(H) <- paste0("factor_", 1:nrow(H))
  colnames(H) <- colnames(matrix)

  W <- model$w
  rownames(W) <- rownames(matrix)
  colnames(W) <- rownames(H)

  return(list(W = W, H = H, model = model))
}

#' Compute clustering statistics for different k values
#'
#' @param M Similarity matrix
#' @param hc Hierarchical clustering object
#' @param k Number of clusters
#' @param w Weight for purity vs separation (default: 0.5)
#' @return Data frame with clustering statistics
#' @export
compute_block_stats <- function(M, hc, k, w = 0.5) {
  # Cut dendrogram
  cl <- cutree(hc, k = k)

  # Group indices by cluster
  blocks <- split(seq_len(nrow(M)), cl)

  # Compute within-block purity
  purity <- sapply(blocks, function(ix) {
    m <- M[ix, ix, drop = FALSE]
    if (length(ix) == 1) {
      mean(m, na.rm = TRUE)
    } else {
      mean(m[upper.tri(m)], na.rm = TRUE)
    }
  })

  # Compute between-block separation
  separation <- sapply(seq_along(blocks), function(i) {
    ix <- blocks[[i]]

    sep_ij <- sapply(seq_along(blocks), function(j) {
      if (i == j) {
        return(NA)
      }
      jx <- blocks[[j]]
      mean(M[ix, jx], na.rm = TRUE)
    })

    max(sep_ij, na.rm = TRUE)
  })

  score <- w * purity + (1 - w) * (-separation)

  # Return results
  data.frame(
    k = k,
    cluster = names(blocks),
    purity = purity,
    separation = -separation,
    score = score,
    stringsAsFactors = FALSE
  )
}

#' Find optimal number of clusters using clustering statistics
#'
#' @param cor.mat Correlation matrix
#' @param k.range Range of k values to test (default: 10:100)
#' @param method Clustering method (default: "complete")
#' @param w Weight for scoring (default: 0.7)
#' @param use.pbapply Use progress bar (default: TRUE)
#' @return List with statistics and optimal k
#' @export
FindOptimalClusters <- function(cor.mat,
                                k.range = 10:100,
                                method = "complete",
                                w = 0.7,
                                use.pbapply = TRUE) {
  # Hierarchical clustering
  column_dend <- hclust(as.dist(1 - cor.mat), method = method)

  # Progress function
  apply_fun <- if (use.pbapply && requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pblapply
  } else {
    lapply
  }

  # Calculate statistics for each k
  clu.stats <- apply_fun(k.range, function(K) {
    compute_block_stats(cor.mat, hc = column_dend, k = K, w = w)
  }) %>% do.call(rbind, .)

  # Find optimal k based on median score
  k.summary <- clu.stats %>%
    group_by(k) %>%
    summarise(
      median_score = median(score, na.rm = TRUE),
      mean_purity = mean(purity, na.rm = TRUE),
      mean_separation = mean(separation, na.rm = TRUE),
      .groups = "drop"
    )

  optimal_k <- k.summary$k[which.max(k.summary$median_score)]

  return(list(
    stats = clu.stats,
    summary = k.summary,
    optimal_k = optimal_k,
    dendrogram = column_dend
  ))
}

#' Calculate spatial neighborhood matrix
#'
#' @param seu Seurat object
#' @param query.cells Query cells (default: all cells)
#' @param ref.cells Reference cells (default: all cells)
#' @param r Radius for neighborhood search (default: 80)
#' @param anno_col Annotation column (default: "annotation")
#' @param use_pbapply Use progress bar (default: TRUE)
#' @param return_neighbors Return neighbor information (default: FALSE)
#' @return List with neighborhood matrices and metadata
#' @export
CalNichMatrix <- function(seu,
                          query.cells = NULL,
                          ref.cells = NULL,
                          r = 80,
                          anno_col = "annotation",
                          use_pbapply = TRUE,
                          return_neighbors = FALSE) {
  # Input validation
  if (!anno_col %in% colnames(seu@meta.data)) {
    stop(sprintf("Column '%s' not found in seu@meta.data.", anno_col))
  }

  # Get coordinates and annotation
  coords <- GetTissueCoordinates(seu)

  if (is.null(ref.cells)) {
    ref.cells <- rownames(coords)
  }

  if (is.null(query.cells)) {
    query.cells <- ref.cells
  }

  # Filter query cells
  query.cells <- intersect(query.cells, ref.cells)
  if (length(query.cells) == 0) {
    stop("No query.cells found in the object.")
  }

  annotation <- seu@meta.data[[anno_col]]
  if (!is.factor(annotation)) annotation <- factor(annotation)
  type.levels <- levels(annotation)

  # Map to integers for efficient counting
  type.int <- setNames(
    match(as.character(annotation), type.levels),
    rownames(seu@meta.data)
  )

  # Neighborhood search
  Q <- coords[query.cells, , drop = FALSE]
  R <- coords[ref.cells, , drop = FALSE]
  nbrs <- dbscan::frNN(R, eps = r, query = Q, sort = FALSE)

  # Progress function
  apply_fun <- if (use_pbapply && requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pblapply
  } else {
    lapply
  }

  # Count cell types in each neighborhood
  niche.list <- apply_fun(nbrs$id, function(nn) {
    if (length(nn) == 0) {
      integer(length(type.levels))
    } else {
      tabulate(type.int[ref.cells[nn]], nbins = length(type.levels))
    }
  })

  niche.mat <- do.call(rbind, niche.list)
  rownames(niche.mat) <- rownames(Q)
  colnames(niche.mat) <- type.levels

  # Normalization
  row_sums <- rowSums(niche.mat)
  niche.mat.norm <- niche.mat
  nz <- row_sums > 0
  niche.mat.norm[nz, ] <- niche.mat[nz, , drop = FALSE] / row_sums[nz]

  # Section-level composition
  niche.mat.section <- colSums(niche.mat)
  total <- sum(niche.mat.section)

  if (total == 0) {
    section.prop <- rep(NA_real_, length(type.levels))
  } else {
    section.prop <- niche.mat.section / total
  }
  names(section.prop) <- type.levels

  # Prepare output
  out <- list(
    niche.mat = niche.mat,
    niche.mat.norm = niche.mat.norm,
    section.prop = section.prop,
    type.levels = type.levels,
    radius = r,
    n.query = length(query.cells),
    n.ref = length(ref.cells)
  )

  if (return_neighbors) {
    nn <- sapply(nbrs$id, function(id) {
      rownames(coords)[id]
    })
    names(nn) <- names(nbrs$id)
    out$neighbors <- nn
  }

  return(out)
}

#' Softmax function for rows
#'
#' @param M Matrix
#' @return Matrix with softmax applied to rows
row_softmax <- function(M) {
  M <- M - apply(M, 1, max) # numerical stability
  E <- exp(M)
  E / rowSums(E)
}

#' Compress similarity matrix to soft assignment matrix using spectral clustering
#'
#' @param S Similarity matrix
#' @param k Number of clusters
#' @param d Number of eigenvectors (default: min(2*k, 50))
#' @param tau Temperature parameter for softmax (default: 1.0)
#' @param make_psd Make matrix positive semi-definite (default: TRUE)
#' @return List with soft assignment matrix and additional information
#' @export
compress_to_A_softmax <- function(S, k, d = min(2 * k, 50), tau = 1.0, make_psd = TRUE) {
  stopifnot(is.matrix(S), nrow(S) == ncol(S))

  # Symmetrize and ensure PSD
  S2 <- (S + t(S)) / 2
  diag(S2) <- 1

  if (make_psd) {
    S2 <- as.matrix(Matrix::nearPD(S2, corr = TRUE)$mat)
  }

  # Spectral decomposition
  eig <- RSpectra::eigs_sym(S2, k = d, which = "LA")
  vals <- Re(eig$values)
  vecs <- Re(eig$vectors)

  # Keep positive eigenvalues
  pos <- vals > 1e-10
  X <- vecs[, pos, drop = FALSE] %*% diag(sqrt(vals[pos]))

  # K-means clustering
  km <- kmeans(X, centers = k, nstart = 20)

  # Compute squared distances to centers
  x2 <- rowSums(X^2)
  c2 <- rowSums(km$centers^2)
  D2 <- outer(x2, c2, "+") - 2 * (X %*% t(km$centers))

  # Apply softmax
  A <- row_softmax(-D2 / tau)

  # Add names
  if (!is.null(rownames(S))) rownames(A) <- rownames(S)
  colnames(A) <- paste0("cluster", seq_len(k))

  return(list(
    A = A,
    embedding = X,
    kmeans = km,
    D2 = D2,
    eigenvalues = vals[pos],
    eigenvectors = vecs[, pos, drop = FALSE]
  ))
}

#' Wrapper function for calculating consensus neighborhood matrix
#'
#' @param seu Seurat object
#' @param A Soft assignment matrix
#' @param section.prefix Section prefix for column names
#' @param r Radius for neighborhood search (default: 12)
#' @param anno_col Annotation column (default: "iter_clusters.r1")
#' @return Consensus neighborhood matrix
#' @export
CNM_wrapper <- function(seu, A, section.prefix = "s1_", r = 12, anno_col = "iter_clusters.r1") {
  # Calculate neighborhood matrix
  niche <- CalNichMatrix(seu,
    query.cells = Cells(seu),
    ref.cells = NULL,
    r = r,
    anno_col = anno_col
  )

  niche.mat <- niche$niche.mat
  colnames(niche.mat) <- paste0(section.prefix, colnames(niche.mat))

  # Matrix multiplication with soft assignment
  available_cols <- intersect(colnames(niche.mat), rownames(A))
  if (length(available_cols) == 0) {
    stop("No matching columns between niche matrix and assignment matrix")
  }

  N <- niche.mat[, available_cols, drop = FALSE] %*% A[available_cols, , drop = FALSE]

  return(N)
}

#' Create consensus heatmap with annotations and cluster boundaries
#'
#' @param cor.mat Correlation matrix
#' @param sections Section labels
#' @param clusters Cluster labels
#' @param K Number of clusters for cutting dendrogram
#' @param method Clustering method (default: "complete")
#' @param colors Color palette for sections
#' @param add_boundaries Whether to add cluster boundaries (default: TRUE)
#' @param heatmap_name Name of the heatmap body for boundaries (default: "PCC")
#' @return ComplexHeatmap object with boundaries added if requested
#' @export
CreateConsensusHeatmap <- function(cor.mat, sections, clusters, K,
                                   method = "complete", colors = NULL,
                                   add_boundaries = TRUE, heatmap_name = "PCC") {
  # Generate colors if not provided
  if (is.null(colors)) {
    uniq_sect <- unique(sections)
    n_sect <- length(uniq_sect)
    pal <- pal_d3("category10")
    colors <- pal(min(n_sect, 10))
    names(colors) <- uniq_sect
  }

  # Hierarchical clustering
  column_dend <- hclust(as.dist(1 - cor.mat), method = method)

  # Column annotation
  ha_col <- HeatmapAnnotation(
    Section = sections,
    Cluster = clusters,
    col = list(Section = colors),
    show_annotation_name = TRUE,
    show_legend = c(Section = TRUE, Cluster = FALSE),
    annotation_name_side = "right",
    annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
    annotation_height = unit(5, "mm")
  )

  # Create heatmap
  ht <- Heatmap(
    matrix = cor.mat,
    name = heatmap_name,
    col = colorRamp2(c(-1, 0, 1), c("#2E349B", "white", "#E53935")),
    border = TRUE,
    cluster_rows = column_dend,
    cluster_columns = column_dend,
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = ha_col
  )

  # If boundaries are not requested, return simple result
  if (!add_boundaries) {
    return(list(heatmap = ht, dendrogram = column_dend))
  }

  # For boundaries, we need to draw and then decorate
  # Draw heatmap first
  ht_drawn <- draw(ht)

  # Add cluster boundaries
  decorate_heatmap_body(heatmap_name, {
    # Use the original dendrogram directly since we know the clustering
    tree <- as.dendrogram(column_dend)
    ind <- cutree(column_dend, k = K)[order.dendrogram(tree)]

    first_index <- function(l) which(l)[1]
    last_index <- function(l) {
      x <- which(l)
      x[length(x)]
    }

    clusters_names <- names(table(ind))
    x1 <- sapply(clusters_names, function(x) first_index(ind == x)) - 1
    x2 <- sapply(clusters_names, function(x) last_index(ind == x))
    x1 <- x1 / length(ind)
    x2 <- x2 / length(ind)

    # Draw rectangles
    grid.rect(
      x = x1, width = (x2 - x1), y = 1 - x1, height = (x1 - x2),
      hjust = 0, vjust = 0, default.units = "npc",
      gp = gpar(fill = NA, col = "black", lwd = 2)
    )

    # Add labels
    x_center <- (x1 + x2) / 2
    y_center <- 1 - x_center
    grid.text(
      label = clusters_names,
      x = x_center,
      y = y_center,
      default.units = "npc",
      hjust = -1.5,
      vjust = 0.5,
      gp = gpar(fontsize = 10, fontface = "bold")
    )
  })

  return(list(heatmap = ht, heatmap_drawn = ht_drawn, dendrogram = column_dend))
}

#' Identify spatial domains using k-means clustering
#'
#' @param niche.mat.norm Normalized neighborhood matrix
#' @param K Number of initial clusters (default: 100)
#' @param Kh Number of final hierarchical clusters (default: 30)
#' @param iter.max Maximum iterations for k-means (default: 20)
#' @param method Hierarchical clustering method (default: "complete")
#' @return List with clustering results and domain assignments
#' @export
IdentifySpatialDomains <- function(niche.mat.norm, K = 100, Kh = 30,
                                   iter.max = 20, method = "complete") {
  # K-means clustering
  kmeans_result <- kmeans(niche.mat.norm, centers = K, iter.max = iter.max)

  # Calculate average profiles
  niche.avg <- sapply(1:K, function(cl) {
    colMeans(niche.mat.norm[kmeans_result$cluster == cl, , drop = FALSE])
  })
  colnames(niche.avg) <- paste0("D", 1:K)

  # Hierarchical clustering of averages
  cor.mat <- cor(niche.avg)
  column_dend <- hclust(as.dist(1 - cor.mat), method = method)
  clusters <- cutree(column_dend, k = Kh)
  names(clusters) <- column_dend$labels

  # Map back to cells
  collapsed.clusters <- paste0("D", clusters[kmeans_result$cluster])
  names(collapsed.clusters) <- names(kmeans_result$cluster)

  return(list(
    kmeans = kmeans_result,
    niche.avg = niche.avg,
    cor.mat = cor.mat,
    dendrogram = column_dend,
    hierarchical_clusters = clusters,
    domain_assignments = collapsed.clusters,
    K = K,
    Kh = Kh
  ))
}

#' Generate spatial plots for clusters or domains
#'
#' @param seu.list List of Seurat objects
#' @param cluster.name Name of the cluster/domain to plot
#' @param ident.name Identity name in Seurat objects
#' @param section.names Names for sections (default: c("Section1", "Section2", "Section3"))
#' @param pt.size.factor Point size factor (default: 5)
#' @param image.alpha Image transparency (default: 0)
#' @param highlight.colors Colors for highlighting (default: c("red", "grey"))
#' @return Combined ggplot object
#' @export
PlotSpatialCluster <- function(seu.list, cluster.name, ident.name,
                               section.names = c("Section1", "Section2", "Section3"),
                               pt.size.factor = 5, image.alpha = 0,
                               highlight.colors = c("red", "grey")) {
  plots <- list()
  n_sections <- length(seu.list)

  for (i in seq_len(n_sections)) {
    seu <- seu.list[[i]]
    Idents(seu) <- ident.name

    cells <- CellsByIdentities(seu)

    if (cluster.name %in% names(cells)) {
      highlight_cells <- cells[[cluster.name]]
    } else {
      highlight_cells <- character(0)
    }

    section_label <- if (i <= length(section.names)) {
      section.names[i]
    } else {
      paste0("Section", i)
    }
    if (length(highlight_cells) > 0) {
      colors.use <- highlight.colors
    } else {
      colors.use <- "grey"
    }
    p <- SpatialDimPlot(seu,
      cells.highlight = highlight_cells,
      pt.size.factor = pt.size.factor,
      image.alpha = image.alpha,
      cols.highlight = colors.use
    ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      ggtitle(paste0(cluster.name, " (", section_label, ")")) +
      NoLegend() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    plots[[i]] <- p
  }

  if (length(plots) == 0L) {
    return(NULL)
  }
  if (length(plots) == 1L) {
    return(plots[[1]])
  }

  # Combine plots horizontally (any number of sections)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(plots, nrow = 1))
  }
  plots
}

#' Complete spatial niche analysis pipeline
#'
#' @param seu.list List of Seurat objects
#' @param section.prefixes Section prefixes (default: c("s1_", "s2_", "s3_"))
#' @param group.by Grouping variable (default: "iter_clusters.r1")
#' @param gini.threshold Gini residual threshold (default: 0.1)
#' @param nmf.k Number of NMF factors (default: 50)
#' @param consensus.k Number of consensus clusters (default: 93)
#' @param softmax.k Number of softmax clusters (default: 39)
#' @param softmax.tau Softmax temperature (default: 0.3)
#' @param niche.radius Neighborhood radius (default: 12)
#' @param domain.K Initial number of domains (default: 100)
#' @param domain.Kh Final number of domains (default: 30)
#' @param cores Number of cores for NMF (default: 5)
#' @return List with all analysis results
#' @export
RunSpatialNicheAnalysis <- function(seu.list,
                                    section.prefixes = c("s1_", "s2_", "s3_"),
                                    group.by = "iter_clusters.r1",
                                    gini.threshold = 0.1,
                                    nmf.k = 50,
                                    consensus.k = 93,
                                    softmax.k = 39,
                                    softmax.tau = 0.3,
                                    niche.radius = 12,
                                    domain.K = 100,
                                    domain.Kh = 30,
                                    cores = 5) {
  message("Starting spatial niche analysis pipeline...")

  # Step 1: Find uneven genes
  message("Step 1: Finding uneven genes...")
  uneven.genes.list <- lapply(seu.list, function(seu) {
    FindUnevenGenes(seu, group.by = group.by)
  })

  # Step 2: Get common features
  message("Step 2: Identifying common features...")
  features.list <- lapply(uneven.genes.list, function(x) {
    x %>%
      filter(gini_residual >= gini.threshold) %>%
      pull(gene)
  })

  # Find intersection
  features.stat <- UpSetR::fromList(features.list)
  rownames(features.stat) <- unique(unlist(features.list))
  features.use <- rownames(features.stat)[rowSums(features.stat) >= length(seu.list)]

  message(sprintf("Found %d common features", length(features.use)))

  # Step 3: Calculate percentage expression matrices
  message("Step 3: Calculating expression matrices...")
  pct.expr.mats <- lapply(seq_along(seu.list), function(i) {
    seu <- seu.list[[i]]
    Idents(seu) <- as.vector(seu[[group.by, drop = TRUE]])
    mat <- CalPctExprMatrix(seu, features.use, 0, 0)
    colnames(mat) <- paste0(section.prefixes[i], colnames(mat))
    return(mat)
  })

  pct.expr.mat <- do.call(cbind, pct.expr.mats)

  # Step 4: Run NMF
  message("Step 4: Running NMF...")
  nmf.result <- RunNMF(pct.expr.mat, k = nmf.k, cores = cores)

  # Step 5: Consensus clustering
  message("Step 5: Performing consensus clustering...")
  cor.mat <- cor(nmf.result$H)

  # Find optimal clusters
  optimal.result <- FindOptimalClusters(cor.mat, w = 0.7)

  # Use specified k or optimal k
  final.k <- ifelse(is.null(consensus.k), optimal.result$optimal_k, consensus.k)

  clusters <- cutree(optimal.result$dendrogram, k = final.k)
  clusters.df <- data.frame(
    sub.id = names(clusters),
    cluster = paste0("C", clusters),
    stringsAsFactors = FALSE
  )

  # Step 6: Softmax compression
  message("Step 6: Creating softmax assignment matrix...")
  softmax.result <- compress_to_A_softmax(cor.mat, k = softmax.k, tau = softmax.tau)
  A <- softmax.result$A
  A[A < 0.05] <- 0
  A <- t(apply(A, 1, function(x) x / sum(x)))

  # Step 7: Calculate consensus neighborhood matrices
  message("Step 7: Calculating neighborhood matrices...")
  niche.mats <- lapply(seq_along(seu.list), function(i) {
    CNM_wrapper(seu.list[[i]], A, section.prefixes[i], r = niche.radius, anno_col = group.by)
  })

  niche.mat <- do.call(rbind, niche.mats)
  niche.mat.norm <- t(apply(niche.mat, 1, function(z) z / sum(z)))

  # Step 8: Identify spatial domains
  message("Step 8: Identifying spatial domains...")
  domain.result <- IdentifySpatialDomains(niche.mat.norm, K = domain.K, Kh = domain.Kh)

  # Step 9: Add annotations to Seurat objects
  message("Step 9: Adding annotations to Seurat objects...")
  for (i in seq_along(seu.list)) {
    # Add consensus clusters
    seu.list[[i]]$consensus_cluster <- plyr::mapvalues(
      x = seu.list[[i]]@meta.data[[group.by]],
      from = sub(paste0("^", section.prefixes[i]), "", clusters.df$sub.id),
      to = clusters.df$cluster,
      warn_missing = FALSE
    )
    seu.list[[i]]$consensus_cluster <- factor(seu.list[[i]]$consensus_cluster,
      levels = paste0("C", 1:final.k)
    )

    # Add spatial domains
    z <- domain.result$domain_assignments[rownames(seu.list[[i]]@meta.data)]
    names(z) <- rownames(seu.list[[i]]@meta.data)
    seu.list[[i]]$spatial_domain <- z
  }

  message("Analysis complete!")

  # Return comprehensive results
  return(list(
    seu.list = seu.list,
    uneven.genes = uneven.genes.list,
    features.use = features.use,
    pct.expr.mat = pct.expr.mat,
    nmf.result = nmf.result,
    cor.mat = cor.mat,
    optimal.clustering = optimal.result,
    consensus.clusters = clusters.df,
    softmax.result = softmax.result,
    A = A,
    niche.mat = niche.mat,
    niche.mat.norm = niche.mat.norm,
    domain.result = domain.result,
    parameters = list(
      gini.threshold = gini.threshold,
      nmf.k = nmf.k,
      consensus.k = final.k,
      softmax.k = softmax.k,
      softmax.tau = softmax.tau,
      niche.radius = niche.radius,
      domain.K = domain.K,
      domain.Kh = domain.Kh
    )
  ))
}
