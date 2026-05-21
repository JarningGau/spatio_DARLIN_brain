find_degs = function(seu, nn.idents, assay=NULL) {
  if (is.null(assay)) {
    assay = DefaultAssay(seu)
  }
  cells.use = CellsByIdentities(seu, idents = nn.idents)
  cells.1 = cells.use[[nn.idents[1]]]
  cells.2 = cells.use[[nn.idents[2]]]
  if(length(cells.1) < 3 || length(cells.2) < 3) {
    return(NULL)
  }
  subsample.cells = function(cells, size=1000) {
    if (length(cells) > size) {
      cells = sample(cells, size = size)
    }
    return(cells)
  }
  cells.1 = subsample.cells(cells.1)
  cells.2 = subsample.cells(cells.2)
  markers = FindMarkers(seu[[assay]], cells.1 = cells.1, cells.2 = cells.2, test.use = "wilcox")
  if (nrow(markers) == 0) {
    return(NULL)
  }
  markers = markers %>% 
    mutate(pct.fg = abs(pct.1 - pct.2) / ifelse(pct.1 > pct.2, pct.1, pct.2)) %>% 
    filter(abs(avg_log2FC) > 1, p_val_adj < 1e-5, pct.fg > 0.5, pct.1 > 0.3 | pct.2 > 0.3)
  return(markers)
}


## Usage:
## seu2 = merge_clusters(seu, cluster.name = "RNA_snn_res.0.8", nn.idents = c("2","1","3"))
## table(seu$RNA_snn_res.0.8)
## table(Idents(seu2))
merge_clusters = function(seu, cluster.name, nn.idents) {
  Idents(seu) = cluster.name
  cells.use = CellsByIdentities(seu, idents = nn.idents)
  renamed.ident = names(sort(sapply(cells.use, length), decreasing = T))[1]
  cells.use = do.call(c, cells.use)
  names(cells.use) = NULL
  idents = as.character(Idents(seu))
  names(idents) = rownames(seu@meta.data)
  idents[cells.use] = renamed.ident
  seu[[cluster.name]] = factor(idents)
  Idents(seu) = cluster.name
  return(seu)
}



drop_small_clusters = function(seu, cluster.name, min.cells = 5) {
  idents.tmp = seu@meta.data[, cluster.name]
  subcluster.cells = table(idents.tmp)
  small.clusters = names(subcluster.cells)[subcluster.cells <= min.cells]
  cells.drop = rownames(seu@meta.data)[idents.tmp %in% small.clusters]
  cells.use = setdiff(rownames(seu@meta.data), cells.drop)
  # print(cells.use)
  if (length(cells.use) == 0) {
    return(NULL)
  }
  seu = subset(seu, cells = cells.use)
  ## drop the levels in factor
  idents.tmp = seu@meta.data[, cluster.name]
  if (class(idents.tmp) == "factor") {
    seu[[cluster.name]] = droplevels(idents.tmp)
  } else {
    seu[[cluster.name]] = factor(idents.tmp)
  }
  Idents(seu) = cluster.name
  return(seu)
}


## Usage:
## seu2 = initial_cluster(seu, cluster.name = "RNA_snn_res.1", num.degs.cutoff = 5)
## table(seu2$RNA_snn_res.1)
## DimPlot(seu, group.by = "RNA_snn_res.1", label=T) + DimPlot(seu2, group.by = "RNA_snn_res.1", label=T)
initial_cluster = function(seu, cluster.name = NULL, num.degs.cutoff = 0) {
  if (is.null(cluster.name) || !cluster.name %in% names(seu@meta.data)) {
    stop(glue("{cluster.name} is not exist, please do cluster using Seurat::FindClusters() first."))
  }
  Idents(seu) = cluster.name
  ## find markers for all paired clusters
  message("Calculating DEGs ...")
  cluster.ids = sort(unique(Idents(seu))) %>% as.character()
  if (length(cluster.ids) < 2) {
    return(seu)
  }
  ## generate a matrix, rows and cols are cluster.ids, values are the number of DEGs
  dist.mat = matrix(0, nrow = length(cluster.ids), ncol = length(cluster.ids))
  rownames(dist.mat) = cluster.ids
  colnames(dist.mat) = cluster.ids
  ## dist.mat is a symmetric matrix, so we only need to calculate the upper triangle
  N = length(cluster.ids)
  for (i in 1:(N-1)) {
    ## progress bar
    message(glue("Calculating DEGs {i}/{N-1} ..."))
    for (j in (i + 1):N) {
      markers = find_degs(seu, c(cluster.ids[i], cluster.ids[j]))
      if (is.null(markers)) {
        dist.mat[i, j] = 0
        dist.mat[j, i] = 0
      } else {
        dist.mat[i, j] = nrow(markers)
        dist.mat[j, i] = nrow(markers)  
      }
    }
  }
  ## transfer the dist.mat to a weight.mat
  ones_like = function(x) {
    y = matrix(1, nrow = nrow(x), ncol = ncol(x))
    rownames(y) = rownames(x)
    colnames(y) = colnames(x)
    y
  }
  weight.mat = ones_like(dist.mat)
  weight.mat[dist.mat >= num.degs.cutoff] = 0
  ## transfer the weight.mat to a graph using igraph
  graph = igraph::graph_from_adjacency_matrix(weight.mat, mode = "undirected", weighted = TRUE)
  ## find the connected components
  components = igraph::components(graph)
  membership = unique(components$membership)
  components = lapply(membership, function(xx) {
    names(components$membership)[components$membership == xx]
  })
  names(components) = membership
  ## merge the clusters in the same connected component if the node number is larger than 1
  for (comp in components) {
    if (length(comp) > 1) {
      info = paste0(comp, collapse = ",")
      message(glue("Merging {info} ..."))
      seu = merge_clusters(seu, cluster.name, comp)
    }
  }
  return(seu)
}


## Usage:
## seu2 = initial_cluster(seu, cluster.name = "RNA_snn_res.0.2", num.degs.cutoff = 0)
## seu3 = iter_cluster(seu, graph.name = "RNA_snn", cluster.name = "RNA_snn_res.1", num.degs.cutoff=5, cells.cutoff = 50, max.round = 4, max.clusters = 200)
iter_cluster = function(seu, graph.name, cluster.name, num.degs.cutoff=1, cells.cutoff=50, max.round=4, max.clusters=200) {
  Idents(seu) = cluster.name
  round = 1
  # the loop terminate when:
  # 1. reach max.round
  # 2. reach max.clusters
  # 3. no new clusters found
  while(TRUE) {
    if (round > max.round) {
      message(glue("Reach max.round: {max.round}, break out."))
      break
    }
    message(glue("========== Round {round} =========="))
    cluster.ids = sort(unique(Idents(seu)))
    subcluster.name = glue("iter_clusters.r{round}")
    num.cl.pre = length(cluster.ids)
    message(glue("Number of clusters: {num.cl.pre}"))
    
    new.idents = as.character(Idents(seu))
    names(new.idents) = rownames(seu@meta.data)
    
    for (cl in cluster.ids) {
      message(glue("Louvain clustering on {cl} ..."))
      
      # cl = "0"
      # graph.name = "Spatial_snn"
      # subcluster.name = "iter_clusters.r1"
      # num.degs.cutoff = 1
      
      cells.use = CellsByIdentities(seu, idents = cl) %>% unlist()
      if (length(cells.use) <= cells.cutoff) {
        message(glue("Too less cells in {cl}, skip."))
        next
      }
      seu = FindSubCluster(seu, cluster = cl, graph.name = graph.name, subcluster.name = subcluster.name, resolution = 0.4)
      
      ## only check whether to merge sub clusters
      seu.sub = subset(seu, idents = cl)
      
      ## deal with small cluters (cells less than cells.cutoff)
      seu.sub = drop_small_clusters(seu.sub, subcluster.name, min.cells = cells.cutoff)
      if(is.null(seu.sub)) {
        next
      }
      seu.sub = initial_cluster(seu.sub, cluster.name = subcluster.name, num.degs.cutoff = num.degs.cutoff)
      new.idents.sub = as.character(seu.sub@meta.data[, subcluster.name])
      names(new.idents.sub) = rownames(seu.sub@meta.data)
      new.idents[names(new.idents.sub)] = new.idents.sub
      seu[[subcluster.name]] = factor(new.idents)
      Idents(seu) = subcluster.name
    }
    ## if now more new clusters, break out
    cluster.cells = table(Idents(seu))
    small.clusters = names(cluster.cells)[cluster.cells <= cells.cutoff]
    valid.clusters = setdiff(unique(Idents(seu)), small.clusters)
    num.cl.post = length(valid.clusters)
    if (num.cl.post < num.cl.pre) {
      message("No more clusters were found, break out.")
      break
    }
    if (num.cl.post > max.clusters) {
      message(glue("Reach max.clusters: {max.clusters}, break out."))
      break
    }
    round = round + 1
  }
  return(seu)
}
