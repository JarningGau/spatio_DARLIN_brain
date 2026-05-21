library(glue)
HOME = here::here()
source(glue("{HOME}/utils/R/helpers.R"))


MakeSegCells = function(seu, K=20) {
  coords = GetTissueCoordinates(seu)
  kmeans = kmeans(coords, centers = K)
  seg.cells = lapply(1:K, function(i) {
    names(kmeans$cluster)[kmeans$cluster == i]
  })
  names(seg.cells) = 1:K
  return(seg.cells)
}


MakeRCTDTask = function(reference, seu, seg.cells, tmp.dir) {
  print_with_time("Preparing coordinates...")
  # Prepare coordinates
  coords = GetTissueCoordinates(seu) %>% as.matrix()
  colnames(coords) = paste0("Spatial_", 1:2)
  
  # Create RCTD task for each segment
  for (seg in names(seg.cells)) {
    print_with_time(glue("Processing segment {seg}/{length(seg.cells)}..."))
    cells = seg.cells[[seg]]
    seu.q = CreateSeuratObject(LayerData(seu, layer = "counts")[, cells], meta.data = seu@meta.data)
    seu.q[["spatial"]] = CreateDimReducObject(coords[cells, ], assay = "RNA")
    
    print_with_time("Creating SpatialRNA object...")
    query = SpatialRNA(
      coords = Embeddings(seu.q, reduction = "spatial") %>% as.data.frame(),
      counts = LayerData(seu.q, layer = "counts"),
      nUMI = colSums(LayerData(seu.q, layer = "counts"))
    )
    
    print_with_time("Creating RCTD object...")
    RCTD = create.RCTD(query, reference, max_cores = 10)
    save_qs(RCTD, glue("{tmp.dir}/RCTD_{seg}.qs"))
  }
}


ParseRCTD = function(seu, res.dir) {
  coords = GetTissueCoordinates(seu) %>% mutate(cell_id = rownames(.), .before = 1)
  input.files = list.files(res.dir, pattern = ".done.qs", recursive = F, full.names = T)
  weights = pbapply::pblapply(input.files, function(input.file) {
    rctd = qread(input.file)
    weights = as.data.frame(rctd@results$weights)
    weights %>%
      mutate(
        cell_id = rownames(.),
        .before = 1
      ) %>%
      left_join(coords, by = "cell_id") %>%
      select(cell_id, segment, imagecol, imagerow, everything())
  }) %>% do.call(rbind, .)
  return(weights)
}

