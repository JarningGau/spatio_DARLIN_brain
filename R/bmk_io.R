library(Seurat)
library(dplyr)


CreateS1000Object = function(
  matrix_path,
  png_path,
  spot_radius = NULL,
  min.cells = 5,
  min.features = 100,
  cell.column = 1, 
  gene.column = 2
  ){
  expr = Seurat::Read10X(matrix_path, cell.column = cell.column, gene.column = gene.column)
  object = Seurat::CreateSeuratObject(counts = expr,
                               assay = 'Spatial',
                               min.cells=min.cells,
                               min.features=min.features)
  #Image zoom rate
  cal_zoom_rate = function(width, height){
    std_width = 1000
    std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
    if (std_width / std_height > width / height){
      scale = width / std_width
    }
    else{
      scale = height / std_height
    }
    return(scale)
  }
  #read png
  png = png::readPNG(png_path)
  zoom_scale =  cal_zoom_rate(dim(png)[2], dim(png)[1])
  #read barcode pos file
  ReadBarcodePos = function(barcode_pos_path){
    barcode_pos = read.table(gzfile(barcode_pos_path),header = F) %>%
      dplyr::rename(Barcode = V1 , pos_w = V2, pos_h = V3)
    return(barcode_pos)
  }
  #get barcode pos file path
  barcode_pos_path = paste0(matrix_path,'/barcodes_pos.tsv.gz')
  barcode_pos = ReadBarcodePos(barcode_pos_path = barcode_pos_path)
  barcode_pos = barcode_pos %>% dplyr::filter(., Barcode %in% rownames(object@meta.data))
  #make spatial coord file for seurat S4 class
  coord = data.frame(tissue = 1,
                      row = barcode_pos$pos_h,
                      col = barcode_pos$pos_w,
                      imagerow = barcode_pos$pos_h,
                      imagecol = barcode_pos$pos_w)
  rownames(coord) = barcode_pos$Barcode
  #spot radius
  spot_radius_lib = c(0.00063, 0.00179, 0.0027, 0.0039, 0.004, 0.0045, 0.005, NA, NA, NA, NA, NA, 0.0120)
  if(is.null(spot_radius)){
    spot_radius = spot_radius_lib[as.numeric(gsub('L', '', strsplit(tail(strsplit(matrix_path, '/')[[1]],1), '_')[[1]][1]))]
  }else{
    spot_radius = spot_radius
  }
  #object
  sample1 =  new(Class = "VisiumV1",
                  image = png,
                  scale.factors = Seurat::scalefactors(zoom_scale, 100, zoom_scale, zoom_scale),
                  coordinates = coord,
                  spot.radius = spot_radius,
                  assay = 'Spatial',
                  key = "sample1_")
  object@images = list(sample1 = sample1)

  return(object)
}


## zoom the coordinates to fit the image (width: 1024, height: 1026)
## only use for cell segmentation data
ZoomImageCoords = function(seu) {
  coords = seu@images[[1]]@coordinates
  scale.factor = 1000/19954
  coords$imagerow = coords$row * scale.factor
  coords$imagecol = coords$col * scale.factor
  seu@images[[1]]@coordinates = coords
  return(seu)
}


Seurat2H5adSpatial = function(seu, output.file) {
  library(reticulate)
  sc = import("scanpy")
  counts = LayerData(seu, layer = "counts")
  counts@x = as.integer(counts@x)
  obs = seu@meta.data
  var = data.frame(row.names = rownames(counts), gene_symbols = rownames(counts))
  adata = sc$AnnData(X = t(counts), obs = obs, var = var)
  barcode_pos = seu@images$sample1@coordinates
  
  # scale.factor = 1026/20000
  # barcode_pos$imagerow = barcode_pos$row * scale.factor
  # barcode_pos$imagecol = barcode_pos$col * scale.factor
  barcode_pos = barcode_pos[, c("tissue", "imagerow", "imagecol")]
  colnames(barcode_pos) = c("in_tissue", "array_col","array_row")
  
  obsm = barcode_pos[, c("array_row", "array_col")]
  obsm = as.matrix(obsm)
  adata$obsm['spatial'] = obsm
  
  # umap = Embeddings(seu, reduction = "umap")
  # adata$obsm["umap"] = umap
  
  he_img = seu@images$sample1@image
  adata$uns$spatial = list()
  adata$uns$spatial[[sample.ID]] = list(
    "images" = list(
      "hires" = he_img
    ),
    "use_quality" = "hires", 
    "scalefactors" = list(
      spot_diameter_fullres = 1,
      tissue_hires_scalef = 1,
      fiducial_diameter_fullres = 1,
      tissue_lowres_scalef = 1
    )
  )
  adata$write_h5ad(output.file)
}
