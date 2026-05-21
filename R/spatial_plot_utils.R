GetZoomInCells = function(seu, zoom.window) {
  coords = GetTissueCoordinates(seu)
  max.row = max(coords$imagerow)
  max.col = max(coords$imagecol)
  min.row = min(coords$imagerow)
  min.col = min(coords$imagecol)
  coords = subset(coords, imagecol > zoom.window[1] & imagecol < zoom.window[3])
  coords = subset(coords, imagerow > max.row-zoom.window[4]+min.row & imagerow < max.row-zoom.window[2]+min.row)
  return(rownames(coords))
}

GetZoomWindow = function(coords, cells.use, half.win=100) {
  xmin = min(coords$imagecol)
  xmax = max(coords$imagecol)
  ymin = min(coords$imagerow)
  ymax = max(coords$imagerow)
  coords.use = coords[rownames(coords) %in% cells.use, ]
  ymid.zoom = ymax - median(coords.use[, "imagerow"])
  xmid.zoom = median(coords.use[, "imagecol"])
  xmin.zoom = max(xmid.zoom - half.win, xmin)
  xmax.zoom = min(xmid.zoom + half.win, xmax)
  ymin.zoom = max(ymid.zoom - half.win, ymin)
  ymax.zoom = min(ymid.zoom + half.win, ymax)
  zoom.window = c(
    xmin.zoom, ymin.zoom, 
    xmax.zoom, ymax.zoom
  )
  return(unname(zoom.window))
}

AddScaleBar = function(seu, scale.bar.pixels, position = "bottom-left", rel.pos = 0.01, color = "black", linewidth = 1) {
  coords = GetTissueCoordinates(seu)
  max.row = max(coords$imagerow)
  max.col = max(coords$imagecol)
  min.row = min(coords$imagerow)
  min.col = min(coords$imagecol)
  range.col = diff(range(coords$imagecol))
  range.row = diff(range(coords$imagerow))
  if (position == "bottom-left") {
    xmin = min.col
    xmax = xmin + scale.bar.pixels
    y = min.row + range.row*rel.pos
  } else if (position == "bottom-right") {
    xmax = max.col
    xmin = xmax - scale.bar.pixels
    xmin = xmin - range.col*rel.pos
    xmax = xmax - range.col*rel.pos
    y = min.row + range.row*rel.pos
  } else if (position == "top-left") {
    xmin = min.col
    xmax = xmin + scale.bar.pixels
    y = max.row - range.row*rel.pos
  } else if (position == "top-right") {
    xmax = max.col
    xmin = xmax - scale.bar.pixels
    y = max.row - range.row*rel.pos
  }
  annotate("segment", x = xmin, y = y, xend = xmax, yend = y, color = color, linewidth = linewidth)
}


AddZoomWindow = function(zoom.window, linetype="solid", color = "black", linewidth=1) {
  annotate("rect", xmin=zoom.window[1], ymin=zoom.window[2], xmax=zoom.window[3], ymax=zoom.window[4],
           color = color, fill="white", alpha=0, linetype=linetype, linewidth=linewidth)
}



FlipYAxis = function(coords, cells.use=NULL) {
  colnames(coords) = c("y", "x")
  O = (min(coords$y) + max(coords$y)) / 2
  coords$y = 2 * O - coords$y  # Mirror flip y-coordinates
  if (!is.null(cells.use)) {
    data.plot = coords[cells.use, ]
  } else {
    data.plot = coords
  }
  return(data.plot)
}


## Set up scale bar function
make_scale_bar_r <- function(x_vals, y_vals, microns_per_pixel = 0.12028, scale_length_um = 500) {
  
  # Adds a scale bar to a ggplot.
  
  #Parameters:
  #x_vals: vector of x coordinates in pixels
  #y_vals: vector of y coordinates in pixels
  #microns_per_pixel: conversion factor from pixels to microns. Default is conversion factor for commercial CosMx instrument
  
  # Example usage:
  # scale_bar = make_scale-bar_r(x_vals = cell_meta$CenterX_global_px, y_vals = cell_meta$CenterY_global_px)
  # ggplot() + scale_bar$bg + scale_bar$rect + scale_bar$label
  
  
  # Calculate x-axis range
  x_range <- range(x_vals, na.rm = TRUE)
  x_length <- diff(x_range)
  x_length_um <- x_length * microns_per_pixel
  
  # Target scale length ~1/4 of the x-axis
  target <- x_length_um / 4
  
  # Compute order of magnitude
  order <- 10^floor(log10(target))
  mantissa <- target / order
  
  # Round mantissa to nearest 1, 2, or 5
  nice_mantissa <- if (mantissa < 1.5) {
    1
  } else if (mantissa < 3.5) {
    2
  } else if (mantissa < 7.5) {
    5
  } else {
    10
  }
  
  # Final scale length in pixels
  if (is.null(scale_length_um)) {
    scale_length_um <- nice_mantissa * order
  }
  
  scale_length_px <- scale_length_um / microns_per_pixel
  
  # Format label
  scale_label <- if (scale_length_um >= 1000) {
    paste0(scale_length_um / 1000, " mm")
  } else {
    paste0(scale_length_um, " µm")
  }
  
  # Set coordinates for the scale bar 
  x_start <- x_range[2] - scale_length_px * 1.1
  x_end <- x_range[2] - scale_length_px * 0.1
  y_pos <- min(y_vals, na.rm = TRUE) + scale_length_px * 0.1
  
  # Generate scale bar background, scale bar and annotation to return
  list(
    bg = annotation_custom(
      grob = rectGrob(gp = gpar(fill = "white", alpha = 0.8, col = NA)),
      xmin = x_start- scale_length_px*0.05, xmax = x_end+ scale_length_px*0.05,
      ymin = y_pos- scale_length_px*0.05, ymax = y_pos + scale_length_px*0.3
    ),
    rect = annotation_custom(
      grob = rectGrob(gp = gpar(fill = "black")),
      xmin = x_start, xmax = x_end,
      ymin = y_pos, ymax = y_pos + scale_length_px * 0.05
    ),
    label = annotation_custom(
      grob = textGrob(scale_label, gp = gpar(col = "black"), just = "center", vjust = 0),
      xmin = (x_start + x_end)/2, xmax = (x_start + x_end)/2,
      ymin = y_pos + scale_length_px * 0.1, ymax = y_pos + scale_length_px * 0.1
    )
  )
}
