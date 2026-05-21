## Log with timestamp
print_with_time <- function(msg) {
  print(glue::glue("[{format(Sys.time(), '%Y-%m-%d %H:%M:%S')}] {msg}"))
}


save_qs <- function(seu, file.name) {
  out.path <- dirname(file.name)
  if (!dir.exists(out.path)) {
    dir.create(out.path, recursive = T)
  }
  print_with_time(glue("Saving to {file.name}"))
  qs::qsave(seu, file.name)
}

## 逆时针旋转90度
rotate_90ccw <- function(img) {
  rotated <- aperm(img, c(2, 1, 3)) # 交换height和width
  rotated <- rotated[nrow(rotated):1, , ] # 上下翻转
  return(rotated)
}

rotate_90ccw_centered <- function(df, H, W) {
  center_row <- (H + 1) / 2
  center_col <- (W + 1) / 2
  x <- df$col - center_col
  y <- -(df$row - center_row)
  x_new <- -y
  y_new <- x
  df$col <- x_new + center_col
  df$row <- -y_new + center_row
  df$row <- df$row
  df$col <- df$col
  return(df)
}

RotateSeurat <- function(seu) {
  image <- seu@images$sample1@image
  image.rot <- rotate_90ccw(image)
  seu@images$sample1@image <- image.rot

  coords <- GetTissueCoordinates(seu)
  df <- data.frame(
    row = coords$imagerow,
    col = coords$imagecol
  )
  df.rot <- rotate_90ccw_centered(df, H = 1024, W = 1026)
  coords$imagerow <- df.rot$row
  coords$imagecol <- df.rot$col

  seu@images$sample1@coordinates <- coords
  seu@images$sample1@scale.factors$lowres <- 1

  return(seu)
}
