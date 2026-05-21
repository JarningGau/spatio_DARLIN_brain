QCVlnPlot <- function(seu, title, save.file, metric = c("UMI", "Gene")) {
  metadata <- seu@meta.data
  if (metric == "Gene") {
    p <- ggplot(metadata, aes(1, nFeature_Spatial)) +
      geom_violin(show.legend = F, fill = "lightblue") +
      stat_summary(fun = mean, geom = "point", color = "black", size = 2, show.legend = F) +
      stat_summary(fun = mean, geom = "text", aes(label = round(after_stat(y), 1)), vjust = -1, size = 5) +
      scale_y_log10() +
      labs(y = "Number of Genes per cell", x = "", title = title)
  } else if (metric == "UMI") {
    p <- ggplot(metadata, aes(1, nCount_Spatial)) +
      geom_violin(show.legend = F, fill = "lightblue") +
      stat_summary(fun = mean, geom = "point", color = "black", size = 2, show.legend = F) +
      stat_summary(fun = mean, geom = "text", aes(label = round(after_stat(y), 1)), vjust = -1, size = 5) +
      scale_y_log10() +
      labs(y = "Number of UMIs per cell", x = "", title = title)
  } else if (metric %in% colnames(seu@meta.data)) {
    p <- ggplot(metadata, aes(1, get(metric))) +
      geom_violin(show.legend = F, fill = "lightblue") +
      stat_summary(fun = mean, geom = "point", color = "black", size = 2, show.legend = F) +
      stat_summary(fun = mean, geom = "text", aes(label = round(after_stat(y), 1)), vjust = -1, size = 5) +
      scale_y_log10() +
      labs(y = metric, x = "", title = title)
  } else {
    warning(glue("Metric: {metric} is not supported."))
  }
  p <- p +
    theme(
      axis.title.y = element_text(size = 16),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = .5, size = 14)
    )
  base_dir <- dirname(save.file)
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
  }
  ggsave(save.file, plot = p, width = 2.5, height = 3.5, units = "in", dpi = 300)
}


QCSpatialPlot <- function(seu, title, save.file, metric = c("UMI", "Gene")) {
  if (metric == "Gene") {
    p <- SpatialFeaturePlot(seu, features = "nFeature_Spatial", max.cutoff = "q99", pt.size.factor = 4) +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      guides(fill = guide_colourbar(title = "Genes", keywidth = 1.5, keyheight = 15)) +
      ggtitle(title)
  } else if (metric == "UMI") {
    p <- SpatialFeaturePlot(seu, features = "nCount_Spatial", max.cutoff = "q99", pt.size.factor = 4) +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      guides(fill = guide_colourbar(title = "UMIs", keywidth = 1.5, keyheight = 15)) +
      ggtitle(title)
  } else if (metric %in% colnames(seu@meta.data)) {
    p <- SpatialFeaturePlot(seu, features = metric, max.cutoff = "q99", pt.size.factor = 4) +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      guides(fill = guide_colourbar(title = metric, keywidth = 1.5, keyheight = 15)) +
      ggtitle(title)
  } else {
    warning(glue("Metric: {metric} is not supported."))
  }
  p <- p +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 15),
      plot.title = element_text(size = 25, hjust = .5)
    )
  base_dir <- dirname(save.file)
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
  }
  ggsave(save.file, plot = p, width = 6, height = 5, units = "in", dpi = 300)
}
