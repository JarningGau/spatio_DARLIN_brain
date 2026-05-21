theme_general <- function(base_size = 15) {
  theme_bw(base_size = base_size) + 
    theme(
      legend.position = "right",
      legend.text = element_text(size = 12),
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 1.5, fill=NA),
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color = "black")
    )
}


theme_polar <- function() {
  theme_void() + 
    theme(
      plot.title = element_text(hjust = .5, face = "bold"),
      strip.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )
}


theme_spatial <- function(spacing = 0) {
  list(
    scale_x_continuous(expand = c(spacing, spacing)),
    scale_y_continuous(expand = c(spacing, spacing)),
    theme_bw(),
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.border = element_rect(linewidth = 2, fill=NA)
    )
  )
}


theme_graphpad_category <- function() {
  list(
    scale_y_continuous(expand = c(0, 0, 0, 0.1)),
    theme_classic(),
    theme(
      plot.title = element_text(size = 15, hjust = 0.5),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(size = 14),
      axis.ticks.x = element_blank()
    )
  )
}

theme_graphpad_continous <- function() {
  list(
    scale_y_continuous(expand = c(0, 0, 0, 0.1)),
    theme_classic(),
    theme(
      plot.title = element_text(size = 15, hjust = 0.5),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(size = 14),
      axis.ticks.length = unit(.2, "cm")
    )
  )
}
