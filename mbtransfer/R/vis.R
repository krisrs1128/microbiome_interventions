
#' @export
subject_order <- function(values_df, taxa, r = 0) {
  values_wide <- values_df |>
    filter(taxon %in% taxa) |>
    select(subject, taxon, value, time) |>
    mutate(time = round(time, r)) |>
    unite("txtime", taxon, time) |>
    pivot_wider(names_from = txtime) |>
    column_to_rownames("subject") |>
    as.matrix()
  
  values_wide[is.na(values_wide)] <- 0
  rownames(values_wide)[hclust(dist(values_wide))$order]
}

#' @export
interaction_hm <- function(values_df, taxa, condition = NULL, r = 0, ...) {
  p <- values_df |>
    filter(taxon %in% taxa) |>
    mutate(subject = factor(subject, subject_order(values_df, taxa, r))) |>
    ggplot() +
    geom_tile(aes(time, subject, fill = value, col = value), ...) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_distiller(direction = 1) +
    scale_color_distiller(direction = 1) +
    theme(
      strip.text.y = element_text(angle = 0),
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 1, fill = NA, color = "#545454"),
      panel.spacing = unit(0, "cm")
    )
  
  if (!is.null(condition)) {
    p <- p + 
      facet_grid(.data[[condition]] ~ reorder(taxon, -value), scales = "free", space = "free")
  }
  
  p
}


#' @importFrom ggh4x facet_nested
#' @export
interaction_barcode <- function(values_df, taxa, condition = NULL, r = 0, ...) {
  p <- values_df |>
    filter(taxon %in% taxa) |>
    mutate(subject = factor(subject, levels = subject_order(values_df, taxa, r))) |>
    ggplot() +
    geom_tile(aes(time, taxon, fill = value, col = value), ...) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_distiller(direction = 1) +
    scale_color_distiller(direction = 1) +
    theme(
      panel.spacing = unit(0, "line"),
      panel.border = element_rect(linewidth = 1, fill = NA, color = "#d3d3d3"),
      strip.text.y = element_text(angle = 0),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  if (!is.null(condition)) {
    p <- p + 
      facet_nested(reorder(.data[[condition]], -value) + subject ~ ., scales = "free", space = "free")
  }
  
  p
}

#' @importFrom ggplot2 theme_minimal theme
#' @export
my_theme <- function() {
  th <- theme_minimal() + 
    theme(
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#f7f7f7"),
      panel.border = element_rect(fill = NA, color = "#0c0c0c", size = 0.6),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.position = "bottom"
    )
}
