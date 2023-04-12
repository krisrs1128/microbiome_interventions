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
