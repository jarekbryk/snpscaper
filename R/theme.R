#' A clean ggplot2 theme for snpscaper plots
#'
#' A minimal theme with light horizontal grid lines, left-aligned titles and
#' no legend title. All snpscaper plots use it, and you can add it to your
#' own plots too, including factoextra plots such as
#' `factoextra::fviz_pca_ind()`.
#'
#' @param base_size Base font size, in points. Other sizes default to this.
#' @param base_family Base font family.
#' @param axis_text_size,axis_title_size Font sizes for axis labels and titles.
#' @param axis_title_x_angle,axis_title_y_angle Rotation of the axis titles.
#' @param strip_text_size Font size of facet labels.
#' @param legend_text_size Font size of legend labels.
#' @param plot_title_size,plot_subtitle_size,plot_caption_size Font sizes of
#'   the plot title, subtitle and caption.
#'
#' @return A ggplot2 theme object, to be added to a plot with `+`.
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   labs(title = "A plot using theme_jb()") +
#'   theme_jb()
theme_jb <- function(
  base_size = 12,
  base_family = "",
  axis_text_size = base_size,
  axis_title_size = base_size,
  axis_title_x_angle = 0,
  axis_title_y_angle = 90,
  strip_text_size = base_size,
  legend_text_size = base_size,
  plot_title_size = base_size + 4,
  plot_subtitle_size = base_size + 1,
  plot_caption_size = base_size - 2
) {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey85",
        linewidth = 0.3
      ),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        size = plot_title_size,
        hjust = 0,
        margin = ggplot2::margin(b = 8)
      ),
      plot.subtitle = ggplot2::element_text(
        size = plot_subtitle_size,
        hjust = 0,
        margin = ggplot2::margin(b = 14)
      ),
      plot.caption = ggplot2::element_text(
        size = plot_caption_size,
        hjust = 1,
        colour = "grey50"
      ),
      axis.text = ggplot2::element_text(size = axis_text_size),
      axis.title.x = ggplot2::element_text(
        size = axis_title_size,
        angle = axis_title_x_angle
      ),
      axis.title.y = ggplot2::element_text(
        size = axis_title_size,
        angle = axis_title_y_angle
      ),
      strip.text = ggplot2::element_text(size = strip_text_size),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend_text_size),
      plot.title.position = "plot"
    )
}
