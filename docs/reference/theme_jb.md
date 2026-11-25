# A clean ggplot2 theme for snpscaper plots

A minimal theme with light horizontal grid lines, left-aligned titles
and no legend title. All snpscaper plots use it, and you can add it to
your own plots too, including factoextra plots such as
[`factoextra::fviz_pca_ind()`](https://rpkgs.datanovia.com/factoextra/reference/fviz_pca.html).

## Usage

``` r
theme_jb(
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
)
```

## Arguments

- base_size:

  Base font size, in points. Other sizes default to this.

- base_family:

  Base font family.

- axis_text_size, axis_title_size:

  Font sizes for axis labels and titles.

- axis_title_x_angle, axis_title_y_angle:

  Rotation of the axis titles.

- strip_text_size:

  Font size of facet labels.

- legend_text_size:

  Font size of legend labels.

- plot_title_size, plot_subtitle_size, plot_caption_size:

  Font sizes of the plot title, subtitle and caption.

## Value

A ggplot2 theme object, to be added to a plot with `+`.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  labs(title = "A plot using theme_jb()") +
  theme_jb()
```
