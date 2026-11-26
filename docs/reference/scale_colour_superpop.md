# Colour and fill scales for populations

ggplot2 scales that apply the fixed snpscaper colours
([superpop_colours](https://jarekbryk.github.io/snpscaper/reference/colours.md),
[pop_colours](https://jarekbryk.github.io/snpscaper/reference/colours.md)),
so the same population always has the same colour.

## Usage

``` r
scale_colour_superpop(...)

scale_fill_superpop(...)

scale_colour_pop(...)

scale_fill_pop(...)
```

## Arguments

- ...:

  Passed on to
  [`ggplot2::scale_colour_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
  or
  [`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html),
  for example `name` or `breaks`.

## Value

A ggplot2 scale.

## Examples

``` r
library(ggplot2)
#> Warning: package ‘ggplot2’ was built under R version 4.6.1
df <- data.frame(super_pop = names(superpop_colours), n = c(661, 347, 504, 503, 489))
ggplot(df, aes(super_pop, n, fill = super_pop)) +
  geom_col() +
  scale_fill_superpop() +
  theme_jb()
```
