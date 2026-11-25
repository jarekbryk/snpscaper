# Plot Fst between pairs of populations

A heatmap of the result of
[`calc_fst()`](https://jarekbryk.github.io/snpscaper/reference/calc_fst.md):
darker squares mean populations are more similar, brighter squares mean
they are more differentiated. Values are printed in the squares when
there are few groups.

## Usage

``` r
plot_fst(fst, title = "Pairwise Fst")
```

## Arguments

- fst:

  The result of
  [`calc_fst()`](https://jarekbryk.github.io/snpscaper/reference/calc_fst.md).

- title:

  Plot title.

## Value

A ggplot object.

## Examples

``` r
example_lct |>
  filter_snps() |>
  keep_pops(pop = c("YRI", "LWK", "GBR", "FIN", "CHB", "JPT")) |>
  calc_fst() |>
  plot_fst()
```
