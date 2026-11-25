# Plot heterozygosity in each population

Shows expected heterozygosity (bars) and observed heterozygosity
(points) for each population, coloured by super-population. See
[`calc_heterozygosity()`](https://jarekbryk.github.io/snpscaper/reference/calc_heterozygosity.md)
for what these measure.

## Usage

``` r
plot_heterozygosity(het, title = "Genetic diversity")
```

## Arguments

- het:

  The result of
  [`calc_heterozygosity()`](https://jarekbryk.github.io/snpscaper/reference/calc_heterozygosity.md).

- title:

  Plot title.

## Value

A ggplot object.

## Examples

``` r
example_lct |>
  filter_snps() |>
  calc_heterozygosity() |>
  plot_heterozygosity()
```
