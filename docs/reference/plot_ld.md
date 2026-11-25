# Plot linkage disequilibrium as a heatmap

Draws the result of
[`calc_ld()`](https://jarekbryk.github.io/snpscaper/reference/calc_ld.md)
as a triangle: each square is a pair of SNPs, coloured from white (r² =
0, no LD) to red (r² = 1, complete LD). SNPs are arranged in order along
the chromosome, so blocks of high LD show up as red triangles along the
diagonal.

## Usage

``` r
plot_ld(ld, title = "Linkage disequilibrium")
```

## Arguments

- ld:

  The result of
  [`calc_ld()`](https://jarekbryk.github.io/snpscaper/reference/calc_ld.md).

- title:

  Plot title.

## Value

A ggplot object.

## Details

SNPs are drawn at equal spacing, in their order along the chromosome,
not at their exact positions; axis labels give the positions in
megabases.

## Examples

``` r
example_lct |>
  filter_snps(min_maf = 0.05) |>
  keep_pops(super_pop = "EUR") |>
  calc_ld() |>
  plot_ld(title = "LD around LCT in Europeans")
#> Dropped 18 SNPs that do not vary among these people.
```
