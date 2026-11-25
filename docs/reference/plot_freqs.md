# Plot the allele frequency of a SNP in each population

A bar chart of how common the alternative allele of one SNP is in each
population, coloured by super-population. Useful for seeing at a glance
whether a variant is shared by everyone, or concentrated in some parts
of the world.

## Usage

``` r
plot_freqs(freqs, snp)
```

## Arguments

- freqs:

  The result of
  [`calc_allele_freqs()`](https://jarekbryk.github.io/snpscaper/reference/calc_allele_freqs.md).

- snp:

  The SNP to plot: an rs number (e.g. `"rs4988235"`, for SNPs in
  [teaching_snps](https://jarekbryk.github.io/snpscaper/reference/teaching_snps.md))
  or an ID from `x$variants$id` (e.g. `"chr2:135851076:G:A"`).

## Value

A ggplot object.

## Examples

``` r
example_lct |>
  calc_allele_freqs() |>
  plot_freqs(snp = "rs4988235")
```
