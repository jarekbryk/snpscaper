# Map the allele frequency of a SNP

Shows on a world map how common the alternative allele of one SNP is in
each population. Each population is a circle, coloured from dark (allele
rare) to yellow (allele common).

## Usage

``` r
map_freqs(freqs, snp)
```

## Arguments

- freqs:

  The result of
  [`calc_allele_freqs()`](https://jarekbryk.github.io/snpscaper/reference/calc_allele_freqs.md)
  with `by = "pop"`.

- snp:

  The SNP to plot: an rs number (e.g. `"rs4988235"`, for SNPs in
  [teaching_snps](https://jarekbryk.github.io/snpscaper/reference/teaching_snps.md))
  or an ID from `x$variants$id` (e.g. `"chr2:135851076:G:A"`).

## Value

A ggplot object.

## Details

Several populations were sampled away from the region their ancestry
comes from (for example Gujarati Indians in Houston). These are drawn at
their region of ancestral origin; see
[pop_locations](https://jarekbryk.github.io/snpscaper/reference/pop_locations.md).
Populations of mixed recent ancestry, such as African Americans (ASW),
are drawn where they were sampled.

## Examples

``` r
example_lct |>
  calc_allele_freqs() |>
  map_freqs(snp = "rs4988235")
```
