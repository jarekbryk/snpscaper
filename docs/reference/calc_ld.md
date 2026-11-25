# Calculate linkage disequilibrium between SNPs

Calculates r², a measure of linkage disequilibrium (LD), between every
pair of SNPs in the data.

## Usage

``` r
calc_ld(x, max_snps = 300)
```

## Arguments

- x:

  A `snpscape` object.

- max_snps:

  Maximum number of SNPs to use.

## Value

A tibble with one row per pair of SNPs: `snp1`, `snp2` (SNP IDs),
`pos1`, `pos2` (positions), `index1`, `index2` (the SNPs' order along
the chromosome) and `r2`.

## Details

### What LD is

Chromosomes are passed from parents to children in long pieces, which
are broken up only occasionally by recombination. SNPs that are close
together therefore tend to be inherited together, so knowing a person's
genotype at one SNP tells you something about their genotype at a nearby
SNP. This non-random association is called linkage disequilibrium.

r² is the squared correlation between the genotypes at two SNPs across
people: 0 means no association, 1 means one SNP perfectly predicts the
other. LD typically decays with distance, but the pattern is patchy:
"blocks" of high LD are separated by recombination hotspots.

LD differs between populations. African populations, which have the
longest history and largest long-term population sizes, generally show
shorter stretches of LD than populations that went through bottlenecks
during migrations out of Africa. A recent selective sweep (as around LCT
in Europeans) creates an unusually long stretch of high LD, because the
favoured allele rose in frequency so fast that recombination had little
time to break up the surrounding chromosome. Compare populations by
using
[`keep_pops()`](https://jarekbryk.github.io/snpscaper/reference/keep_pops.md)
before `calc_ld()`.

### Practicalities

Rare SNPs give unreliable r² values, so filter them out first, e.g. with
`filter_snps(min_maf = 0.05)` (and without LD pruning, which would
remove exactly the pattern you want to see). The number of pairs grows
quickly with the number of SNPs, so if there are more than `max_snps`
SNPs, an evenly spaced subset is used.

## Examples

``` r
example_lct |>
  filter_snps(min_maf = 0.05) |>
  keep_pops(super_pop = "EUR") |>
  calc_ld()
#> Dropped 18 SNPs that do not vary among these people.
#> # A tibble: 14,878 x 7
#>    snp1               snp2                    pos1      pos2 index1 index2    r2
#>    <chr>              <chr>                  <int>     <int>  <int>  <int> <dbl>
#>  1 chr2:135763944:G:A chr2:135764025:G:A 135763944 135764025      1      2 0.120
#>  2 chr2:135763944:G:A chr2:135765105:C:T 135763944 135765105      1      3 0.996
#>  3 chr2:135764025:G:A chr2:135765105:C:T 135764025 135765105      2      3 0.120
#>  4 chr2:135763944:G:A chr2:135765140:A:G 135763944 135765140      1      4 0.859
#>  5 chr2:135764025:G:A chr2:135765140:A:G 135764025 135765140      2      4 0.124
#>  6 chr2:135765105:C:T chr2:135765140:A:G 135765105 135765140      3      4 0.862
#>  7 chr2:135763944:G:A chr2:135765213:T:G 135763944 135765213      1      5 0.299
#>  8 chr2:135764025:G:A chr2:135765213:T:G 135764025 135765213      2      5 0.442
#>  9 chr2:135765105:C:T chr2:135765213:T:G 135765105 135765213      3      5 0.300
#> 10 chr2:135765140:A:G chr2:135765213:T:G 135765140 135765213      4      5 0.343
#> # i 14,868 more rows
```
