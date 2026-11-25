# Measure genetic diversity within populations (heterozygosity)

Calculates two measures of genetic diversity for each population (or
super-population), averaged over all SNPs in the data.

## Usage

``` r
calc_heterozygosity(x, by = c("pop", "super_pop"))
```

## Arguments

- x:

  A `snpscape` object.

- by:

  Calculate per `"pop"` or per `"super_pop"`.

## Value

A tibble with one row per group: the group (`pop` and `super_pop`, or
only `super_pop`), `people`, `snps`, `observed` and `expected`.

## Details

A person is **heterozygous** at a SNP if they carry two different
alleles (genotype count 1). Populations with more genetic diversity have
more heterozygous people.

- **Observed heterozygosity** is the proportion of genotypes that are
  heterozygous.

- **Expected heterozygosity** (also called gene diversity) is the
  proportion we would expect if people paired up alleles at random
  within the population: `2p(1 - p)` for allele frequency `p`, with a
  small correction for sample size.

The two are usually close. Observed heterozygosity clearly below
expected can point to inbreeding or to hidden structure within a
population.

Across the human genome, African populations have the highest diversity,
and diversity decreases with distance from Africa. This is thought to
reflect a series of population bottlenecks as small groups of people
migrated out of Africa and across the world (the "serial founder
effect").

The values depend on which SNPs are included. They are useful for
comparing populations on the same set of SNPs, not as absolute numbers:
for example, removing rare SNPs with
[`filter_snps()`](https://jarekbryk.github.io/snpscaper/reference/filter_snps.md)
raises the values.

## Examples

``` r
example_lct |>
  filter_snps() |>
  calc_heterozygosity()
#> # A tibble: 26 x 6
#>    pop   super_pop people  snps observed expected
#>    <chr> <chr>      <int> <int>    <dbl>    <dbl>
#>  1 ACB   AFR           96   417   0.203    0.200 
#>  2 ASW   AFR           61   417   0.203    0.197 
#>  3 BEB   SAS           86   417   0.154    0.163 
#>  4 CDX   EAS           93   417   0.195    0.178 
#>  5 CEU   EUR           99   417   0.0786   0.0841
#>  6 CHB   EAS          103   417   0.180    0.187 
#>  7 CHS   EAS          105   417   0.201    0.187 
#>  8 CLM   AMR           94   417   0.162    0.158 
#>  9 ESN   AFR           99   417   0.199    0.200 
#> 10 FIN   EUR           99   417   0.134    0.123 
#> # i 16 more rows
```
