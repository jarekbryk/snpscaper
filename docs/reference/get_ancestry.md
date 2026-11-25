# Get ancestry proportions from an admixture analysis

Extracts the estimated ancestry proportions for one or more values of
`K` from the result of
[`run_admixture()`](https://jarekbryk.github.io/snpscaper/reference/run_admixture.md),
as a tidy table.

## Usage

``` r
get_ancestry(adm, k = NULL, wide = FALSE)
```

## Arguments

- adm:

  The result of
  [`run_admixture()`](https://jarekbryk.github.io/snpscaper/reference/run_admixture.md).

- k:

  The `K` values to extract, or `NULL` for all.

- wide:

  If `TRUE`, return one row per person with one column per ancestral
  group (only for a single `K`).

## Value

A tibble. In the default long format: `sample_id`, `pop`, `super_pop`,
`k`, `group` and `proportion`.

## Examples

``` r
# \donttest{
lct <- filter_snps(example_lct, min_maf = 0.05, prune_ld = TRUE)
adm <- run_admixture(lct, k = 2:3, reps = 2)
#> Running admixture for K = 2 and 3 (2 runs each, 12 SNPs)...
get_ancestry(adm, k = 3, wide = TRUE)
#> # A tibble: 2,504 x 7
#>    sample_id pop   super_pop     k   group_1   group_2   group_3
#>    <chr>     <chr> <chr>     <int>     <dbl>     <dbl>     <dbl>
#>  1 HG01879   ACB   AFR           3 0.0001000 1.000     0.0001000
#>  2 HG01880   ACB   AFR           3 0.0001000 0.0001000 1.000    
#>  3 HG01882   ACB   AFR           3 0.0001000 1.000     0.0001000
#>  4 HG01883   ACB   AFR           3 0.0001000 1.000     0.0001000
#>  5 HG01885   ACB   AFR           3 0.503     0.0001000 0.497    
#>  6 HG01886   ACB   AFR           3 1.000     0.0001000 0.0001000
#>  7 HG01889   ACB   AFR           3 0.0001000 1.000     0.0001000
#>  8 HG01890   ACB   AFR           3 0.0001000 1.000     0.0001000
#>  9 HG01894   ACB   AFR           3 0.504     0.0001000 0.496    
#> 10 HG01896   ACB   AFR           3 0.132     0.0001000 0.868    
#> # i 2,494 more rows
# }
```
