# Measure genetic differentiation between populations (Fst)

Calculates Fst for every pair of populations (or super-populations).

## Usage

``` r
calc_fst(x, by = c("pop", "super_pop"))
```

## Arguments

- x:

  A `snpscape` object.

- by:

  Compare populations (`"pop"`) or super-populations (`"super_pop"`).

## Value

A tibble with one row per pair of groups: `group1`, `group2` and `fst`.

## Details

### What Fst measures

Fst (the "fixation index") measures how much of the genetic variation is
due to differences *between* populations rather than variation *within*
them. It is 0 when two populations have identical allele frequencies,
and 1 when they share no alleles at all (each fixed for a different
one). For human populations, genome-wide values are small: roughly 0.01
between neighbouring populations and 0.1 to 0.15 between continents.
Most human genetic variation is found within populations.

Individual genes can deviate strongly from the genome-wide average. A
gene under strong selection in some populations (e.g. SLC24A5 or ACKR1)
can show much higher Fst than a typical region.

### How it is calculated

We use Hudson's estimator, as recommended by Bhatia et al. (2013,
*Genome Research* 23: 1514-1521). For each SNP it compares the allele
frequency difference between the two populations with the diversity
expected if they were one population, correcting for sample size. The
values for all SNPs are then combined as a ratio of sums (the "ratio of
averages"), which is more reliable than averaging Fst values of single
SNPs. Small negative values can occur when populations are nearly
identical; they mean "no detectable difference".

## Examples

``` r
example_lct |>
  filter_snps() |>
  calc_fst(by = "super_pop")
#> # A tibble: 10 × 3
#>    group1 group2    fst
#>    <chr>  <chr>   <dbl>
#>  1 AFR    AMR    0.0880
#>  2 AFR    EAS    0.0797
#>  3 AFR    EUR    0.173 
#>  4 AFR    SAS    0.0692
#>  5 AMR    EAS    0.0243
#>  6 AMR    EUR    0.0372
#>  7 AMR    SAS    0.0116
#>  8 EAS    EUR    0.0686
#>  9 EAS    SAS    0.0388
#> 10 EUR    SAS    0.0824
```
