# Calculate allele frequencies in each population

For every SNP, calculates how common the alternative allele is in each
population (or super-population).

## Usage

``` r
calc_allele_freqs(x, by = c("pop", "super_pop"))
```

## Arguments

- x:

  A `snpscape` object.

- by:

  Calculate frequencies per `"pop"` or per `"super_pop"`.

## Value

A tibble with one row per SNP and group: the SNP details (`id`, `rsid`,
`chrom`, `pos`, `ref`, `alt`), the group (`pop` and `super_pop`, or only
`super_pop`), `people` (number of people with a known genotype) and
`alt_freq` (alternative allele frequency, between 0 and 1).

## Details

Each person carries two copies of every SNP, one from each parent. The
alternative allele frequency in a population is the number of copies of
the alternative allele carried by its people, divided by the total
number of copies (twice the number of people with a known genotype).

Which allele is called "alternative" is a convention: the reference
allele is simply the one in the human reference genome, which was
assembled mostly from a few individuals. It is not necessarily the
ancestral or the most common allele. Check the `ref` and `alt` columns.

## Examples

``` r
freqs <- calc_allele_freqs(example_lct)
freqs[freqs$rsid %in% "rs4988235", ]
#> # A tibble: 26 x 10
#>    id             rsid  chrom    pos ref   alt   pop   super_pop people alt_freq
#>    <chr>          <chr> <chr>  <int> <chr> <chr> <chr> <chr>      <int>    <dbl>
#>  1 chr2:13585107~ rs49~ chr2  1.36e8 G     A     ACB   AFR           96   0.0677
#>  2 chr2:13585107~ rs49~ chr2  1.36e8 G     A     ASW   AFR           61   0.172 
#>  3 chr2:13585107~ rs49~ chr2  1.36e8 G     A     BEB   SAS           86   0.0523
#>  4 chr2:13585107~ rs49~ chr2  1.36e8 G     A     CDX   EAS           93   0     
#>  5 chr2:13585107~ rs49~ chr2  1.36e8 G     A     CEU   EUR           99   0.737 
#>  6 chr2:13585107~ rs49~ chr2  1.36e8 G     A     CHB   EAS          103   0     
#>  7 chr2:13585107~ rs49~ chr2  1.36e8 G     A     CHS   EAS          105   0     
#>  8 chr2:13585107~ rs49~ chr2  1.36e8 G     A     CLM   AMR           94   0.309 
#>  9 chr2:13585107~ rs49~ chr2  1.36e8 G     A     ESN   AFR           99   0     
#> 10 chr2:13585107~ rs49~ chr2  1.36e8 G     A     FIN   EUR           99   0.591 
#> # i 16 more rows
```
