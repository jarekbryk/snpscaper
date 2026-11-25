# Example data: a neutral region of chromosome 18

A `snpscape` object for a 300 kb stretch of chromosome 18
(chr18:1,000,000-1,300,000) that contains no gene known to be under
strong selection. Patterns of variation here are shaped mainly by
population history (migration, population size changes and drift), which
makes it a useful comparison for genes under selection such as
[example_lct](https://jarekbryk.github.io/snpscaper/reference/example_lct.md).
Fetched with `fetch_region("chr18:1000000-1300000")`.

## Usage

``` r
example_neutral
```

## Format

A `snpscape` object (see
[snpscape](https://jarekbryk.github.io/snpscaper/reference/snpscape.md)).

## Source

1000 Genomes Project high-coverage phased panel. Built by
`data-raw/examples.R`.

## Examples

``` r
example_neutral
#> <snpscape> 2,504 people x 6,125 SNPs (GRCh38)
#> Regions:     chr18:1,000,000-1,300,000 
#> Populations: 26 populations in 5 super-populations
#> Steps:       fetched 
```
