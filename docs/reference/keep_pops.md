# Keep only some populations or people

Narrows the data down to the populations you want to compare. You can
choose by population (e.g. `"GBR"`, `"YRI"`), by super-population (e.g.
`"EUR"`, `"AFR"`), by sex, or any combination; people must match all the
choices you make. See
[samples_1kg](https://jarekbryk.github.io/snpscaper/reference/samples_1kg.md)
for the population codes.

## Usage

``` r
keep_pops(x, pop = NULL, super_pop = NULL, sex = NULL)
```

## Arguments

- x:

  A `snpscape` object.

- pop:

  Population codes to keep, or `NULL` to keep all.

- super_pop:

  Super-population codes to keep, or `NULL` to keep all.

- sex:

  `"female"` and/or `"male"`, or `NULL` to keep both.

## Value

The `snpscape` object with only the chosen people.

## Examples

``` r
example_lct |>
  keep_pops(pop = c("YRI", "GBR", "CHB"))
#> <snpscape> 302 people x 1,922 SNPs (GRCh38)
#> Regions:     LCT (chr2:135,762,850-135,862,184) 
#> Populations: 3 populations in 3 super-populations
#> Steps:       fetched -> kept YRI, GBR, CHB 

example_lct |>
  keep_pops(super_pop = c("AFR", "EUR"))
#> <snpscape> 1,164 people x 1,922 SNPs (GRCh38)
#> Regions:     LCT (chr2:135,762,850-135,862,184) 
#> Populations: 12 populations in 2 super-populations
#> Steps:       fetched -> kept AFR, EUR 
```
