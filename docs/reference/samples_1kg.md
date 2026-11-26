# People in the 1000 Genomes Project

The 2,504 unrelated people from 26 populations whose genomes were
sequenced by the 1000 Genomes Project. These are the people included in
every `snpscape` object.

## Usage

``` r
samples_1kg
```

## Format

A tibble with 2,504 rows and 6 columns:

- sample_id:

  Sample identifier, e.g. `"HG00096"`.

- pop:

  Three-letter population code, e.g. `"GBR"`.

- super_pop:

  Super-population code, e.g. `"EUR"`.

- sex:

  `"female"` or `"male"`.

- pop_name:

  Full population name.

- super_pop_name:

  Full super-population name.

## Source

1000 Genomes Project, International Genome Sample Resource:
<https://www.internationalgenome.org>. Built by
`data-raw/samples_1kg.R`.

## Details

The 26 populations are grouped into five continental
**super-populations**: African (AFR), Admixed American (AMR), East Asian
(EAS), European (EUR) and South Asian (SAS). "Admixed American"
populations have recent ancestry from several continents (Indigenous
American, European and African). Some populations were sampled outside
their region of ancestral origin, for example Gujarati Indians in
Houston (GIH) and Sri Lankan Tamils in the UK (STU).

## Examples

``` r
samples_1kg
#> # A tibble: 2,504 × 6
#>    sample_id pop   super_pop sex    pop_name                      super_pop_name
#>    <chr>     <chr> <chr>     <chr>  <chr>                         <chr>         
#>  1 HG01879   ACB   AFR       male   African Caribbean in Barbados African       
#>  2 HG01880   ACB   AFR       female African Caribbean in Barbados African       
#>  3 HG01882   ACB   AFR       male   African Caribbean in Barbados African       
#>  4 HG01883   ACB   AFR       female African Caribbean in Barbados African       
#>  5 HG01885   ACB   AFR       male   African Caribbean in Barbados African       
#>  6 HG01886   ACB   AFR       female African Caribbean in Barbados African       
#>  7 HG01889   ACB   AFR       female African Caribbean in Barbados African       
#>  8 HG01890   ACB   AFR       male   African Caribbean in Barbados African       
#>  9 HG01894   ACB   AFR       female African Caribbean in Barbados African       
#> 10 HG01896   ACB   AFR       female African Caribbean in Barbados African       
#> # ℹ 2,494 more rows
table(samples_1kg$super_pop)
#> 
#> AFR AMR EAS EUR SAS 
#> 661 347 504 503 489 
```
