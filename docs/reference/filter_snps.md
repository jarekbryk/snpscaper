# Remove rare, poorly genotyped or redundant SNPs

Most analyses of population structure work best on a cleaned-up set of
SNPs. `filter_snps()` applies up to three filters, in this order, and
records how many SNPs remain after each one (see
[`summary.snpscape()`](https://jarekbryk.github.io/snpscaper/reference/summary.snpscape.md)).

## Usage

``` r
filter_snps(
  x,
  min_maf = 0.01,
  max_missing = 0.05,
  prune_ld = FALSE,
  ld_r2 = 0.2,
  ld_window = 50
)
```

## Arguments

- x:

  A `snpscape` object.

- min_maf:

  Minimum minor allele frequency, between 0 and 0.5.

- max_missing:

  Maximum proportion of people with a missing genotype.

- prune_ld:

  Remove SNPs in strong linkage disequilibrium with a neighbour (`TRUE`
  or `FALSE`).

- ld_r2:

  Maximum r² allowed between SNPs within a window when
  `prune_ld = TRUE`.

- ld_window:

  Window size, in number of SNPs, for LD pruning.

## Value

The filtered `snpscape` object.

## Details

### 1. Minor allele frequency (`min_maf`)

The minor allele frequency (MAF) is the frequency of the less common of
a SNP's two alleles. Most SNPs in the human genome are very rare: they
are carried by only one or a handful of people. Rare SNPs say little
about differences between populations, add noise, and (when nobody in
the data carries the minor allele at all) cannot be analysed by PCA.
Removing SNPs with a MAF below 1% (`min_maf = 0.01`) is a common choice.

MAF is calculated from the people currently in the data, so if you later
remove populations with
[`keep_pops()`](https://jarekbryk.github.io/snpscaper/reference/keep_pops.md),
some SNPs may become rare again.

### 2. Missing genotypes (`max_missing`)

SNPs for which genotypes are unknown for many people are removed. The
1000 Genomes high-coverage data have very few missing genotypes, so this
filter rarely removes anything, but it is good practice.

### 3. Linkage disequilibrium pruning (`prune_ld`)

Nearby SNPs tend to be inherited together, because chromosomes are
passed down in long chunks that are only occasionally broken up by
recombination. This correlation between neighbouring SNPs is called
linkage disequilibrium (LD). In a small region, many SNPs may be
strongly correlated and so carry almost the same information. Methods
like PCA and admixture analysis assume each SNP is independent evidence;
a block of correlated SNPs effectively gets counted many times and can
dominate the result.

LD pruning keeps a set of SNPs that are only weakly correlated with each
other. It moves a window of `ld_window` SNPs along the region and,
within each window, whenever two SNPs have a squared correlation (r²)
above `ld_r2`, removes the one with the lower minor allele frequency.
This is a simplified version of the widely used `--indep-pairwise`
method of the PLINK program.

Pruning is off by default. Turn it on (`prune_ld = TRUE`) for PCA and
admixture analysis. Leave it off when you are interested in LD itself or
in particular SNPs.

## Examples

``` r
lct <- example_lct |>
  filter_snps(min_maf = 0.05, prune_ld = TRUE)
summary(lct)
#> <snpscape> 2,504 people x 12 SNPs (GRCh38)
#> Regions:     LCT (chr2:135,762,850-135,862,184) 
#> Populations: 26 populations in 5 super-populations
#> Steps:       fetched -> MAF >= 0.05 -> missing <= 0.05 -> LD pruned (r2 <= 0.2) 
#> 
#> People per population:
#> # A tibble: 26 × 4
#>    super_pop pop   pop_name                                               people
#>    <chr>     <chr> <chr>                                                   <int>
#>  1 AFR       ACB   African Caribbean in Barbados                              96
#>  2 AFR       ASW   African Ancestry in Southwest US                           61
#>  3 AFR       ESN   Esan in Nigeria                                            99
#>  4 AFR       GWD   Gambian in Western Division, The Gambia                   113
#>  5 AFR       LWK   Luhya in Webuye, Kenya                                     99
#>  6 AFR       MSL   Mende in Sierra Leone                                      85
#>  7 AFR       YRI   Yoruba in Ibadan, Nigeria                                 108
#>  8 AMR       CLM   Colombian in Medellin, Colombia                            94
#>  9 AMR       MXL   Mexican Ancestry in Los Angeles, CA, USA                   64
#> 10 AMR       PEL   Peruvian in Lima, Peru                                     85
#> 11 AMR       PUR   Puerto Rican in Puerto Rico                               104
#> 12 EAS       CDX   Chinese Dai in Xishuangbanna, China                        93
#> 13 EAS       CHB   Han Chinese in Beijing, China                             103
#> 14 EAS       CHS   Han Chinese South, China                                  105
#> 15 EAS       JPT   Japanese in Tokyo, Japan                                  104
#> 16 EAS       KHV   Kinh in Ho Chi Minh City, Vietnam                          99
#> 17 EUR       CEU   Utah residents with Northern and Western European anc…     99
#> 18 EUR       FIN   Finnish in Finland                                         99
#> 19 EUR       GBR   British from England and Scotland                          91
#> 20 EUR       IBS   Iberian populations in Spain                              107
#> 21 EUR       TSI   Toscani in Italy                                          107
#> 22 SAS       BEB   Bengali in Bangladesh                                      86
#> 23 SAS       GIH   Gujarati Indians in Houston, TX, USA                      103
#> 24 SAS       ITU   Indian Telugu in the UK                                   102
#> 25 SAS       PJL   Punjabi in Lahore, Pakistan                                96
#> 26 SAS       STU   Sri Lankan Tamil in the UK                                102
#> 
#> Processing steps:
#> # A tibble: 4 × 3
#>   step                  people  snps
#>   <chr>                  <int> <int>
#> 1 fetched                 2504  1922
#> 2 MAF >= 0.05             2504   191
#> 3 missing <= 0.05         2504   191
#> 4 LD pruned (r2 <= 0.2)   2504    12
```
