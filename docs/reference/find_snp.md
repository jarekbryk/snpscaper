# Find the position of a SNP from its rs number

Most well-studied SNPs have an "rs number" (reference SNP ID) given by
the dbSNP database, for example `rs4988235`, the main lactase
persistence variant in Europeans. The 1000 Genomes files do not use rs
numbers, so snpscaper looks up where each rs number sits in the genome
(GRCh38) and finds it by position.

## Usage

``` r
find_snp(rsid)
```

## Arguments

- rsid:

  An rs number such as `"rs4988235"`, or a vector of them.

## Value

A tibble with one row per SNP: `rsid`, `chrom`, `pos` and `alleles` (as
reported by Ensembl, reference allele first, e.g. `"G/A"`).

## Examples

``` r
if (FALSE) { # \dontrun{
find_snp("rs4988235")
} # }
```
