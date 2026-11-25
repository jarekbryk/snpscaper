# Key SNPs in the teaching genes

One or two well-studied SNPs in each of the
[teaching_genes](https://jarekbryk.github.io/snpscaper/reference/teaching_genes.md),
usually the variant responsible for the gene's striking pattern of
variation.

## Usage

``` r
teaching_snps
```

## Format

A tibble with one row per SNP:

- gene:

  Gene symbol.

- rsid:

  rs number (dbSNP identifier).

- chrom, pos:

  Position on GRCh38.

- alleles:

  Alleles reported by Ensembl, reference allele first.

- note:

  What the variant is known for.

## Source

SNP positions from Ensembl (GRCh38). Built by
`data-raw/teaching_genes.R`.
