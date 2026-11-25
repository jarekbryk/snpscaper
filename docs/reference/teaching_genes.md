# Genes with well-known patterns of variation between populations

A curated set of genes that make good teaching examples, from genes
under strong selection in particular populations (e.g. LCT, SLC24A5,
ACKR1) to genes with more ordinary patterns of variation (CYP1A2, HTT).
Any of them can be passed to
[`fetch_region()`](https://jarekbryk.github.io/snpscaper/reference/fetch_region.md)
by name.

## Usage

``` r
teaching_genes
```

## Format

A tibble with one row per gene:

- gene:

  Gene symbol.

- chrom, gene_start, gene_end, strand:

  Gene position on GRCh38.

- region_start, region_end:

  The region fetched by default: the gene plus 25 kb on each side, which
  includes nearby regulatory variants.

- description:

  A short, neutral description of what the gene does.

- staff_notes:

  The expected pattern of variation and its explanation, for teaching
  staff.

## Source

Gene positions from Ensembl (GRCh38). Built by
`data-raw/teaching_genes.R`.

## See also

[teaching_snps](https://jarekbryk.github.io/snpscaper/reference/teaching_snps.md)
for the key SNPs in each gene.
