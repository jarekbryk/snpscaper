# Fetch genotypes for a region of the genome

Downloads the genotypes of the 2,504 unrelated people in the 1000
Genomes Project for one or more regions of the genome, and returns them
as a `snpscape` object ready for analysis.

## Usage

``` r
fetch_region(regions, flank = 25000, verbose = TRUE)
```

## Arguments

- regions:

  A character vector of gene symbols and/or coordinates.

- flank:

  Number of base pairs to add on each side of a gene. Ignored for
  coordinates.

- verbose:

  Print progress messages.

## Value

A `snpscape` object (see
[snpscape](https://jarekbryk.github.io/snpscaper/reference/snpscape.md)).

## Details

### What you can ask for

Each element of `regions` can be:

- a **gene symbol**, such as `"LCT"` or `"MC1R"`. The gene is looked up
  with
  [`find_gene()`](https://jarekbryk.github.io/snpscaper/reference/find_gene.md)
  and `flank` base pairs are added on each side, because variants that
  regulate a gene often lie just outside it.

- a set of **coordinates** on the GRCh38 assembly, such as
  `"chr2:135787850-135862184"`. The `chr` prefix and commas are
  optional, so `"2:135,787,850-135,862,184"` also works.

Only the autosomes (chromosomes 1 to 22) are supported. Large regions
take longer and need more memory: a 1 Mb region typically holds around
25,000 SNPs.

### What you get

Only **biallelic SNPs** are kept: positions where exactly two DNA
letters are seen in the population (for example `G` and `A`).
Insertions, deletions, larger structural variants, and positions with
more than two alleles are dropped, because the analyses in this package
assume each variant has exactly two forms.

Each person's genotype at a SNP is stored as the number of copies of the
alternative allele (`alt`) they carry: 0, 1 or 2. The reference allele
(`ref`) is simply the letter found in the human reference genome; it is
not necessarily the ancestral or the most common allele.

The 698 relatives (mostly children of parent pairs) who are also in the
high-coverage panel are not included, because close relatives share long
stretches of DNA and would distort population analyses.

### Caching

Fetched regions are saved on disk, so fetching the same region again is
instant and works offline. See
[`clear_cache()`](https://jarekbryk.github.io/snpscaper/reference/clear_cache.md).

## Examples

``` r
if (FALSE) { # \dontrun{
lct <- fetch_region("LCT")
lct

# several regions at once
pigmentation <- fetch_region(c("MC1R", "SLC24A5"))

# a stretch of chromosome 18 with no well-known genes under selection
neutral <- fetch_region("chr18:1000000-1500000")
} # }
```
