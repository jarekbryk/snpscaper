# The snpscape object

A `snpscape` object holds everything snpscaper knows about a set of
genotypes. It is created by
[`fetch_region()`](https://jarekbryk.github.io/snpscaper/reference/fetch_region.md)
and passed along to the other functions. Printing it gives a short
summary.

## Details

It is a list with four parts:

- `genotypes`: a matrix with one row per person and one column per SNP.
  Each value is the number of copies of the alternative allele that the
  person carries (0, 1 or 2), or `NA` if unknown.

- `variants`: a tibble describing each SNP (one row per column of
  `genotypes`): `id`, `chrom`, `pos` (position on GRCh38), `ref`
  (reference allele), `alt` (alternative allele), `region` (which
  requested region it came from) and `rsid` (rs number, for well-known
  SNPs only).

- `samples`: a tibble describing each person (one row per row of
  `genotypes`); see
  [samples_1kg](https://jarekbryk.github.io/snpscaper/reference/samples_1kg.md).

- `provenance`: where the data came from and what has been done to them:
  the regions fetched, genome build, data source, date, package version
  and a log of filtering steps.

You rarely need to touch these parts directly, but you can, e.g.
`x$samples` or `x$variants`.
