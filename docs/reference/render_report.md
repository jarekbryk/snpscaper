# Create an HTML report for a region of the genome

Produces a handout-style HTML report for a gene or region: a summary of
the data, PCA, admixture, Fst, heterozygosity, allele frequencies of key
SNPs (with maps) and LD, plus a methods section recording exactly how
the results were produced. The report is a single self-contained HTML
file that can be emailed, uploaded to a virtual learning environment, or
printed.

## Usage

``` r
render_report(
  x,
  pop = NULL,
  super_pop = NULL,
  focal_snps = NULL,
  k = 2:5,
  reps = 5,
  min_maf = 0.01,
  seed = 1,
  version = c("staff", "student"),
  title = NULL,
  output_dir = ".",
  verbose = TRUE
)
```

## Arguments

- x:

  What to report on. Either a character vector of gene symbols or
  coordinates (one report is made for each; see
  [`fetch_region()`](https://jarekbryk.github.io/snpscaper/reference/fetch_region.md)),
  or a `snpscape` object you have prepared yourself (e.g. combining
  several regions).

- pop, super_pop:

  Populations or super-populations to include (see
  [`keep_pops()`](https://jarekbryk.github.io/snpscaper/reference/keep_pops.md)).
  `NULL` includes everyone.

- focal_snps:

  rs numbers or SNP IDs to show allele frequencies for.

- k:

  The numbers of ancestral groups for admixture analysis.

- reps:

  Number of admixture runs for each `K`.

- min_maf:

  Minimum minor allele frequency for SNPs (see
  [`filter_snps()`](https://jarekbryk.github.io/snpscaper/reference/filter_snps.md)).

- seed:

  Random seed, so that the same settings always give the same report.

- version:

  `"staff"` (with notes for staff) or `"student"` (without).

- title:

  Report title. By default, made from the region name.

- output_dir:

  Folder to save reports in. Created if needed.

- verbose:

  Print progress messages.

## Value

The paths of the created HTML files, invisibly.

## Details

### Staff and student versions

The staff version (the default) adds a "Notes for staff" section with
the expected pattern and its explanation, for genes in
[teaching_genes](https://jarekbryk.github.io/snpscaper/reference/teaching_genes.md).
The student version leaves this out, so students can interpret the
results themselves.

### Which data are used for what

The report applies
[`filter_snps()`](https://jarekbryk.github.io/snpscaper/reference/filter_snps.md)
itself, using `min_maf`:

- PCA and admixture use LD-pruned SNPs, because these methods assume
  SNPs are independent.

- Allele frequencies, Fst and heterozygosity use all SNPs that pass the
  frequency filter.

- LD heatmaps use SNPs with a minor allele frequency of at least 5%, one
  heatmap per super-population, without pruning.

### Key SNPs

The allele frequency section shows the SNPs given in `focal_snps`. If
none are given, it uses the SNPs listed for the gene in
[teaching_snps](https://jarekbryk.github.io/snpscaper/reference/teaching_snps.md),
or, for other regions, the SNP whose frequency differs most between
super-populations.

### Requirements and run time

Reports are made with Quarto, which is installed together with RStudio
(or from <https://quarto.org>). Making a report takes from under a
minute to a few minutes, mostly for the admixture analysis; use fewer
`k` values or `reps` to speed it up.

## Examples

``` r
if (FALSE) { # \dontrun{
render_report("LCT")

# one report per gene, for a chosen set of populations
render_report(
  c("LCT", "SLC24A5", "ACKR1"),
  pop = c("YRI", "LWK", "GBR", "FIN", "CHB", "JPT", "GIH", "PEL"),
  output_dir = "handouts"
)

# student version, for a region with no well-known genes
render_report("chr18:1000000-1500000", version = "student")
} # }
```
