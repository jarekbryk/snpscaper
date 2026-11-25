# Making reports and preparing classes

[`render_report()`](https://jarekbryk.github.io/snpscaper/reference/render_report.md)
produces a handout-style HTML report for any gene or region, with all
the analyses from the other vignettes in one document. Reports can be
printed, emailed or uploaded to a virtual learning environment as a
single file. They are intended both as teaching material (results for
students to interpret) and as a quick way for staff to explore a region.

## Making a report

``` r

library(snpscaper)

render_report("LCT")
```

This fetches the region, runs every analysis and saves
`snpscaper_LCT_staff.html` in the current folder. It usually takes well
under a minute.

The report contains:

1.  **The region**: its coordinates, and a short description for genes
    in `teaching_genes`.
2.  **The data**: the number of people per population, and the number of
    SNPs left after each filtering step.
3.  **PCA**: the first principal components, coloured by
    super-population and by population, and the share of variation each
    component captures.
4.  **Admixture analysis**: cross-entropy for each *K*, and ancestry
    proportions grouped by super-population and by population.
5.  **Fst**: a heatmap for pairs of populations, and a table for pairs
    of super-populations.
6.  **Heterozygosity**: a plot and a table for each population.
7.  **Allele frequencies of key SNPs**: a bar chart and a world map for
    each key SNP.
8.  **Linkage disequilibrium**: a heatmap for each super-population.
9.  **Methods**: data source, filters, settings, random seed, package
    version and references, so that results can be reproduced and cited.
10. **Notes for staff** (staff version only): the expected pattern and
    its explanation.

Figures and tables are numbered, so students can refer to them (e.g.
“Figure 6”) in their answers.

## Staff and student versions

The staff version is the default. It adds notes on the expected patterns
for genes in `teaching_genes`. The student version leaves these out:

``` r

render_report("LCT", version = "student")
```

## Choosing what goes in

``` r

render_report(
  "SLC24A5",
  pop = c("YRI", "LWK", "GBR", "FIN", "CHB", "JPT", "GIH", "PEL"),
  focal_snps = "rs1426654",
  k = 2:4,
  title = "Skin pigmentation: SLC24A5",
  output_dir = "handouts"
)
```

- `pop` or `super_pop` choose the populations. Fewer populations make
  the plots easier to read.
- `focal_snps` chooses the SNPs for the allele frequency section, as rs
  numbers or SNP IDs. By default the report uses the key SNPs in
  `teaching_snps`, or, for any other region, the SNP whose frequency
  differs most between super-populations.
- `k` and `reps` control the admixture analysis, which is the slowest
  part.
- `min_maf` sets the minor allele frequency filter.
- `seed` makes the admixture results reproducible: the same settings
  always give the same report.

### Several reports at once

Give several genes or regions to make one report for each:

``` r

render_report(
  c("LCT", "SLC24A5", "EDAR", "ALDH2", "ACKR1", "HBB"),
  version = "student",
  output_dir = "handouts"
)
```

### Combining regions in one report

To analyse several regions together, fetch them first and pass the
result:

``` r

pigmentation <- fetch_region(c("SLC24A5", "MC1R"))
render_report(pigmentation, title = "Pigmentation genes")
```

## Preparing for a class

### Fetching data in advance

Every region is downloaded once and then kept in a cache folder, so
students working on the same regions do not all download them during the
class. By default, each user has their own cache. To share one, set the
cache folder before fetching, for example to a shared network drive:

``` r

options(snpscaper.cache_dir = "S:/teaching/snpscaper_cache")

for (gene in teaching_genes$gene) {
  fetch_region(gene)
}
```

Students then set the same option (for example at the top of a script
you provide), and every
[`fetch_region()`](https://jarekbryk.github.io/snpscaper/reference/fetch_region.md)
call for those regions is instant and works offline.

### Requirements

Reports need [Quarto](https://quarto.org), which is installed with
recent versions of RStudio, and the FactoMineR and factoextra packages.
Fetching data needs internet access to `ftp.1000genomes.ebi.ac.uk` and,
for genes not in `teaching_genes`, `rest.ensembl.org`.
