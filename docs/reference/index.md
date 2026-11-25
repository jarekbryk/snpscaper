# Package index

## Fetching data

- [`fetch_region()`](https://jarekbryk.github.io/snpscaper/reference/fetch_region.md)
  : Fetch genotypes for a region of the genome
- [`find_gene()`](https://jarekbryk.github.io/snpscaper/reference/find_gene.md)
  : Find the position of a gene
- [`find_snp()`](https://jarekbryk.github.io/snpscaper/reference/find_snp.md)
  : Find the position of a SNP from its rs number
- [`clear_cache()`](https://jarekbryk.github.io/snpscaper/reference/clear_cache.md)
  : Delete downloaded data
- [`snpscape`](https://jarekbryk.github.io/snpscaper/reference/snpscape.md)
  : The snpscape object
- [`summary(`*`<snpscape>`*`)`](https://jarekbryk.github.io/snpscaper/reference/summary.snpscape.md)
  : Summarise a snpscape object

## Preparing data

- [`filter_snps()`](https://jarekbryk.github.io/snpscaper/reference/filter_snps.md)
  : Remove rare, poorly genotyped or redundant SNPs
- [`keep_pops()`](https://jarekbryk.github.io/snpscaper/reference/keep_pops.md)
  : Keep only some populations or people
- [`as_pca_table()`](https://jarekbryk.github.io/snpscaper/reference/as_pca_table.md)
  : Prepare genotypes for principal component analysis

## Admixture

- [`run_admixture()`](https://jarekbryk.github.io/snpscaper/reference/run_admixture.md)
  : Estimate ancestry proportions (admixture analysis)
- [`get_ancestry()`](https://jarekbryk.github.io/snpscaper/reference/get_ancestry.md)
  : Get ancestry proportions from an admixture analysis
- [`plot_admixture()`](https://jarekbryk.github.io/snpscaper/reference/plot_admixture.md)
  : Plot ancestry proportions
- [`plot_cross_entropy()`](https://jarekbryk.github.io/snpscaper/reference/plot_cross_entropy.md)
  : Plot cross-entropy to help choose K

## Allele frequencies, diversity and LD

- [`calc_allele_freqs()`](https://jarekbryk.github.io/snpscaper/reference/calc_allele_freqs.md)
  : Calculate allele frequencies in each population
- [`plot_freqs()`](https://jarekbryk.github.io/snpscaper/reference/plot_freqs.md)
  : Plot the allele frequency of a SNP in each population
- [`map_freqs()`](https://jarekbryk.github.io/snpscaper/reference/map_freqs.md)
  : Map the allele frequency of a SNP
- [`calc_fst()`](https://jarekbryk.github.io/snpscaper/reference/calc_fst.md)
  : Measure genetic differentiation between populations (Fst)
- [`plot_fst()`](https://jarekbryk.github.io/snpscaper/reference/plot_fst.md)
  : Plot Fst between pairs of populations
- [`calc_heterozygosity()`](https://jarekbryk.github.io/snpscaper/reference/calc_heterozygosity.md)
  : Measure genetic diversity within populations (heterozygosity)
- [`plot_heterozygosity()`](https://jarekbryk.github.io/snpscaper/reference/plot_heterozygosity.md)
  : Plot heterozygosity in each population
- [`calc_ld()`](https://jarekbryk.github.io/snpscaper/reference/calc_ld.md)
  : Calculate linkage disequilibrium between SNPs
- [`plot_ld()`](https://jarekbryk.github.io/snpscaper/reference/plot_ld.md)
  : Plot linkage disequilibrium as a heatmap

## Reports

- [`render_report()`](https://jarekbryk.github.io/snpscaper/reference/render_report.md)
  : Create an HTML report for a region of the genome

## Plot styling

- [`theme_jb()`](https://jarekbryk.github.io/snpscaper/reference/theme_jb.md)
  : A clean ggplot2 theme for snpscaper plots
- [`superpop_colours`](https://jarekbryk.github.io/snpscaper/reference/colours.md)
  [`pop_colours`](https://jarekbryk.github.io/snpscaper/reference/colours.md)
  : Colours for populations and super-populations
- [`scale_colour_superpop()`](https://jarekbryk.github.io/snpscaper/reference/scale_colour_superpop.md)
  [`scale_fill_superpop()`](https://jarekbryk.github.io/snpscaper/reference/scale_colour_superpop.md)
  [`scale_colour_pop()`](https://jarekbryk.github.io/snpscaper/reference/scale_colour_superpop.md)
  [`scale_fill_pop()`](https://jarekbryk.github.io/snpscaper/reference/scale_colour_superpop.md)
  : Colour and fill scales for populations

## Data

- [`samples_1kg`](https://jarekbryk.github.io/snpscaper/reference/samples_1kg.md)
  : People in the 1000 Genomes Project
- [`pop_locations`](https://jarekbryk.github.io/snpscaper/reference/pop_locations.md)
  : Map locations of the 1000 Genomes populations
- [`teaching_genes`](https://jarekbryk.github.io/snpscaper/reference/teaching_genes.md)
  : Genes with well-known patterns of variation between populations
- [`teaching_snps`](https://jarekbryk.github.io/snpscaper/reference/teaching_snps.md)
  : Key SNPs in the teaching genes
- [`example_lct`](https://jarekbryk.github.io/snpscaper/reference/example_lct.md)
  : Example data: the LCT region
- [`example_neutral`](https://jarekbryk.github.io/snpscaper/reference/example_neutral.md)
  : Example data: a neutral region of chromosome 18
