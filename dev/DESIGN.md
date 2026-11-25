# snpscaper — design notes

Living document recording the agreed architecture and decisions for turning
the `snpscaper-alpha` notebooks into the `snpscaper` R package.

Last updated: 2026-11-25

## Purpose

Teach variation in human populations (final-year biology undergraduates) using
the 1000 Genomes Project data. A user names a genomic region or gene; the
package fetches the genotypes, prepares them, and produces PCA, admixture,
allele-frequency, Fst, heterozygosity and LD results, plus an HTML handout.

**Key principle: everything happens in R.** No command line, no browser-based
data slicer, no external binaries (PLINK, ADMIXTURE). Users are students and
staff who know R functions but not bioinformatics tooling.

## Code and documentation principles

- **Simple, straightforward code.** Short functions doing one thing, readable
  over clever, no unnecessary abstraction or dependencies. A biology lecturer
  reading a function should be able to follow it. The alpha code is a
  reference for *what* to do, not *how* — rewrite rather than port.
- **Generous documentation.** Assume a reader who does not know what a VCF is,
  what the 1000 Genomes Project is, or why we filter on allele frequency or
  prune SNPs in high LD. Roxygen docs, vignettes and in-code comments explain
  the *why* (the biology and the statistics), not just the *what*.
- Comments in code explain reasoning and concepts; they do not narrate each
  line.

## Decisions

| Topic | Decision |
|---|---|
| Package name | `snpscaper` |
| Repository | `github.com/jarekbryk/snpscaper` (GitHub only, no CRAN) |
| Licence | MIT (`MIT + file LICENSE`) |
| Code style | tidyverse style guide; snake_case; plain-verb function names |
| Pipe | native `\|>`; `Depends: R (>= 4.1)` |
| Plot theme | `theme_jb()` (house style), exported by the package |
| Genome build | GRCh38 only |
| Data source | 1000 Genomes 30x high-coverage phased panel (EBI FTP), remote tabix access; related samples removed by default (2,504 unrelated) |
| Chromosomes | autosomes first; chrX later (male hemizygosity) |
| Platform | mostly Windows lab machines, imaged with full R stack + Rtools; install burden not a concern |
| Installation | `pak::pak("jarekbryk/snpscaper")` (handles Bioconductor deps) |
| Interface | functions only, no Shiny app |
| Wrapping policy | wrap what users never see (fetching, transforming, admixture, report); use standard packages directly where students learn them elsewhere (PCA via FactoMineR/factoextra) |
| LD pruning | parameter: `filter_snps(prune_ld = TRUE)` |
| Report | Quarto, HTML only (self-contained, print-friendly); **staff version is the default**; student version later |
| Report users | anyone can call `render_report()` → needs friendly input validation and messages |

## Typical workflow

```r
pacman::p_load(snpscaper, FactoMineR, factoextra)

lct <- fetch_region("LCT", flank = 50e3) |>
  filter_snps(min_maf = 0.01, prune_ld = TRUE) |>
  keep_pops(super_pop = c("AFR", "EUR", "EAS"))

lct   # readable summary

# PCA with standard tools
pca <- PCA(as_pca_table(lct), quali.sup = 1:3, graph = FALSE)
fviz_pca_ind(pca, habillage = "super_pop", geom = "point",
             palette = superpop_colours) + theme_jb()

# Admixture (wrapped LEA::snmf)
adm <- run_admixture(lct, k = 2:6)
plot_cross_entropy(adm)
plot_admixture(adm, k = 3, group = pop)

# Other analyses
calc_allele_freqs(lct) |> plot_freqs(snp = "rs4988235")
map_freqs(lct, snp = "rs4988235")
calc_fst(lct, by = pop) |> plot_fst()
calc_heterozygosity(lct) |> plot_heterozygosity()
plot_ld(lct)

# Handout
render_report("LCT", pops = c("YRI", "GBR", "CHB", "PEL"),
              focal_snps = "rs4988235", output_dir = "handouts")
```

## The `snpscape` object

S3 list, created by `fetch_region()`, carried through every step:

- `genotypes` — integer matrix, samples × SNPs, values 0/1/2/NA (alt allele count)
- `variants` — tibble: chrom, pos, rsid, ref, alt, alt allele frequency, gene/region
- `samples` — tibble: sample_id, pop, super_pop, sex, population name
- `provenance` — regions requested, build, data source, filters applied, fetch date, package version

Samples are always matched **by ID**, never by row position (fixes the label
misalignment risk in the alpha code). Methods: `print()`, `summary()`, `dim()`.

## Functions

| Area | Functions |
|---|---|
| Fetch | `find_gene()` (symbol → GRCh38 coordinates via Ensembl REST), `fetch_region()` (genes and/or coordinates, multiple regions, unrelated samples by default, cached), `clear_cache()` |
| Transform | `filter_snps(min_maf, max_missing, prune_ld = FALSE, ld_r2 = 0.2, ld_window)`, `keep_pops(pop, super_pop, sex)`, `as_pca_table()` |
| Analyse | `run_admixture(k, reps, seed)`, `get_ancestry(adm, k, run = "best")`, `calc_allele_freqs(by)`, `calc_fst(by)` (pairwise Hudson), `calc_heterozygosity(by)`, `calc_ld()` |
| Plot (return ggplot objects, `theme_jb()`) | `plot_admixture()`, `plot_cross_entropy()`, `plot_freqs()`, `map_freqs()`, `plot_fst()`, `plot_heterozygosity()`, `plot_ld()` |
| Style | `theme_jb()`, `pop_colours`, `superpop_colours`, `scale_colour_superpop()`, `scale_fill_superpop()` |
| Report | `render_report()` |

Notes:

- `as_pca_table()` returns sample_id, pop, super_pop first, then SNP columns, ready for `FactoMineR::PCA(quali.sup = 1:3)`. Monomorphic SNPs are removed by filtering, so `scale.unit = TRUE` is safe.
- LD pruning: greedy windowed r² pruner in R (no SNPRelate/bigsnpr).
- Fst and heterozygosity computed in-house from allele counts.
- `run_admixture()` runs `LEA::snmf` in a temporary directory; no `setwd()`; seeded.
- `map_freqs()` uses bundled population coordinates placed at ancestral origin, with a caption noting diaspora sampling (CEU, GIH, ITU, STU, ASW, ACB).
- Caching: `tools::R_user_dir("snpscaper", "cache")`, overridable with `options(snpscaper.cache_dir = ...)` so staff can pre-populate a shared cache before a practical.

## Report (`render_report()`)

- Input: a region/gene (or vector of them, one report each), or a prepared `snpscape` object.
- Parameters: populations, filters, K range, admixture reps, focal SNPs, seed, title, output_dir.
- Template in `inst/quarto/report.qmd`, copied to a temp directory and rendered via the `quarto` R package; requires the Quarto CLI (bundled with RStudio).
- Output: single self-contained HTML with print CSS; numbered figures and tables.
- Default = **staff version**: includes teaching notes / expected patterns from `teaching_genes`. A student version (results only, no interpretation, optional questions) comes later.
- Contents:
  1. Header: gene, GRCh38 coordinates, populations, short gene description
  2. Data summary: samples per population; SNPs before/after each filter
  3. PCA: PC1–PC2, PC1–PC3, scree plot
  4. Admixture: cross-entropy plot; bar plots for each K
  5. Fst: heatmap and table
  6. Heterozygosity: table and plot
  7. Focal SNP(s): frequency plot and map
  8. LD heatmap
  9. Methods and provenance (source, build, filters, seed, package version, date)
  10. Staff notes (staff version only)
- Admixture is the slow step: lighter defaults in the report (e.g. 5 reps, LD pruning on); cache per-region results for batch rendering.

## Bundled data

- `samples_1kg` — tidy sample metadata with population names
- `teaching_genes` — curated regions (MC1R, LCT, ALDH2, HBB, CYP1A2, HTT, plus neutral blocks) with GRCh38 coordinates, a neutral description, and staff notes
- `pop_locations` — coordinates for mapping
- `example_*` — small pre-fetched `snpscape` objects so examples, tests and vignettes run offline
- `data-raw/` — scripts that build all of the above

## Package layout

```
snpscaper/
├── DESCRIPTION, NAMESPACE (roxygen2), LICENSE, README.md, NEWS.md
├── R/
│   fetch.R lookup.R cache.R snpscape.R filter.R ld.R admixture.R
│   freqs.R fst.R heterozygosity.R maps.R report.R theme.R palettes.R
│   data.R globals.R helpers.R
├── data/  data-raw/  inst/extdata/  inst/quarto/
├── tests/testthat/
├── vignettes/   getting-started, pca-with-factominer, admixture, teaching-genes
├── man/figures/
├── _pkgdown.yml   docs/
└── .github/workflows/   R-CMD-check (incl. windows-latest)
```

## Conventions (following maudr)

- `DESCRIPTION`: `Authors@R` with ORCID, `URL` and `BugReports`, `Roxygen: list(markdown = TRUE)`, knitr/rmarkdown vignettes.
- README: short intro, emoji section headers, GitHub-only install, a step-by-step description of each main function, worked example, acknowledgements.
- pkgdown: Bootstrap 5, `simplex` bootswatch, site built into `docs/` and served via GitHub Pages; `man/figures` as assets.
- Explicit namespacing (`dplyr::mutate()`) inside package functions; `utils::globalVariables()` in `R/globals.R`.
- `# --- Section ---` comment headers inside longer functions.
- `verbose = TRUE` argument with `message()` progress output for staff-facing functions.
- Package-level settings stored in `options(snpscaper.*)`.
- Demo/offline mode via bundled data.
- `.gitignore` / `.Rbuildignore` as in maudr (incl. `.claude`, `.positai`, `.quarto`).

## Risks / first tasks

1. **Remote tabix on Windows** (Rsamtools/Rhtslib with network access). Spike first on a lab-type Windows machine. Fallback: shared pre-fetched cache via `options(snpscaper.cache_dir)`.
2. **Data source details**: confirm exact GRCh38 panel file URLs, chromosome naming (`chr2`), and the list of 2,504 unrelated samples.
3. Ensembl REST availability during class → cache gene lookups; bundle coordinates for `teaching_genes`.

## Data source findings (2026-11-25, tested on macOS)

- Genotypes: `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr{N}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz` (+ `.tbi`). chrX file has a `.v2` suffix.
- Samples: `.../1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt` (3,202 samples, sex, pop, super-pop). The 2,504 unrelated are listed in `1000G_2504_high_coverage.sequence.index`; same IDs as the alpha `populations.txt`.
- Chromosomes are named `chr2` etc. All alpha sample files are GRCh38 too.
- Remote tabix via `Rsamtools::scanTabix()` works: LCT region (~70 kb, ~1,700 records) in ~15 s on first call, <1 s once the index is local.
- **htslib saves the `.tbi` index into the working directory** unless given one. The package downloads indexes into the cache dir itself and passes `TabixFile(url, index = ...)`.
- **No rsIDs**: the ID column is `chr:pos:ref:alt` (or `HGSV_*` for SVs). rsIDs (e.g. for `focal_snps`) must be resolved to positions, e.g. via Ensembl REST; bundle positions for `teaching_genes` focal SNPs.
- **Structural variants** are in the same files (symbolic ALT like `<INS:ME:SVA>`), and long SVs starting far upstream are returned for overlapping regions. Keep only biallelic SNVs.
- Text size: ~15 MB per 50 kb region for 3,202 samples, so a 1.5 Mb region is ~0.5 GB of text. Parse efficiently (e.g. `VariantAnnotation::readVcf()` with sample subsetting and GT only) and warn on large regions.
- Teaching note: the lactase persistence variant rs4988235 (chr2:135,851,076) lies ~14 kb upstream of LCT, in MCM6. `teaching_genes` coordinates must include such regulatory regions, not just the gene body.
- Windows: untested. Test script in `spikes/windows_fetch_test.R`.

## Known issues in the alpha code (for reference)

- Functions duplicated across notebooks with drift; global `population_data`; metadata assumed in first 4 columns.
- `cbind()` of admixture Q-matrix with metadata assumes matching sample order.
- `order("pop")` always returns 1 (sorting never happens).
- `rename()` arguments appear reversed.
- `get_maf_data()` uses undefined `mafdata_id`; deprecated `count_()`.
- `setwd()` inside `create_admixture_data()`; hard-coded file names.
- No missing-data handling, MAF filtering or LD pruning.
- Many unused package dependencies.

## Later

- Student version of the report (no interpretation, optional questions file).
- Per-student region allocation (`assign_regions()`), personalised reports.
- `export = TRUE`: write figures (PNG) and tables (CSV) alongside the report.
- chrX support.
- PDF output (Typst).

## Progress

- 2026-11-25: package created at `software/snpscaper` (git initialised, nothing committed). Done: bundled data (`samples_1kg`, `pop_locations`, `teaching_genes` (9 genes, added SLC24A5, EDAR, ACKR1), `teaching_snps`, `example_lct`), `find_gene()`, `find_snp()`, `fetch_region()` with caching, `snpscape` object, `filter_snps()` (incl. LD pruning), `keep_pops()`, `as_pca_table()`, `theme_jb()`, colours and scales, `clear_cache()`. Tests pass; R CMD check clean.
- 2026-11-25: added `run_admixture()`, `get_ancestry()`, `plot_admixture()`, `plot_cross_entropy()`, `calc_allele_freqs()`, `plot_freqs()`, `map_freqs()`, `calc_fst()`, `plot_fst()`, `calc_heterozygosity()`, `plot_heterozygosity()`, `calc_ld()`, `plot_ld()`. Tests pass (incl. live fetch); R CMD check clean.
  - Interface decisions made while building: grouping arguments are strings (`by = "pop"`, `group = "super_pop"`), not bare names; every plot function takes the matching `calc_*()` result (`calc_ld(x) |> plot_ld()`, `calc_allele_freqs(x) |> map_freqs(snp = ...)`); admixture groups are matched across K by correlation so colours stay comparable.
  - Teaching observations: LD around LCT is a long high-r² block in EUR, broken up in AFR (sweep). Single gene regions need not show AFR as most diverse (LCT: CHB > YRI after MAF filter); neutral chr18 block does.
- 2026-11-25: added `render_report()` with `inst/quarto/report.qmd` + `report.css`. A character `x` gives one report per element; a prepared `snpscape` gives one report. `version = "student"` simply omits the staff notes (fuller student version still later). Focal SNPs: user-given, else `teaching_snps` for the gene, else the SNP with the largest frequency range across super-populations. List parameters are passed to Quarto as comma-separated text (Quarto rejects empty lists). Reports take ~15-25 s and are 4-6 MB. Fixed: empty fetch windows crashed `fetch_region()`.
- 2026-11-25: added `example_neutral` (chr18:1,000,000-1,300,000), four vignettes (getting started, population structure, diversity and selection, reports), README with figures (`data-raw/readme_figures.R`), `_pkgdown.yml` (site built into `docs/`), GitHub Actions R-CMD-check on Windows/macOS/Linux (runs the live-fetch test, so it will test Windows fetching). Moved this file to `dev/` so pkgdown does not publish it. Full R CMD check incl. vignettes clean.
- Next: create GitHub repo, first commit and push; confirm Windows CI; enable GitHub Pages from `docs/`. Later items below.
