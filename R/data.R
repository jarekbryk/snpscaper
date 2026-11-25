#' People in the 1000 Genomes Project
#'
#' The 2,504 unrelated people from 26 populations whose genomes were
#' sequenced by the 1000 Genomes Project. These are the people included in
#' every `snpscape` object.
#'
#' The 26 populations are grouped into five continental
#' **super-populations**: African (AFR), Admixed American (AMR), East Asian
#' (EAS), European (EUR) and South Asian (SAS). "Admixed American" populations
#' have recent ancestry from several continents (Indigenous American,
#' European and African). Some populations were sampled outside their region
#' of ancestral origin, for example Gujarati Indians in Houston (GIH) and Sri
#' Lankan Tamils in the UK (STU).
#'
#' @format A tibble with 2,504 rows and 6 columns:
#' \describe{
#'   \item{sample_id}{Sample identifier, e.g. `"HG00096"`.}
#'   \item{pop}{Three-letter population code, e.g. `"GBR"`.}
#'   \item{super_pop}{Super-population code, e.g. `"EUR"`.}
#'   \item{sex}{`"female"` or `"male"`.}
#'   \item{pop_name}{Full population name.}
#'   \item{super_pop_name}{Full super-population name.}
#' }
#' @source 1000 Genomes Project, International Genome Sample Resource:
#'   <https://www.internationalgenome.org>. Built by
#'   `data-raw/samples_1kg.R`.
#' @examples
#' samples_1kg
#' table(samples_1kg$super_pop)
"samples_1kg"

#' Map locations of the 1000 Genomes populations
#'
#' Approximate coordinates used to draw maps. Populations sampled away from
#' their region of ancestral origin (CEU, GIH, ITU, STU, MXL) are placed at
#' that region, which is more informative on a map of allele frequencies;
#' this is recorded in `placed_at`. Populations with recent ancestry from
#' several continents (ASW, ACB, and the admixed American populations other
#' than MXL) are placed where they were sampled.
#'
#' @format A tibble with 26 rows and 4 columns:
#' \describe{
#'   \item{pop}{Population code.}
#'   \item{latitude, longitude}{Approximate location, in decimal degrees.}
#'   \item{placed_at}{Whether the point marks the sampling location or the
#'     ancestral origin.}
#' }
#' @source Built by `data-raw/pop_locations.R`.
"pop_locations"

#' Genes with well-known patterns of variation between populations
#'
#' A curated set of genes that make good teaching examples, from genes under
#' strong selection in particular populations (e.g. LCT, SLC24A5, ACKR1) to
#' genes with more ordinary patterns of variation (CYP1A2, HTT). Any of them
#' can be passed to [fetch_region()] by name.
#'
#' @format A tibble with one row per gene:
#' \describe{
#'   \item{gene}{Gene symbol.}
#'   \item{chrom, gene_start, gene_end, strand}{Gene position on GRCh38.}
#'   \item{region_start, region_end}{The region fetched by default: the gene
#'     plus 25 kb on each side, which includes nearby regulatory variants.}
#'   \item{description}{A short, neutral description of what the gene does.}
#'   \item{staff_notes}{The expected pattern of variation and its
#'     explanation, for teaching staff.}
#' }
#' @source Gene positions from Ensembl (GRCh38). Built by
#'   `data-raw/teaching_genes.R`.
#' @seealso [teaching_snps] for the key SNPs in each gene.
"teaching_genes"

#' Key SNPs in the teaching genes
#'
#' One or two well-studied SNPs in each of the [teaching_genes], usually the
#' variant responsible for the gene's striking pattern of variation.
#'
#' @format A tibble with one row per SNP:
#' \describe{
#'   \item{gene}{Gene symbol.}
#'   \item{rsid}{rs number (dbSNP identifier).}
#'   \item{chrom, pos}{Position on GRCh38.}
#'   \item{alleles}{Alleles reported by Ensembl, reference allele first.}
#'   \item{note}{What the variant is known for.}
#' }
#' @source SNP positions from Ensembl (GRCh38). Built by
#'   `data-raw/teaching_genes.R`.
"teaching_snps"

#' Example data: the LCT region
#'
#' A `snpscape` object for the lactase gene (LCT) plus 25 kb on each side,
#' fetched with `fetch_region("LCT")`. It is included so that examples and
#' tutorials work without an internet connection.
#'
#' @format A `snpscape` object (see [snpscape]).
#' @source 1000 Genomes Project high-coverage phased panel. Built by
#'   `data-raw/examples.R`.
#' @examples
#' example_lct
"example_lct"

#' Example data: a neutral region of chromosome 18
#'
#' A `snpscape` object for a 300 kb stretch of chromosome 18
#' (chr18:1,000,000-1,300,000) that contains no gene known to be under strong
#' selection. Patterns of variation here are shaped mainly by population
#' history (migration, population size changes and drift), which makes it a
#' useful comparison for genes under selection such as [example_lct].
#' Fetched with `fetch_region("chr18:1000000-1300000")`.
#'
#' @format A `snpscape` object (see [snpscape]).
#' @source 1000 Genomes Project high-coverage phased panel. Built by
#'   `data-raw/examples.R`.
#' @examples
#' example_neutral
"example_neutral"
