#' Prepare genotypes for principal component analysis
#'
#' Turns a `snpscape` object into a data frame that can be passed straight to
#' [FactoMineR::PCA()], with population labels included so that factoextra
#' plots can be coloured by population.
#'
#' @details
#' ## What PCA does with genetic data
#'
#' Each person is described by hundreds or thousands of numbers (their
#' genotype at each SNP). Principal component analysis (PCA) finds the
#' directions in which people differ most and summarises them as a few new
#' variables, the principal components. When people are plotted on the first
#' two components, those with similar genotypes appear close together. Because
#' people from the same population share more of their genetic variation,
#' populations often form clusters, and the arrangement of the clusters
#' reflects human population history.
#'
#' ## The table
#'
#' * Row names are the people's sample IDs.
#' * The first three columns are `pop`, `super_pop` and `sex` (as factors).
#'   Tell `PCA()` they are descriptive, not part of the analysis, with
#'   `quali.sup = 1:3`.
#' * The remaining columns are SNPs (allele counts 0, 1 or 2).
#'
#' SNPs that do not vary among the people in the data are dropped, because
#' they contain no information and would cause errors in PCA when scaling.
#' Missing genotypes are left as `NA`; `PCA()` replaces them with the SNP's
#' average.
#'
#' It is usually best to run [filter_snps()] with `prune_ld = TRUE` first.
#'
#' @param x A `snpscape` object.
#'
#' @return A data frame with one row per person.
#' @export
#'
#' @examples
#' lct <- example_lct |>
#'   filter_snps(min_maf = 0.05, prune_ld = TRUE)
#'
#' pca_table <- as_pca_table(lct)
#' pca_table[1:5, 1:6]
#'
#' if (requireNamespace("FactoMineR", quietly = TRUE) &&
#'     requireNamespace("factoextra", quietly = TRUE)) {
#'   pca <- FactoMineR::PCA(pca_table, quali.sup = 1:3, graph = FALSE)
#'   factoextra::fviz_pca_ind(
#'     pca,
#'     habillage = "super_pop",
#'     palette = superpop_colours,
#'     geom = "point",
#'     title = "LCT region"
#'   ) +
#'     theme_jb()
#' }
as_pca_table <- function(x) {
  check_snpscape(x)

  varies <- apply(x$genotypes, 2, stats::var, na.rm = TRUE) > 0
  if (any(!varies)) {
    cli::cli_inform(
      "Dropped {sum(!varies)} SNP{?s} that do not vary among these people."
    )
  }
  if (sum(varies) < 2) {
    cli::cli_abort("At least two SNPs that vary are needed for PCA.")
  }

  labels <- data.frame(
    pop = factor(x$samples$pop, levels = names(pop_colours)),
    super_pop = factor(x$samples$super_pop, levels = names(superpop_colours)),
    sex = factor(x$samples$sex),
    row.names = x$samples$sample_id
  ) |>
    droplevels()

  snps <- as.data.frame(x$genotypes[, varies, drop = FALSE])
  cbind(labels, snps)
}
