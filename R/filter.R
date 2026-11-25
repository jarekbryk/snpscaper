#' Remove rare, poorly genotyped or redundant SNPs
#'
#' Most analyses of population structure work best on a cleaned-up set of
#' SNPs. `filter_snps()` applies up to three filters, in this order, and
#' records how many SNPs remain after each one (see [summary.snpscape()]).
#'
#' @details
#' ## 1. Minor allele frequency (`min_maf`)
#'
#' The minor allele frequency (MAF) is the frequency of the less common of a
#' SNP's two alleles. Most SNPs in the human genome are very rare: they are
#' carried by only one or a handful of people. Rare SNPs say little about
#' differences between populations, add noise, and (when nobody in the data
#' carries the minor allele at all) cannot be analysed by PCA. Removing SNPs
#' with a MAF below 1% (`min_maf = 0.01`) is a common choice.
#'
#' MAF is calculated from the people currently in the data, so if you later
#' remove populations with [keep_pops()], some SNPs may become rare again.
#'
#' ## 2. Missing genotypes (`max_missing`)
#'
#' SNPs for which genotypes are unknown for many people are removed. The 1000
#' Genomes high-coverage data have very few missing genotypes, so this filter
#' rarely removes anything, but it is good practice.
#'
#' ## 3. Linkage disequilibrium pruning (`prune_ld`)
#'
#' Nearby SNPs tend to be inherited together, because chromosomes are passed
#' down in long chunks that are only occasionally broken up by recombination.
#' This correlation between neighbouring SNPs is called linkage
#' disequilibrium (LD). In a small region, many SNPs may be strongly
#' correlated and so carry almost the same information. Methods like PCA and
#' admixture analysis assume each SNP is independent evidence; a block of
#' correlated SNPs effectively gets counted many times and can dominate the
#' result.
#'
#' LD pruning keeps a set of SNPs that are only weakly correlated with each
#' other. It moves a window of `ld_window` SNPs along the region and, within
#' each window, whenever two SNPs have a squared correlation (r²) above
#' `ld_r2`, removes the one with the lower minor allele frequency. This is a
#' simplified version of the widely used `--indep-pairwise` method of the
#' PLINK program.
#'
#' Pruning is off by default. Turn it on (`prune_ld = TRUE`) for PCA and
#' admixture analysis. Leave it off when you are interested in LD itself
#' or in particular SNPs.
#'
#' @param x A `snpscape` object.
#' @param min_maf Minimum minor allele frequency, between 0 and 0.5.
#' @param max_missing Maximum proportion of people with a missing genotype.
#' @param prune_ld Remove SNPs in strong linkage disequilibrium with a
#'   neighbour (`TRUE` or `FALSE`).
#' @param ld_r2 Maximum r² allowed between SNPs within a window when
#'   `prune_ld = TRUE`.
#' @param ld_window Window size, in number of SNPs, for LD pruning.
#'
#' @return The filtered `snpscape` object.
#' @export
#'
#' @examples
#' lct <- example_lct |>
#'   filter_snps(min_maf = 0.05, prune_ld = TRUE)
#' summary(lct)
filter_snps <- function(
  x,
  min_maf = 0.01,
  max_missing = 0.05,
  prune_ld = FALSE,
  ld_r2 = 0.2,
  ld_window = 50
) {
  check_snpscape(x)

  # --- 1. Minor allele frequency ---------------------------------------------
  maf <- minor_allele_freq(x$genotypes)
  x <- subset_snpscape(x, snps = maf >= min_maf, step = paste0("MAF >= ", min_maf))

  # --- 2. Missing genotypes --------------------------------------------------
  missing <- colMeans(is.na(x$genotypes))
  x <- subset_snpscape(
    x,
    snps = missing <= max_missing,
    step = paste0("missing <= ", max_missing)
  )

  # --- 3. LD pruning ------------------------------------------------------------
  if (prune_ld) {
    keep <- prune_by_ld(x$genotypes, x$variants, r2 = ld_r2, window = ld_window)
    x <- subset_snpscape(x, snps = keep, step = paste0("LD pruned (r2 <= ", ld_r2, ")"))
  }

  if (ncol(x$genotypes) == 0) {
    cli::cli_warn("No SNPs are left after filtering. Try a lower {.arg min_maf}.")
  }
  x
}

#' Keep only some populations or people
#'
#' Narrows the data down to the populations you want to compare. You can
#' choose by population (e.g. `"GBR"`, `"YRI"`), by super-population (e.g.
#' `"EUR"`, `"AFR"`), by sex, or any combination; people must match all the
#' choices you make. See [samples_1kg] for the population codes.
#'
#' @param x A `snpscape` object.
#' @param pop Population codes to keep, or `NULL` to keep all.
#' @param super_pop Super-population codes to keep, or `NULL` to keep all.
#' @param sex `"female"` and/or `"male"`, or `NULL` to keep both.
#'
#' @return The `snpscape` object with only the chosen people.
#' @export
#'
#' @examples
#' example_lct |>
#'   keep_pops(pop = c("YRI", "GBR", "CHB"))
#'
#' example_lct |>
#'   keep_pops(super_pop = c("AFR", "EUR"))
keep_pops <- function(x, pop = NULL, super_pop = NULL, sex = NULL) {
  check_snpscape(x)
  check_codes(pop, x$samples$pop, "pop")
  check_codes(super_pop, x$samples$super_pop, "super_pop")
  check_codes(sex, c("female", "male"), "sex")

  keep <- rep(TRUE, nrow(x$samples))
  if (!is.null(pop)) keep <- keep & x$samples$pop %in% pop
  if (!is.null(super_pop)) keep <- keep & x$samples$super_pop %in% super_pop
  if (!is.null(sex)) keep <- keep & x$samples$sex %in% sex

  if (!any(keep)) {
    cli::cli_abort("No people match all of these choices.")
  }

  chosen <- c(pop, super_pop, sex)
  subset_snpscape(x, people = keep, step = paste0("kept ", paste(chosen, collapse = ", ")))
}

# --- Helpers ------------------------------------------------------------------------

# Frequency of the less common allele at each SNP (column).
minor_allele_freq <- function(genotypes) {
  alt_freq <- colMeans(genotypes, na.rm = TRUE) / 2
  pmin(alt_freq, 1 - alt_freq)
}

check_snpscape <- function(x) {
  if (!inherits(x, "snpscape")) {
    cli::cli_abort(c(
      "Expected a {.cls snpscape} object, but got {.cls {class(x)}}.",
      "i" = "Create one with {.fn fetch_region}."
    ))
  }
}

check_codes <- function(chosen, available, arg) {
  unknown <- setdiff(chosen, available)
  if (length(unknown) > 0) {
    cli::cli_abort(c(
      "Unknown {.arg {arg}} code{?s}: {.val {unknown}}.",
      "i" = "Available: {.val {sort(unique(available))}}."
    ))
  }
}
