# Linkage disequilibrium (LD) pruning. See ?filter_snps for the explanation.
#
# The window moves along the region `step` SNPs at a time. Within each window
# we calculate r² (the squared correlation of genotype counts) between every
# pair of SNPs still kept, and while any pair has r² above the threshold we
# remove the SNP with the lower minor allele frequency from that pair.
#
# Returns a logical vector: TRUE for SNPs to keep.
prune_by_ld <- function(genotypes, variants, r2 = 0.2, window = 50, step = 5) {
  n_snps <- ncol(genotypes)
  keep <- rep(TRUE, n_snps)
  if (n_snps < 2) {
    return(keep)
  }

  # Work through each chromosome separately, with SNPs in position order.
  # (SNPs on different chromosomes are never in LD in this sense.)
  maf <- minor_allele_freq(genotypes)
  genotypes <- fill_missing_with_mean(genotypes)

  for (chrom in unique(variants$chrom)) {
    on_chrom <- which(variants$chrom == chrom)
    on_chrom <- on_chrom[order(variants$pos[on_chrom])]

    for (first in seq(1, length(on_chrom), by = step)) {
      in_window <- on_chrom[first:min(first + window - 1, length(on_chrom))]
      in_window <- in_window[keep[in_window]]
      if (length(in_window) < 2) next

      r2_matrix <- suppressWarnings(stats::cor(genotypes[, in_window]))^2
      r2_matrix[is.na(r2_matrix)] <- 0
      r2_matrix[lower.tri(r2_matrix, diag = TRUE)] <- 0 # look at each pair once

      # Remove SNPs one at a time until no pair in the window is above r2.
      repeat {
        pairs <- which(r2_matrix > r2, arr.ind = TRUE)
        if (nrow(pairs) == 0) break
        a <- pairs[1, 1]
        b <- pairs[1, 2]
        drop <- if (maf[in_window[a]] < maf[in_window[b]]) a else b
        keep[in_window[drop]] <- FALSE
        r2_matrix[drop, ] <- 0
        r2_matrix[, drop] <- 0
      }
    }
  }
  keep
}

# Replace missing genotypes with the average for that SNP, so correlations
# can be calculated quickly. (The 1000 Genomes data have very few missing
# genotypes, so this has a negligible effect on the results.)
fill_missing_with_mean <- function(genotypes) {
  if (!anyNA(genotypes)) {
    return(genotypes)
  }
  means <- colMeans(genotypes, na.rm = TRUE)
  missing <- which(is.na(genotypes), arr.ind = TRUE)
  genotypes[missing] <- means[missing[, 2]]
  genotypes
}

#' Calculate linkage disequilibrium between SNPs
#'
#' Calculates r², a measure of linkage disequilibrium (LD), between every
#' pair of SNPs in the data.
#'
#' @details
#' ## What LD is
#'
#' Chromosomes are passed from parents to children in long pieces, which are
#' broken up only occasionally by recombination. SNPs that are close together
#' therefore tend to be inherited together, so knowing a person's genotype at
#' one SNP tells you something about their genotype at a nearby SNP. This
#' non-random association is called linkage disequilibrium.
#'
#' r² is the squared correlation between the genotypes at two SNPs across
#' people: 0 means no association, 1 means one SNP perfectly predicts the
#' other. LD typically decays with distance, but the pattern is patchy:
#' "blocks" of high LD are separated by recombination hotspots.
#'
#' LD differs between populations. African populations, which have the
#' longest history and largest long-term population sizes, generally show
#' shorter stretches of LD than populations that went through bottlenecks
#' during migrations out of Africa. A recent selective sweep (as around LCT
#' in Europeans) creates an unusually long stretch of high LD, because the
#' favoured allele rose in frequency so fast that recombination had little
#' time to break up the surrounding chromosome. Compare populations by using
#' [keep_pops()] before `calc_ld()`.
#'
#' ## Practicalities
#'
#' Rare SNPs give unreliable r² values, so filter them out first, e.g. with
#' `filter_snps(min_maf = 0.05)` (and without LD pruning, which would remove
#' exactly the pattern you want to see). The number of pairs grows quickly
#' with the number of SNPs, so if there are more than `max_snps` SNPs, an
#' evenly spaced subset is used.
#'
#' @param x A `snpscape` object.
#' @param max_snps Maximum number of SNPs to use.
#'
#' @return A tibble with one row per pair of SNPs: `snp1`, `snp2` (SNP IDs),
#'   `pos1`, `pos2` (positions), `index1`, `index2` (the SNPs' order along the
#'   chromosome) and `r2`.
#' @export
#'
#' @examples
#' example_lct |>
#'   filter_snps(min_maf = 0.05) |>
#'   keep_pops(super_pop = "EUR") |>
#'   calc_ld()
calc_ld <- function(x, max_snps = 300) {
  check_snpscape(x)
  if (length(unique(x$variants$chrom)) > 1) {
    cli::cli_abort("LD can only be calculated within one chromosome; fetch one region at a time.")
  }

  # SNPs that do not vary among these people (for example, after keeping
  # only some populations) have no defined correlation with other SNPs.
  varies <- apply(x$genotypes, 2, stats::var, na.rm = TRUE) > 0
  if (any(!varies)) {
    cli::cli_inform("Dropped {sum(!varies)} SNP{?s} that do not vary among these people.")
  }
  snps <- which(varies)[order(x$variants$pos[varies])]
  if (length(snps) > max_snps) {
    cli::cli_inform(
      "Using {max_snps} evenly spaced SNPs out of {comma(length(snps))} (see {.arg max_snps})."
    )
    snps <- snps[round(seq(1, length(snps), length.out = max_snps))]
  }
  if (length(snps) < 2) {
    cli::cli_abort("At least two SNPs are needed to calculate LD.")
  }

  genotypes <- fill_missing_with_mean(x$genotypes[, snps, drop = FALSE])
  r2 <- suppressWarnings(stats::cor(genotypes))^2

  pairs <- which(upper.tri(r2), arr.ind = TRUE)
  variants <- x$variants[snps, ]
  tibble::tibble(
    snp1 = variants$id[pairs[, 1]],
    snp2 = variants$id[pairs[, 2]],
    pos1 = variants$pos[pairs[, 1]],
    pos2 = variants$pos[pairs[, 2]],
    index1 = pairs[, 1],
    index2 = pairs[, 2],
    r2 = r2[pairs]
  )
}

#' Plot linkage disequilibrium as a heatmap
#'
#' Draws the result of [calc_ld()] as a triangle: each square is a pair of
#' SNPs, coloured from white (r² = 0, no LD) to red (r² = 1, complete LD).
#' SNPs are arranged in order along the chromosome, so blocks of high LD show
#' up as red triangles along the diagonal.
#'
#' SNPs are drawn at equal spacing, in their order along the chromosome, not
#' at their exact positions; axis labels give the positions in megabases.
#'
#' @param ld The result of [calc_ld()].
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' example_lct |>
#'   filter_snps(min_maf = 0.05) |>
#'   keep_pops(super_pop = "EUR") |>
#'   calc_ld() |>
#'   plot_ld(title = "LD around LCT in Europeans")
plot_ld <- function(ld, title = "Linkage disequilibrium") {
  n_snps <- max(ld$index2)
  positions <- numeric(n_snps)
  positions[ld$index1] <- ld$pos1
  positions[ld$index2] <- ld$pos2
  breaks <- unique(round(seq(1, n_snps, length.out = 6)))
  labels <- sprintf("%.2f", positions[breaks] / 1e6)

  ggplot2::ggplot(ld, ggplot2::aes(.data$index1, .data$index2, fill = .data$r2)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1), name = expression(r^2)) +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = breaks, labels = labels, expand = c(0, 0)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = "Position (Mb)", y = "Position (Mb)") +
    theme_jb() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_text()
    )
}
