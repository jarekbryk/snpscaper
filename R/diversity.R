#' Measure genetic differentiation between populations (Fst)
#'
#' Calculates Fst for every pair of populations (or super-populations).
#'
#' @details
#' ## What Fst measures
#'
#' Fst (the "fixation index") measures how much of the genetic variation is
#' due to differences *between* populations rather than variation *within*
#' them. It is 0 when two populations have identical allele frequencies, and
#' 1 when they share no alleles at all (each fixed for a different one). For
#' human populations, genome-wide values are small: roughly 0.01 between
#' neighbouring populations and 0.1 to 0.15 between continents. Most human
#' genetic variation is found within populations.
#'
#' Individual genes can deviate strongly from the genome-wide average. A gene
#' under strong selection in some populations (e.g. SLC24A5 or ACKR1) can show
#' much higher Fst than a typical region.
#'
#' ## How it is calculated
#'
#' We use Hudson's estimator, as recommended by Bhatia et al. (2013,
#' *Genome Research* 23: 1514-1521). For each SNP it compares the allele
#' frequency difference between the two populations with the diversity
#' expected if they were one population, correcting for sample size. The
#' values for all SNPs are then combined as a ratio of sums (the "ratio of
#' averages"), which is more reliable than averaging Fst values of single
#' SNPs. Small negative values can occur when populations are nearly
#' identical; they mean "no detectable difference".
#'
#' @param x A `snpscape` object.
#' @param by Compare populations (`"pop"`) or super-populations
#'   (`"super_pop"`).
#'
#' @return A tibble with one row per pair of groups: `group1`, `group2` and
#'   `fst`.
#' @export
#'
#' @examples
#' example_lct |>
#'   filter_snps() |>
#'   calc_fst(by = "super_pop")
calc_fst <- function(x, by = c("pop", "super_pop")) {
  check_snpscape(x)
  by <- match.arg(by)
  counts <- group_allele_counts(x, by)
  groups <- rownames(counts$alt_freq)
  if (length(groups) < 2) {
    cli::cli_abort("At least two groups are needed to calculate Fst.")
  }

  pairs <- utils::combn(groups, 2)
  fst <- apply(pairs, 2, function(pair) {
    hudson_fst(
      p1 = counts$alt_freq[pair[1], ],
      p2 = counts$alt_freq[pair[2], ],
      n1 = counts$chromosomes[pair[1], ],
      n2 = counts$chromosomes[pair[2], ]
    )
  })

  tibble::tibble(group1 = pairs[1, ], group2 = pairs[2, ], fst = fst)
}

# Hudson's Fst estimator (Bhatia et al. 2013, equation 10), combined across
# SNPs as a ratio of averages.
#   p1, p2: alternative allele frequency in each population, per SNP
#   n1, n2: number of chromosomes sampled (2 x people), per SNP
hudson_fst <- function(p1, p2, n1, n2) {
  numerator <- (p1 - p2)^2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
  denominator <- p1 * (1 - p2) + p2 * (1 - p1)
  sum(numerator, na.rm = TRUE) / sum(denominator, na.rm = TRUE)
}

#' Plot Fst between pairs of populations
#'
#' A heatmap of the result of [calc_fst()]: darker squares mean populations
#' are more similar, brighter squares mean they are more differentiated.
#' Values are printed in the squares when there are few groups.
#'
#' @param fst The result of [calc_fst()].
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' example_lct |>
#'   filter_snps() |>
#'   keep_pops(pop = c("YRI", "LWK", "GBR", "FIN", "CHB", "JPT")) |>
#'   calc_fst() |>
#'   plot_fst()
plot_fst <- function(fst, title = "Pairwise Fst") {
  levels <- c(names(superpop_colours), names(pop_colours))
  groups <- intersect(levels, c(fst$group1, fst$group2))

  # Show each pair once, in the lower triangle of the square.
  tiles <- fst |>
    dplyr::mutate(
      a = factor(.data$group1, levels = groups),
      b = factor(.data$group2, levels = groups),
      x = pmin(as.integer(.data$a), as.integer(.data$b)),
      y = pmax(as.integer(.data$a), as.integer(.data$b))
    )

  plot <- ggplot2::ggplot(tiles, ggplot2::aes(.data$x, .data$y, fill = .data$fst)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_viridis_c(name = "Fst") +
    ggplot2::scale_x_continuous(breaks = seq_along(groups), labels = groups, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_along(groups), labels = groups, expand = c(0, 0)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    theme_jb() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)
    )

  if (length(groups) <= 10) {
    plot <- plot +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data$fst)), colour = "white", size = 3)
  }
  plot
}

#' Measure genetic diversity within populations (heterozygosity)
#'
#' Calculates two measures of genetic diversity for each population (or
#' super-population), averaged over all SNPs in the data.
#'
#' @details
#' A person is **heterozygous** at a SNP if they carry two different alleles
#' (genotype count 1). Populations with more genetic diversity have more
#' heterozygous people.
#'
#' * **Observed heterozygosity** is the proportion of genotypes that are
#'   heterozygous.
#' * **Expected heterozygosity** (also called gene diversity) is the
#'   proportion we would expect if people paired up alleles at random within
#'   the population: `2p(1 - p)` for allele frequency `p`, with a small
#'   correction for sample size.
#'
#' The two are usually close. Observed heterozygosity clearly below expected
#' can point to inbreeding or to hidden structure within a population.
#'
#' Across the human genome, African populations have the highest diversity,
#' and diversity decreases with distance from Africa. This is thought to
#' reflect a series of population bottlenecks as small groups of people
#' migrated out of Africa and across the world (the "serial founder effect").
#'
#' The values depend on which SNPs are included. They are useful for
#' comparing populations on the same set of SNPs, not as absolute numbers:
#' for example, removing rare SNPs with [filter_snps()] raises the values.
#'
#' @param x A `snpscape` object.
#' @param by Calculate per `"pop"` or per `"super_pop"`.
#'
#' @return A tibble with one row per group: the group (`pop` and
#'   `super_pop`, or only `super_pop`), `people`, `snps`, `observed` and
#'   `expected`.
#' @export
#'
#' @examples
#' example_lct |>
#'   filter_snps() |>
#'   calc_heterozygosity()
calc_heterozygosity <- function(x, by = c("pop", "super_pop")) {
  check_snpscape(x)
  by <- match.arg(by)
  counts <- group_allele_counts(x, by)
  groups <- x$samples[[by]]

  p <- counts$alt_freq
  n <- counts$chromosomes
  expected_per_snp <- n / (n - 1) * 2 * p * (1 - p)

  heterozygous <- rowsum((x$genotypes == 1L) * 1L, groups, na.rm = TRUE)
  known <- rowsum((!is.na(x$genotypes)) * 1L, groups)

  out <- tibble::tibble(
    group = rownames(p),
    people = as.integer(table(groups)[rownames(p)]),
    snps = ncol(x$genotypes),
    observed = unname(rowSums(heterozygous) / rowSums(known)),
    expected = unname(rowMeans(expected_per_snp, na.rm = TRUE))
  )
  names(out)[1] <- by

  if (by == "pop") {
    out <- out |>
      dplyr::left_join(dplyr::distinct(x$samples, .data$pop, .data$super_pop), by = "pop") |>
      dplyr::relocate("super_pop", .after = "pop")
  }
  out
}

#' Plot heterozygosity in each population
#'
#' Shows expected heterozygosity (bars) and observed heterozygosity (points)
#' for each population, coloured by super-population. See
#' [calc_heterozygosity()] for what these measure.
#'
#' @param het The result of [calc_heterozygosity()].
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' example_lct |>
#'   filter_snps() |>
#'   calc_heterozygosity() |>
#'   plot_heterozygosity()
plot_heterozygosity <- function(het, title = "Genetic diversity") {
  group <- if ("pop" %in% names(het)) "pop" else "super_pop"
  levels <- if (group == "pop") names(pop_colours) else names(superpop_colours)
  het[[group]] <- factor(het[[group]], levels = intersect(levels, het[[group]]))

  ggplot2::ggplot(het, ggplot2::aes(.data[[group]], .data$expected, fill = .data$super_pop)) +
    ggplot2::geom_col() +
    ggplot2::geom_point(ggplot2::aes(y = .data$observed), shape = 21, fill = "white", size = 2) +
    scale_fill_superpop() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      title = title,
      subtitle = "Bars: expected heterozygosity. Points: observed heterozygosity.",
      x = NULL,
      y = "Heterozygosity"
    ) +
    theme_jb() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = if (group == "pop") 90 else 0, vjust = 0.5)
    )
}

# Alternative allele frequencies and numbers of sampled chromosomes, as
# matrices with one row per group and one column per SNP. rowsum() adds up
# the rows (people) that belong to the same group.
group_allele_counts <- function(x, by) {
  groups <- x$samples[[by]]
  known <- !is.na(x$genotypes)
  alt_copies <- rowsum(ifelse(known, x$genotypes, 0L), groups)
  chromosomes <- 2 * rowsum(known * 1L, groups)
  list(alt_freq = alt_copies / chromosomes, chromosomes = chromosomes)
}
