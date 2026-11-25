#' Calculate allele frequencies in each population
#'
#' For every SNP, calculates how common the alternative allele is in each
#' population (or super-population).
#'
#' @details
#' Each person carries two copies of every SNP, one from each parent. The
#' alternative allele frequency in a population is the number of copies of
#' the alternative allele carried by its people, divided by the total number
#' of copies (twice the number of people with a known genotype).
#'
#' Which allele is called "alternative" is a convention: the reference allele
#' is simply the one in the human reference genome, which was assembled
#' mostly from a few individuals. It is not necessarily the ancestral or the
#' most common allele. Check the `ref` and `alt` columns.
#'
#' @param x A `snpscape` object.
#' @param by Calculate frequencies per `"pop"` or per `"super_pop"`.
#'
#' @return A tibble with one row per SNP and group: the SNP details (`id`,
#'   `rsid`, `chrom`, `pos`, `ref`, `alt`), the group (`pop` and `super_pop`,
#'   or only `super_pop`), `people` (number of people with a known genotype)
#'   and `alt_freq` (alternative allele frequency, between 0 and 1).
#' @export
#'
#' @examples
#' freqs <- calc_allele_freqs(example_lct)
#' freqs[freqs$rsid %in% "rs4988235", ]
calc_allele_freqs <- function(x, by = c("pop", "super_pop")) {
  check_snpscape(x)
  by <- match.arg(by)
  counts <- group_allele_counts(x, by)
  n_groups <- nrow(counts$alt_freq)

  # Turn the groups x SNPs matrices into one row per group and SNP.
  freqs <- tibble::tibble(
    group = rep(rownames(counts$alt_freq), times = ncol(x$genotypes)),
    snp = rep(seq_len(ncol(x$genotypes)), each = n_groups),
    people = as.integer(counts$chromosomes / 2),
    alt_freq = as.numeric(counts$alt_freq)
  )

  group_info <- if (by == "pop") {
    dplyr::distinct(x$samples, .data$pop, .data$super_pop)
  } else {
    dplyr::distinct(x$samples, .data$super_pop)
  }

  out <- x$variants[freqs$snp, c("id", "rsid", "chrom", "pos", "ref", "alt")]
  out[[by]] <- freqs$group
  out |>
    dplyr::left_join(group_info, by = by) |>
    dplyr::mutate(people = freqs$people, alt_freq = freqs$alt_freq) |>
    dplyr::relocate(dplyr::any_of(c("pop", "super_pop")), .after = "alt")
}

#' Plot the allele frequency of a SNP in each population
#'
#' A bar chart of how common the alternative allele of one SNP is in each
#' population, coloured by super-population. Useful for seeing at a glance
#' whether a variant is shared by everyone, or concentrated in some parts of
#' the world.
#'
#' @param freqs The result of [calc_allele_freqs()].
#' @param snp The SNP to plot: an rs number (e.g. `"rs4988235"`, for SNPs in
#'   [teaching_snps]) or an ID from `x$variants$id` (e.g.
#'   `"chr2:135851076:G:A"`).
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' example_lct |>
#'   calc_allele_freqs() |>
#'   plot_freqs(snp = "rs4988235")
plot_freqs <- function(freqs, snp) {
  one <- pick_snp(freqs, snp)
  group <- if ("pop" %in% names(one)) "pop" else "super_pop"
  levels <- if (group == "pop") names(pop_colours) else names(superpop_colours)
  one[[group]] <- factor(one[[group]], levels = intersect(levels, one[[group]]))

  ggplot2::ggplot(one, ggplot2::aes(.data[[group]], .data$alt_freq, fill = .data$super_pop)) +
    ggplot2::geom_col() +
    scale_fill_superpop() +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(
      title = snp_title(one),
      x = NULL,
      y = paste0("Frequency of allele ", one$alt[1])
    ) +
    theme_jb() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = if (group == "pop") 90 else 0, vjust = 0.5)
    )
}

#' Map the allele frequency of a SNP
#'
#' Shows on a world map how common the alternative allele of one SNP is in
#' each population. Each population is a circle, coloured from dark (allele
#' rare) to yellow (allele common).
#'
#' Several populations were sampled away from the region their ancestry comes
#' from (for example Gujarati Indians in Houston). These are drawn at their
#' region of ancestral origin; see [pop_locations]. Populations of mixed
#' recent ancestry, such as African Americans (ASW), are drawn where they were
#' sampled.
#'
#' @inheritParams plot_freqs
#' @param freqs The result of [calc_allele_freqs()] with `by = "pop"`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' example_lct |>
#'   calc_allele_freqs() |>
#'   map_freqs(snp = "rs4988235")
map_freqs <- function(freqs, snp) {
  if (!"pop" %in% names(freqs)) {
    cli::cli_abort("Maps need frequencies per population: use {.code calc_allele_freqs(x, by = \"pop\")}.")
  }
  one <- pick_snp(freqs, snp) |>
    dplyr::left_join(pop_locations, by = "pop")

  world <- ggplot2::map_data("world")

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(.data$long, .data$lat, group = .data$group),
      fill = "grey90",
      colour = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_point(
      data = one,
      ggplot2::aes(.data$longitude, .data$latitude, fill = .data$alt_freq),
      shape = 21,
      size = 5,
      colour = "grey20"
    ) +
    ggplot2::scale_fill_viridis_c(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      name = paste0("Frequency of ", one$alt[1]),
      guide = ggplot2::guide_colourbar(barwidth = 12, title.vjust = 0.8)
    ) +
    ggplot2::coord_quickmap(xlim = c(-125, 150), ylim = c(-35, 65)) +
    ggplot2::labs(
      title = snp_title(one),
      caption = "Populations sampled outside their region of ancestry are shown at that region.",
      x = NULL,
      y = NULL
    ) +
    theme_jb() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(),
      legend.position = "bottom"
    )
}

# Select the rows for one SNP from an allele frequency table, by rs number
# or by ID.
pick_snp <- function(freqs, snp) {
  if (length(snp) != 1) {
    cli::cli_abort("{.arg snp} must be a single SNP.")
  }
  one <- freqs[freqs$id %in% snp | freqs$rsid %in% snp, ]
  if (nrow(one) == 0) {
    cli::cli_abort(c(
      "SNP {.val {snp}} is not in these data.",
      "i" = "It may be outside the region, or removed by {.fn filter_snps}.",
      "i" = "rs numbers are only known for SNPs in {.code teaching_snps}; for others use the ID, e.g. {.val {freqs$id[1]}}."
    ))
  }
  one
}

snp_title <- function(one) {
  name <- if (!is.na(one$rsid[1])) paste0(one$rsid[1], " (", one$id[1], ")") else one$id[1]
  paste0(name, ", alleles ", one$ref[1], "/", one$alt[1])
}
