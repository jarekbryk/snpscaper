#' Estimate ancestry proportions (admixture analysis)
#'
#' Runs an admixture analysis: it assumes that the people in the data descend
#' from `K` ancestral groups, and estimates what proportion of each person's
#' genome comes from each group. The calculation is done by the `snmf()`
#' function of the LEA package.
#'
#' @details
#' ## What admixture analysis does
#'
#' Imagine there were `K` ancestral populations, each with its own allele
#' frequencies at every SNP. Admixture analysis looks for the set of
#' ancestral frequencies, and the mixture of ancestral groups in each person,
#' that best explains the genotypes we see. The result for each person is `K`
#' proportions that add up to 1, usually drawn as a stacked bar (see
#' [plot_admixture()]).
#'
#' The ancestral groups are statistical constructs, not real historical
#' populations, and the method finds `K` groups whether or not they exist.
#' Results also depend on which populations are included and how many people
#' each contributes. Treat them as a description of the structure in these
#' data, not a literal history.
#'
#' ## Choosing K
#'
#' We usually try a range of `K` values. For each, snmf leaves out a small
#' fraction of genotypes, then measures how well the model predicts them (the
#' "cross-entropy"; lower is better). A good `K` is one where the
#' cross-entropy stops dropping noticeably. See [plot_cross_entropy()].
#'
#' ## Repeated runs
#'
#' The method starts from random guesses and can get stuck in different
#' solutions, so each `K` is run `reps` times and the run with the lowest
#' cross-entropy is kept. Setting `seed` makes results reproducible.
#'
#' ## Before you run it
#'
#' Use [filter_snps()] with `prune_ld = TRUE` first: admixture analysis
#' assumes SNPs are independent. Larger numbers of SNPs, people and `K`
#' values take longer; a few hundred SNPs for all 2,504 people take seconds.
#'
#' The ancestral groups are numbered arbitrarily by snmf. To make plots for
#' different `K` easier to compare, the groups are renumbered so that each
#' matches, as closely as possible, a group at the previous `K`.
#'
#' @param x A `snpscape` object.
#' @param k The numbers of ancestral groups to try, e.g. `2:6`.
#' @param reps Number of repeated runs for each `K`.
#' @param seed Random seed, for reproducible results.
#' @param verbose Print progress messages.
#'
#' @return An `snpscape_admixture` object: a list with
#'   * `ancestry`: a tibble with one row per person, `K` and group, giving the
#'     estimated `proportion` of ancestry (best run for each `K` only). Use
#'     [get_ancestry()] to extract it.
#'   * `cross_entropy`: a tibble with the cross-entropy of every run.
#'   * `samples`: the people analysed.
#' @export
#'
#' @examples
#' \donttest{
#' lct <- filter_snps(example_lct, min_maf = 0.05, prune_ld = TRUE)
#' adm <- run_admixture(lct, k = 2:4, reps = 3)
#' plot_cross_entropy(adm)
#' plot_admixture(adm, k = 3)
#' }
run_admixture <- function(x, k = 2:6, reps = 5, seed = 1, verbose = TRUE) {
  check_snpscape(x)
  if (any(k < 2) || any(k >= nrow(x$samples))) {
    cli::cli_abort("{.arg k} must be at least 2 and smaller than the number of people.")
  }
  genotypes <- x$genotypes[, apply(x$genotypes, 2, stats::var, na.rm = TRUE) > 0, drop = FALSE]
  if (ncol(genotypes) < 10) {
    cli::cli_abort(c(
      "Only {ncol(genotypes)} SNP{?s} vary among these people; admixture needs more.",
      "i" = "Fetch a larger region, or lower {.arg min_maf} in {.fn filter_snps}."
    ))
  }

  # --- Run snmf in a temporary folder -----------------------------------------
  # LEA works with files: it reads genotypes from a ".geno" file and writes
  # its results next to it. We do all this in a temporary folder that is
  # deleted afterwards, so nothing is left in the user's working folder.
  work_dir <- tempfile("snpscaper_admixture_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  geno_file <- file.path(work_dir, "genotypes.geno")
  genotypes[is.na(genotypes)] <- 9L # LEA's code for a missing genotype
  LEA::write.geno(genotypes, geno_file)

  if (verbose) {
    cli::cli_inform(
      "Running admixture for K = {k} ({reps} run{?s} each, {comma(ncol(genotypes))} SNPs)..."
    )
  }
  # snmf prints a lot of technical output; capture it to keep the console tidy.
  utils::capture.output(
    project <- LEA::snmf(
      geno_file,
      K = k,
      repetitions = reps,
      entropy = TRUE,
      project = "new",
      seed = seed
    )
  )

  # --- Collect the results ----------------------------------------------------
  cross_entropy <- purrr::map(k, function(this_k) {
    tibble::tibble(
      k = this_k,
      run = seq_len(reps),
      cross_entropy = as.numeric(LEA::cross.entropy(project, K = this_k))
    )
  }) |>
    dplyr::bind_rows()

  # For each K, keep the run with the lowest cross-entropy.
  best_q <- purrr::map(k, function(this_k) {
    runs <- cross_entropy[cross_entropy$k == this_k, ]
    best_run <- runs$run[which.min(runs$cross_entropy)]
    LEA::Q(project, K = this_k, run = best_run)
  })
  best_q <- match_groups_across_k(best_q)

  ancestry <- purrr::map2(best_q, k, function(q, this_k) {
    tibble::tibble(
      sample_id = rep(x$samples$sample_id, times = this_k),
      pop = rep(x$samples$pop, times = this_k),
      super_pop = rep(x$samples$super_pop, times = this_k),
      k = this_k,
      group = rep(seq_len(this_k), each = nrow(q)),
      proportion = as.numeric(q)
    )
  }) |>
    dplyr::bind_rows()

  structure(
    list(ancestry = ancestry, cross_entropy = cross_entropy, samples = x$samples),
    class = "snpscape_admixture"
  )
}

# snmf numbers ancestral groups arbitrarily, so "group 1" at K = 3 may be
# unrelated to "group 1" at K = 4. Going up through the K values, reorder the
# columns of each ancestry matrix so that each group lines up with the most
# similar group (highest correlation across people) at the previous K. Any
# new group goes last.
match_groups_across_k <- function(q_list) {
  if (length(q_list) < 2) {
    return(q_list)
  }
  for (i in 2:length(q_list)) {
    previous <- q_list[[i - 1]]
    current <- q_list[[i]]
    similarity <- suppressWarnings(stats::cor(previous, current))
    similarity[is.na(similarity)] <- -Inf

    order <- integer(0)
    for (group in seq_len(ncol(previous))) {
      candidates <- setdiff(seq_len(ncol(current)), order)
      order <- c(order, candidates[which.max(similarity[group, candidates])])
    }
    order <- c(order, setdiff(seq_len(ncol(current)), order))
    q_list[[i]] <- current[, order, drop = FALSE]
  }
  q_list
}

#' @export
print.snpscape_admixture <- function(x, ...) {
  k <- unique(x$cross_entropy$k)
  cat(
    "<snpscape_admixture> ", comma(nrow(x$samples)), " people, K = ",
    paste(range(k), collapse = " to "), ", ",
    max(x$cross_entropy$run), " runs each\n",
    sep = ""
  )
  invisible(x)
}

#' Get ancestry proportions from an admixture analysis
#'
#' Extracts the estimated ancestry proportions for one or more values of `K`
#' from the result of [run_admixture()], as a tidy table.
#'
#' @param adm The result of [run_admixture()].
#' @param k The `K` values to extract, or `NULL` for all.
#' @param wide If `TRUE`, return one row per person with one column per
#'   ancestral group (only for a single `K`).
#'
#' @return A tibble. In the default long format: `sample_id`, `pop`,
#'   `super_pop`, `k`, `group` and `proportion`.
#' @export
#'
#' @examples
#' \donttest{
#' lct <- filter_snps(example_lct, min_maf = 0.05, prune_ld = TRUE)
#' adm <- run_admixture(lct, k = 2:3, reps = 2)
#' get_ancestry(adm, k = 3, wide = TRUE)
#' }
get_ancestry <- function(adm, k = NULL, wide = FALSE) {
  check_admixture(adm, k)
  ancestry <- adm$ancestry
  if (!is.null(k)) {
    ancestry <- ancestry[ancestry$k %in% k, ]
  }
  if (wide) {
    if (length(unique(ancestry$k)) != 1) {
      cli::cli_abort("{.arg wide = TRUE} needs a single value of {.arg k}.")
    }
    ancestry <- tidyr::pivot_wider(
      ancestry,
      names_from = "group",
      names_prefix = "group_",
      values_from = "proportion"
    )
  }
  ancestry
}

#' Plot ancestry proportions
#'
#' Draws the classic admixture bar plot: each person is a thin vertical bar,
#' divided into coloured segments showing the estimated proportion of their
#' genome from each ancestral group. People are grouped by population or
#' super-population, and sorted within each group by their main ancestry.
#'
#' Colours represent the ancestral groups found by the analysis, not
#' populations, so they are deliberately different from the population
#' colours used in other plots. Several `K` values can be drawn at once, one
#' row per `K`; groups are matched across `K` values (see [run_admixture()]),
#' so the same colour usually means a similar group.
#'
#' @param adm The result of [run_admixture()].
#' @param k One or more `K` values to plot.
#' @param group Group people by `"super_pop"` or `"pop"`.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \donttest{
#' lct <- filter_snps(example_lct, min_maf = 0.05, prune_ld = TRUE)
#' adm <- run_admixture(lct, k = 2:4, reps = 2)
#' plot_admixture(adm, k = 3, group = "pop")
#' plot_admixture(adm, k = 2:4)
#' }
plot_admixture <- function(adm, k, group = c("super_pop", "pop"), title = NULL) {
  group <- match.arg(group)
  check_admixture(adm, k)
  if (max(k) > length(ancestry_colours)) {
    cli::cli_abort("At most {length(ancestry_colours)} ancestral groups can be plotted.")
  }
  ancestry <- get_ancestry(adm, k = k)

  # Sort people within each population by their main ancestral group at the
  # smallest K plotted, so that similar people sit next to each other.
  first_k <- ancestry[ancestry$k == min(k), ]
  main <- first_k |>
    dplyr::group_by(.data$sample_id) |>
    dplyr::slice_max(.data$proportion, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  person_order <- main$sample_id[order(main$group, -main$proportion)]

  group_levels <- if (group == "pop") names(pop_colours) else names(superpop_colours)
  k_levels <- paste0("K = ", sort(k))

  ancestry <- ancestry |>
    dplyr::mutate(
      sample_id = factor(.data$sample_id, levels = person_order),
      facet = factor(.data[[group]], levels = intersect(group_levels, .data[[group]])),
      group = factor(.data$group),
      k = factor(paste0("K = ", .data$k), levels = k_levels)
    )

  ggplot2::ggplot(ancestry, ggplot2::aes(.data$sample_id, .data$proportion, fill = .data$group)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$k),
      cols = ggplot2::vars(.data$facet),
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::scale_fill_manual(values = ancestry_colours, guide = "none") +
    ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    ggplot2::labs(title = title, x = NULL, y = "Ancestry proportion") +
    theme_jb() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.spacing.x = ggplot2::unit(1, "pt"),
      strip.text.x = ggplot2::element_text(angle = if (group == "pop") 90 else 0)
    )
}

#' Plot cross-entropy to help choose K
#'
#' Shows the cross-entropy of every admixture run for each `K` (grey points)
#' and the best run for each `K` (line). Lower cross-entropy means the model
#' predicts left-out genotypes better. A sensible choice of `K` is where the
#' line stops falling steeply or reaches its lowest point. In small genomic
#' regions the curve is often flat or noisy, which itself tells you that the
#' region carries little information about population structure.
#'
#' @param adm The result of [run_admixture()].
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \donttest{
#' lct <- filter_snps(example_lct, min_maf = 0.05, prune_ld = TRUE)
#' adm <- run_admixture(lct, k = 2:5, reps = 3)
#' plot_cross_entropy(adm)
#' }
plot_cross_entropy <- function(adm) {
  check_admixture(adm)
  best <- adm$cross_entropy |>
    dplyr::group_by(.data$k) |>
    dplyr::summarise(cross_entropy = min(.data$cross_entropy))

  ggplot2::ggplot(adm$cross_entropy, ggplot2::aes(.data$k, .data$cross_entropy)) +
    ggplot2::geom_point(colour = "grey60") +
    ggplot2::geom_line(data = best, colour = "black") +
    ggplot2::geom_point(data = best, colour = "black", size = 2) +
    ggplot2::scale_x_continuous(breaks = unique(adm$cross_entropy$k)) +
    ggplot2::labs(
      title = "Choosing the number of ancestral groups",
      x = "K (number of ancestral groups)",
      y = "Cross-entropy (lower is better)"
    ) +
    theme_jb()
}

# Colours for ancestral groups (Tableau 10 palette).
ancestry_colours <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
)

check_admixture <- function(adm, k = NULL) {
  if (!inherits(adm, "snpscape_admixture")) {
    cli::cli_abort("Expected the result of {.fn run_admixture}.")
  }
  available <- unique(adm$cross_entropy$k)
  if (!is.null(k) && !all(k %in% available)) {
    cli::cli_abort("{.arg k} must be one of the values analysed: {available}.")
  }
}
