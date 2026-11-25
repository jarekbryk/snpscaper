#' The snpscape object
#'
#' A `snpscape` object holds everything snpscaper knows about a set of
#' genotypes. It is created by [fetch_region()] and passed along to the other
#' functions. Printing it gives a short summary.
#'
#' It is a list with four parts:
#'
#' * `genotypes`: a matrix with one row per person and one column per SNP.
#'   Each value is the number of copies of the alternative allele that the
#'   person carries (0, 1 or 2), or `NA` if unknown.
#' * `variants`: a tibble describing each SNP (one row per column of
#'   `genotypes`): `id`, `chrom`, `pos` (position on GRCh38), `ref`
#'   (reference allele), `alt` (alternative allele), `region` (which requested
#'   region it came from) and `rsid` (rs number, for well-known SNPs only).
#' * `samples`: a tibble describing each person (one row per row of
#'   `genotypes`); see [samples_1kg].
#' * `provenance`: where the data came from and what has been done to them:
#'   the regions fetched, genome build, data source, date, package version and
#'   a log of filtering steps.
#'
#' You rarely need to touch these parts directly, but you can, e.g.
#' `x$samples` or `x$variants`.
#'
#' @name snpscape
NULL

# Create a snpscape object. Used by fetch_region().
new_snpscape <- function(genotypes, variants, samples, regions) {
  stopifnot(
    identical(rownames(genotypes), samples$sample_id),
    identical(colnames(genotypes), variants$id)
  )
  structure(
    list(
      genotypes = genotypes,
      variants = variants,
      samples = samples,
      provenance = list(
        regions = regions,
        build = "GRCh38",
        source = "1000 Genomes Project high-coverage phased panel (NYGC, 2022)",
        fetched = Sys.Date(),
        package_version = as.character(utils::packageVersion("snpscaper")),
        steps = tibble::tibble(
          step = "fetched",
          people = nrow(samples),
          snps = nrow(variants)
        )
      )
    ),
    class = "snpscape"
  )
}

# Keep some people (rows) and/or SNPs (columns), keeping all parts in step,
# and record what was done in the provenance log.
subset_snpscape <- function(x, people = TRUE, snps = TRUE, step) {
  x$genotypes <- x$genotypes[people, snps, drop = FALSE]
  x$samples <- x$samples[people, ]
  x$variants <- x$variants[snps, ]
  x$provenance$steps <- dplyr::bind_rows(
    x$provenance$steps,
    tibble::tibble(step = step, people = nrow(x$samples), snps = nrow(x$variants))
  )
  x
}

#' @export
print.snpscape <- function(x, ...) {
  regions <- x$provenance$regions
  n_pops <- length(unique(x$samples$pop))
  n_superpops <- length(unique(x$samples$super_pop))

  cat("<snpscape>", comma(nrow(x$samples)), "people x", comma(nrow(x$variants)), "SNPs (GRCh38)\n")
  cat(
    "Regions:    ",
    paste0(
      regions$label,
      ifelse(
        grepl(":", regions$label),
        "",
        paste0(" (", format_region(regions$chrom, regions$start, regions$end), ")")
      ),
      collapse = ", "
    ),
    "\n"
  )
  cat("Populations:", n_pops, "populations in", n_superpops, "super-populations\n")
  cat("Steps:      ", paste0(x$provenance$steps$step, collapse = " -> "), "\n")
  invisible(x)
}

#' @export
dim.snpscape <- function(x) {
  dim(x$genotypes)
}

#' Summarise a snpscape object
#'
#' Counts how many people from each population are in the data, and shows
#' each processing step with the number of people and SNPs remaining.
#'
#' @param object A `snpscape` object.
#' @param ... Not used.
#'
#' @return A list with `people` (a tibble of counts per population) and
#'   `steps` (the processing log), invisibly. Both are also printed.
#' @export
summary.snpscape <- function(object, ...) {
  people <- object$samples |>
    dplyr::count(.data$super_pop, .data$pop, .data$pop_name, name = "people")
  steps <- object$provenance$steps

  print(object)
  cat("\nPeople per population:\n")
  print(people, n = Inf)
  cat("\nProcessing steps:\n")
  print(steps, n = Inf)

  invisible(list(people = people, steps = steps))
}
