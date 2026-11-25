#' Create an HTML report for a region of the genome
#'
#' Produces a handout-style HTML report for a gene or region: a summary of the
#' data, PCA, admixture, Fst, heterozygosity, allele frequencies of key SNPs
#' (with maps) and LD, plus a methods section recording exactly how the
#' results were produced. The report is a single self-contained HTML file that
#' can be emailed, uploaded to a virtual learning environment, or printed.
#'
#' @details
#' ## Staff and student versions
#'
#' The staff version (the default) adds a "Notes for staff" section with the
#' expected pattern and its explanation, for genes in [teaching_genes]. The
#' student version leaves this out, so students can interpret the results
#' themselves.
#'
#' ## Which data are used for what
#'
#' The report applies [filter_snps()] itself, using `min_maf`:
#'
#' * PCA and admixture use LD-pruned SNPs, because these methods assume SNPs
#'   are independent.
#' * Allele frequencies, Fst and heterozygosity use all SNPs that pass the
#'   frequency filter.
#' * LD heatmaps use SNPs with a minor allele frequency of at least 5%, one
#'   heatmap per super-population, without pruning.
#'
#' ## Key SNPs
#'
#' The allele frequency section shows the SNPs given in `focal_snps`. If none
#' are given, it uses the SNPs listed for the gene in [teaching_snps], or,
#' for other regions, the SNP whose frequency differs most between
#' super-populations.
#'
#' ## Requirements and run time
#'
#' Reports are made with Quarto, which is installed together with RStudio
#' (or from <https://quarto.org>). Making a report takes from under a minute
#' to a few minutes, mostly for the admixture analysis; use fewer `k` values
#' or `reps` to speed it up.
#'
#' @param x What to report on. Either a character vector of gene symbols or
#'   coordinates (one report is made for each; see [fetch_region()]), or a
#'   `snpscape` object you have prepared yourself (e.g. combining several
#'   regions).
#' @param pop,super_pop Populations or super-populations to include (see
#'   [keep_pops()]). `NULL` includes everyone.
#' @param focal_snps rs numbers or SNP IDs to show allele frequencies for.
#' @param k The numbers of ancestral groups for admixture analysis.
#' @param reps Number of admixture runs for each `K`.
#' @param min_maf Minimum minor allele frequency for SNPs (see
#'   [filter_snps()]).
#' @param seed Random seed, so that the same settings always give the same
#'   report.
#' @param version `"staff"` (with notes for staff) or `"student"` (without).
#' @param title Report title. By default, made from the region name.
#' @param output_dir Folder to save reports in. Created if needed.
#' @param verbose Print progress messages.
#'
#' @return The paths of the created HTML files, invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' render_report("LCT")
#'
#' # one report per gene, for a chosen set of populations
#' render_report(
#'   c("LCT", "SLC24A5", "ACKR1"),
#'   pop = c("YRI", "LWK", "GBR", "FIN", "CHB", "JPT", "GIH", "PEL"),
#'   output_dir = "handouts"
#' )
#'
#' # student version, for a region with no well-known genes
#' render_report("chr18:1000000-1500000", version = "student")
#' }
render_report <- function(
  x,
  pop = NULL,
  super_pop = NULL,
  focal_snps = NULL,
  k = 2:5,
  reps = 5,
  min_maf = 0.01,
  seed = 1,
  version = c("staff", "student"),
  title = NULL,
  output_dir = ".",
  verbose = TRUE
) {
  version <- match.arg(version)
  check_report_tools()

  # A character vector gives one report per element.
  if (is.character(x)) {
    paths <- vapply(
      x,
      function(region) {
        render_report(
          fetch_region(region, verbose = verbose),
          pop = pop, super_pop = super_pop, focal_snps = focal_snps,
          k = k, reps = reps, min_maf = min_maf, seed = seed,
          version = version, title = title, output_dir = output_dir,
          verbose = verbose
        )
      },
      character(1)
    )
    return(invisible(unname(paths)))
  }
  check_snpscape(x)

  # --- Prepare the data ----------------------------------------------------
  if (!is.null(pop) || !is.null(super_pop)) {
    x <- keep_pops(x, pop = pop, super_pop = super_pop)
  }
  if (length(focal_snps) > 0) {
    x <- add_rsids(x, focal_snps)
  }
  label <- paste(x$provenance$regions$label, collapse = " + ")
  title <- title %||% paste("Genetic variation in", label)

  # --- Render the Quarto template in a temporary folder --------------------
  # The template reads the prepared data from a file, because Quarto
  # parameters can only be simple values (text and numbers).
  work_dir <- tempfile("snpscaper_report_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  template_dir <- system.file("quarto", package = "snpscaper")
  file.copy(list.files(template_dir, full.names = TRUE), work_dir)
  data_file <- file.path(work_dir, "data.rds")
  saveRDS(x, data_file)

  if (verbose) {
    cli::cli_inform("Making the {version} report for {label} (this can take a few minutes)...")
  }
  quarto::quarto_render(
    input = file.path(work_dir, "report.qmd"),
    execute_params = list(
      data = data_file,
      focal_snps = paste(focal_snps, collapse = ","),
      k = paste(k, collapse = ","),
      reps = reps,
      min_maf = min_maf,
      seed = seed,
      version = version
    ),
    metadata = list(title = title),
    quiet = TRUE
  )

  # --- Copy the finished report to the output folder -----------------------
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  file_label <- gsub("[^A-Za-z0-9]+", "_", label)
  output_file <- file.path(output_dir, paste0("snpscaper_", file_label, "_", version, ".html"))
  file.copy(file.path(work_dir, "report.html"), output_file, overwrite = TRUE)

  if (verbose) {
    cli::cli_inform("Saved {.file {output_file}}.")
  }
  invisible(output_file)
}

# Check that everything needed to make a report is installed.
check_report_tools <- function() {
  rlang::check_installed(
    c("quarto", "knitr", "rmarkdown", "FactoMineR", "factoextra"),
    reason = "to make reports."
  )
  if (is.null(quarto::quarto_path())) {
    cli::cli_abort(c(
      "The Quarto program is needed to make reports, but was not found.",
      "i" = "Quarto comes with recent versions of RStudio, or can be installed from {.url https://quarto.org}."
    ))
  }
}

# Look up the positions of rs numbers that are not already labelled in the
# data, and label the matching SNPs. SNP IDs (chr:pos:ref:alt) are left as
# they are.
add_rsids <- function(x, snps) {
  rs_numbers <- setdiff(snps[grepl("^rs[0-9]+$", snps)], x$variants$rsid)
  if (length(rs_numbers) == 0) {
    return(x)
  }
  found <- find_snp(rs_numbers)
  match_row <- match(
    paste(found$chrom, found$pos),
    paste(x$variants$chrom, x$variants$pos)
  )
  x$variants$rsid[match_row[!is.na(match_row)]] <- found$rsid[!is.na(match_row)]
  x
}
