# Fetch genotypes from the 1000 Genomes Project.
#
# The data: the 1000 Genomes "high-coverage phased panel" (GRCh38), produced
# by the New York Genome Center by sequencing every genome about 30 times
# over. It is stored on the EBI FTP server as one VCF file per chromosome.
# Each file is compressed in blocks ("bgzip") and comes with a small index
# (".tbi", made by the tabix program). The index lets us download just the
# blocks that cover the region we want, instead of the whole multi-gigabyte
# file. The Rsamtools and VariantAnnotation packages (from Bioconductor) do
# the reading.

panel_url <- paste0(
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/",
  "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/",
  "1kGP_high_coverage_Illumina.%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)

# Regions are read in windows of this size, and each window is converted to
# compact allele counts before reading the next, which keeps memory use low.
window_size <- 250000

#' Fetch genotypes for a region of the genome
#'
#' Downloads the genotypes of the 2,504 unrelated people in the 1000 Genomes
#' Project for one or more regions of the genome, and returns them as a
#' `snpscape` object ready for analysis.
#'
#' @details
#' ## What you can ask for
#'
#' Each element of `regions` can be:
#'
#' * a **gene symbol**, such as `"LCT"` or `"MC1R"`. The gene is looked up
#'   with [find_gene()] and `flank` base pairs are added on each side, because
#'   variants that regulate a gene often lie just outside it.
#' * a set of **coordinates** on the GRCh38 assembly, such as
#'   `"chr2:135787850-135862184"`. The `chr` prefix and commas are optional,
#'   so `"2:135,787,850-135,862,184"` also works.
#'
#' Only the autosomes (chromosomes 1 to 22) are supported. Large regions take
#' longer and need more memory: a 1 Mb region typically holds around 25,000
#' SNPs.
#'
#' ## What you get
#'
#' Only **biallelic SNPs** are kept: positions where exactly two DNA letters
#' are seen in the population (for example `G` and `A`). Insertions,
#' deletions, larger structural variants, and positions with more than two
#' alleles are dropped, because the analyses in this package assume each
#' variant has exactly two forms.
#'
#' Each person's genotype at a SNP is stored as the number of copies of the
#' alternative allele (`alt`) they carry: 0, 1 or 2. The reference allele
#' (`ref`) is simply the letter found in the human reference genome; it is
#' not necessarily the ancestral or the most common allele.
#'
#' The 698 relatives (mostly children of parent pairs) who are also in the
#' high-coverage panel are not included, because close relatives share long
#' stretches of DNA and would distort population analyses.
#'
#' ## Caching
#'
#' Fetched regions are saved on disk, so fetching the same region again is
#' instant and works offline. See [clear_cache()].
#'
#' @param regions A character vector of gene symbols and/or coordinates.
#' @param flank Number of base pairs to add on each side of a gene. Ignored
#'   for coordinates.
#' @param verbose Print progress messages.
#'
#' @return A `snpscape` object (see [snpscape]).
#' @export
#'
#' @examples
#' \dontrun{
#' lct <- fetch_region("LCT")
#' lct
#'
#' # several regions at once
#' pigmentation <- fetch_region(c("MC1R", "SLC24A5"))
#'
#' # a stretch of chromosome 18 with no well-known genes under selection
#' neutral <- fetch_region("chr18:1000000-1500000")
#' }
fetch_region <- function(regions, flank = 25000, verbose = TRUE) {
  region_table <- resolve_regions(regions, flank)
  check_region_size(region_table)

  pieces <- purrr::pmap(
    region_table,
    function(label, chrom, start, end) {
      if (verbose) {
        where <- format_region(chrom, start, end)
        shown <- if (label == where) label else paste0(label, " (", where, ")")
        cli::cli_inform("Fetching {shown}...")
      }
      piece <- read_region(chrom, start, end)
      piece$variants$region <- rep(label, nrow(piece$variants))
      piece
    }
  )

  genotypes <- do.call(cbind, purrr::map(pieces, "genotypes"))
  variants <- dplyr::bind_rows(purrr::map(pieces, "variants"))

  # Overlapping regions would give the same SNP twice; keep it once.
  keep <- !duplicated(variants$id)

  # The 1000 Genomes files do not include rs numbers. Add them for the
  # well-known SNPs in `teaching_snps`, so they can be referred to by name.
  variants$rsid <- teaching_snps$rsid[match(
    paste(variants$chrom, variants$pos),
    paste(teaching_snps$chrom, teaching_snps$pos)
  )]

  x <- new_snpscape(
    genotypes = genotypes[, keep, drop = FALSE],
    variants = variants[keep, ],
    samples = samples_1kg,
    regions = region_table
  )
  if (verbose) {
    cli::cli_inform("Done: {comma(nrow(x$samples))} people, {comma(nrow(x$variants))} SNPs.")
  }
  x
}

# --- Working out which regions to fetch ---------------------------------------

# Turn the user's input into a tibble with one row per region:
# label, chrom, start, end.
resolve_regions <- function(regions, flank) {
  if (!is.character(regions) || length(regions) == 0) {
    cli::cli_abort(
      "{.arg regions} must be gene symbols or coordinates, e.g. {.val LCT} or {.val chr2:135787850-135862184}."
    )
  }
  # lapply() rather than purrr::map(), so that error messages stay short.
  lapply(regions, resolve_one_region, flank = flank) |>
    dplyr::bind_rows()
}

resolve_one_region <- function(region, flank) {
  coords <- parse_coordinates(region)
  if (!is.null(coords)) {
    return(coords)
  }
  gene <- find_gene(region)
  tibble::tibble(
    label = gene$gene,
    chrom = gene$chrom,
    start = max(1L, gene$start - as.integer(flank)),
    end = gene$end + as.integer(flank)
  )
}

# Read "chr2:135,787,850-135,862,184" (or without "chr" or commas).
# Returns NULL if the text does not look like coordinates.
parse_coordinates <- function(text) {
  clean <- gsub("[,[:space:]]", "", text)
  pattern <- "^(chr)?([0-9XYMT]+):([0-9]+)-([0-9]+)$"
  if (!grepl(pattern, clean, ignore.case = TRUE)) {
    return(NULL)
  }
  parts <- regmatches(clean, regexec(pattern, clean, ignore.case = TRUE))[[1]]
  chrom <- paste0("chr", toupper(parts[3]))
  start <- as.integer(parts[4])
  end <- as.integer(parts[5])

  if (!chrom %in% paste0("chr", 1:22)) {
    cli::cli_abort(c(
      "Only chromosomes 1 to 22 are supported, not {.val {chrom}}.",
      "i" = "Sex chromosomes and mitochondrial DNA are inherited differently and need separate handling."
    ))
  }
  if (start >= end) {
    cli::cli_abort("In {.val {text}}, the start must be smaller than the end.")
  }
  tibble::tibble(
    label = format_region(chrom, start, end),
    chrom = chrom,
    start = start,
    end = end
  )
}

check_region_size <- function(region_table) {
  total <- sum(region_table$end - region_table$start + 1)
  if (total > 5e6) {
    cli::cli_abort(c(
      "The requested regions add up to {round(total / 1e6, 1)} Mb, which is too large.",
      "i" = "Please ask for 5 Mb or less in total."
    ))
  }
  if (total > 1e6) {
    cli::cli_warn(c(
      "The requested regions add up to {round(total / 1e6, 1)} Mb.",
      "i" = "Fetching may take a few minutes and use a lot of memory."
    ))
  }
}

# --- Reading genotypes ----------------------------------------------------------

# Read one region, from the cache if we have fetched it before.
# Returns a list with `genotypes` (people x SNPs) and `variants` (a tibble).
read_region <- function(chrom, start, end) {
  cache_file <- file.path(
    cache_dir("regions"),
    paste0(chrom, "_", start, "_", end, ".rds")
  )
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  url <- sprintf(panel_url, chrom)
  index <- get_index(url, chrom)

  window_starts <- seq(start, end, by = window_size)
  windows <- purrr::map(window_starts, function(window_start) {
    window_end <- min(window_start + window_size - 1, end)
    read_window(url, index, chrom, window_start, window_end)
  })

  piece <- list(
    genotypes = do.call(cbind, purrr::map(windows, "genotypes")),
    variants = dplyr::bind_rows(purrr::map(windows, "variants"))
  )

  # Positions with more than two alleles appear as several records at the
  # same position; drop them all to keep only biallelic SNPs.
  multiallelic <- piece$variants$pos %in%
    piece$variants$pos[duplicated(piece$variants$pos)]
  piece$genotypes <- piece$genotypes[, !multiallelic, drop = FALSE]
  piece$variants <- piece$variants[!multiallelic, ]

  saveRDS(piece, cache_file)
  piece
}

# Download the index file for a chromosome into the cache, if needed.
# (If we let Rsamtools do this, it saves the index in the user's current
# working folder.)
get_index <- function(url, chrom) {
  index <- file.path(cache_dir("index"), paste0(chrom, ".vcf.gz.tbi"))
  if (!file.exists(index)) {
    tryCatch(
      utils::download.file(
        paste0(url, ".tbi"),
        index,
        mode = "wb", # binary mode, needed on Windows
        quiet = TRUE
      ),
      error = function(e) {
        unlink(index)
        cli::cli_abort(c(
          "Could not download data from the 1000 Genomes server.",
          "i" = "Check your internet connection and try again.",
          "x" = conditionMessage(e)
        ))
      }
    )
  }
  index
}

# Read the SNPs in one window and convert genotypes to allele counts.
read_window <- function(url, index, chrom, start, end) {
  param <- VariantAnnotation::ScanVcfParam(
    which = GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end)),
    samples = samples_1kg$sample_id,
    fixed = "ALT",
    info = NA, # skip the INFO column, which we do not need
    geno = "GT" # genotypes only
  )
  vcf <- VariantAnnotation::readVcf(
    Rsamtools::TabixFile(url, index = index),
    genome = "hg38",
    param = param
  )

  # Keep single-letter SNPs that start inside this window. (Long structural
  # variants that start earlier but overlap the window are also returned by
  # the query; they are removed here.)
  pos <- BiocGenerics::start(SummarizedExperiment::rowRanges(vcf))
  keep <- VariantAnnotation::isSNV(vcf, singleAltOnly = TRUE) &
    pos >= start & pos <= end
  vcf <- vcf[keep]
  pos <- pos[keep]

  if (length(pos) == 0) {
    return(empty_window())
  }

  ref <- as.character(VariantAnnotation::ref(vcf))
  alt <- as.character(unlist(VariantAnnotation::alt(vcf)))
  ids <- paste(chrom, pos, ref, alt, sep = ":")

  genotypes <- genotype_counts(VariantAnnotation::geno(vcf)$GT)
  genotypes <- genotypes[samples_1kg$sample_id, , drop = FALSE]
  colnames(genotypes) <- ids

  variants <- tibble::tibble(
    id = ids,
    chrom = rep(chrom, length(pos)),
    pos = as.integer(pos),
    ref = ref,
    alt = alt
  )

  list(genotypes = genotypes, variants = variants)
}

# What read_window() returns for a window with no SNPs.
empty_window <- function() {
  list(
    genotypes = matrix(
      integer(0),
      nrow = nrow(samples_1kg),
      ncol = 0,
      dimnames = list(samples_1kg$sample_id, NULL)
    ),
    variants = tibble::tibble(
      id = character(0),
      chrom = character(0),
      pos = integer(0),
      ref = character(0),
      alt = character(0)
    )
  )
}

# Convert VCF genotype text to allele counts.
#
# In a VCF, a genotype such as "0|1" lists the two alleles a person carries:
# 0 is the reference allele, 1 the alternative allele. "|" means the data are
# phased (we know which allele came from which parent) and "/" means unphased.
# We only need the number of alternative alleles: "0|0" -> 0, "0|1" or
# "1|0" -> 1, "1|1" -> 2. Anything else (e.g. missing, ".") becomes NA.
#
# Input: a matrix of SNPs x people. Output: an integer matrix of people x SNPs.
genotype_counts <- function(gt) {
  counts <- c(
    "0|0" = 0L, "0|1" = 1L, "1|0" = 1L, "1|1" = 2L,
    "0/0" = 0L, "0/1" = 1L, "1/0" = 1L, "1/1" = 2L
  )
  out <- matrix(
    unname(counts[gt]),
    nrow = nrow(gt),
    ncol = ncol(gt), # stated explicitly so that a window with no SNPs works
    dimnames = dimnames(gt)
  )
  t(out)
}
