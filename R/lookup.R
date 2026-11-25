# Look up genes and SNPs in Ensembl.
#
# Ensembl (https://www.ensembl.org) is a genome database. Its REST API lets
# us ask simple questions over the internet, such as "where is the LCT gene?"
# or "where is SNP rs4988235?", and get the answer back as JSON. All answers
# use the GRCh38 genome assembly, the same one as the 1000 Genomes data that
# snpscaper uses.

ensembl_url <- "https://rest.ensembl.org"

# Send one request to the Ensembl REST API and return the parsed JSON, or
# NULL if Ensembl says the thing we asked for does not exist (HTTP 400).
ensembl_get <- function(...) {
  response <- httr2::request(ensembl_url) |>
    httr2::req_url_path_append(...) |>
    httr2::req_headers(Accept = "application/json") |>
    httr2::req_user_agent("snpscaper R package") |>
    httr2::req_retry(max_tries = 3) |>
    httr2::req_error(is_error = function(resp) FALSE) |>
    httr2::req_perform()

  status <- httr2::resp_status(response)
  if (status == 400) {
    return(NULL)
  }
  if (status >= 300) {
    cli::cli_abort(c(
      "Could not reach the Ensembl database (HTTP status {status}).",
      "i" = "Check your internet connection, or try again in a few minutes."
    ))
  }
  httr2::resp_body_json(response)
}

#' Find the position of a gene
#'
#' Looks up a human gene by its symbol (for example `"LCT"` or `"MC1R"`) in
#' the Ensembl database and returns where it sits in the genome (GRCh38
#' assembly).
#'
#' Gene symbols are the short official names approved by the HUGO Gene
#' Nomenclature Committee. If you are unsure of a symbol, search for the gene
#' at <https://www.genenames.org> or <https://www.ensembl.org>.
#'
#' Genes listed in [teaching_genes] are looked up there first, so they work
#' without an internet connection.
#'
#' @param symbol A gene symbol, or a vector of them.
#'
#' @return A tibble with one row per gene: `gene`, `chrom` (chromosome, e.g.
#'   `"chr2"`), `start`, `end`, `strand` and `description`.
#' @export
#'
#' @examples
#' \dontrun{
#' find_gene("LCT")
#' find_gene(c("MC1R", "SLC24A5"))
#' }
find_gene <- function(symbol) {
  lapply(symbol, find_one_gene) |>
    dplyr::bind_rows()
}

find_one_gene <- function(symbol) {
  known <- teaching_genes[teaching_genes$gene == toupper(symbol), ]
  if (nrow(known) == 1) {
    return(dplyr::select(
      known,
      "gene", "chrom", "start" = "gene_start", "end" = "gene_end",
      "strand", "description"
    ))
  }

  gene <- ensembl_get("lookup", "symbol", "homo_sapiens", symbol)
  if (is.null(gene)) {
    cli::cli_abort(c(
      "Gene {.val {symbol}} was not found in Ensembl.",
      "i" = "Check the spelling, or search for the official gene symbol at
             {.url https://www.genenames.org}."
    ))
  }

  tibble::tibble(
    gene = gene$display_name,
    chrom = paste0("chr", gene$seq_region_name),
    start = as.integer(gene$start),
    end = as.integer(gene$end),
    strand = as.integer(gene$strand),
    description = sub(" \\[.*\\]$", "", gene$description %||% NA_character_)
  )
}

#' Find the position of a SNP from its rs number
#'
#' Most well-studied SNPs have an "rs number" (reference SNP ID) given by the
#' dbSNP database, for example `rs4988235`, the main lactase persistence
#' variant in Europeans. The 1000 Genomes files do not use rs numbers, so
#' snpscaper looks up where each rs number sits in the genome (GRCh38) and
#' finds it by position.
#'
#' @param rsid An rs number such as `"rs4988235"`, or a vector of them.
#'
#' @return A tibble with one row per SNP: `rsid`, `chrom`, `pos` and
#'   `alleles` (as reported by Ensembl, reference allele first, e.g. `"G/A"`).
#' @export
#'
#' @examples
#' \dontrun{
#' find_snp("rs4988235")
#' }
find_snp <- function(rsid) {
  lapply(rsid, find_one_snp) |>
    dplyr::bind_rows()
}

find_one_snp <- function(rsid) {
  known <- teaching_snps[teaching_snps$rsid == rsid, ]
  if (nrow(known) == 1) {
    return(dplyr::select(known, "rsid", "chrom", "pos", "alleles"))
  }

  snp <- ensembl_get("variation", "human", rsid)
  if (is.null(snp)) {
    cli::cli_abort("SNP {.val {rsid}} was not found in Ensembl.")
  }

  # A SNP can map to several places (e.g. alternative versions of a region);
  # keep the mapping on a main chromosome.
  on_main_chromosome <- purrr::keep(
    snp$mappings,
    function(m) m$seq_region_name %in% c(1:22, "X", "Y")
  )
  if (length(on_main_chromosome) == 0) {
    cli::cli_abort("SNP {.val {rsid}} does not map to a main chromosome.")
  }
  mapping <- on_main_chromosome[[1]]

  tibble::tibble(
    rsid = rsid,
    chrom = paste0("chr", mapping$seq_region_name),
    pos = as.integer(mapping$start),
    alleles = mapping$allele_string
  )
}
