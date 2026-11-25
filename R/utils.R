#' @importFrom rlang %||% .data
NULL

# ggplot2::map_data() draws its world map from the maps package.
#' @importFrom maps map
NULL

# Format a number with thousands separators, e.g. 2504 -> "2,504".
comma <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

# Format a genomic region as "chr2:135,787,850-135,862,184".
format_region <- function(chrom, start, end) {
  paste0(chrom, ":", comma(start), "-", comma(end))
}
