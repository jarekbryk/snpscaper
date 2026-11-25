# snpscaper keeps a cache folder on disk with two kinds of files:
#
# * index files (.tbi), one per chromosome, each about 200 KB. An index lists
#   where in the big chromosome file each stretch of the genome is stored, so
#   that we can download only the part we need.
# * fetched regions (.rds), so that asking for the same region again is
#   instant and works offline.
#
# The folder is chosen by R (see tools::R_user_dir()) unless the user sets
# options(snpscaper.cache_dir = "some/folder"), which is useful for staff who
# want to prepare data for a class in a shared folder in advance.

cache_dir <- function(...) {
  dir <- getOption(
    "snpscaper.cache_dir",
    tools::R_user_dir("snpscaper", which = "cache")
  )
  path <- file.path(dir, ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

#' Delete downloaded data
#'
#' snpscaper saves every region it fetches, plus small index files for each
#' chromosome, in a cache folder, so that fetching the same region again is
#' instant and works offline. `clear_cache()` deletes these files, for example
#' to free disk space or to force a fresh download.
#'
#' By default the cache lives in a folder chosen by R for your user account.
#' To use a different folder, for example a shared folder that staff fill
#' before a practical, set `options(snpscaper.cache_dir = "path/to/folder")`.
#'
#' @param verbose Print the folder that was cleared.
#'
#' @return The path of the cache folder, invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' clear_cache()
#' }
clear_cache <- function(verbose = TRUE) {
  dir <- cache_dir()
  unlink(list.files(dir, full.names = TRUE), recursive = TRUE)
  if (verbose) {
    cli::cli_inform("Cleared the snpscaper cache at {.path {dir}}.")
  }
  invisible(dir)
}
