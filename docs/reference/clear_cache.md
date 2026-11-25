# Delete downloaded data

snpscaper saves every region it fetches, plus small index files for each
chromosome, in a cache folder, so that fetching the same region again is
instant and works offline. `clear_cache()` deletes these files, for
example to free disk space or to force a fresh download.

## Usage

``` r
clear_cache(verbose = TRUE)
```

## Arguments

- verbose:

  Print the folder that was cleared.

## Value

The path of the cache folder, invisibly.

## Details

By default the cache lives in a folder chosen by R for your user
account. To use a different folder, for example a shared folder that
staff fill before a practical, set
`options(snpscaper.cache_dir = "path/to/folder")`.

## Examples

``` r
if (FALSE) { # \dontrun{
clear_cache()
} # }
```
