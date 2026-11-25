# Find the position of a gene

Looks up a human gene by its symbol (for example `"LCT"` or `"MC1R"`) in
the Ensembl database and returns where it sits in the genome (GRCh38
assembly).

## Usage

``` r
find_gene(symbol)
```

## Arguments

- symbol:

  A gene symbol, or a vector of them.

## Value

A tibble with one row per gene: `gene`, `chrom` (chromosome, e.g.
`"chr2"`), `start`, `end`, `strand` and `description`.

## Details

Gene symbols are the short official names approved by the HUGO Gene
Nomenclature Committee. If you are unsure of a symbol, search for the
gene at <https://www.genenames.org> or <https://www.ensembl.org>.

Genes listed in
[teaching_genes](https://jarekbryk.github.io/snpscaper/reference/teaching_genes.md)
are looked up there first, so they work without an internet connection.

## Examples

``` r
if (FALSE) { # \dontrun{
find_gene("LCT")
find_gene(c("MC1R", "SLC24A5"))
} # }
```
