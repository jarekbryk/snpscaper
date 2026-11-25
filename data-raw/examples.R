# Build the example snpscape objects shipped with the package, so that
# examples, tests and vignettes work offline.
#
# * example_lct: a gene under strong recent selection in Europeans.
# * example_neutral: a 300 kb stretch of chromosome 18 with no gene known to
#   be under strong selection, as a comparison.

pkgload::load_all()

example_lct <- fetch_region("LCT")
example_neutral <- fetch_region("chr18:1000000-1300000")

usethis::use_data(example_lct, example_neutral, overwrite = TRUE, compress = "xz")
