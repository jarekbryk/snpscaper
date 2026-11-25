test_that("a report is rendered to a self-contained HTML file", {
  skip_on_cran()
  skip_if_not_installed("quarto")
  skip_if_not_installed("FactoMineR")
  skip_if_not_installed("factoextra")
  skip_if(is.null(quarto::quarto_path()), "Quarto is not installed")

  output_dir <- withr::local_tempdir()
  x <- keep_pops(example_lct, pop = c("YRI", "GBR", "CHB"))
  path <- render_report(
    x,
    k = 2:3,
    reps = 1,
    title = "Test report",
    output_dir = output_dir,
    verbose = FALSE
  )

  expect_true(file.exists(path))
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(html, "Test report")
  expect_match(html, "Notes for staff")
  expect_match(html, "rs4988235")
})

test_that("rs numbers not in teaching_snps are labelled for the report", {
  skip_on_cran()
  skip_if_offline("rest.ensembl.org")
  x <- add_rsids(example_lct, "rs182549") # another SNP in the LCT region
  expect_true("rs182549" %in% x$variants$rsid)
})
