test_that("coordinates are parsed in several formats", {
  expected <- tibble::tibble(
    label = "chr2:135,787,850-135,862,184",
    chrom = "chr2",
    start = 135787850L,
    end = 135862184L
  )
  expect_equal(parse_coordinates("chr2:135787850-135862184"), expected)
  expect_equal(parse_coordinates("2:135,787,850-135,862,184"), expected)
  expect_null(parse_coordinates("LCT"))
})

test_that("unsupported or malformed coordinates give clear errors", {
  expect_error(parse_coordinates("chrX:1-1000"), "chromosomes 1 to 22")
  expect_error(parse_coordinates("chr2:2000-1000"), "start must be smaller")
})

test_that("teaching genes are resolved without the internet", {
  regions <- resolve_regions("LCT", flank = 25000)
  gene <- teaching_genes[teaching_genes$gene == "LCT", ]
  expect_equal(regions$start, gene$region_start)
  expect_equal(regions$end, gene$region_end)
})

test_that("VCF genotypes become allele counts", {
  gt <- matrix(
    c("0|0", "0|1", "1|0", "1|1", ".", "0/1"),
    nrow = 2,
    dimnames = list(c("snp1", "snp2"), c("p1", "p2", "p3"))
  )
  counts <- genotype_counts(gt)
  expect_equal(dim(counts), c(3, 2))
  expect_equal(counts["p1", ], c(snp1 = 0L, snp2 = 1L))
  expect_equal(counts["p2", ], c(snp1 = 1L, snp2 = 2L))
  expect_equal(counts["p3", ], c(snp1 = NA, snp2 = 1L))

  # a window with no SNPs
  empty <- matrix(character(0), nrow = 0, ncol = 3, dimnames = list(NULL, c("p1", "p2", "p3")))
  expect_equal(dim(genotype_counts(empty)), c(3, 0))
})

test_that("fetching a small region from the 1000 Genomes server works", {
  skip_on_cran()
  skip_if_offline("ftp.1000genomes.ebi.ac.uk")
  withr::local_options(snpscaper.cache_dir = withr::local_tempdir())

  x <- fetch_region("chr2:135851000-135851200", verbose = FALSE)
  expect_s3_class(x, "snpscape")
  expect_equal(nrow(x$samples), 2504)
  expect_true("rs4988235" %in% x$variants$rsid)
  expect_true(all(x$genotypes %in% c(0L, 1L, 2L, NA)))
})
