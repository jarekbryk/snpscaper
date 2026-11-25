small <- example_lct |>
  filter_snps(min_maf = 0.05) |>
  keep_pops(pop = c("YRI", "GBR", "CHB"))

test_that("allele frequencies match a hand calculation", {
  freqs <- calc_allele_freqs(small)
  snp <- small$variants$id[1]
  gbr <- small$samples$pop == "GBR"
  by_hand <- sum(small$genotypes[gbr, snp]) / (2 * sum(gbr))
  expect_equal(freqs$alt_freq[freqs$id == snp & freqs$pop == "GBR"], by_hand)
  expect_equal(nrow(freqs), 3 * ncol(small$genotypes))

  super <- calc_allele_freqs(small, by = "super_pop")
  expect_false("pop" %in% names(super))
})

test_that("Hudson Fst behaves as expected", {
  expect_equal(hudson_fst(p1 = c(0, 1), p2 = c(1, 0), n1 = c(1e6, 1e6), n2 = c(1e6, 1e6)), 1, tolerance = 1e-4)
  expect_lt(abs(hudson_fst(p1 = c(0.3, 0.6), p2 = c(0.3, 0.6), n1 = c(200, 200), n2 = c(200, 200))), 0.01)

  fst <- calc_fst(small)
  expect_equal(nrow(fst), 3)
  expect_true(all(fst$fst > 0 & fst$fst < 1))
})

test_that("heterozygosity matches a hand calculation", {
  het <- calc_heterozygosity(small)
  gbr <- small$genotypes[small$samples$pop == "GBR", ]
  expect_equal(het$observed[het$pop == "GBR"], mean(gbr == 1))

  p <- colMeans(gbr) / 2
  n <- 2 * nrow(gbr)
  expect_equal(het$expected[het$pop == "GBR"], mean(n / (n - 1) * 2 * p * (1 - p)))
})

test_that("LD values are between 0 and 1 and capped in number", {
  ld <- suppressMessages(calc_ld(small, max_snps = 50))
  expect_equal(nrow(ld), 50 * 49 / 2)
  expect_true(all(ld$r2 >= 0 & ld$r2 <= 1 + 1e-12))
})

test_that("admixture runs and proportions add up to 1", {
  skip_if_not_installed("LEA")
  x <- filter_snps(small, prune_ld = TRUE, min_maf = 0.01)
  adm <- run_admixture(x, k = 2:3, reps = 2, verbose = FALSE)
  expect_s3_class(adm, "snpscape_admixture")
  totals <- tapply(adm$ancestry$proportion, paste(adm$ancestry$sample_id, adm$ancestry$k), sum)
  expect_true(all(abs(totals - 1) < 1e-3))
  expect_equal(nrow(get_ancestry(adm, k = 3, wide = TRUE)), nrow(x$samples))
  expect_s3_class(plot_admixture(adm, k = 2:3), "ggplot")
  expect_s3_class(plot_cross_entropy(adm), "ggplot")
})

test_that("plots are created", {
  freqs <- calc_allele_freqs(example_lct)
  expect_s3_class(plot_freqs(freqs, "rs4988235"), "ggplot")
  expect_s3_class(map_freqs(freqs, "rs4988235"), "ggplot")
  expect_error(plot_freqs(freqs, "rs123"), "not in these data")
  expect_s3_class(plot_fst(calc_fst(small)), "ggplot")
  expect_s3_class(plot_heterozygosity(calc_heterozygosity(small)), "ggplot")
  expect_s3_class(plot_ld(suppressMessages(calc_ld(small, max_snps = 30))), "ggplot")
})
