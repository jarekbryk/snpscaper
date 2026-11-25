test_that("the example object is consistent", {
  x <- example_lct
  expect_s3_class(x, "snpscape")
  expect_equal(rownames(x$genotypes), x$samples$sample_id)
  expect_equal(colnames(x$genotypes), x$variants$id)
  expect_equal(dim(x), c(2504L, nrow(x$variants)))
  expect_output(print(x), "2,504 people")
})

test_that("filter_snps removes rare SNPs and logs each step", {
  x <- filter_snps(example_lct, min_maf = 0.05)
  maf <- minor_allele_freq(x$genotypes)
  expect_true(all(maf >= 0.05))
  expect_equal(x$provenance$steps$step, c("fetched", "MAF >= 0.05", "missing <= 0.05"))
})

test_that("LD pruning leaves no strongly correlated SNPs within a window", {
  # Pruning compares SNPs within a sliding window of 50, moving 5 at a time,
  # so any two SNPs up to 45 positions apart are guaranteed to have been
  # compared. SNPs further apart may still be correlated.
  x <- filter_snps(example_lct, min_maf = 0.05)
  keep <- prune_by_ld(x$genotypes, x$variants, r2 = 0.2, window = 50, step = 5)
  expect_lt(sum(keep), length(keep))

  kept <- which(keep)
  r2 <- stats::cor(x$genotypes[, kept])^2
  close_pairs <- abs(outer(kept, kept, "-")) <= 45 & upper.tri(r2)
  expect_true(all(r2[close_pairs] <= 0.2))
})

test_that("keep_pops selects people and rejects unknown codes", {
  x <- keep_pops(example_lct, pop = c("GBR", "YRI"))
  expect_setequal(unique(x$samples$pop), c("GBR", "YRI"))
  expect_equal(nrow(x$genotypes), nrow(x$samples))

  x <- keep_pops(example_lct, super_pop = "EUR", sex = "female")
  expect_true(all(x$samples$super_pop == "EUR" & x$samples$sex == "female"))

  expect_error(keep_pops(example_lct, pop = "XYZ"), "Unknown")
})

test_that("as_pca_table works with FactoMineR", {
  skip_if_not_installed("FactoMineR")
  x <- example_lct |>
    filter_snps(min_maf = 0.05, prune_ld = TRUE) |>
    keep_pops(pop = c("GBR", "YRI", "CHB"))
  tab <- suppressMessages(as_pca_table(x))
  expect_equal(names(tab)[1:3], c("pop", "super_pop", "sex"))
  expect_equal(rownames(tab), x$samples$sample_id)

  pca <- FactoMineR::PCA(tab, quali.sup = 1:3, graph = FALSE)
  expect_equal(nrow(pca$ind$coord), nrow(x$samples))
})
