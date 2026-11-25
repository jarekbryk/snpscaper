# Make the figures shown in README.md (saved in man/figures).

pkgload::load_all()

save_figure <- function(plot, name, width = 7, height = 4.5) {
  ggplot2::ggsave(
    file.path("man", "figures", paste0(name, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 110,
    bg = "white"
  )
}

# Allele frequency map: SLC24A5, a pigmentation gene
fetch_region("SLC24A5") |>
  calc_allele_freqs() |>
  map_freqs(snp = "rs1426654") |>
  save_figure("readme_map", height = 4.2)

# PCA of the neutral region
pca_table <- example_neutral |>
  filter_snps(prune_ld = TRUE) |>
  as_pca_table()
pca <- FactoMineR::PCA(pca_table, quali.sup = 1:3, graph = FALSE)
factoextra::fviz_pca_ind(
  pca,
  habillage = "super_pop",
  palette = unname(superpop_colours),
  geom = "point",
  pointshape = 19,
  mean.point = FALSE,
  title = "PCA, 300 kb of chromosome 18"
) +
  theme_jb()
ggplot2::last_plot() |> save_figure("readme_pca")

# LD around LCT in Europeans
example_lct |>
  keep_pops(super_pop = "EUR") |>
  filter_snps(min_maf = 0.05) |>
  calc_ld() |>
  plot_ld(title = "LD around LCT in Europeans") |>
  save_figure("readme_ld", width = 5.5, height = 5)

# Admixture, neutral region
example_neutral |>
  filter_snps(prune_ld = TRUE) |>
  run_admixture(k = 2:3, reps = 3) |>
  plot_admixture(k = 2:3) |>
  save_figure("readme_admixture", height = 3.5)
