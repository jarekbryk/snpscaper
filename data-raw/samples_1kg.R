# Build `samples_1kg`: one row per person in the 1000 Genomes high-coverage
# panel who is not closely related to anyone else in the panel (2,504 people).
#
# Sources (1000 Genomes FTP at EBI):
# * 20130606_g1k_3202_samples_ped_population.txt: all 3,202 samples with sex,
#   population and super-population. The extra 698 are relatives (mostly
#   children in parent-child trios) sequenced later.
# * 1000G_2504_high_coverage.sequence.index: lists the original 2,504
#   unrelated samples. We keep only these, because close relatives share long
#   stretches of DNA and would distort PCA and admixture results.

base_url <- paste0(
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/",
  "1000G_2504_high_coverage/"
)

ped <- readr::read_delim(
  paste0(base_url, "20130606_g1k_3202_samples_ped_population.txt"),
  delim = " ",
  show_col_types = FALSE
)

unrelated_ids <- readr::read_tsv(
  paste0(base_url, "1000G_2504_high_coverage.sequence.index"),
  comment = "##",
  show_col_types = FALSE
) |>
  dplyr::pull(SAMPLE_NAME) |>
  trimws() |>
  unique()

stopifnot(length(unrelated_ids) == 2504)

pop_names <- tibble::tribble(
  ~pop,  ~pop_name,
  "ACB", "African Caribbean in Barbados",
  "ASW", "African Ancestry in Southwest US",
  "ESN", "Esan in Nigeria",
  "GWD", "Gambian in Western Division, The Gambia",
  "LWK", "Luhya in Webuye, Kenya",
  "MSL", "Mende in Sierra Leone",
  "YRI", "Yoruba in Ibadan, Nigeria",
  "CLM", "Colombian in Medellin, Colombia",
  "MXL", "Mexican Ancestry in Los Angeles, CA, USA",
  "PEL", "Peruvian in Lima, Peru",
  "PUR", "Puerto Rican in Puerto Rico",
  "CDX", "Chinese Dai in Xishuangbanna, China",
  "CHB", "Han Chinese in Beijing, China",
  "CHS", "Han Chinese South, China",
  "JPT", "Japanese in Tokyo, Japan",
  "KHV", "Kinh in Ho Chi Minh City, Vietnam",
  "CEU", "Utah residents with Northern and Western European ancestry",
  "FIN", "Finnish in Finland",
  "GBR", "British from England and Scotland",
  "IBS", "Iberian populations in Spain",
  "TSI", "Toscani in Italy",
  "BEB", "Bengali in Bangladesh",
  "GIH", "Gujarati Indians in Houston, TX, USA",
  "ITU", "Indian Telugu in the UK",
  "PJL", "Punjabi in Lahore, Pakistan",
  "STU", "Sri Lankan Tamil in the UK"
)

superpop_names <- tibble::tribble(
  ~super_pop, ~super_pop_name,
  "AFR", "African",
  "AMR", "Admixed American",
  "EAS", "East Asian",
  "EUR", "European",
  "SAS", "South Asian"
)

samples_1kg <- ped |>
  dplyr::filter(SampleID %in% unrelated_ids) |>
  dplyr::transmute(
    sample_id = SampleID,
    pop = Population,
    super_pop = Superpopulation,
    sex = dplyr::if_else(Sex == 1, "male", "female")
  ) |>
  dplyr::left_join(pop_names, by = "pop") |>
  dplyr::left_join(superpop_names, by = "super_pop") |>
  dplyr::arrange(super_pop, pop, sample_id)

stopifnot(
  nrow(samples_1kg) == 2504,
  !anyNA(samples_1kg)
)

usethis::use_data(samples_1kg, overwrite = TRUE)
