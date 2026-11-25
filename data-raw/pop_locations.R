# Build `pop_locations`: approximate map coordinates for the 26 populations.
#
# Several 1000 Genomes populations were sampled away from the region their
# ancestry comes from (for example Gujarati Indians in Houston, Sri Lankan
# Tamils in the UK, and Utah residents of Northern and Western European
# ancestry). For maps of allele frequencies it is more informative to place
# these populations at their region of ancestral origin, so we do that and
# record it in `placed_at`. The two African-diaspora populations (ASW, ACB)
# and the admixed American populations are placed where they were sampled,
# because their ancestry comes from several continents.
#
# Coordinates are approximate (city or region level) and are meant only for
# drawing maps.

pop_locations <- tibble::tribble(
  ~pop,  ~latitude, ~longitude, ~placed_at,
  "ACB",  13.10,  -59.62, "sampling location",
  "ASW",  33.50,  -97.50, "sampling location",
  "ESN",   6.70,    6.30, "sampling location",
  "GWD",  13.45,  -16.58, "sampling location",
  "LWK",   0.60,   34.77, "sampling location",
  "MSL",   8.48,  -13.23, "sampling location",
  "YRI",   7.38,    3.93, "sampling location",
  "CLM",   6.25,  -75.56, "sampling location",
  "MXL",  19.43,  -99.13, "ancestral origin (sampled in Los Angeles)",
  "PEL", -12.05,  -77.04, "sampling location",
  "PUR",  18.22,  -66.59, "sampling location",
  "CDX",  22.01,  100.80, "sampling location",
  "CHB",  39.90,  116.40, "sampling location",
  "CHS",  26.07,  119.30, "sampling location",
  "JPT",  35.68,  139.69, "sampling location",
  "KHV",  10.82,  106.63, "sampling location",
  "CEU",  50.50,    5.00, "ancestral origin (sampled in Utah)",
  "FIN",  60.17,   24.94, "sampling location",
  "GBR",  53.00,   -2.00, "sampling location",
  "IBS",  40.42,   -3.70, "sampling location",
  "TSI",  43.77,   11.25, "sampling location",
  "BEB",  23.81,   90.41, "sampling location",
  "GIH",  22.26,   71.19, "ancestral origin (sampled in Houston)",
  "ITU",  16.50,   80.64, "ancestral origin (sampled in the UK)",
  "PJL",  31.55,   74.34, "sampling location",
  "STU",   9.66,   80.02, "ancestral origin (sampled in the UK)"
)

stopifnot(nrow(pop_locations) == 26)

usethis::use_data(pop_locations, overwrite = TRUE)
