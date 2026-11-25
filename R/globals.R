# Datasets shipped with the package and used inside its functions. Declaring
# them here stops R CMD check from reporting them as undefined variables.
utils::globalVariables(c(
  "pop_locations",
  "samples_1kg",
  "teaching_genes",
  "teaching_snps"
))
