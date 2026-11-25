#' snpscaper: explore human genetic variation with 1000 Genomes data
#'
#' @description
#' snpscaper is a teaching package. It lets you pick any region of the human
#' genome, fetch the genotypes of 2,504 people from 26 populations sampled by
#' the 1000 Genomes Project, and explore how genetic variation is distributed
#' within and between those populations. Everything happens inside R: there
#' is no need for the command line or for downloading files by hand.
#'
#' ## Key ideas, briefly
#'
#' * **SNP** (single nucleotide polymorphism): a position in the genome where
#'   people differ by a single DNA letter, for example some carry `G` and
#'   others `A`. Each person carries two copies of every autosome, so at a SNP
#'   they carry 0, 1 or 2 copies of the alternative letter (allele). snpscaper
#'   stores genotypes as these counts.
#' * **VCF** (variant call format): the standard text file format for storing
#'   genotypes. The 1000 Genomes Project publishes one compressed VCF per
#'   chromosome, each tens of gigabytes. snpscaper reads only the small slice
#'   you ask for, using an index file that tells it where in the big file a
#'   region lives.
#' * **1000 Genomes Project**: an international project that sequenced the
#'   genomes of people from 26 populations, grouped into five continental
#'   "super-populations": African (AFR), Admixed American (AMR), East Asian
#'   (EAS), European (EUR) and South Asian (SAS). See [samples_1kg].
#'
#' ## Typical workflow
#'
#' 1. [fetch_region()] gets genotypes for a gene or a set of coordinates.
#' 2. [filter_snps()] removes rare and poorly genotyped SNPs, and optionally
#'    SNPs that are strongly correlated with their neighbours (LD pruning).
#' 3. [keep_pops()] narrows the data to the populations you want to compare.
#' 4. [as_pca_table()] prepares the data for principal component analysis
#'    with the FactoMineR and factoextra packages.
#'
#' @keywords internal
"_PACKAGE"
