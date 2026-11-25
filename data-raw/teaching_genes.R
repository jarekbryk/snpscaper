# Build `teaching_genes` and `teaching_snps`: a curated set of genes with
# well-known patterns of variation between human populations, and one or two
# "focal" SNPs in each that drive the pattern.
#
# Positions are looked up in Ensembl (GRCh38) rather than typed in by hand,
# so they are always consistent with the 1000 Genomes data.
#
# The `description` is neutral and suitable for student handouts.
# `staff_notes` describe the expected pattern and are shown only in the staff
# version of the report.

pkgload::load_all()

# Fetch a little extra sequence either side of each gene, so that nearby
# regulatory variants are included. For LCT this matters: the main European
# lactase persistence variant (rs4988235) lies about 14 kb upstream of LCT,
# inside the neighbouring MCM6 gene.
flank <- 25000

genes <- tibble::tribble(
  ~gene, ~description, ~staff_notes,

  "LCT",
  "Lactase, the enzyme that digests the milk sugar lactose. In most mammals its production stops after weaning.",
  "Lactase persistence (continued lactase production in adults) is common in Northern Europe and parts of Africa and the Middle East, and rare in East Asia. The European persistence allele at rs4988235 (in an MCM6 intron, acting as an LCT enhancer) is at high frequency in FIN, GBR and CEU and near absent in EAS. It sits on a long, little-broken haplotype: a classic signal of recent positive selection linked to dairying. African persistence alleles (e.g. rs145946881) are mostly in pastoralist populations not sampled by 1000 Genomes.",

  "MC1R",
  "Melanocortin 1 receptor, which controls the switch between dark (eumelanin) and red/yellow (pheomelanin) pigment in skin and hair.",
  "Several amino-acid changing variants (e.g. rs1805007 R151C, rs1805008 R160W) cause red hair and fair skin and are found mainly in Europeans. African populations show very few amino-acid changing MC1R variants, consistent with purifying selection for dark pigmentation under high UV, while outside Africa constraint appears relaxed.",

  "SLC24A5",
  "A potassium-dependent sodium/calcium exchanger found in melanosomes, involved in skin pigmentation.",
  "The derived A allele at rs1426654 (Ala111Thr) is close to fixation in Europeans, common in South Asians, and rare in East Asians and most Africans. One of the strongest single-SNP differences between populations; a textbook example of selection on light skin pigmentation.",

  "EDAR",
  "Ectodysplasin A receptor, involved in the development of hair, teeth and sweat glands.",
  "The derived allele at rs3827760 (Val370Ala) is common in East Asians and Native American-ancestry populations (high in PEL, also MXL) and essentially absent in Europeans and Africans. Associated with thicker hair shafts and shovel-shaped incisors; thought to have risen in frequency by selection in East Asia.",

  "ALDH2",
  "Aldehyde dehydrogenase 2, a mitochondrial enzyme that breaks down acetaldehyde produced when the body metabolises alcohol.",
  "The ALDH2*2 allele at rs671 (Glu504Lys) inactivates the enzyme and causes the alcohol flushing reaction. It is found almost only in East Asians (common in CHB, CHS, JPT, KHV) and is essentially absent elsewhere. Carriers have a higher risk of oesophageal cancer if they drink.",

  "ACKR1",
  "Atypical chemokine receptor 1 (Duffy antigen), a red blood cell surface protein that the malaria parasite Plasmodium vivax uses to enter cells.",
  "The Duffy-null allele at rs2814778 (a promoter variant that switches off expression on red cells) is near fixation in sub-Saharan Africa and rare elsewhere, giving resistance to P. vivax malaria. One of the most extreme allele frequency differences between populations in the human genome; very high Fst. ASW and ACB show intermediate frequencies due to admixture.",

  "HBB",
  "Haemoglobin subunit beta, part of the oxygen-carrying protein in red blood cells.",
  "The sickle cell allele (HbS, rs334, Glu6Val) is found mainly in African populations (and the African-ancestry ASW and ACB) and is a classic case of balancing selection: heterozygotes are protected against severe falciparum malaria, homozygotes have sickle cell disease. Allele frequencies are moderate (a few percent to ~10%), not high.",

  "CYP1A2",
  "Cytochrome P450 1A2, a liver enzyme that metabolises caffeine and many medicines.",
  "rs762551 (CYP1A2*1F) affects how inducible the enzyme is (faster caffeine metabolism). Allele frequencies vary moderately between populations without an extreme pattern. A useful contrast to the strongly differentiated genes: most genes look like this.",

  "HTT",
  "Huntingtin, a large protein needed for nerve cell function. Expansions of a CAG repeat in this gene cause Huntington's disease.",
  "The disease-causing CAG repeat expansion is not captured by SNP data. Population structure here mostly reflects shared ancestry rather than selection, so HTT serves as a large, relatively 'ordinary' region: PCA and admixture patterns resemble genome-wide structure (Africa most diverse, separation by continent)."
)

focal <- tibble::tribble(
  ~gene,      ~rsid,         ~note,
  "LCT",      "rs4988235",   "European lactase persistence variant, known as -13910 C>T (T = persistence). LCT lies on the reverse strand, so on the genome's forward strand used here it is G>A (A = persistence)",
  "MC1R",     "rs1805007",   "R151C, red hair variant",
  "MC1R",     "rs1805008",   "R160W, red hair variant",
  "SLC24A5",  "rs1426654",   "Ala111Thr, light skin pigmentation",
  "EDAR",     "rs3827760",   "Val370Ala, hair thickness",
  "ALDH2",    "rs671",       "Glu504Lys (ALDH2*2), alcohol flush",
  "ACKR1",    "rs2814778",   "Duffy-null promoter variant",
  "HBB",      "rs334",       "Glu6Val, sickle cell (HbS)",
  "CYP1A2",   "rs762551",    "CYP1A2*1F, caffeine metabolism"
)

# --- Look up gene positions --------------------------------------------------
gene_positions <- purrr::map(genes$gene, function(symbol) {
  g <- ensembl_get("lookup", "symbol", "homo_sapiens", symbol)
  tibble::tibble(
    gene = symbol,
    chrom = paste0("chr", g$seq_region_name),
    gene_start = as.integer(g$start),
    gene_end = as.integer(g$end),
    strand = as.integer(g$strand)
  )
}) |>
  dplyr::bind_rows()

teaching_genes <- genes |>
  dplyr::left_join(gene_positions, by = "gene") |>
  dplyr::mutate(
    region_start = gene_start - flank,
    region_end = gene_end + flank
  ) |>
  dplyr::select(
    "gene", "chrom", "gene_start", "gene_end", "strand",
    "region_start", "region_end", "description", "staff_notes"
  )

# --- Look up focal SNP positions ---------------------------------------------
snp_positions <- purrr::map(focal$rsid, function(rsid) {
  s <- ensembl_get("variation", "human", rsid)
  m <- purrr::keep(s$mappings, function(m) m$seq_region_name %in% c(1:22, "X"))[[1]]
  tibble::tibble(
    rsid = rsid,
    chrom = paste0("chr", m$seq_region_name),
    pos = as.integer(m$start),
    alleles = m$allele_string
  )
}) |>
  dplyr::bind_rows()

teaching_snps <- focal |>
  dplyr::left_join(snp_positions, by = "rsid") |>
  dplyr::select("gene", "rsid", "chrom", "pos", "alleles", "note")

# Every focal SNP must fall inside the region fetched for its gene.
check <- teaching_snps |>
  dplyr::left_join(teaching_genes, by = c("gene", "chrom")) |>
  dplyr::mutate(inside = pos >= region_start & pos <= region_end)
stopifnot(all(check$inside))

usethis::use_data(teaching_genes, teaching_snps, overwrite = TRUE)
