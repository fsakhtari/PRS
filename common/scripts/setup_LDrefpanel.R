### Setup the LD reference panel data for the PRS pipeline ###


### init() ###

# Load packages and variables

# The common PRS directory for all diseases
PRS_COMMON_DIR <- Sys.getenv("PRS_COMMON_DIR")
# The PRS directory of the specific phenotype for which PRS is being computed
PRS_PHENO_DIR <- Sys.getenv("PRS_PHENO_DIR")

source(paste0(PRS_COMMON_DIR, "/scripts/prs_common.R"))


### main() ###

# 1000G Phase3 data downloaded via bigsnpr, for 2490 mostly unrelated
# individuals and ~1.7M SNPs in common with either HapMap3 or the UK Biobank.
ldref_in_dir <- paste0(PRS_COMMON_DIR, "/raw_data/LDref_panel/1000G_phase3_common/")
ldref_out_dir <- paste0(PRS_PHENO_DIR, "/derived_data/LDref_panel/1000G_phase3_common/")

# QC the PLINK data files and create new QCed PLINK files
bedfile <- snp_plinkQC(
  plink.path = "/ddn/gs1/home/akhtarifs/software/plink2",
  prefix.in = paste0(ldref_in_dir, "1000G_phase3_common_norel"),
  prefix.out = paste0(ldref_out_dir, "1000g_phase3_geno"),
  maf = 0.01,
  geno = 0.01,
  mind = 0.1,
  hwe = 1e-10,
  autosome.only = TRUE,
  extra.options = " --snps-only"
)

# Read from QCed bed/bim/fam PLINK files. This generates .bk and .rds files.
rdsfile <- snp_readBed2(bedfile, ncores = NCORES)

# Attach the "bigSNP" object from the .rds file in R session
ldref.bigSNP <- snp_attach(rdsfile)

# See how the file looks like
print("LD reference panel bigsnpr object:")
str(ldref.bigSNP, max.level = 2, strict.width = "cut")

# Get MAF for QCed data
maf <- snp_MAF(ldref.bigSNP$genotypes, ncores = NCORES)

# Get population relationship info for each sample
popn_info <- ldref.bigSNP$fam
ped <- bigreadr::fread2(paste0(ldref_in_dir, "20130606_g1k.ped"))
popn_info <- dplyr::left_join(popn_info, ped,
  by = c("sample.ID" = "Individual ID")
)
popns <- bigreadr::fread2(paste0(ldref_in_dir, "20131219.populations.tsv"))
popn_info <- dplyr::left_join(popn_info, popns,
  by = c("Population" = "Population Code")
)
popn_info$popn_code <- paste(
  popn_info$`Super Population`, popn_info$Population,
  sep = "_"
)
write.table(
  popn_info, paste0(ldref_out_dir, "1000g_phase3_geno_popinfo.txt"),
  row.names = FALSE
)

# Relatedness
# https://cran.r-project.org/web/packages/bigsnpr/bigsnpr.pdf states that these
# are mostly unrelated. Upon examination, this subset does seem mostly unrelated,
# i.e. one of a pair of close relatives has been removed.
# However, if you only want to keep individuals who have Relationship == 'unrel',
#  then the following code will keep 907 samples from the 2490.
# ind_norel <- which(popn_info$Relationship == "unrel")
# # Get MAF using the subset of unrel individuals
# maf <- snp_MAF(ldref.bigSNP$genotypes, ind_unrel, ncores = NCORES)
# # Only keep the 'unrel' individuals and SNPs which pass MAF for this sample subset
# rdsfile_norel <- snp_subset(ldref.bigSNP, ind.row = ind_norel, ind.col = which(maf > 0.01), backingfile = '_norel')
# # note: we do lose more snps here when filtering for maf after removing relatives.
# # Attach the "bigSNP" object from the .rds file in R session
# ldref.bigSNP <- snp_attach(rdsfile_norel)
# print("LD reference panel bigsnpr object with 'unrel' only:")
# str(ldref.bigSNP, max.level = 2, strict.width = "cut")

G <- ldref.bigSNP$genotypes
print("LD reference panel summary (after QC):")
print("LD Reference panel used: 1000 Genome Phase 3 common")
cat("Number of SNPs = ", ncol(G), "\n")
cat("Number of individuals/samples = ", nrow(G), "\n")
cat("Minor allele frequency (MAF) range = ", range(maf), "\n")
print("Number of individuals per super-population:")
print(table(popn_info$`Super Population`))
print("Number of individuals per sub-population:")
print(table(popn_info$Population))
