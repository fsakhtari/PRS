### Setup the summary statistics dataset for the PRS pipeline ###

# TODO:
# 1. From https://www.nature.com/articles/ng.3406#Sec8 under "We then applied the following filters (implemented in the script munge_summstats.py included with ldsc):"
# 2. QC = remove snps with OR/beta/CI/pvalue = NA, etc....


### init() ###

# Load packages and variables

# The common PRS directory for all diseases
PRS_COMMON_DIR <- Sys.getenv("PRS_COMMON_DIR")
# The PRS directory of the specific phenotype for which PRS is being computed
PRS_PHENO_DIR <- Sys.getenv("PRS_PHENO_DIR")

source(paste0(PRS_COMMON_DIR, "/scripts/prs_common.R"))


### main() ###

## Setup

# Read summary statistics data from meta-analysis
orig_file <- paste0(PRS_PHENO_DIR, "/raw_data/summary_stats/DIAGRAMv3.2012DEC17.txt")
summstats <- bigreadr::fread2(orig_file)
cat("Summary statistics data from file", orig_file, ":\n")
str(summstats)


## Modify summary statistics columns as required by LDpred2

# Convert OR and OR CIs to beta and beta_se
summstats$beta <- log(summstats$OR)
summstats$beta_se <- (log(summstats$OR_95U) - log(summstats$OR_95L)) / (2 * 1.96)

# Compute effective sample size for binary traits
summstats$n_eff <- 4 / ((1 / summstats$N_CASES) + (1 / summstats$N_CONTROLS))

# Remove unwanted columns
summstats$OR <- summstats$OR_95L <- summstats$OR_95U <- summstats$N_CASES <- summstats$N_CONTROLS <- NULL
print("Summary statistics after removing unwanted columns:")
str(summstats)

# Rename columns
names(summstats) <- c("rsid", "chr", "pos", "a1", "a0", "p", "beta", "beta_se", "n_eff")
print("Summary statistics after renaming columns:")
str(summstats)

# write formatted summstats file
fmtd_file <- paste0(PRS_PHENO_DIR, '/derived_data/summary_stats/summstats_fmtd.txt')
bigreadr::fwrite2(summstats, file = fmtd_file)

cat("\n=== Summary statistics information ===\n")
print("Summary statistics data used: An expanded genome-wide association study of type 2 diabetes in Europeans. Diabetes, 2017, Scott et al.")
cat("Summary statistics file (original/input):", orig_file, "\n")
cat("Summary statistics file (formatted/output):", fmtd_file, "\n")
cat("Number of SNPs =", length(summstats$rsid), "\n")