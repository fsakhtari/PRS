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
orig_file <- paste0(PRS_PHENO_DIR, "/raw_data/summary_stats/29273806-GCST005212-EFO_0000270.h.tsv")
summstats <- bigreadr::fread2(orig_file)
cat("Summary statistics data from file", orig_file, ":\n")
str(summstats)


## Modify summary statistics columns as required by LDpred2

summstats <- summstats %>% 
  select(hm_rsid, hm_chrom, hm_pos, hm_other_allele, hm_effect_allele, hm_beta, standard_error, p_value) %>%
  drop_na()
names(summstats) <- c("rsid", "chr", "pos", "a0", "a1", "beta", "beta_se", "p")

# Compute effective sample size for binary traits
# From $PRS_PHENO_DIR/raw_data/summary_stats/41588_2017_14_MOESM3_ESM.xlsx, sum of 'Number of cases' and 'Number of controls'
n_cases <- 23948
n_controls <- 118538
summstats$n_eff <- 4 / ((1 / n_cases) + (1 / n_controls))

print("Formatted summary statistics:")
str(summstats)

# write formatted summstats file
fmtd_file <- paste0(PRS_PHENO_DIR, '/derived_data/summary_stats/summstats_fmtd.txt')
bigreadr::fwrite2(summstats, file = fmtd_file)

cat("\n=== Summary statistics information ===\n")
print("Summary statistics data used: Demenais, Florence, et al. Multiancestry association study identifies new asthma risk loci that colocalize with immune-cell enhancer marks. Nature genetics (2018).")
cat("Summary statistics file (original/input):", orig_file, "\n")
cat("Summary statistics file (formatted/output):", fmtd_file, "\n")
cat("Number of SNPs =", length(summstats$rsid), "\n")