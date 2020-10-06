### Match SNPs and alleles across all datasets + QC across datasets ###

# TODO: Change genotype data
# snp_match(summstats, ldref.map) - 843,276 variants matched. With full phase1 data (12M SNPs after QC), 1.8M snps matched. Get full 1KG phase3 data.
# snp_match(info_snp.ld, validn.map) - 146,335 variants matched. The UKBB genotyped-only data has the fewest #snps of all 4 datasets. Get the genotyped+imputed UKBB data.


### init() ###

# Load packages and variables

# The common PRS directory for all diseases
PRS_COMMON_DIR <- Sys.getenv("PRS_COMMON_DIR")
# The PRS directory of the specific phenotype for which PRS is being computed
PRS_PHENO_DIR <- Sys.getenv("PRS_PHENO_DIR")

source(paste0(PRS_COMMON_DIR, "/scripts/prs_common.R"))

# Retrieve command line arguments/options
opts_spec <- list(
  make_option(c("--summstats_fmtd"), type="character", default=paste0(PRS_PHENO_DIR, '/derived_data/summary_stats/summstats_fmtd.txt'), 
      help="Formatted summary statistics file [default \"%default\"]", metavar="character"),

  make_option(c("--ldref_geno"), type="character", default=paste0(PRS_PHENO_DIR, '/derived_data/LDref_panel/1000G_phase3_common/1000g_phase3_geno.rds'), 
    help="LD reference panel genotype data .rds file containing the bigSNP object [default \"%default\"]", metavar="character"),

  make_option(c("--validn_geno"), type="character", default=paste0(PRS_PHENO_DIR, '/derived_data/validation_data/ukb_geno.rds'), 
    help="Validation genotype data .rds file containing the bigSNP object [default \"%default\"]", metavar="character"),

  make_option(c("--test_geno"), type="character", default=paste0(PRS_PHENO_DIR,'/derived_data/test_data/epr_geno.rds'),
    help="Test genotype data .rds file containing the bigSNP object [default \"%default\"]", metavar="character")
);
opts <- get_opts(opts_spec)


### main() ###

## Setup

# Summary stats output directory to save files to
ss_dir <-  paste0(PRS_PHENO_DIR, '/derived_data/summary_stats/')

# Summary statistics data
summstats <- bigreadr::fread2(opts$summstats_fmtd)
print("Formatted Summary statistics data:")
str(summstats)

# LD reference panel data
ldref.bigSNP <- snp_attach(opts$ldref_geno)
print("LD reference panel bigsnpr object:")
str(ldref.bigSNP, max.level = 2, strict.width = "cut")

# Validation data
validn.bigSNP <- snp_attach(opts$validn_geno)
print("Validation data bigsnpr object:")
str(validn.bigSNP, max.level = 2, strict.width = "cut")

# Test data
test.bigSNP <- snp_attach(opts$test_geno)
print("Test data bigsnpr object:")
str(test.bigSNP, max.level = 2, strict.width = "cut")

# Get aliases for the .map data and remove unwanted columns
ldref.map <- ldref.bigSNP$map %>% select(-genetic.dist)
validn.map <- validn.bigSNP$map %>% select(-genetic.dist)
test.map <- test.bigSNP$map %>% select(-genetic.dist)
summstats$pos <- NULL

# The fol. column names are used to match variants between the dataframes.
# Note: In plink format map file allele1 = minor allele, allele2= major allele
names(ldref.map) <- names(validn.map) <- names(test.map) <- c("chr", "rsid", "pos", "a1", "a0")

## Summary statistics QC as per https://github.com/privefl/paper-ldpred2/blob/master/paper/paper-ldpred2-supp.pdf - Drop variants whose standard deviations are out of range.

# Get the intersection of SNPs present in all datasets
common_snps <- Reduce(intersect, list(ldref.map$rsid, validn.map$rsid, test.map$rsid, summstats$rsid))
cat("Number of common snps across 4 datasets = ", length(common_snps), "\n")
# Only keep the common intersection of SNPs for the analysis
summstats <- summstats %>% filter(rsid %in% common_snps)
cat("Summary statistics dimensions with common snps only: ", dim(summstats), "\n")

# Compute standard deviations for summary statistics, validation & test datasets
summstats$SDss <- 2 / (summstats$beta_se * sqrt(summstats$n_eff))
SD.ss <- summstats %>% select(rsid, SDss)
summstats$SDss <- NULL

SD.validn <- data.frame(rsid = validn.bigSNP$map$marker.ID, SDval = sqrt(big_colstats(validn.bigSNP$genotypes)$var))

SD.test <- data.frame(rsid = test.bigSNP$map$marker.ID, SDtest = sqrt(big_colstats(test.bigSNP$genotypes)$var))

SD.all <- Reduce(function(x, y) merge(x, y, by = 'rsid'), list(SD.ss, SD.validn, SD.test))
print("Standard deviations of SNPs in summary statistics, validation and test datasets:")
print(str(SD.all))

# For each SNP (row), compare standard deviations across summary statistics, validation & test datasets and flag violations
SD.all$keep <- apply(SD.all %>% select(-rsid), 1, check_sd)
write.table(SD.all, paste0(ss_dir, 'SD_all.txt'), row.names = FALSE, quote = FALSE)

# Plots to check for outliers/violations
pdf(paste0(PRS_PHENO_DIR, "/plots/sdplots.pdf"), onefile = TRUE)

ggplot(SD.all, aes(x = SDval, y = SDss, colour = keep)) +
geom_point() + 
geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'Red')  +
labs(x = 'Standard deviations in the validation dataset', y = 'Standard deviations derived from the summary statistics') +
theme_bw()

ggplot(SD.all, aes(x = SDtest, y = SDss, colour = keep)) +
geom_point() + 
geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'Red')  +
labs(x = 'Standard deviations in the test dataset', y = 'Standard deviations derived from the summary statistics') +
theme_bw()

dev.off()

# Remove the flagged SNPs from the summary statistics. Only the SNPs remaining in the summary statistics will be used for further analysis.
snps.keep <- SD.all %>% filter(keep == TRUE) %>% pull(rsid)
cat("Number of SNPs to keep after SD check = ", length(snps.keep), "\n")
summstats <- summstats %>% filter(rsid %in% snps.keep)
cat("Summary statistics dimensions after filtering variants that violated standard deviations check: ", dim(summstats), "\n")


## Match SNPs & alleles across all 4 datasets
print("Match variants between (A) summary statistics and (B) LDrefpanel = (AB)")
info_snp.ld <- snp_match(summstats, ldref.map, join_by_pos = FALSE)
# rename fields because same fields will be created in the next snp_match()'es
info_snp.ld <- rename(info_snp.ld, c("_NUM_ID_.summstats" = "_NUM_ID_.ss",
                                     "_NUM_ID_.ldref" = "_NUM_ID_",
                                     "beta.summstats" = "beta", 
                                     "pos.ldref" = "pos"))
info_snp.ld$beta <- 1
str(info_snp.ld)

print("Match variants between (AB) & (C) validation genotype data = (ABC)")
info_snp.validn <- snp_match(info_snp.ld, validn.map, join_by_pos = FALSE, 
                             match.min.prop = 0)
str(info_snp.validn)

print("Match variants between (AB) & (D) test genotype data = (ABD)")
info_snp.test <- snp_match(info_snp.ld, test.map, join_by_pos = FALSE)
str(info_snp.test)

print("Merge (ABC) and (ABD) to get the final matched info across A,B,C & D")
final_info <- inner_join(info_snp.validn, info_snp.test, 
                         by = setdiff(c(names(info_snp.ld), "_NUM_ID_.ss"),
                                      c("beta", "_NUM_ID_", "pos")),
                         suffix = c(".validn", ".test"))
final_info$`_NUM_ID_.ss` <- NULL
str(final_info)
head(final_info)

# Use genetic maps to interpolate physical positions (in bp) to genetic positions (in cM).
# TODO: GRCh/hg versions need to be same for the LDrefpanel (we are using $pos from there) and for whatever snp snp_asGeneticPos is using.
final_info$genetic.pos <- snp_asGeneticPos(final_info$chr, final_info$pos.ldref, dir = paste0(PRS_COMMON_DIR, '/raw_data/1000G_genetic_maps'), ncores = NCORES)

# Save the matched SNPs & info across all 4 datasets
write.table(final_info, file = paste0(ss_dir, 'summstats_matched.txt'), quote=FALSE, row.names=FALSE)

print("Final matched SNP info/matched summary statistics:")
cat("Number of SNPs = ", nrow(final_info), "\n")
str(final_info)
head(final_info)