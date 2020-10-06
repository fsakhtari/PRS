### Setup the Test data for the PRS pipeline ###

# This script assumes that genotype QC similar to that done prior to a GWAS has already been done for the Test genotype data.

## TODO:
# Remove related individuals using kinship matrix
# Add make_option for covariates
# Add PCs to covariates


### init() ###

# Enable memory profiling
# Rprof(memory.profiling = TRUE, gc.profiling = TRUE, line.profiling = TRUE)

# Load packages and variables

# The common PRS directory for all diseases
PRS_COMMON_DIR <- Sys.getenv("PRS_COMMON_DIR")
# The PRS directory of the specific phenotype for which PRS is being computed
PRS_PHENO_DIR <- Sys.getenv("PRS_PHENO_DIR")

source(paste0(PRS_COMMON_DIR, "/scripts/prs_common.R"))

# Retrieve command line arguments/options
opts_spec <- list(
  make_option(c("--pheno"), type = "character", default = NULL, 
              help = "Variable name in the EPR H&E Survey data of the phenotype for which to compute PRS", metavar = "character")
);
opts <- get_opts(opts_spec)


### main() ###

## Setup

# Directory to save files to
test_out_dir <-  paste0(PRS_PHENO_DIR, "/derived_data/test_data/")

# EPR H&E Survey data - loads epr.he & epr.he.meta
load("/ddn/gs1/project/controlled/EPR/data_S3/Surveys/Health and Exposure/healthexposure_02jan20_v04.RData")

# EPR demographic data - loads epr.bcbb.map & epr.bcbb.map.meta
load("/ddn/gs1/project/controlled/EPR/data_S3/MAP Data/bcbb_map_02jan20_v05.RData")

# Link file linking the EPR participant numbers to the WGS plink IIDs
#epr_broad_map <- read_excel('/ddn/gs1/project/nextgen/controlled/wgsepr/Contolled/plink/Burkholder_WGSDemographicsV2IF_03172020.xlsx')
epr_broad_map <- read.csv("/ddn/gs1/project/controlled/EPR/data_S3/WGS Data/WGS_ID_Key.csv")
epr_broad_map <- epr_broad_map %>% select(-control_flag)


## Phenotype and covariate data from EPR Surveys

# Merge all data
epr_survey_cc <- merge(epr.he, epr.bcbb.map, by = 'epr_number')
cat("dimensions (epr + map) =", nrow(epr_survey_cc), "\n")
epr_survey_cc <- merge(epr_survey_cc, epr_broad_map, by = 'epr_number')
cat("dimensions (epr_survey_cc + broad map) =", dim(epr_survey_cc), "\n")
print(paste("number of unique epr ids", length(unique(epr_survey_cc$epr_number))))
print(paste("number of unique broad ids", length(unique(epr_survey_cc$Sample_ID))))

# Format the data as required + keep complete cases only
epr_survey_cc <- epr_survey_cc %>%
  select(c(epr_number, Sample_ID, age_derived, gender, he_bmi_derived, race, opts$pheno)) %>%
  rename(Y = !!opts$pheno) %>%
  #mutate(Y = .data[[opts$pheno]]) %>%
  replace_with_na_all(condition = ~.x %in% EPR_NA_STRINGS) %>%
  drop_na() %>%
  mutate(Sample_ID = as.character(Sample_ID)) %>%
  mutate_at(vars(epr_number, age_derived, he_bmi_derived), as.numeric) %>%
  mutate_at(vars(gender, Y, race), as.factor)

print("EPR survey data (complete cases):")
str(epr_survey_cc)

# Save the WGS IDs for complete cases in plink's required format
indivs_keep <- paste0(test_out_dir, 'indivs.keep')
write.table(data.frame(FID = epr_survey_cc$Sample_ID, IID = epr_survey_cc$Sample_ID), file = indivs_keep, quote = FALSE, row.names = FALSE, col.names = FALSE)


## Genotype data

# Filter genotype data + only keep genotype samples for complete cases
bedfile <- snp_plinkQC(plink.path = '/ddn/gs1/home/akhtarifs/software/plink2',
  prefix.in = "/ddn/gs1/project/nextgen/controlled/wgsepr/Contolled/plink/all_confirmed_filtered",
  prefix.out = paste0(test_out_dir, 'epr_geno'),
  maf = 0.01, geno = 0.01, mind = 0.1, hwe = 1e-10, autosome.only = TRUE,extra.options = paste(" --snps-only --keep", indivs_keep)
  )

# Read from filtered bed/bim/fam PLINK files. This generates .bk and .rds files.
rdsfile <- snp_readBed2(bedfile, ncores = NCORES)

# Attach the "bigSNP" object from the .rds file in R session
test.bigSNP <- snp_attach(rdsfile)
print("Test genotype data (bigsnpr object):")
str(test.bigSNP, max.level = 2, strict.width = "cut")

# Impute missing genotypes, since missing values are not allowed when computing PRS
print("Imputing missing genotypes")
test.bigSNP$genotypes <- snp_fastImputeSimple(test.bigSNP$genotypes, method = "random")
# Save the modified bigSNP object to the RDS file
rdsfile <- snp_save(test.bigSNP)

# TODO: rm the bed/bim/fam files created above and do snp_writebed() to create new ones with the imputed genotypes.

# Complete cases with Survey + WGS data
cc.iids <- test.bigSNP$fam$sample.ID
# Only keep complete cases in the survey data
epr_survey_cc <- epr_survey_cc %>% filter(Sample_ID %in% cc.iids)
print("Test survey data:")
str(epr_survey_cc)

saveRDS(epr_survey_cc, file = paste0(test_out_dir, 'epr_pheno.rds'))


## Summarize data

cat("\n=== Test data summary (after QC & filtering) ===\n")

print("Phenotype and covariate data from EPR H&E Survey")
num_cases <- sum(epr_survey_cc$Y == 1)
num_controls = sum(epr_survey_cc$Y == 0)
num_total = nrow(epr_survey_cc)

cat("Phenotype = ", opts$pheno, "\n")
cat("Cases = ", num_cases, "\n")
cat("Controls = ", num_controls, "\n")
cat("Total = ", num_total, "\n")
cat("%Cases = ", num_cases/num_total, "\n")

print("Summary of 'age_derived':")
print(summary(epr_survey_cc$age_derived))

print("Summary of 'he_bmi_derived' :")
print(summary(epr_survey_cc$he_bmi_derived))

print("Cross-tabulation of 'race' :")
print(table(epr_survey_cc$race))
print(xtabs(~ Y + race, data = epr_survey_cc))

print("Cross-tabulation of 'gender' :")
print(table(epr_survey_cc$gender))
print(xtabs(~ Y + gender, data = epr_survey_cc))

print("Genotype data from EPR WGS")
G <- test.bigSNP$genotypes
maf <- snp_MAF(G, ncores = NCORES)
cat("Number of SNPs = ", ncol(G), "\n")
cat("Number of individuals/samples = ", nrow(G), "\n")
cat("Minor allele frequency (MAF) range = ", range(maf, na.rm = TRUE), "\n")


### cleanup() ###

# Turn off memory profiling
# Rprof(NULL)
