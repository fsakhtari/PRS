### Setup the validation data for the PRS pipeline ###


### init() ###

# Enable memory profiling
Rprof(memory.profiling = TRUE, gc.profiling = TRUE, line.profiling = TRUE)

# Load packages and variables

# The common PRS directory for all diseases
PRS_COMMON_DIR <- Sys.getenv("PRS_COMMON_DIR")
# The PRS directory of the specific phenotype for which PRS is being computed
PRS_PHENO_DIR <- Sys.getenv("PRS_PHENO_DIR")

source(paste0(PRS_COMMON_DIR, "/scripts/prs_common.R"))

# The first 10 genetic principal components in the UK Biobank data
pcs <- setNames(as.list(paste0("f.22009.0.", 1:10)), paste0("PC", 1:10))

# Retrieve command line arguments/options
opts_spec <- list(
  make_option(c("--pheno"),
    type = "character",
    default = NULL,
    help = "Phenotype for which to compute PRS",
    metavar = "character"
  ),

  make_option(c("--covariates"),
    type = "list",
    default = c(list(
      age = "f.21022.0.0", sex = "f.31.0.0",
      ethnicity = "f.21000.0.0", bmi = "f.21001.0.0"
    ), pcs),
    help = "covariates to use in the association analysis specified in a named list of covariate names and their corresponding UIDs in the UK Biobank data [default \"%default\"]",
    metavar = "named list"
  )
)
opts <- get_opts(opts_spec)


### main() ###

## Setup

# Directory to save files to
valdn_out_dir <- paste0(PRS_PHENO_DIR, "/derived_data/validation_data/")

# UK Biobank data
ukb_data <- fread(
  "/ddn/gs1/shared/ukbiobank/data_versions/data_v1/dataformat_R/ukb41671.tab"
)
print("UK Biobank (UKB) survey data :")
cat("dimensions =", dim(ukb_data), "\n")
print("first 50 column names :")
print(colnames(ukb_data)[1:50])


## Phenotype and covariate data

# The values/special codes to turn into missing = NA
na_strings <- c(-1, -3)

# covariate column names to extract from the dataset
covars <- as.character(unlist(opts$covariates))

# Identify cases and controls for the specified phenotype/disease
ukb_pheno_prepd <- prepare_ukb_phenotype(
  ukb_data = ukb_data, phenotype = opts$pheno
)

ukb_data <- merge(ukb_data, ukb_pheno_prepd, by = "f.eid")

ukb_survey_cc <- ukb_data %>%
  select(all_of(c("f.eid", "Y", covars))) %>%
  rename_at(vars(covars), function(x) names(opts$covariates)) %>%
  replace_with_na_at(
    .vars = c("ethnicity"),
    condition = ~ .x %in% na_strings
  ) %>%
  drop_na() %>%
  # Group the ethnicities into the top-level ethnic group which is represented
  # by the first number in the string
  mutate(ethnicity = substr(ethnicity, 1, 1)) %>%
  mutate(f.eid = as.integer(f.eid)) %>%
  mutate_at(vars(age, bmi), as.numeric) %>%
  mutate_at(vars(sex, Y, ethnicity), as.factor)

print("UKB survey data (complete cases):")
str(ukb_survey_cc)

# Identify related samples with the phenotype of interest, to remove
reldness <- read.table("/ddn/gs1/shared/ukbiobank/data_versions/data_v1/genotype_data/common_data/ukb57849_rel_s488264.dat",
  header = TRUE
)
indivs.remove <- ukb_gen_samples_to_remove(
  data = reldness,
  ukb_with_data = ukb_survey_cc$f.eid,
  cutoff = 0.0884
)
cat("Number of related individuals to remove: ", length(indivs.remove), "\n")

ukb_survey_cc <- ukb_survey_cc %>%
  anti_join(data.frame(f.eid = indivs.remove), by = "f.eid")
cat("UKB dimensions after removing related individuals:", dim(ukb_survey_cc), "\n")

# Save the WGS IDs for complete cases in plink's required format
indivs_keep <- paste0(valdn_out_dir, "indivs.keep")
write.table(data.frame(FID = ukb_survey_cc$f.eid, IID = ukb_survey_cc$f.eid),
  file = indivs_keep, quote = FALSE, row.names = FALSE, col.names = FALSE
)

# Free-up memory
rm(ukb_data, ukb_pheno_prepd, reldness, indivs.remove)
gc()


## Genotype data

# Filter genotype data + only keep genotype samples for complete cases
print("starting snp_plinkqc")
bedfile <- snp_plinkQC(
  plink.path = "/ddn/gs1/home/akhtarifs/software/plink2",
  prefix.in = "/ddn/gs1/shared/ukbiobank/data_versions/data_v1/genotype_data/genotypes/all_chromosomes/ukb_cal_all_v2",
  prefix.out = paste0(valdn_out_dir, "ukb_geno"),
  maf = 0.01, geno = 0.01, mind = 0.1, hwe = 1e-10, autosome.only = TRUE,
  extra.options = paste0(" --snps-only --keep ", indivs_keep)
)
print("finished snp_plinkqc")

# Read from filtered bed/bim/fam PLINK files. This generates .bk and .rds files.
rdsfile <- snp_readBed2(bedfile, ncores = NCORES)

# Attach the "bigSNP" object from the .rds file in R session
validn.bigSNP <- snp_attach(rdsfile)
print("Validation data bigsnpr object:")
str(validn.bigSNP, max.level = 2, strict.width = "cut")

# Impute missing genotypes, since missing values are not allowed when computing PRS
print("Imputing missing genotypes")
validn.bigSNP$genotypes <- snp_fastImputeSimple(validn.bigSNP$genotypes,
  method = "random", ncores = NCORES
)
# Save the modified bigSNP object to the RDS file
rdsfile <- snp_save(validn.bigSNP)

# Complete cases with genotype + phenotype data
cc.iids <- validn.bigSNP$fam$sample.ID
# Only keep complete cases in the survey data
ukb_survey_cc <- ukb_survey_cc %>% filter(f.eid %in% cc.iids)
cat("Final UKB survey data dimensions:", dim(ukb_survey_cc), "\n")

saveRDS(ukb_survey_cc, file = paste0(valdn_out_dir, "ukb_pheno.rds"))


## Summarize data

cat("\n=== Validation data summary (after QC & filtering) ===\n")

print("Phenotype and covariate data from UK Biobank")
num_cases <- sum(ukb_survey_cc$Y == 1)
num_controls <- sum(ukb_survey_cc$Y == 0)
num_total <- nrow(ukb_survey_cc)

cat("Phenotype = ", opts$pheno, "\n")
cat("Cases = ", num_cases, "\n")
cat("Controls = ", num_controls, "\n")
cat("Total = ", nrow(ukb_survey_cc), "\n")
cat("%Cases = ", num_cases / num_total, "\n")

print("Summary of 'age':")
print(summary(ukb_survey_cc$age))

print("Summary of 'bmi' :")
print(summary(ukb_survey_cc$bmi))

print("Cross-tabulation of 'ethnicity' :")
print(table(ukb_survey_cc$ethnicity))
print(xtabs(~ Y + ethnicity, data = ukb_survey_cc))

print("Cross-tabulation of 'sex' :")
print(table(ukb_survey_cc$sex))
print(xtabs(~ Y + sex, data = ukb_survey_cc))

print("Genotype data from UK Biobank")
maf <- snp_MAF(validn.bigSNP$genotypes, ncores = NCORES)
cat("Number of SNPs = ", ncol(validn.bigSNP$genotypes), "\n")
cat("Number of individuals/samples = ", nrow(validn.bigSNP$genotypes), "\n")
cat("Minor allele frequency (MAF) range = ", range(maf, na.rm = TRUE), "\n")


### cleanup() ###

# Turn off memory profiling
Rprof(NULL)
