### Compute a grid of betas/effect sizes and polygenic risk scores (PRS) for the range of hyper-parameters. Use the validation data to find the best PRS. ###


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
  make_option(c("--chr"),
    type = "integer",
    default = NULL,
    help = "Chromosome # for which to compute PRS",
    metavar = "integer"
  ),

  make_option(c("--summstats_matched"),
    type = "character",
    default = paste0(PRS_PHENO_DIR, "/derived_data/summary_stats/summstats_matched.txt"),
    help = "Matched summary statistics file",
    metavar = "character"
  ),

  make_option(c("--ldref_geno"),
    type = "character",
    default = paste0(PRS_PHENO_DIR, "/derived_data/LDref_panel/1000G_phase3_common/1000g_phase3_geno.rds"),
    help = "LD reference panel genotype data .rds file containing the bigSNP object [default \"%default\"]",
    metavar = "character"
  ),

  make_option(c("--validn_geno"),
    type = "character",
    default = paste0(PRS_PHENO_DIR, "/derived_data/validation_data/ukb_geno.rds"),
    help = "Validation genotype data .rds file containing the bigSNP object [default \"%default\"]",
    metavar = "character"
  ),

  make_option(c("--validn_pheno"),
    type = "character",
    default = paste0(PRS_PHENO_DIR, "/derived_data/validation_data/ukb_pheno.rds"),
    help = "Validation data .rds file containing the phenotype & covariates [default \"%default\"]",
    metavar = "character"
  ),

  make_option(c("--test_geno"),
    type = "character",
    default = paste0(PRS_PHENO_DIR, "/derived_data/test_data/epr_geno.rds"),
    help = "Test genotype data .rds file containing the bigSNP object [default \"%default\"]",
    metavar = "character"
  )
)
opts <- get_opts(opts_spec)


### main() ###

## Setup

# Directory to save output files to
out_dir <- paste0(PRS_PHENO_DIR, "/output/")

# Matched summary statistics data
matched_summstats <- bigreadr::fread2(opts$summstats_matched)
print("Summary statistics data matched across all datasets:")
str(matched_summstats)

# LD reference panel data
ldref.bigSNP <- snp_attach(opts$ldref_geno)
print("LD reference panel bigsnpr object:")
str(ldref.bigSNP, max.level = 2, strict.width = "cut")

# Validation genotype data
validn.bigSNP <- snp_attach(opts$validn_geno)
print("Validation data bigsnpr object:")
str(validn.bigSNP, max.level = 2, strict.width = "cut")

# Validation survey data
validn.survey <- readRDS(opts$validn_pheno)
print("Validation survey data:")
str(validn.survey, max.level = 2, strict.width = "cut")

# Test genotype data
test.bigSNP <- snp_attach(opts$test_geno)
print("Test data bigsnpr object:")
str(test.bigSNP, max.level = 2, strict.width = "cut")


## Compute correlation between variants
## Recommmended window size = 3cM (see ldpred2 paper).

# indices in the matched summary statistics data
ind.chr.ss <- which(matched_summstats$chr == opts$chr)
print(sprintf("Chromosome %d : Number of SNPs in matched summary statistics = %d",
  opts$chr, length(ind.chr.ss)))

df_beta <- matched_summstats[ind.chr.ss, c("beta.summstats", "beta_se", "n_eff")]
df_beta <- df_beta %>% rename(beta = beta.summstats)

# indices in the LDref panel genotype data
ind.chr.ldref <- matched_summstats$`_NUM_ID_.ldref`[ind.chr.ss]
corr_mat <- snp_cor(
  ldref.bigSNP$genotypes, ind.col = ind.chr.ldref,
  ncores = NCORES,
  infos.pos = as.vector(unlist(matched_summstats[ind.chr.ss, "genetic.pos"])),
  size = 3 / 1000, alpha = 0.9
  )
corr_sfbm <- bigsparser::as_SFBM(as(corr_mat, "dgCMatrix"))

cat("Chromosome", opts$chr, ": Dimensions of SNP correlation matrix = ",
    dim(corr_mat), "\n")
if (anyNA(corr_mat)) {
  stop("corr_mat contains NAs")
} else {
  print("corrmat OK")
}

## Build parameter grid with hyper-parameters - h2, p(rho), sparse

# Estimate values for h2
ldsc <- snp_ldsc2(corr_mat, df_beta, ncores = NCORES)
h2_est <- ldsc[["h2"]]
print(paste("Estimated h2 =", h2_est))

# Grid of models
print("Values used for h2 (heritability):")
(h2_seq <- signif(seq_log(from = h2_est / 100, to = h2_est * 10, length.out = 10), 4))

print("Values used for rho/p (proportion of SNPs with non-zero effect sizes):")
(p_seq <- signif(seq_log(from = 1e-4, to = 1, length.out = 10), 2))

print("Parameter grid:")
(params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))
print(paste("Number of parameter combinations =", nrow(params)))


## LDpred2 - Infer posterior mean effect size for each SNP by using a prior on
# effect sizes and LD information from the LD reference panel

# Adjust betas for each parameter combination
ldref.beta_grid <- snp_ldpred2_grid(
  corr = corr_sfbm, df_beta = df_beta, grid_param = params, ncores = NCORES
  )

write.table(
  ldref.beta_grid, paste0(out_dir, "ldref_adjbetas_na_chr", opts$chr, ".txt"),
  quote = FALSE, row.names = FALSE
)

# ldpred2 fails for some parameter combinations and returns all NA values.
# Remove these invalid values.
invalid_param.rows <- (1:ncol(ldref.beta_grid))[colSums(is.na(ldref.beta_grid)) == nrow(ldref.beta_grid)]

print("All adjusted beta values are 'NA' for the following parameter grid combinations:")
print(invalid_param.rows)

# Only keep the valid columns in beta_grid and valid rows in the param grid
valid_param.rows <- setdiff(1:nrow(params), invalid_param.rows)
ldref.beta_grid <- ldref.beta_grid[, valid_param.rows]
params <- params[valid_param.rows, ]

# Check if there are other NAs in the beta_grid
if (anyNA(ldref.beta_grid)) stop("Adjusted betas matrix still contains NAs.")

print(paste("Number of 'non-NA' parameter combinations =", nrow(params)))
print("'Non-NA' parameter grid:")
print(params)

# ldpred2 returns absurd/extreme beta values for certain parameter combinations.
# Remove these absurd, extreme values and keep the sane values.
std_devs <- apply(ldref.beta_grid, 2, sd)
sane <- (1:ncol(ldref.beta_grid))[abs(std_devs - median(std_devs)) < 3 * mad(std_devs)]

print(
  "Number of parameter combinations that gave insane beta values =",
  length(setdiff(1:nrow(params), sane)), "\n"
)

# Only keep the non-extreme/sane values
ldref.beta_grid <- ldref.beta_grid[, sane]
params <- params[sane, ]

print(paste("Number of sane parameter combinations =", nrow(params)))
print("sane parameter grid:")
print(params)

write.table(
  ldref.beta_grid, paste0(out_dir, "ldref_adjbetas_chr", opts$chr, ".txt"),
  quote = FALSE, row.names = FALSE
)


## Compute multiple PRSes in the validation data for each parameter combination

# Change signs for the betas for the validation genotypes
validn.beta_grid <- apply(ldref.beta_grid, 2, function(beta_vec) {
  return(beta_vec * matched_summstats$beta.validn[ind.chr.ss])
})

# Compute PRS
pred_grid <- big_prodMat(
  validn.bigSNP$genotypes,
  validn.beta_grid,
  ind.col = matched_summstats$`_NUM_ID_.validn`[ind.chr.ss]
  )

write.table(
  pred_grid, paste0(out_dir, "validn_prsgrid_chr", opts$chr, ".txt"),
  quote = FALSE, row.names = FALSE
)


## Determine best PRS based on AUC

tic("GLM logistic regression")
logreg_op <- NULL
for (i in 1:ncol(pred_grid)) {
  print(sprintf("Testing pred_grid col %d built with params row %d : ", i, i))
  print(params[i, ])
  prs_df <- cbind("f.eid" = validn.bigSNP$fam$sample.ID, PRS = pred_grid[, i])
  df <- merge(validn.survey %>% select(-ethnicity), prs_df, by = "f.eid")
  logreg_op <- rbind(logreg_op, do_glm(glm_df = df %>% select(-f.eid)))
}
toc()

# TODO: see /Downloads/prs_todos.R
# QC: required QC on logreg op - invalid pval, zval, beta, se...
# checks: if highest AUC is the one with the min p, min h or max h, then print a
# warning + suggest changing hyperparameters
# info: save whole param table with auc. save best param values for chr.
# merge for all chrs and make plots, print ranges, etc.

write.table(
  logreg_op, paste0(out_dir, "validn_glm_chr", opts$chr, ".txt"),
  quote = FALSE, row.names = FALSE
)

# Plots and tables to check for convergence of hyperparameters

# The rows in params and logreg are in the same order as the columns in
# *.beta_grid and pred_grid
params$auc <- logreg_op$auc
params %>%
  mutate(sparsity = colMeans(validn.beta_grid == 0), id = row_number()) %>%
  arrange(desc(auc)) %>%
  mutate_at(c("auc", "sparsity"), round, digits = 3) %>%
  slice(1:10)

pdf(file = paste0(PRS_PHENO_DIR, "/plots/validn_auc_chr", opts$chr, ".pdf"))
ggplot(params, aes(x = p, y = auc, color = as.factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(minor_breaks = params$p) +
  facet_wrap(~sparse, labeller = label_both) +
  labs(y = "AUC", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))
dev.off()

# Extract the id for the params that corresponds to the PRS with the highest AUC
best_id <- logreg_op %>%
  mutate(id = row_number()) %>%
  arrange(desc(auc)) %>%
  slice(1) %>%
  pull(id)

cat("Param with highest AUC = row", best_id, "\n")
print(params[best_id, ])

# Save the best PRS for the validation data
validn.best_prs <- pred_grid[, best_id]
write.table(
  set_names(data.frame(validn.best_prs), paste0("chr", opts$chr)),
  paste0(out_dir, "validn_prs_chr", opts$chr, ".txt"),
  quote = FALSE, row.names = FALSE
)


## Compute PRS in the test data using the best adjusted beta vector

# Extract the best adjusted beta vector that corresponds to the best params,
# i.e. the PRS with the highest AUC. Here, we must extract the appropriate
# column from ldref.beta_grid because the signs for the betas in the test data
# (beta.test from snp_match()) are relative to the LD reference panel.
ldref.best_betas <- ldref.beta_grid[, best_id]

# Change signs for the betas for the test genotypes
test.beta_vec <- ldref.best_betas * matched_summstats$beta.test[ind.chr.ss]

write.table(
  set_names(data.frame(test.beta_vec), paste0("chr", opts$chr)),
  paste0(out_dir, "test_bestbetas_chr", opts$chr, ".txt"),
  quote = FALSE, row.names = FALSE
)

# Compute PRS
prs_test <- big_prodVec(
  test.bigSNP$genotypes, test.beta_vec,
  ind.col = matched_summstats$`_NUM_ID_.test`[ind.chr.ss]
  )

write.table(
  set_names(data.frame(prs_test), paste0("chr", opts$chr)),
  paste0(out_dir, "test_prs_chr", opts$chr, ".txt"),
  quote = FALSE, row.names = FALSE
)


### cleanup() ###

# Turn off memory profiling
# Rprof(NULL)



### Unused code ###

# Compare AUC with auto vs grid vs average of all autos

# ## Automatic model
# # Adjust betas using the Automatic model which automatically estimates
# hyperparameters
# multi_auto <- snp_ldpred2_auto(corr_sfbm, df_beta, h2_init = h2_est,
# vec_p_init = seq_log(1e-2, 1, length.out = NCORES), sparse = TRUE, ncores = NCORES)
# str(multi_auto)

# saveRDS(multi_auto, paste0('./logs/ldpred_auto_chr', opts$chr, '.rds'))

# # TODO: add error checking + convergence checks of p_est, h2_est values for
# each list in multi-auto - see ldpred2 tutorial, compute_prs_auto.R & Downloads/random_ldpred.R

# # Extract the beta values for each estimated parameter combination in
# multi_auto. These betas are appropriate for the genotypes/alleles in the LD
# reference panel. Number of columns in the beta grid is 2xlength(multi_auto).
# The first n columns are for sparse=FALSE and the remaining n columns are for
# sparse=TRUE.
# ldref.beta_grid <- sapply(multi_auto, function(auto) return(auto$beta_est))
# ldref.beta_grid <- cbind(ldref.beta_grid, sapply(multi_auto, function(auto) return(auto$beta_est_sparse)))
