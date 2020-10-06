### Combine PRS scores across chromosomes and use the combined PRS on the Test dataset. ###


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
  ),

  make_option(c("--test_pheno"),
    type = "character",
    default = paste0(PRS_PHENO_DIR, "/derived_data/test_data/epr_pheno.rds"),
    help = "Test data .rds file containing the phenotype & covariates[default \"%default\"]",
    metavar = "character"
  )
)
opts <- get_opts(opts_spec)


### Functions ###

# Worker function to compute & test total PRS
# dataset_type = c("validn", "test") indicates which dataset to compute total PRS for.
compute_total_prs <- function(dataset_type) {
  id_var <- geno_data <- pheno_data <- NULL

  if (dataset_type == "validn") {
    id_var <- "f.eid"
    geno_data <- validn.bigSNP
    pheno_data <- validn.survey
  }
  else if (dataset_type == "test") {
    id_var <- "Sample_ID"
    geno_data <- test.bigSNP
    pheno_data <- test.survey
  }
  else {
    stop("compute_total_prs: dataset_type must be one of 'validn' or 'test' to indicate which dataset you want to compute total PRS for.")
  }

  ## Compute total PRS

  # Files containing per-chromosome PRS
  prs_files <- list.files(
    path = paste0(PRS_PHENO_DIR, "/output"),
    pattern = paste0(dataset_type, "_prs_chr\\d*.txt"),
    full.names = TRUE
  )

  cat("\n==Computing and testing total PRS for", dataset_type,
      "dataset from the following files ==\n")
  print(prs_files)

  # Put the per-chromosome PRS into one dataframe
  prs_per_chr <- NULL
  for (f in prs_files) {
    df <- read.table(f, header = TRUE)
    prs_per_chr <- c(prs_per_chr, df)
  }
  prs_all <- bind_cols(prs_per_chr)

  # Total PRS = sum of per-chromosome PRS
  prs_all$PRS <- rowSums(prs_all)


  ## Test the association of total PRS with phenotype

  # Create dataframe for logistic regression
  prs_all[, id_var] <- geno_data$fam$sample.ID
  prs_total <- prs_all[, c(id_var, "PRS")]
  glm_df <- merge(pheno_data, prs_total, by = id_var)

  cat("PRS dataframe for", dataset_type, "dataset:\n")
  str(glm_df)
  write.table(glm_df,
    paste0(PRS_PHENO_DIR, "/output/", dataset_type, "_prs_df.txt"),
    row.names = FALSE, quote = FALSE
  )

  # TODO: temporarily using race in lieu of PCs. Grouping race into 3 races.
  if (dataset_type == "test") {
    glm_df <- glm_df %>%
      mutate(race = recode(race, `3` = "black", `5` = "white", .default = "other"))
    print("Test survey data after recoding race:")
    xtabs(~ Y + race, data = glm_df)
    str(glm_df)
  }

  # Remove columns not used in glm
  remove_cols <- c("f.eid", "ethnicity", "epr_number", "Sample_ID")
  glm_df <- glm_df %>% select(-any_of(remove_cols))

  # Logistic regression of phenotype (Y) ~ covariates + PRS
  glm_out <- do_glm(
    glm_df = glm_df,
    plot_file = paste0(PRS_PHENO_DIR, "/plots/", dataset_type, "_roc.pdf")
  )
  cat("GLM output for", dataset_type, "dataset:\n")
  print(glm_out)

  # Estimate the variance explained by PRS
  # TODO: replace glm in do_glm() with lrm to get both glm output & r2
  print("Calculating pseudo-R2:")
  print("Full model (with PRS):")
  lrm_full <- lrm(Y ~ ., data = glm_df)
  print(lrm_full)

  print("Reduced model (without PRS):")
  lrm_reduced <- lrm(Y ~ ., data = glm_df %>% select(-PRS))
  print(lrm_reduced)

  varexp_prs <- lrm_full$stats["R2"] - lrm_reduced$stats["R2"]
  cat(
    "Proportion of variance explained by PRS using Nagelkerke’s pseudo-R2 metric =",
    varexp_prs, "\n")


  ## Flipping the independent and dependent variables
  ## anova of PRS ~ Y to see if the mean PRS differs between cases and controls
  cat("Proportion of cases =", sum(glm_df$Y == 1) / nrow(glm_df), "\n")
  cat("mean PRS for controls =", mean(glm_df[glm_df$Y == 0, "PRS"]), "\n")
  cat("mean PRS for cases =", mean(glm_df[glm_df$Y == 1, "PRS"]), "\n")

  print("Anova of PRS ~ Y :")
  aov_obj <- aov(PRS ~ Y, data = glm_df)
  print(summary(aov_obj))


  ## Make plots

  pdf(paste0(PRS_PHENO_DIR, "/plots/", dataset_type, "_prs_plots.pdf"))

  plot_df <- glm_df
  # Scale/standardize the PRS
  plot_df$PRS_scaled <- scale(plot_df$PRS, scale = TRUE, center = TRUE)
  # PRS percentile
  plot_df$PRS_ptile <- ntile(plot_df$PRS, 100)

  # Boxplots of PRS ~ Y (phenotype)
  boxplot_prs(plot_data = plot_df, prs = "PRS")
  boxplot_prs(plot_data = plot_df, prs = "PRS_scaled")
  boxplot_prs(plot_data = plot_df, prs = "PRS_ptile")

  # disease prevalence ~ prs percentile plot
  prev_ptile <- plot_df %>%
    group_by(PRS_ptile) %>%
    summarise(prev = mean(as.numeric(levels(Y))[Y]))

  ggplot(
    data = prev_ptile,
    aes(x = PRS_ptile, y = prev * 100, colour = PRS_ptile)
  ) +
    geom_point() +
    scale_color_gradientn(colours = brewer.pal(n = 9, name = "Blues")) +
    ylab("Disease prevalence %") +
    xlab("PRS percentile") +
    theme_classic()

  dev.off()
}

# Wrapper function to make boxplots of PRS ~ Y (phenotype)
# prs = character(); column name of the PRS value to plot
boxplot_prs <- function(plot_data, prs) {
  ggplot(data = plot_data, aes_string(x = "Y", y = prs, fill = "Y")) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Paired") +
    ylab(prs) +
    xlab("Disease phenotype") +
    theme_classic() +
    theme(legend.position = "none") +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold", size = 14)
    )
}


### main() ###

## Setup

# Validation genotype data
validn.bigSNP <- snp_attach(opts$validn_geno)
print("Validation genotype data:")
str(validn.bigSNP, max.level = 2, strict.width = "cut")

# Validation phenotype/survey data
validn.survey <- readRDS(opts$validn_pheno)
print("Validation survey data:")
str(validn.survey)

# Test genotype data
test.bigSNP <- snp_attach(opts$test_geno)
print("Test genotype data:")
str(test.bigSNP, max.level = 2, strict.width = "cut")

# Test phenotype/survey data
test.survey <- readRDS(opts$test_pheno)
print("Test survey data:")
str(test.survey)


## Compute and test total PRS

# Validation data
compute_total_prs(dataset_type = "validn")

# Test data
compute_total_prs(dataset_type = "test")


### cleanup() ###

# Turn off memory profiling
# Rprof(NULL)
