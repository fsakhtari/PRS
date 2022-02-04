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

  cat(
    "\n==Computing and testing total PRS for", dataset_type,
    "dataset from the following files ==\n"
  )
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

  # Create PRS dataframe
  prs_all[, id_var] <- geno_data$fam$sample.ID
  prs_total <- merge(pheno_data, prs_all[, c(id_var, "PRS")], by = id_var)

  cat("PRS dataframe for", dataset_type, "dataset:\n")
  str(prs_total)
  write.table(prs_total,
    paste0(PRS_PHENO_DIR, "/output/", dataset_type, "_prs_df.txt"),
    row.names = FALSE, quote = FALSE
  )

  # Remove columns not used in glm
  remove_cols <- c("f.eid", "ethnicity", "epr_number", "Sample_ID", "race")
  glm_df <- prs_total %>% select(-any_of(remove_cols))
  glm_df <- glm_df %>% select(-any_of(paste0('PC', 11:40)))

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
    varexp_prs, "\n"
  )


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

  plot_df <- prs_total
  # Scale/standardize the PRS
  plot_df$PRS_scaled <- scale(plot_df$PRS, scale = TRUE, center = TRUE)
  # PRS percentile
  plot_df$PRS_ptile <- ntile(plot_df$PRS, 100)
  # PRS levels - low, med, high
  plot_df <- plot_df %>% mutate(PRS_level = case_when(
    PRS_ptile < 25                    ~ "low",
    PRS_ptile >= 25 & PRS_ptile <= 75 ~ "med",
    PRS_ptile > 75                    ~ "high",
    TRUE ~ NA_character_
  ))

  write.table(plot_df,
    paste0(PRS_PHENO_DIR, "/output/", dataset_type, "_prs_plot_df.txt"),
    row.names = FALSE, quote = FALSE
  )
  save.image(file = paste0(PRS_PHENO_DIR, "/output/", dataset_type,
   "_data.RData"))

  print("Number of ppts & cases in each PRS level:")
  plot_df %>%
    group_by(PRS_level) %>%
    summarise(
      N = n(),
      N_percent = (N/nrow(.))*100,
      N_cases = sum(as.numeric(levels(Y))[Y]),
      percent_cases = mean(as.numeric(levels(Y))[Y])
    ) %>%
  print()

  # Density plot to show PRS distribution
  print(
    ggplot(data = plot_df, aes(x = PRS_scaled)) +
      geom_density(outline.type = "both") +
      ylab("Density") +
      xlab("Scaled PRS") +
      theme_classic()
  )

  # Boxplots of PRS ~ Y (disease/phenotype)
  boxplot_prs(plot_data = plot_df, prs = "PRS")
  boxplot_prs(plot_data = plot_df, prs = "PRS_scaled")
  boxplot_prs(plot_data = plot_df, prs = "PRS_ptile")

  # disease prevalence ~ prs percentile plot
  prev_ptile <- plot_df %>%
    group_by(PRS_ptile) %>%
    summarise(prev = mean(as.numeric(levels(Y))[Y]))

  print(
    ggplot(
      data = prev_ptile,
      aes(x = PRS_ptile, y = prev * 100, colour = PRS_ptile)
    ) +
      geom_point() +
      scale_color_gradientn(colours = brewer.pal(n = 9, name = "Blues")) +
      ylab("Disease prevalence %") +
      xlab("PRS percentile") +
      theme_classic()
  )

  # PRS ~ race/ethnicity plot to see racial differences in PRS
  race_var <- NULL
  if (dataset_type == "validn") {
    race_var <- "ethnicity"
  }
  else if (dataset_type == "test") {
    race_var <- "race"
  }

  print(
    ggplot(
      data = plot_df,
      aes_string(x = race_var, y = "PRS_scaled", fill = race_var)
    ) +
      geom_boxplot() +
      scale_fill_brewer(palette = "Blues") +
      ylab("PRS_scaled") +
      xlab(race_var) +
      theme_classic() +
      theme(legend.position = "none") +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14)
      )
  )

  dev.off()
}

# Wrapper function to make boxplots of PRS ~ Y (disease/phenotype)
# prs = character(); column name of the PRS value to plot
boxplot_prs <- function(plot_data, prs) {
  print(
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


# Compute the proportion of SNPs used in PRS computation across the genome

# Files containing the best betas per chromosome
beta_files <- list.files(
  path = paste0(PRS_PHENO_DIR, "/output"),
  pattern = "test_bestbetas_chr\\d*.txt",
  full.names = TRUE
)

all_betas <- NULL
for (f in beta_files) {
  df <- read.table(f, header = TRUE)
  all_betas <- c(all_betas, df[, 1])
}
sparsity <- mean(all_betas == 0)
cat("Total number of SNPs input to LDpred = ", length(all_betas), "\n")
cat(
  "Proportion of SNPs used in PRS computation (i.e. SNPs with non-zero beta
  values) =",
  (1 - sparsity), "\n"
)


## Compute and test total PRS for:

# Validation data
compute_total_prs(dataset_type = "validn")

# Test data
compute_total_prs(dataset_type = "test")


### cleanup() ###

# Turn off memory profiling
# Rprof(NULL)
