### Case-Only Analysis of Gene-Environment Interactions Using Polygenic Risk
# Scores.

### TODOs:
# Re-add non-binary exposures
# *QC* the data before modeling - Make plots, tables, etc. Remove variables
# with large missing% or N < 30 ?
# put in separate txt files: vars.dropped, vars.kept
# put betas in .out
# Go through dropped and kept variables and see if this is what we want.
# Prune more variables, since several variables don't make sense to keep


### init() ###

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

# Computed PRS scores for the EPR ppts
test_prs <- read.table(paste0(PRS_PHENO_DIR, '/output/test_prs_df.txt'), header
    = TRUE)
test_prs$PRS <- scale(test_prs$PRS, scale = TRUE, center = TRUE)

# EPR Survey data
# loads epr.he & epr.he.meta
load("/ddn/gs1/project/controlled/EPR/data_S3/Surveys/Health and Exposure/healthexposure_02jan20_v05.RData")
# loads epr.ea & epr.ea.meta
load("/ddn/gs1/project/controlled/EPR/data_S3/Surveys/Exposome/exposomea_02jan20_v02.RData")
# loads epr.eb & epr.eb.meta
load("/ddn/gs1/project/controlled/EPR/data_S3/Surveys/Exposome/exposomeb_02jan20_v02.RData")


### Functions - move to prs_common.R

# Convert columns in EPR dataframes to appropriate type
epr_convert_type <- function(epr_df, epr_df_meta) {
  for (col in colnames(epr_df)) {
    print(col)
    #true_class <- attributes(epr_df[[col]])$true_class
    true_class <- epr_df_meta %>% filter(long_variable_name == col) %>% pull(true_class)
    print(paste("true class = ",true_class))

   epr_df[[col]] <- switch(true_class,
                   "numeric" = as.numeric(epr_df[[col]]),
                   # coding binary as factor so that we can check the minimum n per level
                   "binary" = as.factor(epr_df[[col]]), #factor or numeric?
                   "character" = as.character(epr_df[[col]]),
                  "factor" = as.factor(epr_df[[col]]),
                  # coding ordered factor as numeric since we are doing lm() with these
                  "ordered factor" = as.numeric(epr_df[[col]]), #factor or numeric?
                   epr_df[[col]]
    )
  }
  return(epr_df)
}

# Create dataframe with selected variables with appropriate class/type
create_exp_df <- function(df, df_meta, vars) {
  exp_df <- df %>%
  select(all_of(c('epr_number', vars))) %>%
  replace_with_na_all(condition = ~.x %in% EPR_NA_STRINGS)

  exp_df <- epr_convert_type(exp_df, df_meta)
  return(exp_df)
}

# temporary wrapper
get_all_exps <- function() {
# Get all exposures from the EPR Survey data

# Drop height, weight, bmi variables from H&E
he_drop_vars <- qc(he_a001a_height_ft, he_a001b_height_in,
 he_a001c_height_in_derived, he_a002_weight,he_bmi_derived, he_bmi_cat_derived)

he_exp_vars <- epr.he.meta %>%
  filter(survey_section %in% c("About Your General Health", "Occupation",
   "Exposures", "Lifestyle", "About Your Home Life", "About Your Mood")) %>%
  filter(! long_variable_name %in% he_drop_vars) %>%
  # TODO: tmp only using binary exposures
  filter(true_class != "character") %>%
  # TODO: tmp only using 1st level questions
  filter(is_child == "N") %>%
  pull(long_variable_name) %>%
  c(opts$pheno)

expa_exp_vars <- epr.ea.meta %>%
  filter(! survey_section %in% c("Not Applicable", "Characteristics of Current and Past Residences")) %>%
  # TODO: tmp only using binary exposures
  filter(true_class != "character") %>%
  # TODO: tmp only using 1st level questions
  filter(is_child == "N") %>%
  pull(long_variable_name)

expb_exp_vars <- epr.eb.meta %>%
  filter(! survey_section %in% c("Not Applicable", "Genetic History")) %>%
  # TODO: tmp only using binary exposures
  filter(true_class != "character") %>%
  # TODO: tmp only using 1st level questions
  filter(is_child == "N") %>%
  pull(long_variable_name)

# Create exposures dataframe
he_exp <- create_exp_df(epr.he, epr.he.meta, he_exp_vars)
expa <- create_exp_df(epr.ea, epr.ea.meta, expa_exp_vars)
expb <- create_exp_df(epr.eb, epr.eb.meta, expb_exp_vars)
all_exps <- merge(he_exp, expa, by = "epr_number", all = TRUE)
all_exps <- merge(all_exps, expb, by = "epr_number", all = TRUE)
cat("dims of all_exps = ", dim(all_exps), "\n")

saveRDS(all_exps, paste0(PRS_PHENO_DIR, "/results/tmp/all_exps.rds"))
return(all_exps)
}

all_exps <- get_all_exps()
#all_exps <- readRDS(paste0(PRS_PHENO_DIR, "/results/tmp/all_exps.rds"))

# Merge with PRS data
prs_e <- merge(test_prs[qc(epr_number, PRS)], all_exps, by = "epr_number")
cat("dims of prs_e = ", dim(prs_e), "\n")
saveRDS(prs_e, paste0(PRS_PHENO_DIR, "/results/tmp/prs_e.rds"))

# only keep phenotype cases for this analysis
prs_e <- prs_e %>% filter(.data[[opts$pheno]] == 1)
cat("dims of prs_e for cases only = ", dim(prs_e), "\n")
print(unique(prs_e[opts$pheno]))

prs_e$epr_number <- prs_e[opts$pheno] <- NULL
print("colnames of prs_e:")
print(colnames(prs_e))

# Drop variables with low sample sizes
vars.drop <- NULL

# Drop numeric variables with n < 30
num_vars <- prs_e %>% select(where(is.numeric)) %>% select(-PRS) %>% colnames()
print("Numeric variables:")
print(num_vars)
for (col in num_vars) {
  cc <- sum(!is.na(prs_e[, col]))
  if (cc < 30) {
    vars.drop <- c(vars.drop, col)
  }
}
print("Dropping following numeric variables:")
print(vars.drop)

# Drop categorical/factor variable with levels < 2 or n < 30 for any level
cat_vars <- prs_e %>% select(where(is.factor)) %>% colnames()
print("Factor variables:")
print(cat_vars)
for (col in cat_vars) {
  nlevels <- prs_e %>% count(.data[[col]]) %>% drop_na() %>% nrow()
  if (nlevels < 2) {
    print(paste("variable", col, "has < 2 levels"))
    vars.drop <- c(vars.drop, col)
    next
  }
  min_n <- prs_e %>%
    count(.data[[col]]) %>%
    drop_na() %>%
    slice_min(n, n=1, with_ties=FALSE) %>%
    pull(n)
  print(paste("variable:", col, "; min_n =", min_n))

  if (min_n < 30) {
    vars.drop <- c(vars.drop, col)
  }
}
print(paste("Dropping following", length(vars.drop), "variables"))
print(vars.drop)

prs_e <- prs_e %>% select(-all_of(vars.drop))

print("prs_e after QC:")
cat("dims of prs_e = ", dim(prs_e), "\n")
print("colnames of prs_e:")
print(colnames(prs_e))

pdf(paste0(PRS_PHENO_DIR, "/results/tmp/co_prs_e.pdf"))
all_res <- NULL

for (col in colnames(prs_e)) {
  if (col == 'PRS') next
  cat("Testing for exposure:", col, "\n")
  test_df <- prs_e %>% select(PRS, !!col)

  # TODO: reconsider model/test used: not right to use anova with severely
  # unbalanced data
  aov_obj <- aov(as.formula(paste('PRS ~', col)), data = test_df, singular.ok = FALSE)
  print(aov_obj)
  summ <- summary(aov_obj)
  print(summ)
  print(aov_obj$coefficients)
  pval <- summ[[1]][col, 'Pr(>F)']

  # lm_obj <- lm(as.formula(paste('PRS ~', col)), data = test_df, singular.ok = FALSE)
  # print(lm_obj)
  # summ <- summary(lm_obj)
  # print(summ)
  # print(coefficients(lm_obj))
  # pval <- coefficients(summ)[col, 'Pr(>|t|)']

  all_res <- rbind(all_res, data.frame(exposure = col, pvalue = pval, n = nobs(aov_obj)))

  if (class(test_df[, col]) == 'factor') {
    test_df %>%
    count(.data[[col]], sort = TRUE) %>% print()

    print(ggplot(data = test_df, aes_string(x = col, y = "PRS", fill = col)) +
    geom_boxplot() +
    scale_fill_brewer(palette="Paired") +
    ylab("PRS") +
    xlab(col) +
    theme_classic() +
    theme(legend.position = "none") +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text = element_text(face = "bold", size = 14),
          axis.title = element_text(face = "bold", size = 14)))

  }
}
# TODO: else for numeric variables, do lm plot to see direction

# MTC
all_res$p.adjust <- p.adjust(all_res$pvalue, method = 'fdr')

write.table(all_res, paste0(PRS_PHENO_DIR, "/results/tmp/co_prs_e.out"),
 row.names = FALSE, quote = FALSE)
dev.off()



### Old code ###

# Since all exposures are binary for now, just sum to see #exposed
# prs_e <- prs_e[which(colSums(prs_e) > 30)]
# TODO : some variables have very low #unexposed - e.g. do you have ac, heat,
# job,etc. So need to check all levels of variable.


# # Drop factor variables with low sample count in any level
# vars.drop <- NULL
# for (col in colnames(prs_e)) {
#   nlevels <- prs_e %>% count(.data[[col]]) %>% drop_na() %>% nrow()

#   # for now assuming that variables with < 10 levels are factors/categorical. could switch to only doing this for factor variables, if variable types are set appropriately before.
#   if (nlevels > 0 & nlevels <= 10) {
#     min_n <- prs_e %>%
#       count(.data[[col]]) %>%
#       drop_na() %>%
#       slice_min(n, n=1, with_ties=FALSE) %>%
#       pull(n)
#     print(paste("variable:", col, "; min_n =", min_n))

#     if (min_n < 30) {
#       vars.drop <- c(vars.drop, col)
#     }
#   }
# }
# print(paste("Variables to drop = ", length(vars.drop)))
# print(vars.drop)
