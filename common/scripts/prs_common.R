### Source file to load common packages, functions and variables ###


options(warn = 1)


### Load required packages ###

# A package to load packages :)
if (!"librarian" %in% rownames(installed.packages())) {
  install.packages("librarian")
}
library(librarian)
librarian::shelf(
  privefl/bigsnpr, data.table, naniar, optparse, pROC, rms, RColorBrewer,
  tictoc, tidyverse, ukbtools, wrapr
)



### Common variables ###

# The main PEGS directory
PEGS_DIR <- "/ddn/gs1/project/controlled/PEGS"

# special codes
EPR_NA_STRINGS <- qc(.M, .S, -444444, -555555, -666666, -777777, -888888, -999999)

# TODO: pass this through the slurm file as cmdline arg so that you dont exceed
# number of cores per node
# Number of cores
NCORES <- nb_cores()
print(paste("num cores available = ", NCORES))
# Leave one core for the OS/UI
NCORES <- NCORES - 10
print(paste("using", NCORES, "cores"))
# Set number of cores for BLAS (matrix computations)
bigparallelr::set_blas_ncores(NCORES)



### Common functions ###

# Retrieve command line arguments/options
get_opts <- function(opts_spec) {
  opt_parser <- OptionParser(option_list = opts_spec)
  opts <- parse_args(opt_parser)
  print("Command-line options:")
  print(opts)

  # Check and complain if required options are missing
  missing_opts <- 0
  sapply(opts_spec, function(x) {
    if (is.null(opts[[x@dest]])) missing_opts <<- missing_opts + 1
  })

  if (missing_opts) {
    stop("Required command-line options to script are missing.")
    print_help(opt_parser)
  }
  else {
    return(opts)
  }
}


# Convert columns in the EPR dataframe to appropriate variable type as
# specified in the metadata file.
epr_convert_type <- function(epr_df, epr_df_meta) {
  # All special codes must be replaced with NA before type conversion
  epr_df <- epr_df %>% replace_with_na_all(condition = ~ .x %in% EPR_NA_STRINGS)

  for (col in colnames(epr_df)) {
    print(col)
    true_class <- epr_df_meta %>%
      filter(long_variable_name == col) %>%
      pull(true_class)
    print(paste("true class = ", true_class))

    epr_df[[col]] <- switch(true_class,
      "numeric" = as.numeric(epr_df[[col]]),
      # coding binary as factor so that we can check the minimum n per level
      "binary" = as.factor(epr_df[[col]]),
      "character" = as.character(epr_df[[col]]),
      "factor" = as.factor(epr_df[[col]]),
      # coding ordered factor as numeric since these will be treated as such in
      # regression models
      "ordered factor" = as.numeric(epr_df[[col]]),
      epr_df[[col]]
    )
  }
  return(epr_df)
}


# Wrapper function to do logistic regression and return formatted output in a
# dataframe.
# Arguments:-
# glm_df : A dataframe containing the variables in the model.
# Required columns are:
# Y = phenotype,
# PRS = Polygenic Risk Score (covariate of interest).
# All other columns are treated as covariates in the model.
do_glm <- function(glm_df, plot_file = NULL) {
  print("In function do_glm() :")
  print("data frame for glm model:")
  print(str(glm_df))

  # Identify covariate columns.
  # Required: Y = phenotype, PRS = Polygenic Risk Score
  covars <- glm_df %>%
    select(-c(Y, PRS)) %>%
    colnames()

  # order of covars does not matter, just put all covars
  # Model formula
  glm_formula <- as.formula(paste("Y ~", paste(c(covars, "PRS"), collapse = "+")))
  print("Fitting GLM model:")
  print(glm_formula)

  # Fit a generalized linear model - logistic regression model
  glm_model <- glm(glm_formula,
    family = "binomial",
    data = glm_df,
    singular.ok = FALSE
  )
  print("GLM model summary:")
  print(summary(glm_model))

  # Output for covariate of interest
  glm_op <- as.data.frame(t(coef(summary(glm_model))["PRS", ]))
  colnames(glm_op) <- c("estimate", "std_error", "zvalue", "pvalue")

  # Compute AUC
  predicted_values <- predict.glm(glm_model, type = "response")
  roc_obj <- roc(response = glm_df$Y, predictor = predicted_values)
  glm_op$auc <- as.numeric(auc(roc_obj))

  # Plot ROC curve
  if (!is.null(plot_file)) {
    pdf(plot_file)
    plot(roc_obj, legacy.axes = TRUE, print.thres = "best", print.auc = TRUE)
    plot(smooth(roc_obj, add = TRUE),
      legacy.axes = TRUE, print.thres = "best", print.auc = TRUE
    )
    dev.off()
  }

  return(glm_op)
}


# Summary statistics standard deviation check as per
# https://github.com/privefl/paper-ldpred2/blob/master/paper/paper-ldpred2-supp.pdf.
check_sd <- function(x) {
  if ((x["SDss"] < 0.5 * x["SDval"])
  | (x["SDss"] < 0.5 * x["SDtest"])
  | (x["SDss"] > 0.1 + x["SDval"])
  | (x["SDss"] > 0.1 + x["SDtest"])
  | (x["SDss"] < 0.1)
  | (x["SDval"] < 0.05)
  | (x["SDtest"] < 0.05)) {
    keep <- FALSE
  }
  else {
    keep <- TRUE
  }
  return(keep)
}


# Define the phenotype in the validation data (UK Biobank data). Identify cases
# and controls for the specified disease/phenotype based on survey data and
# disease-specific inclusion/exclusion criteria.
prepare_ukb_phenotype <- function(ukb_data, phenotype) {

  ## Extract disease columns
  # TODO: include cancer columns

  # Non-cancer illness code, self-reported
  ukb_illness <- ukb_data[, c(
    paste0("f.20002.0.", 0:33), paste0("f.20002.1.", 0:33),
    paste0("f.20002.2.", 0:33), paste0("f.20002.3.", 0:33)
  )]

  # ICD9 codes
  ukb_ICD9 <- ukb_data[, c(
    paste0("f.41203.0.", 0:27), paste0("f.41205.0.", 0:29),
    paste0("f.41271.0.", 0:46)
  )]

  # ICD10 codes
  ukb_ICD10 <- ukb_data[, c(
    paste0("f.40001.", 0:1, ".0"), paste0("f.40002.0.", 0:13),
    paste0("f.40002.1.", 0:13), paste0("f.41202.0.", 0:65),
    paste0("f.41204.0.", 0:183), paste0("f.41270.0.", 0:212)
  )]


  ## Identify cases and controls for specified phenotype:

  # Type 2 diabetes (T2D)
  if (phenotype == "T2D") {
    # Identify individuals with different subtypes:
    # (gets row numbers in the UK Biobank data for the matches)

    # any diabetes
    ind_diabetes <- sort(unique(unlist(c(
      lapply(ukb_illness, function(x) which(x %in% 1220:1223)),
      lapply(ukb_ICD9, function(x) which(substr(x, 1, 3) == 250)),
      lapply(ukb_ICD10, function(x) which(substr(x, 1, 3) %in% paste0("E", 10:14)))
    ))))

    # Type 1 diabetes (T1D)
    ind_T1D <- sort(unique(unlist(c(
      lapply(ukb_illness, function(x) which(x == 1222)),
      lapply(ukb_ICD9, function(x) which(x %in% c(25001, 25011, 25021, 25091))),
      lapply(ukb_ICD10, function(x) which(substr(x, 1, 3) == "E10"))
    ))))

    # Type 2 diabetes
    ind_T2D <- sort(unique(unlist(c(
      lapply(ukb_illness, function(x) which(x == 1223)),
      lapply(ukb_ICD9, function(x) which(x %in% c(25000, 25010, 25020, 25090))),
      lapply(ukb_ICD10, function(x) which(substr(x, 1, 3) == "E11"))
    ))))

    print(paste("#ppts with any type of diabetes =", length(ind_diabetes)))
    print(paste("#ppts with T1D =", length(ind_T1D)))
    print(paste("#ppts with T2D =", length(ind_T2D)))

    ## Mark cases and controls

    # initialize phenotype vector
    Y <- rep(0, nrow(ukb_data))
    # exclude individuals with any diabetes subtype (from controls)
    Y[ind_diabetes] <- NA
    # mark individuals with T2D as cases
    Y[ind_T2D] <- 1
    # exclude individuals with T1D (from cases and controls)
    Y[ind_T1D] <- NA

    cat("Phenotype = ", phenotype, "\n")
    print(table(Y, exclude = NULL))

    # Create cases & controls dataframe
    ukb_pheno_prepd <- data.frame(f.eid = ukb_data$f.eid, Y = Y)

    return(ukb_pheno_prepd)
  }
  else {
    stop("Unrecognized phenotype string : 'phenotype'")
  }
}
