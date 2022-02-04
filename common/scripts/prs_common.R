### Source file to load common packages, functions and variables ###


options(warn = 1)


### Load required packages ###

# A package to load packages :)
if (!"librarian" %in% rownames(installed.packages())) {
  install.packages("librarian")
}
library(librarian)
librarian::shelf(
  privefl / bigsnpr, data.table, naniar, optparse, pROC, rms, RColorBrewer,
  tictoc, tidyverse, ukbtools, wrapr
)



### Common variables ###

# The main PEGS directory
PEGS_DIR <- "/ddn/gs1/project/controlled/PEGS"
# The PEGS freeze 1 directory
PEGS_F1_DIR <- "/ddn/gs1/project/controlled/PEGS/Data_Staging_Area/Freeze 1 Surveys with Age at Completion"

# special codes
EPR_NA_STRINGS <- qc(.M, .S, -444444, -555555, -666666, -777777, -888888, -999999)

# TODO: pass this through the slurm file as cmdline arg so that you dont exceed
# the number of cores per node
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
  } else {
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
      # may be better as logical? but how to convert?
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

  # Order of covars does not matter, just put all covars
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
  } else {
    keep <- TRUE
  }
  return(keep)
}


# Define the phenotype in the UK Biobank data (validation data). Identify cases
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

    # Any diabetes
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

    # Initialize phenotype vector
    Y <- rep(0, nrow(ukb_data))
    # Exclude individuals with any diabetes subtype (from controls)
    Y[ind_diabetes] <- NA
    # Mark individuals with T2D as cases
    Y[ind_T2D] <- 1
    # Exclude individuals with T1D (from cases and controls)
    Y[ind_T1D] <- NA

    cat("Phenotype = ", phenotype, "\n")
    print(table(Y, exclude = NULL))

    # Create cases & controls dataframe
    ukb_pheno_prepd <- data.frame(f.eid = ukb_data$f.eid, Y = Y)

    return(ukb_pheno_prepd)
  } else if (phenotype == "asthma") {
    # TODO
    stop("Asthma phenotype code TODO")
    asthma <- "f.22127.0.0"
  } else {
    stop("Unrecognized phenotype string : 'phenotype'")
  }
}


# Define the phenotype in the PEGS/EPR data (test data). Identify cases
# and controls for the specified disease/phenotype based on survey data and
# disease-specific inclusion/exclusion criteria.
# Missing vs skipped responses for child questions are correctly accounted for
# by examining the response for the parent question in the exclusion criteria.
prepare_epr_phenotype <- function(epr_data, phenotype) {

  print("in function prepare_epr_phenotype")
  # Type 2 diabetes (T2D)
  if (phenotype == "T2D") {
    # Exclude participants unlikely to have Type 2 diabetes - i.e. those
    # with gestational diabetes only or those who were diagnosed with diabetes
    # at age < 20 (these are likely Type 1 diabetes ppts)
    # epr_data <- epr_data %>%
    #   mutate(Y = replace(
    #     he_c022_diabetes_PARQ,
    #     (he_c022e_diabetes_age_CHILDQ < 20 | he_c022a_diabetes_preg_CHILDQ == 1),
    #     NA
    #   ))

    diabetes_exclusions <- epr_data %>%
      filter(he_c022_diabetes_PARQ == 1) %>%
      filter(he_c022e_diabetes_age_CHILDQ < 20
             | is.na(he_c022e_diabetes_age_CHILDQ)) %>%
      pull(epr_number)

    diabetes_exclusions <- sort(unique(c(
      diabetes_exclusions,
      (epr_data %>%
        filter(he_c022_diabetes_PARQ == 1
               & .data[["_he_gender_"]] %in% c(1, NA)) %>%
        filter(he_c022a_diabetes_preg_CHILDQ == 1
               | is.na(he_c022a_diabetes_preg_CHILDQ)) %>%
        pull(epr_number)))))

    epr_data <- epr_data %>%
      mutate(Y = replace(
        he_c022_diabetes_PARQ,
        epr_number %in% diabetes_exclusions,
        NA
      ))
  } else if (phenotype == "prediabetes") {
    # Exclude participants with likely different disease etiology, i.e. those
    # who were diagnosed at age < 20
    prediabetes_exclusions <- epr_data %>%
      filter(he_c021_pre_diabetes_PARQ == 1) %>%
      filter(he_c021a_pre_diabetes_age_CHILDQ < 20
             | is.na(he_c021a_pre_diabetes_age_CHILDQ)) %>%
      pull(epr_number)

    epr_data <- epr_data %>%
      mutate(Y = replace(
        he_c021_pre_diabetes_PARQ,
        epr_number %in% prediabetes_exclusions,
        NA
      ))
  } else if (phenotype == "hypertension") {
    # Female ppts with gestational hypertension only are treated as controls
    hypertension_controls <- epr_data %>%
      filter(he_b007_hypertension_PARQ == 1
             & .data[["_he_gender_"]] %in% c(1, NA)) %>%
      filter(he_b007a_hypertension_preg_CHILDQ == 1) %>%
      pull(epr_number)

    epr_data <- epr_data %>%
      mutate(Y = replace(
        he_b007_hypertension_PARQ,
        epr_number %in% hypertension_controls,
        0
      ))
  } else if (phenotype == "hypercholesterolemia") {
    epr_data <- epr_data %>%
      mutate(Y = he_b008_high_cholesterol)
  } else if (phenotype == "asthma") {
    # Exclude participants with COPD
    epr_data <- epr_data %>%
      mutate(Y = replace(
        he_d030_asthma_PARQ,
        he_d025_copd %in% c(1, NA),
        NA
      ))
  } else {
    print("did not match any phenotype")
    stop(paste("Unrecognized phenotype string :"), phenotype)
  }

  # Create cases & controls dataframe
  epr_pheno_prepd <- epr_data %>% select(epr_number, Y)
  cat("Phenotype = ", phenotype, "\n")
  print(table(epr_pheno_prepd$Y, exclude = NULL))

  return(epr_pheno_prepd)
}


# Convert the specified label string from any PEGS metadata file to a named
# character vector containing key-value pairs.
create_label_vector <- function(labels) {

  # Separate the multiple labels in the label string (at ';') and then separate each label into key and value pairs (at '=').
  kv_pairs <- strsplit(unlist(str_split(labels, ";")), "=")

  # Remove unwanted characters from the keys and values
  kv_pairs <- lapply(kv_pairs, function(x) {x <- trimws(x); x <- gsub("(^')|('$)", "", x); return(x)})

  # Convert to a two-column data frame with the keys in one column and the
  # values in the other column.
  kv.df <- do.call(rbind.data.frame, kv_pairs)
  names(kv.df) <- c("key", "value")

  # Convert to a named character vector containing key-value pairs
  label_vector <- kv.df %>%
    pmap(~set_names(..2, ..1)) %>%
    unlist()

  return(label_vector)
}
