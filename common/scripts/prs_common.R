### Source file to load common packages, functions and variables ###


options(warn=1)

### Load required packages ###

# A package to load packages :)
if (! 'librarian' %in% rownames(installed.packages())) {
  install.packages('librarian')
}
library(librarian)
librarian::shelf(privefl/bigsnpr, dplyr, naniar, optparse, pROC, rms, RColorBrewer, tictoc, tidyverse, ukbtools, wrapr)



### Common variables ###

# special codes
EPR_NA_STRINGS <- qc(.M, .S, -444444, -555555, -666666, -777777, -888888, -999999)

# pass this through the slurm file as cmdline arg so that you dont exceed number of cores per node
# Number of cores
NCORES <- nb_cores()
print(paste("num cores available = ", NCORES))
# Don't use all the cores, since we may have multiple tasks per node
NCORES <- min(10, NCORES)
print(paste("using", NCORES, "cores"))
# Set number of cores for BLAS (matrix computations)
#bigparallelr::set_blas_ncores(NCORES)



### Common functions ###

# Function to retrieve command line arguments/options
get_opts <- function(opts_spec) {
  opt_parser = OptionParser(option_list=opts_spec);
  opts = parse_args(opt_parser);
  print("Command-line options:")
  print(opts)

  # Check and complain if required options are missing
  missing_opts <- 0
  sapply(opts_spec, function(x) { if (is.null(opts[[x@dest]])) missing_opts <<- missing_opts + 1})

  if (missing_opts) {
    stop("Required command-line options to script are missing.")
    print_help(opt_parser)
  }
  else {
    return(opts)
  }
}


# Wrapper function to do logistic regression and return formatted output in a dataframe.
# Arguments:-
# glm_df : A dataframe containing the variables in the model. Required columns are: Y = phenotype, PRS = Polygenic Risk Score (covariate of interest). All other columns are treated as covariates in the model.
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
  glm_formula <- as.formula(paste('Y ~', paste(c(covars, 'PRS'), collapse ='+')))
  print("Fitting GLM model:")
  print(glm_formula)

  # Fit a generalized linear model - logistic regression model
  glm_model <- glm(glm_formula, family = "binomial", data = glm_df, singular.ok = FALSE)
  print("GLM model summary:")
  print(summary(glm_model))

  # Output for covariate of interest
  glm_op <- as.data.frame(t(coef(summary(glm_model))["PRS",]))
  colnames(glm_op) <- c('estimate', 'std_error', 'zvalue', 'pvalue')

  # Compute AUC
  predicted_values <- predict.glm(glm_model, type = "response")
  roc_obj <- roc(response = glm_df$Y, predictor = predicted_values)
  glm_op$auc <- as.numeric(auc(roc_obj))

  # Plot ROC curve
  if (!is.null(plot_file)) {
    pdf(plot_file)
    plot(roc_obj, legacy.axes = TRUE, print.thres = "best", print.auc = TRUE)
    plot(smooth(roc_obj, add = TRUE), legacy.axes = TRUE, print.thres = "best", print.auc = TRUE)
    dev.off()
  }

  return(glm_op)
}


# Summary statistics standard deviation check as per https://github.com/privefl/paper-ldpred2/blob/master/paper/paper-ldpred2-supp.pdf.
check_sd <- function(x) {

  if((x['SDss'] < 0.5 * x['SDval']) | (x['SDss'] < 0.5 * x['SDtest']) | (x['SDss'] > 0.1 + x['SDval']) | (x['SDss'] > 0.1 + x['SDtest']) | (x['SDss'] < 0.1) | (x['SDval'] < 0.05) | (x['SDtest'] < 0.05)) {
    keep <- FALSE
  }
  else {
    keep <- TRUE
  }
  return(keep)
}
