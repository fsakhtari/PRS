### Risk score analyses ###


### init ###

#
# Load packages and variables
#

# The common PRS directory for all diseases
PRS_COMMON_DIR <- Sys.getenv("PRS_COMMON_DIR")
# The PRS directory of the specific phenotype for which PRS is being computed
PRS_PHENO_DIR <- Sys.getenv("PRS_PHENO_DIR")
source(paste0(PRS_COMMON_DIR, "/scripts/prs_common.R"))

librarian::shelf(
  car, caret, coefplot, corrplot, DataExplorer, e1071, glmnet,
  interactions, nricens, stringr
)

# WGS data
plink_bfileset <- "/ddn/gs1/project/nextgen/controlled/wgsepr/Contolled/plink/final_4768_confirmed_filtered_nocontrols_nooutliers_unrelated"


# Retrieve command line arguments/options
opts_spec <- list(
  make_option(c("--prediabetes"),
    type = "logical",
    default = TRUE,
    help = "Include prediabetes in CRS computation",
    metavar = "character"
  ),
  make_option(c("--reuse_pcs"),
    type = "logical",
    default = FALSE,
    help = "Reuse existing PCs",
    metavar = "character"
  )
)
opts <- get_opts(opts_spec)

# Output dir
out_dir <- paste0(PRS_PHENO_DIR, "/risk_scores_tt50")

script_info <- c(
  "Running risk_scores_tt50.R :", "\n",
  "Uses significant variables from ExWAS on nonWGS ppts for H&E survey only", "\n",
  "Computes weights for risk scores on training dataset", "\n",
  "Predicts risk on testing dataset", "\n"
)

if (!opts$prediabetes) {
  out_dir <- paste0(out_dir, "_nopd")
  script_info <- c(script_info, "CRS computed without prediabetes", "\n")
}

cat(script_info)
cat("Output dir = ", out_dir, "\n")

#
# Function defns
#

# Compute genomic principal components for the specified individuals/samples.
# This function also does the QC & filtering of the SNPs to ensure that both the
# SNPs and the samples used for PCA are the same as those used for analysis. The
# PCA output is saved in the specified output directory with plink's default
# names, i.e. - plink.eigenval and plink.eigenvec
# Arguments:
# bfile = plink binary fileset
# pca_dir = directory for  output files
# indivs_keep = file containing the samples for whom the PCs need to be computed
# (in the same format as required by plink's --keep option)
compute_pcs <- function(bfile, pca_dir, indivs_keep) {
  plink_cmd <- "plink"
  plink_eigenvec_out <- paste0(pca_dir, "/plink.eigenvec")

  if (file.exists(plink_eigenvec_out) & opts$reuse_pcs) {
    print(sprintf(
      "Skipping PCA because %s exists and 'reuse_pcs' flag is set",
      plink_eigenvec_out
    ))
    return()
  } else {
    print("Computing PCs")
    # Filter on --maf and --geno before filtering on --mind (samples) since some
    # sample level QC was done during the WGS QC.
    system2(plink_cmd,
      args = paste(
        "--bfile", bfile,
        "--keep", indivs_keep,
        "--maf 0.01",
        "--geno 0.05",
        "--make-bed",
        "--out", paste0(pca_dir, "/rs_1")
      )
    )
    system2(plink_cmd,
      args = paste(
        "--bfile", paste0(pca_dir, "/rs_1"),
        "--mind 0.05",
        "--hwe 1e-10", # ???
        "--make-bed",
        "--out", paste0(pca_dir, "/rs_2")
      )
    )

    # LD pruning
    system2(plink_cmd,
      args = paste(
        "--bfile", paste0(pca_dir, "/rs_2"),
        "--chr 1-22",
        "--indep-pairwise 50 5 0.5",
        "--threads 20",
        "--out", paste0(pca_dir, "/plink")
      )
    )

    # Compute PCs for the specified samples using the LD-pruned SNPs
    system2(plink_cmd,
      args = paste(
        "--bfile", paste0(pca_dir, "/rs_2"),
        "--extract", paste0(pca_dir, "/plink.prune.in"),
        "--pca 10",
        "--out", paste0(pca_dir, "/plink")
      )
    )

    # Write column names to plink's eigenvector output
    pcs <- read.table(plink_eigenvec_out, header = FALSE, quote = "")
    colnames(pcs) <- c("epr_number", "IID", paste0("PC", 1:(ncol(pcs) - 2)))
    write.table(pcs, plink_eigenvec_out, row.names = FALSE, quote = FALSE)

    return()
  }
}

# Data QC & filtering
filter_data <- function(unfiltered_data) {
  print("unfiltered data:")
  cat("dim:", dim(unfiltered_data), "\n")
  print("columns:")
  print(colnames(unfiltered_data))
  str(unfiltered_data)

  # Drop variables with low sample sizes
  vars.drop <- NULL
  # Drop any variable with n < 30
  for (col in colnames(unfiltered_data)) {
    cc <- sum(!is.na(unfiltered_data[, col]))
    if (cc < 30) {
      vars.drop <- c(vars.drop, col)
    }
  }
  print("Following variables have n < 30:")
  print(vars.drop)

  # For categorical/factor variables:
  cat_vars <- unfiltered_data %>%
    select(where(is.factor)) %>%
    select(-Y) %>%
    colnames()
  print("Factor variables:")
  print(cat_vars)

  for (col in cat_vars) {
    # Drop levels with n < 10, i.e. nullify (set to NA) values with n < 10 at
    # any level
    low_cnt_lvls <- unfiltered_data %>%
      count(.data[[col]]) %>%
      drop_na() %>%
      filter(n < 10) %>%
      pull(col)
    if (!is_empty(low_cnt_lvls)) {
      print(paste("variable", col, "has n < 10 for levels:", low_cnt_lvls))
      cat("Dropping/nullifying levels:", low_cnt_lvls, "\n")
      # factor values/levels to nullify
      na_levels <- setNames(
        rep(NA_character_, length(low_cnt_lvls)),
        low_cnt_lvls
      )
      print(na_levels)
      unfiltered_data[, col] <- recode_factor(
        unlist(unfiltered_data[, col]),
        !!!na_levels
      )
    }

    # Drop variables with low cell counts in crosstabs with phenotype (Y)
    # because cells with small counts make the logistic regression model
    # unstable.
    low_cnt_lvls <- unfiltered_data %>%
      group_by(Y, .data[[col]], .drop = FALSE) %>%
      tally() %>%
      filter(n < 1) %>%  # should do n < 5
      pull(col)
    if (!is_empty(low_cnt_lvls)) {
      print(paste(
        "variable", col, "has n < 5 in crosstab with Y for levels:",
        low_cnt_lvls
      ))
      cat("Dropping/nullifying levels:", low_cnt_lvls, "\n")
      # factor values/levels to nullify
      na_levels <- setNames(
        rep(NA_character_, length(low_cnt_lvls)),
        low_cnt_lvls
      )
      print(na_levels)
      unfiltered_data[, col] <- recode_factor(
        unlist(unfiltered_data[, col]),
        !!!na_levels
      )
    }

    # Now that we have dropped levels for some variables, we need to drop any
    # categorical/factor variables with levels < 2
    nlevels <- unfiltered_data %>%
      count(.data[[col]]) %>%
      drop_na() %>%
      nrow()
    if (nlevels < 2) {
      print(paste("variable", col, "has < 2 levels"))
      vars.drop <- c(vars.drop, col)
      next
    }
  }

  print(paste("Dropping following", length(vars.drop), "variables"))
  print(vars.drop)
  write(vars.drop, paste0(out_dir, "/vars.drop"))

  filtered_data <- unfiltered_data %>% select(-all_of(vars.drop))
  write(colnames(filtered_data), paste0(out_dir, "/vars.keep"))

  print("filtered data:")
  cat("dim:", dim(filtered_data), "\n")
  print("columns:")
  print(colnames(filtered_data))

  return(filtered_data)
}

# Make PCA biplots
pca_biplot <- function(pc_data, pca, pcb) {
  print(
    ggplot(
      pc_data,
      aes_string(x = pca, y = pcb, color = "race", shape = "ethnicity")
    ) +
      geom_point() +
      xlab(pca) +
      ylab(pcb)
  )
}

# Create model/design matrix
create_model_matrix <- function(model_df) {
  model_df <- model_df %>% select(-epr_number)

  # Remove the intercept column since the glmnet model fits an intercept.
  X <- model.matrix(Y ~ ., model_df)[, -1]
  Y <- model_df %>% pull(Y)

  print("model matrix:")
  cat("dimensions:", dim(X), "\n")
  print(head(X))

  return(list(X = X, Y = Y))
}

# Variable selection
do_variable_selection <- function(model_df, desc, plot = FALSE) {
  mm <- create_model_matrix(model_df)
  X <- mm$X
  Y <- mm$Y

  cat("Fitting glmnet for", desc, "\n")
  set.seed(123)
  cvfit <- cv.glmnet(X, Y,
    type.measure = "auc",
    nfolds = 10,
    family = "binomial",
    alpha = 1,
    standardize = TRUE,
  )

  print(sprintf("cvfit for %s:", desc))
  print(cvfit)
  plot(cvfit, main = desc)

  lambda_min <- cvfit$lambda.min
  lambda_1se <- cvfit$lambda.1se
  cat("lambda.min = ", lambda_min, "\n")
  cat("lambda.1se = ", lambda_1se, "\n")
  cat("coefficients at lambda = :", lambda_min, lambda_1se, "\n")
  print(coef(cvfit, s = c(lambda_min, lambda_1se)))

  # Obtain variables & their coefficients for non-zero, non-control variables
  coeffs <- coef(cvfit, s = "lambda.min")
  coeffs_df <- data.frame(variable = rownames(coeffs), beta = coeffs[, 1])
  rownames(coeffs_df) <- NULL

  exclude_vars <- c("(Intercept)", "gender1", control_vars, pc_vars)
  nonzero_coeffs <- coeffs_df %>%
    filter(!variable %in% exclude_vars) %>%
    filter(beta != 0)

  print("Variables of interest with non-zero coefficients:")
  print(nonzero_coeffs)

  if (plot == TRUE) {
    # TODO: These plots are with the coefficients on the original scale and
    # hence the coefficients are not comparable to each other. Need to convert
    # the coefficients to the standardized scale to plot and compare.
    print(coefplot(cvfit,
      lambda = "lambda.min",
      sort = "magnitude",
      coefficients = nonzero_coeffs$variable,
      ylab = "Estimate",
      xlab = "Variable",
      title = desc
    ))

    print(ggplot(
      data = nonzero_coeffs,
      aes(x = variable, y = beta)
    ) +
      geom_col() +
      coord_flip() +
      ylab("Estimate") +
      xlab("Variable") +
      ggtitle(desc) +
      theme_classic())
  }

  return(nonzero_coeffs)
}

# Calculate risk scores from precomputed betas/weights
calc_risk_scores <- function(model_df, coeffs_df) {
  mm <- create_model_matrix(model_df)
  X <- mm$X

  nonzero_coeffs <- coeffs_df

  # Compute risk scores using the coefficients for the variables
  risk_scores <- X[, nonzero_coeffs$variable] %*% nonzero_coeffs$beta
  risk_score_df <- data.frame(
    epr_number = model_df$epr_number,
    risk_score = risk_scores
  )

  return(risk_score_df)
}

# Compute risk scores
compute_risk_scores <- function(model_df, desc, plot = FALSE) {
  cat("redoing model matrix for:", desc, "\n")
  mm <- create_model_matrix(model_df)
  X <- mm$X

  # Get coefficients for each variable
  nonzero_coeffs <- do_variable_selection(
    model_df = model_df, desc = desc, plot = plot
  )
  write.table(nonzero_coeffs, paste0(out_dir, "/coeffs_", make.names(desc)),
    quote = FALSE, row.names = FALSE
  )

  # Compute risk scores using the coefficients for the variables
  risk_scores <- X[, nonzero_coeffs$variable] %*% nonzero_coeffs$beta
  risk_score_df <- data.frame(
    epr_number = model_df$epr_number,
    risk_score = risk_scores
  )

  return(risk_score_df)
}

# Summarize and plot data characteristics
summarize_data <- function(df, plotfile, info_string = NULL) {
  pdf(plotfile)

  print(sprintf("== Data summaries for %s ==", info_string))

  cat("dimensions:", dim(df), "\n")
  print("columns:")
  print(colnames(df))
  print("variable type:")
  str(df)

  df %>%
    group_by(Y) %>%
    tally() %>%
    mutate(percent = prop.table(n) * 100) %>%
    print()

  df %>%
    group_by(gender) %>%
    tally() %>%
    print()

  cont_vars <- df %>%
    select_if(is.numeric) %>%
    colnames()
  print("Continuous variables:")
  print(cont_vars)
  cat_vars <- df %>%
    select_if(is.factor) %>%
    colnames()
  print("Categorical variables:")
  print(cat_vars)

  print("Contingency tables for all categorical variables by Y :")
  for (cvar in cat_vars) {
    df %>%
      select(all_of(c("Y", cvar))) %>%
      table() %>%
      addmargins() %>%
      print()
  }

  # Histograms for key continuous variables:
  # Age
  age_mean <- round(mean(df$he_age_derived), 1)
  print("Summary statistics for age at H&E survey completion:")
  cat("mean = ", age_mean, "\n")
  cat("median = ", median(df$he_age_derived), "\n")
  cat("std.dev = ", sd(df$he_age_derived), "\n")
  cat("range = ", range(df$he_age_derived), "\n")

  ggplot(df, aes(x = he_age_derived)) +
    geom_histogram(
      binwidth = 10, fill = "lightblue", color = "lightblue", alpha = 0.7
    ) +
    geom_vline(xintercept = age_mean, color = "steelblue") +
    geom_text(aes(label = paste("mean =", age_mean), y = -3, x = age_mean),
      vjust = "outward", col = "steelblue", size = 3
    ) +
    theme_classic() +
    xlab("Age (at H&E completion) in years") +
    ylab("Count")

  # BMI
  print("Summary statistics for BMI:")
  bmi_mean <- round(mean(df$he_bmi_derived), 1)
  cat("mean = ", bmi_mean, "\n")
  cat("median = ", median(df$he_bmi_derived), "\n")
  cat("std.dev = ", sd(df$he_bmi_derived), "\n")
  cat("range = ", range(df$he_bmi_derived), "\n")

  ggplot(df, aes(x = he_bmi_derived)) +
    geom_histogram(
      binwidth = 7, fill = "lightblue", color = "lightblue", alpha = 0.7
    ) +
    geom_vline(xintercept = bmi_mean, color = "steelblue") +
    geom_text(aes(label = paste("mean =", bmi_mean), y = -3, x = bmi_mean),
      vjust = "outward", col = "steelblue", size = 3
    ) +
    theme_classic() +
    xlab("Body mass index (BMI)") +
    ylab("Count")

  # Correlation plots:
  # Convert to numeric to compute correlation
  all_factor_cols <- df %>%
    select_if(is.factor) %>%
    colnames()
  all_data_numeric <- df %>%
    select(-c(epr_number, broad_wgs_sample_id_CHILDQ)) %>%
    mutate_at(all_of(all_factor_cols), ~ as.numeric(levels(.))[.])

  # correlation amongst all data
  all_cor <- cor(all_data_numeric, method = "spearman")
  corrplot.mixed(all_cor,
    lower.col = "black", tl.pos = "lt", tl.cex = 0.6, number.cex = 0.6,
    main = info_string
  )
  corrplot(all_cor, tl.cex = 0.6, order = "hclust", addrect = 3, main = info_string)

  cortest <- cor.mtest(all_data_numeric, conf.level = 0.95)
  corrplot.mixed(all_cor,
    p.mat = cortest$p, sig.level = 0.05, number.cex = 0.4,
    insig = "blank", tl.cex = 0.6, lower.col = "black", tl.pos = "lt", main =
      paste("Significant correlations in", info_string)
  )

  # PCA plots

  # PC correlations
  pc_cor <- cor(df %>% select(starts_with("PC")))
  corrplot.mixed(pc_cor, lower.col = "black", tl.pos = "lt", main = "PCs")

  # PCA biplots
  pca_plot_data <- inner_join(df, epr.bcbb.map.conv, by = "epr_number") %>%
    select(starts_with("PC"), race, ethnicity, gender.x)
  for (i in 1:4) {
    pca_biplot(pca_plot_data, pca = paste0("PC", i), pcb = paste0("PC", i + 1))
  }
  print("Race and ethnicity in final complete cases data:")
  print(table(pca_plot_data$race))
  print(table(pca_plot_data$ethnicity))

  dev.off()
}


# Association & prediction of risk scores with phenotype for the given data
# train_data: data frame; training dataset
# test_data : data frame; testing dataset
# rs_outdir: character; The output directory to store risk score results in
test_rs_assoc <- function(train_data, test_data, rs_outdir) {
  cat("== Association & prediction of risk scores with phenotype ==", "\n")

  # Models specified as: term of interest = expression in glm formula
  rs_models <- c(
    "CRS" = "CRS",
    "ERS" = "ERS",
    "PRS" = "PRS",
    "ORS" = "ORS",
    "ERS" = "CRS+ERS",
    "PRS" = "CRS+PRS",
    "ERS" = "CRS+ERS+PRS",
    "CRS:ERS" = "CRS*ERS",
    "CRS:PRS" = "CRS*PRS",
    "ERS:PRS" = "ERS*PRS"
  )

  strata_control_vars <- control_vars

  # When stratifying by gender, remove the gender covariate from the control
  # variables
  if (str_detect(rs_outdir, "(?i)male$")) {
    strata_control_vars <- setdiff(control_vars, "gender")
  }

  # Results for each risk score
  rs_results <- data.frame()

  # Predicted risk of disease for all individuals from each risk score model
  indiv_pred_risk <- data.frame(
    matrix(
      data = NA,
      nrow = nrow(test_data),
      ncol = length(rs_models),
      dimnames = list(c(), rs_models)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  pred_risks_clin <- NA

  for (i in seq_along(rs_models)) {
    model_term <- unname(rs_models[i])
    cat("Modeling", model_term, "\n")
    toi <- names(rs_models)[i]
    cat("Model term of interest is", toi, "\n")
    cat("Risk score results dir is", rs_outdir, "\n")

    print("Fitting GLM model:")
    glm_formula <- as.formula(paste(
      "Y ~",
      paste(c(strata_control_vars, paste(pc_vars, collapse = "+"), model_term),
        collapse = "+"
      )
    ))
    print(glm_formula)

    print("train data GLM model summary")
    train_glm_model <- glm(
      glm_formula,
      family = "binomial",
      data = train_data,
      singular.ok = FALSE
    )
    print(summary(train_glm_model))
    train_glm_op <- coef(summary(train_glm_model)) %>% as.data.frame()
    train_or_ci <- exp(confint(train_glm_model)[toi, ])

    print("test data GLM model summary")
    test_glm_model <- glm(
      glm_formula,
      family = "binomial",
      data = test_data,
      singular.ok = FALSE
    )
    print(summary(test_glm_model))
    test_glm_op <- coef(summary(test_glm_model)) %>% as.data.frame()
    test_or_ci <- exp(confint(test_glm_model)[toi, ])

    print("Prediction on test data")
    predicted_values <- predict.glm(
      train_glm_model,
      type = "response",
      newdata = test_data
    )
    indiv_pred_risk[, model_term] <- predicted_values
    roc_obj <- roc(response = test_data$Y, predictor = predicted_values)
    auc_ci <- as.numeric(ci(roc_obj, of = "auc"))
    cat("AUC = ", auc_ci[2], "\n")

    print("Calculating pseudo-R2:")
    print(sprintf("Full model (with %s):", model_term))
    lrm_full <- lrm(glm_formula, data = test_data)
    print(lrm_full)

    print(sprintf("Reduced model (without %s):", model_term))
    lrm_reduced <- lrm(
      as.formula(paste(
        "Y ~",
        paste(c(strata_control_vars, paste(pc_vars, collapse = "+")),
          collapse = "+"
        )
      )),
      data = test_data
    )
    print(lrm_reduced)

    varexp_rs <- lrm_full$stats[["R2"]] - lrm_reduced$stats[["R2"]]
    print(paste(
      "Proportion of variance explained by", model_term,
      "using Nagelkerke’s pseudo-R2 metric =", varexp_rs
    ))

    if (model_term == "CRS") {
      # Use the individual predicted risk from the CRS model as baseline to
      # compute NRI
      pred_risks_clin <- indiv_pred_risk$CRS
      names(pred_risks_clin) <- seq_len(length(pred_risks_clin))
      nri_cont_val <- NA
    } else {
      print("Calculating NRI")
      # Continuous NRI
      nri_cont <- nribin(
        event = as.numeric(levels(test_data$Y))[test_data$Y],
        p.std = pred_risks_clin,
        p.new = predicted_values,
        updown = "diff",
        cut = 0,
        niter = 0
      )
      nri_cont_val <- nri_cont$nri["NRI", "Estimate"]
    }

    rs_res <- data.frame(
      risk_score = model_term,
      toi = toi,
      train_OR = exp(train_glm_op[toi, "Estimate"]),
      train_OR_L95 = train_or_ci[["2.5 %"]],
      train_OR_U95 = train_or_ci[["97.5 %"]],
      train_pvalue = train_glm_op[toi, "Pr(>|z|)"],
      test_OR = exp(test_glm_op[toi, "Estimate"]),
      test_OR_L95 = test_or_ci[["2.5 %"]],
      test_OR_U95 = test_or_ci[["97.5 %"]],
      test_pvalue = test_glm_op[toi, "Pr(>|z|)"],
      test_AUC = auc_ci[2],
      test_AUC_L95 = auc_ci[1],
      test_AUC_U95 = auc_ci[3],
      test_pseudo_R2 = varexp_rs,
      test_NRI_cont = nri_cont_val
    )
    rs_results <- rbind(rs_results, rs_res)
  }

  write.table(rs_results, paste0(rs_outdir, "/rs_results.txt"),
    quote = FALSE, row.names = FALSE
  )
  write.table(indiv_pred_risk, paste0(rs_outdir, "/indiv_predicted_risk.txt"),
    quote = FALSE, row.names = FALSE
  )
}


# Make risk score plots
# test_data : data frame; data to plot
# rs_outdir: character; The output directory to store risk score plots in
plot_rs_results <- function(train_data, test_data, rs_outdir) {
  cat("== Plotting risk score results ==", "\n")

  # Output file for plots
  pdf(paste0(rs_outdir, "/risk_scores.pdf"))

  rs_results <- read.table(paste0(rs_outdir, "/rs_results.txt"), header = TRUE) %>%
    mutate(model_term = paste0(toi, " (", risk_score, ")")) %>%
    arrange(train_OR) %>%
    mutate(model_term = factor(model_term, levels = model_term))

  # Forest plots for risk score odds ratios
  # For training data
  print(ggplot(
    data = rs_results,
    aes(x = model_term, y = train_OR, ymin = train_OR_L95, ymax = train_OR_U95)
  ) +
    geom_pointrange() +
    geom_hline(yintercept = 1, lty = 2) +
    coord_flip() +
    xlab("Term of interest (Model)") +
    ylab("95% CI for Odds ratio") +
    ggtitle("OR CI for Risk scores (Training data)") +
    theme_bw())

  # Combined plot for training & test data
  train_or <- rs_results %>%
    select(model_term, train_OR, train_OR_L95, train_OR_U95) %>%
    rename_with(~ gsub("^train_", "", .x))
  test_or <- rs_results %>%
    select(model_term, test_OR, test_OR_L95, test_OR_U95) %>%
    rename_with(~ gsub("^test_", "", .x))

  print(ggplot(
    data = train_or,
    aes(x = model_term, y = OR, ymin = OR_L95, ymax = OR_U95, colour = "black")
  ) +
    geom_pointrange() +
    geom_pointrange(data = test_or, position = "jitter", aes(colour = "red")) +
    geom_hline(yintercept = 1, lty = 2) +
    scale_color_identity(name = "Dataset",
                         breaks = c("black", "red"),
                         labels = c("Training", "Test"),
                         guide = "legend") +
    coord_flip() +
    xlab("Term of interest (Model)") +
    ylab("95% confidence interval for risk score odds ratios") +
    theme_bw()
  )

  # Percentile plots
  for (rs in rs_names) {
    prev_ptile <- test_data %>%
      select(Y, !!rs) %>%
      mutate(percentile = ntile(as.numeric(.data[[rs]]), 100)) %>%
      group_by(percentile) %>%
      summarise(prev = mean(as.numeric(levels(Y))[Y]))

    print(
      ggplot(
        data = prev_ptile,
        aes(x = percentile, y = prev * 100, colour = percentile)
      ) +
        geom_point() +
        scale_color_gradientn(colours = brewer.pal(n = 9, name = "Blues")) +
        ylab("T2D prevalence %") +
        xlab(paste(rs, "percentile")) +
        theme_classic() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          axis.text = element_text(face = "bold", size = 14),
          axis.title = element_text(face = "bold", size = 14)
        )
    )
  }

  # Decile plots
  # Compute prevalence in every decile for each risk score
  prev_decile <- data.frame(decile = 1:10)
  for (rs in rs_names) {
    pd <- test_data %>%
      select(Y, !!rs) %>%
      mutate(decile = ntile(as.numeric(.data[[rs]]), 10)) %>%
      group_by(decile) %>%
      summarise(!!rs := mean(as.numeric(levels(Y))[Y]))

    prev_decile <- merge(prev_decile, pd, by = "decile")
  }

  # Convert to long format for ggplot
  prev_decile_long <- melt(prev_decile, id.vars = "decile", value.name = "prev")

  print(ggplot(
    data = prev_decile_long,
    aes(x = decile, y = prev * 100, colour = variable)
  ) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = c(1:10)) +
    scale_colour_brewer(palette = "Dark2") +
    ylab("T2D prevalence %") +
    xlab("Risk score decile") +
    theme_classic())

  # Correlation plots

  # TODO: For the plots containing factor or ordinal variables, instead of
  # converting to numeric, consider doing polychoric, etc. types of
  # correlations - See GWAS pipeline code.

  # with Y
  cor_y <- test_data %>%
    select(Y, CRS, ERS, PRS, ORS) %>%
    mutate(Y = as.numeric(levels(Y))[Y]) %>%
    cor(method = "spearman")
  corrplot.mixed(cor_y, lower.col = "black")

  # with Y, PCs
  cor_y_pcs <- test_data %>%
    select(Y, ends_with("RS", ignore.case = FALSE), starts_with("PC")) %>%
    mutate(Y = as.numeric(levels(Y))[Y]) %>%
    cor(method = "spearman")
  corrplot.mixed(cor_y_pcs, tl.pos = "lt", tl.cex = 0.7, lower.col = "black")

  # Plots relevant for the entire dataset only
  # if (str_detect(rs_outdir, "skipping for now")) {
  if (str_detect(rs_outdir, "/all")) {

    # correlations with Y, control variables
    cor_y_cntl <- test_data %>%
      # Removing this race variable which contains strings
      select(-race) %>%
      inner_join(epr.bcbb.map.conv[, qc(epr_number, race, ethnicity)],
        by = "epr_number"
      ) %>%
      select(!!reqd_vars, !!control_vars, race, ethnicity, CRS, ERS, PRS, ORS) %>%
      select(-epr_number) %>%
      mutate_if(is.factor, ~ as.numeric(levels(.))[.]) %>%
      cor(method = "spearman", use = "pairwise.complete.obs")
    corrplot.mixed(cor_y_cntl, tl.pos = "lt", lower.col = "black")

    # Risk scores stratified by race, for all races
    rs_race <- test_data %>%
      select(ends_with("RS", ignore.case = FALSE), race) %>%
      group_by(race) %>%
      mutate(n = n()) %>%
      mutate(label = paste(race, "\nN =", n))

    for (rs in rs_names) {
      print(
        ggplot(
          data = rs_race,
          aes(
            x = reorder(label, .data[[rs]], FUN = median),
            y = .data[[rs]],
            fill = label
          )
        ) +
          geom_boxplot() +
          scale_fill_brewer(palette = "Dark2") +
          ylab(rs) +
          xlab("Race") +
          coord_flip() +
          theme_classic() +
          theme(legend.position = "none")
      )
    }

    # Risk scores stratified by race for Black and White races
    rs_race_bw <- test_data %>%
      filter(race %in% c("Black or African American", "White")) %>%
      mutate(race = recode(race, "Black or African American" = "Black")) %>%
      select(ends_with("RS", ignore.case = FALSE), race) %>%
      melt(id.vars = "race")

    print(
      ggplot(
      data = rs_race_bw, aes(x = race, y = value, fill = race)
      ) +
      facet_grid(. ~ variable) +
      geom_boxplot() +
      scale_fill_brewer(palette = "Paired", direction = -1) +
      ylab("Risk score") +
      xlab("Race") +
      theme_classic() +
      guides(fill = guide_legend(title = "Race"))
    )

    # Risk scores stratified by gender
    rs_gender <- test_data %>%
      select(ends_with("RS", ignore.case = FALSE), gender) %>%
      melt(id.vars = "gender")

    print(
      ggplot(
        data = rs_gender, aes(x = gender, y = value, fill = gender)
      ) +
        facet_grid(. ~ variable) +
        geom_boxplot() +
        scale_fill_brewer(palette = "Paired", direction = -1) +
        ylab("Risk score") +
        xlab("Gender") +
        theme_classic() +
        guides(fill = guide_legend(title = "Gender"))
    )
  }

  dev.off()
}

# Do "RS"+E analyses for the given data
# test_data : data frame; testing dataset
# rs_term : character; The risk score to use
# rs_outdir: character; The output directory to store risk score results in
test_rs_e <- function(test_data, rs_term, rs_outdir) {
  cat("== Testing", rs_term, "+ E ==", "\n")

  strata_control_vars <- control_vars
  # When stratifying by gender, remove the gender covariate from the control
  # variables
  if (str_detect(rs_outdir, "(?i)male$")) {
    strata_control_vars <- setdiff(control_vars, "gender")
  }

  all_model_results <- all_exp_results <- NULL

  for (sev in ers_vars) {
    print("Fitting GLM model:")
    glm_formula <- as.formula(paste(
      "Y ~",
      paste(c(
        strata_control_vars,
        paste(pc_vars, collapse = "+"),
        rs_term,
        sev
      ),
      collapse = "+"
      )
    ))
    print(glm_formula)

    glm_model <- glm(
      glm_formula,
      family = "binomial",
      data = test_data,
      singular.ok = FALSE
    )
    print("GLM model summary")
    print(summary(glm_model))

    glm_op <- coef(summary(glm_model)) %>% as.data.frame()

    print("Anova table")
    anova_out <- Anova(glm_model, singular.ok = FALSE)
    print(anova_out)

    print("Get model statistics")
    lrm_model <- lrm(glm_formula, data = test_data)
    print(lrm_model)

    model_res <- data.frame(
      exposure = sev,
      R2 = lrm_model$stats["R2"],
      N = lrm_model$stats["Obs"],
      N_controls = lrm_model$freq[["0"]],
      N_cases = lrm_model$freq[["1"]]
      # TODO
      # AUC = auc_ci[2],
      # AUC_L95 = auc_ci[1],
      # AUC_U95 = auc_ci[3],
    )
    all_model_results <- rbind(all_model_results, model_res)

    exp_res <- data.frame(
      variable = sev,
      pvalue = anova_out[sev, "Pr(>Chisq)"]
    )
    all_exp_results <- rbind(all_exp_results, exp_res)
  }

  write.table(
    all_model_results,
    paste0(rs_outdir, "/", rs_term, "_e_models.txt"),
    row.names = FALSE,
    quote = TRUE
  )

  all_exp_results$pvalue_adj <- p.adjust(all_exp_results$pvalue, method = "fdr")
  all_exp_results <- all_exp_results %>% arrange(pvalue_adj)
  write.table(all_exp_results, paste0(rs_outdir, "/", rs_term, "_e_exp.txt"),
    row.names = FALSE, quote = TRUE
  )
}

# Association & prediction of risk scores with prediabetes for the given data
# train_data: data frame; training dataset
# test_data : data frame; testing dataset
# rs_outdir: character; The output directory to store risk score results in
test_rs_pd <- function(train_data, test_data, rs_outdir) {
  cat("== Association & prediction of risk scores with prediabetes ==", "\n")

  # Models specified as: term of interest = expression in glm formula
  rs_models <- c(
    # "CRS" = "CRS",
    "ERS" = "ERS"
    # "ERS" = "CRS+ERS",
    # "CRS:ERS" = "CRS*ERS",
  )

  strata_control_vars <- control_vars

  # When stratifying by gender, remove the gender covariate from the control
  # variables
  if (str_detect(rs_outdir, "(?i)male$")) {
    strata_control_vars <- setdiff(control_vars, "gender")
  }

  # Association testing results for each risk score
  rs_results <- data.frame()

  # Predicted risk of disease for all individuals from each risk score model
  indiv_pred_risk <- data.frame(
    matrix(
      data = NA,
      nrow = nrow(test_data),
      ncol = length(rs_models),
      dimnames = list(c(), rs_models)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  pred_risks_clin <- NA

  for (i in seq_along(rs_models)) {
    model_term <- unname(rs_models[i])
    cat("Modeling", model_term, "\n")
    toi <- names(rs_models)[i]
    cat("Model term of interest is", toi, "\n")
    cat("Risk score results dir is", rs_outdir, "\n")

    glm_formula <- as.formula(paste(
      "he_c021_pre_diabetes_PARQ ~",
      paste(c(strata_control_vars, paste(pc_vars, collapse = "+"), model_term),
        collapse = "+"
      )
    ))
    print("Fitting GLM model:")
    print(glm_formula)

    print("train data GLM model summary")
    train_glm_model <- glm(
      glm_formula,
      family = "binomial",
      data = train_data,
      singular.ok = FALSE
    )
    print(summary(train_glm_model))
    train_glm_op <- coef(summary(train_glm_model)) %>% as.data.frame()
    train_or_ci <- exp(confint(train_glm_model)[toi, ])

    print("test data GLM model summary")
    test_glm_model <- glm(
      glm_formula,
      family = "binomial",
      data = test_data,
      singular.ok = FALSE
    )
    print(summary(test_glm_model))
    test_glm_op <- coef(summary(test_glm_model)) %>% as.data.frame()
    test_or_ci <- exp(confint(test_glm_model)[toi, ])

    print("Prediction on test data")
    predicted_values <- predict.glm(
      train_glm_model,
      type = "response",
      newdata = test_data
    )
    indiv_pred_risk[, model_term] <- predicted_values
    roc_obj <- roc(response = test_data$Y, predictor = predicted_values)
    auc_ci <- as.numeric(ci(roc_obj, of = "auc"))
    cat("AUC = ", auc_ci[2], "\n")

    print("Calculating pseudo-R2:")
    print(sprintf("Full model (with %s):", model_term))
    lrm_full <- lrm(glm_formula, data = test_data)
    print(lrm_full)

    print(sprintf("Reduced model (without %s):", model_term))
    lrm_reduced <- lrm(
      as.formula(paste(
        "he_c021_pre_diabetes_PARQ ~",
        paste(c(strata_control_vars, paste(pc_vars, collapse = "+")),
          collapse = "+"
        )
      )),
      data = test_data
    )
    print(lrm_reduced)

    varexp_rs <- lrm_full$stats[["R2"]] - lrm_reduced$stats[["R2"]]
    print(paste(
      "Proportion of variance explained by", model_term,
      "using Nagelkerke’s pseudo-R2 metric =", varexp_rs
    ))

    if (model_term == "CRS") {
      # Use the individual predicted risk from the CRS model as baseline to
      # compute NRI
      pred_risks_clin <- indiv_pred_risk$CRS
      names(pred_risks_clin) <- seq_len(length(pred_risks_clin))
      nri_cont_val <- NA
    } else {
      print("Calculating NRI")
      # Continuous NRI
      nri_cont <- nribin(
        event = as.numeric(levels(test_data$Y))[test_data$Y],
        p.std = pred_risks_clin,
        p.new = predicted_values,
        updown = "diff",
        cut = 0,
        niter = 0
      )
      nri_cont_val <- nri_cont$nri["NRI", "Estimate"]
    }

    rs_res <- data.frame(
      risk_score = model_term,
      toi = toi,
      train_OR = exp(train_glm_op[toi, "Estimate"]),
      train_OR_L95 = train_or_ci[["2.5 %"]],
      train_OR_U95 = train_or_ci[["97.5 %"]],
      train_pvalue = train_glm_op[toi, "Pr(>|z|)"],
      test_pseudo_R2 = varexp_rs,
      test_OR = exp(test_glm_op[toi, "Estimate"]),
      test_OR_L95 = test_or_ci[["2.5 %"]],
      test_OR_U95 = test_or_ci[["97.5 %"]],
      test_pvalue = test_glm_op[toi, "Pr(>|z|)"],
      test_AUC = auc_ci[2],
      test_AUC_L95 = auc_ci[1],
      test_AUC_U95 = auc_ci[3],
      test_NRI_cont = nri_cont_val
    )
    rs_results <- rbind(rs_results, rs_res)
  }

  write.table(rs_results, paste0(rs_outdir, "/rs_pd_results.txt"),
    quote = FALSE, row.names = FALSE
  )
}

### main ###

#
# Setup
#

# Define each set of variables for the risk scores computation
reqd_vars <- qc(epr_number, Y)

control_vars <- qc(he_age_derived, gender)

# Principal components
pc_vars <- paste0("PC", 1:10)

clinical_vars <- qc(he_bmi_derived, hypertension_derived, he_b008_high_cholesterol)
if (opts$prediabetes) {
  clinical_vars <- c(clinical_vars, "he_c021_pre_diabetes_PARQ")
}

# ExWAS results
# loads anova_list
load("/ddn/gs1/project/controlled/PEGS/StatGen/lloyddt/Diabetes/ExWAS_Final/Output/Farida_Files/anova_list_HE_no_WGS.Rdata")

# Use the significant exposures from ExWAS results for ERS
exwas_results <- NULL
for (evar in names(anova_list)) {
  exwas_results <- rbind(
    exwas_results,
    data.frame(variable = evar, pvalue = anova_list[[evar]][evar, "Pr(>Chisq)"])
  )
}
# In this analysis, we will use BMI as a clinical variable, hence remove it and
# other variables highly correlated with it from the exposure variables
exwas_results <- exwas_results %>%
  filter(!variable %in% qc(
    he_bmi_derived, he_a001b_height_in, he_a002_weight,
    he_a005_physical_health, he_a006_health_comparison
  ))

# pvalue histogram
pdf(paste0(out_dir, "/overall.pdf"))
hist(exwas_results$pvalue,
  breaks = 10, main = paste("Histogram of ExWAS pvalues")
)
dev.off()

# Multiple testing correction of ExWAS results
exwas_results$pvalue_adj <- p.adjust(exwas_results$pvalue, method = "fdr")
exwas_results <- exwas_results %>% arrange(pvalue_adj)
print("ExWAS results:")
print(exwas_results)
exp_vars <- exwas_results %>%
  filter(pvalue_adj < 0.1) %>%
  pull(variable)
print("Significant exposure variables from ExWAS at FDR < 0.1 :")
print(exp_vars)

# Load data
print("Loading data")
# loads epr.he & epr.he.meta
load("/ddn/gs1/project/controlled/PEGS/Data_Freezes/freeze_v1.1/Surveys/Health_and_Exposure/healthexposure_02jan20_v1.1.RData")
# epr.he.conv <- epr_convert_type(epr.he, epr.he.meta)
# saveRDS(epr.he.conv, paste0(out_dir, "/epr.he.conv.rds"))
epr.he.conv <- readRDS(paste0(out_dir, "/epr.he.conv.rds"))
# loads epr.bcbb.map & epr.bcbb.map.meta
load("/ddn/gs1/project/controlled/PEGS/Data_Freezes/freeze_v1.1/Map/bcbb_map_02jan20_v1.1.RData")
# epr.bcbb.map.conv <- epr_convert_type(epr.bcbb.map, epr.bcbb.map.meta)
# saveRDS(epr.bcbb.map.conv, paste0(out_dir, "/epr.bcbb.map.conv.rds"))
epr.bcbb.map.conv <- readRDS(paste0(out_dir, "/epr.bcbb.map.conv.rds"))

# WGSed QCed ppts - no controls, no outliers, unrelated ppts
wgsed_ppts <- fread(paste0(plink_bfileset, ".fam"), select = c(1)) %>%
  setNames(., c("epr_number"))

# ppts whose self-reported race does not match their PCA-inferred race are
# removed from the analyses. These ppts were identified in
# /ddn/gs1/project/controlled/PEGS/StatGen/akhtarifs/projects/epr_explore/PC_strata_debug/strata_pc_outliers.txt
race_outliers <- c(
  100836, 300048, 303703, 306123, 306258, 303694, 306038,
  310624, 314324
)

# Filter H&E data
epr.he.conv <- epr.he.conv %>% filter(he_bmi_derived >= 15)
# Only keep the WGSed QC ppts for this analysis
epr.he.conv <- epr.he.conv %>%
  filter(epr_number %in% wgsed_ppts$epr_number) %>%
  filter(!epr_number %in% race_outliers)
cat("epr.he.conv dim with only wgs ppts: ", dim(epr.he.conv), "\n")

# PRS data
prs_data <- read.table(paste0(PRS_PHENO_DIR, "/output/test_prs_df.txt"),
  header = TRUE
) %>%
  select(epr_number, PRS)


# Create data
reqd_data <- prepare_epr_phenotype(epr_data = epr.he.conv, phenotype = "T2D")

control_data <- inner_join(epr.he.conv, epr.bcbb.map.conv, by = "epr_number") %>%
  select(all_of(c("epr_number", control_vars)))

ht <- prepare_epr_phenotype(epr_data = epr.he.conv, phenotype = "hypertension") %>%
  rename(hypertension_derived = Y)

clin_data <- inner_join(epr.he.conv, ht, by = "epr_number") %>%
  select(all_of(c("epr_number", clinical_vars)))

# Split data into train & test
he_y <- merge(epr.he.conv, reqd_data, by = "epr_number")
set.seed(123)
train_indxs <- createDataPartition(y = he_y$Y, p = .5, list = FALSE)
training_ids <- he_y[train_indxs, ] %>% pull(epr_number)
testing_ids <- he_y[-train_indxs, ] %>% pull(epr_number)
cat("num training_ids = ", length(training_ids), "\n")
cat("num testing_ids = ", length(testing_ids), "\n")

#
# Phase I - select exposure variables that are jointly associated with phenotype
#

P1_data <- inner_join(epr.he.conv, reqd_data, by = "epr_number") %>%
  inner_join(epr.bcbb.map.conv, by = "epr_number") %>%
  select(all_of(c(
    reqd_vars, control_vars, exp_vars,
    "broad_wgs_sample_id_CHILDQ"
  ))) %>%
  filter(epr_number %in% training_ids) %>%
  drop_na()

print("Phase 1 data / training data:")
cat("dimensions:", dim(P1_data), "\n")
print("Y:")
table(P1_data$Y, exclude = NULL)
P1_data %>%
  group_by(Y) %>%
  tally() %>%
  mutate(percent = prop.table(n) * 100) %>%
  print()

rs_ppts_P1 <- P1_data %>% select(epr_number, broad_wgs_sample_id_CHILDQ)
pca_P1_dir <- paste0(out_dir, "/pca_P1")
if (!dir.exists(pca_P1_dir)) dir.create(pca_P1_dir)
write.table(rs_ppts_P1, paste0(pca_P1_dir, "/rs_ppts"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

# Compute PCs for this sample
compute_pcs(
  bfile = plink_bfileset,
  pca_dir = pca_P1_dir,
  indivs_keep = paste0(pca_P1_dir, "/rs_ppts")
)
pcs_P1 <- read.table(paste0(pca_P1_dir, "/plink.eigenvec"),
  header = TRUE, quote = ""
)

P1_data_pc <- inner_join(P1_data, pcs_P1, by = "epr_number") %>%
  select(all_of(c(reqd_vars, control_vars, pc_vars, exp_vars))) %>%
  drop_na()

print("Phase 1 data with PCs:")
cat("dimensions:", dim(P1_data_pc), "\n")
print("Y:")
table(P1_data_pc$Y, exclude = NULL)

# Data QC & filtering
P1_data_filtered <- filter_data(P1_data_pc) %>% drop_na()
# Do variable selection
vars_seld <- do_variable_selection(
  model_df = P1_data_filtered,
  desc = "Exposure variables"
)
sel_exp_vars <- vars_seld$variable
sel_exp_vars <- sel_exp_vars %>%  str_replace(., "\\d$", "") %>% unique()


#
# Phase II - create the final complete cases data for all Risk scores computation
#

# Create data
created_data <- Reduce(
  function(x, y) merge(x, y, by = "epr_number"),
  list(reqd_data, control_data, clin_data, prs_data)
)
P2_data <- inner_join(epr.bcbb.map.conv, epr.he.conv, by = "epr_number") %>%
  select(epr_number, broad_wgs_sample_id_CHILDQ, !!sel_exp_vars) %>%
  inner_join(created_data, by = "epr_number") %>%
  filter(epr_number %in% training_ids) %>%
  select(-PRS) %>%
  drop_na()

# If this is the CRS_without_prediabetes analyses, add the prediabetes variable
# to the final dataset:
# (1) so that this sample subset is the same as that of the
# CRS_with_prediabetes analyses for comparison
# (2) because this variable is required to test its association with risk scores.
if (!opts$prediabetes) {
  P2_data <- merge(P2_data, epr.he.conv[, qc(epr_number, he_c021_pre_diabetes_PARQ)])
}

print("Phase 2 data / training data:")
cat("dimensions:", dim(P2_data), "\n")
print("Y:")
table(P2_data$Y, exclude = NULL)
P2_data %>%
  group_by(Y) %>%
  tally() %>%
  mutate(percent = prop.table(n) * 100) %>%
  print()

# Data QC & filtering
P2_data_filtered <- filter_data(P2_data) %>% drop_na()

print("Phase 2 data after filtering:")
cat("dimensions:", dim(P2_data_filtered), "\n")
print("Y:")
table(P2_data_filtered$Y, exclude = NULL)

rs_ppts_P2 <- P2_data_filtered %>% select(epr_number, broad_wgs_sample_id_CHILDQ)
pca_P2_dir <- paste0(out_dir, "/pca_P2")
if (!dir.exists(pca_P2_dir)) dir.create(pca_P2_dir)
write.table(rs_ppts_P2, paste0(pca_P2_dir, "/rs_ppts"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

# Compute PCs for this sample
compute_pcs(
  bfile = plink_bfileset,
  pca_dir = pca_P2_dir,
  indivs_keep = paste0(pca_P2_dir, "/rs_ppts")
)
pcs_P2 <- read.table(paste0(pca_P2_dir, "/plink.eigenvec"),
  header = TRUE, quote = ""
)
pc_data <- pcs_P2 %>% select(epr_number, !!pc_vars)

#
# Training data
#

training_data <- inner_join(P2_data_filtered, pc_data, by = "epr_number")

# Pre-analysis data visualization & summary
summarize_data(
  df = training_data,
  plotfile = paste0(out_dir, "/training_data_summary.pdf"),
  info_string = "training data"
)

# Risk score computation for training data

# CRS
clin_cc <- training_data %>%
  select(all_of(c(reqd_vars, control_vars, pc_vars, clinical_vars)))
crs <- compute_risk_scores(
  model_df = clin_cc, desc = "Clinical variables", plot = TRUE
)
colnames(crs)[2] <- c("CRS")
crs_scaled <- crs %>% mutate(CRS = as.numeric(scale(CRS)))

# ERS
exp_cc <- training_data %>%
  select(all_of(c(reqd_vars, control_vars, pc_vars, sel_exp_vars)))
ers <- compute_risk_scores(
  model_df = exp_cc, desc = "Exposure variables", plot = TRUE
)
colnames(ers)[2] <- c("ERS")
ers_scaled <- ers %>% mutate(ERS = as.numeric(scale(ERS)))

# PRS
prs <- inner_join(training_data, prs_data, by = "epr_number") %>%
  select(epr_number, PRS)
prs_scaled <- prs %>% mutate(PRS = as.numeric(scale(PRS)))

# ORS
all_rs_train_scaled <- Reduce(
  function(x, y) merge(x, y, by = "epr_number"),
  list(crs_scaled, ers_scaled, prs_scaled)
)
all_rs_train_scaled <- all_rs_train_scaled %>%
  mutate(ORS = CRS + ERS + PRS) %>%
  mutate(ORS = as.numeric(scale(ORS)))

write.table(all_rs_train_scaled, paste0(out_dir, "/risk_scores_scaled_training.txt"),
  quote = FALSE, row.names = FALSE
)

rs_names <- all_rs_train_scaled %>%
  select(ends_with("RS", ignore.case = FALSE)) %>%
  colnames()

# Save the unscaled risk scores
all_rs_train <- Reduce(
  function(x, y) merge(x, y, by = "epr_number"),
  list(crs, ers, prs)
)
write.table(all_rs_train, paste0(out_dir, "/risk_scores_training.txt"),
  quote = FALSE, row.names = FALSE
)

# Recode race and gender values
race_labels <- create_label_vector(
  epr.bcbb.map.meta %>% filter(variable_name == "race") %>% pull(label)
)
gender_labels <- create_label_vector(
  epr.bcbb.map.meta %>% filter(variable_name == "gender") %>% pull(label)
)

training_data <- epr.bcbb.map.conv %>%
  select(epr_number, race) %>%
  inner_join(training_data, by = "epr_number") %>%
  inner_join(all_rs_train_scaled, by = "epr_number") %>%
  mutate(race = recode_factor(race, !!!race_labels)) %>%
  mutate(gender = recode_factor(gender, !!!gender_labels))

saveRDS(training_data, paste0(out_dir, "/training_data.rds"))

#
# Testing data
#

ers_coeffs <- read.table(
  paste0(out_dir, "/coeffs_", make.names("Exposure variables")),
  header = TRUE, quote = ""
)
ers_vars <- ers_coeffs$variable %>%
  str_replace(., "\\d$", "") %>%
  unique()

testing_data <- inner_join(epr.bcbb.map.conv, epr.he.conv, by = "epr_number") %>%
  select(epr_number, broad_wgs_sample_id_CHILDQ, !!ers_vars) %>%
  inner_join(created_data, by = "epr_number") %>%
  filter(epr_number %in% testing_ids) %>%
  select(-PRS)

# If this is the CRS_without_prediabetes analyses, add the prediabetes variable
# to the final dataset:
# (1) so that this sample subset is the same as that of the
# CRS_with_prediabetes analyses for comparison
# (2) because this variable is required to test its association with risk scores.
if (!opts$prediabetes) {
  testing_data <- merge(
    testing_data,
    epr.he.conv[, qc(epr_number, he_c021_pre_diabetes_PARQ)],
    by = "epr_number"
  )
}

if (length(intersect(
  testing_data$epr_number,
  unique(c(P1_data$epr_number, P2_data$epr_number, training_data$epr_number))
))) {
  stop("Overlap detected between training and test sets")
}

testing_data <- testing_data %>% drop_na()

# Compute PCs
testing_ppts <- testing_data %>% select(epr_number, broad_wgs_sample_id_CHILDQ)
testing_pca_dir <- paste0(out_dir, "/pca_test")
write.table(testing_ppts, paste0(testing_pca_dir, "/testing_ppts"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

compute_pcs(
  bfile = plink_bfileset,
  pca_dir = testing_pca_dir,
  indivs_keep = paste0(testing_pca_dir, "/testing_ppts")
)
pcs_testing <- read.table(
  paste0(testing_pca_dir, "/plink.eigenvec"),
  header = TRUE, quote = ""
) %>%
  select(epr_number, !!pc_vars)
testing_data <- inner_join(testing_data, pcs_testing, by = "epr_number")

# Pre-analysis data visualization & summary
summarize_data(
  df = testing_data,
  plotfile = paste0(out_dir, "/testing_data_summary.pdf"),
  info_string = "testing data"
)

# Calculate risk scores for testing data

# CRS
crs_coeffs <- read.table(
  paste0(out_dir, "/coeffs_", make.names("Clinical variables")),
  header = TRUE, quote = ""
)
clin_cc <- testing_data %>%
  select(all_of(c(reqd_vars, control_vars, pc_vars, clinical_vars)))
crs_test <- calc_risk_scores(model_df = clin_cc, coeffs_df = crs_coeffs)
colnames(crs_test)[2] <- c("CRS")
crs_test_scaled <- crs_test %>% mutate(CRS = as.numeric(scale(CRS)))

# ERS
exp_cc <- testing_data %>%
  select(all_of(c(reqd_vars, control_vars, pc_vars, ers_vars)))
ers_test <- calc_risk_scores(model_df = exp_cc, coeffs_df = ers_coeffs)
colnames(ers_test)[2] <- c("ERS")
ers_test_scaled <- ers_test %>% mutate(ERS = as.numeric(scale(ERS)))

# PRS
prs_test <- inner_join(testing_data, prs_data, by = "epr_number") %>%
  select(epr_number, PRS)
prs_test_scaled <- prs_test %>% mutate(PRS = as.numeric(scale(PRS)))

# ORS
all_rs_test_scaled <- Reduce(
  function(x, y) merge(x, y, by = "epr_number"),
  list(crs_test_scaled, ers_test_scaled, prs_test_scaled)
)
all_rs_test_scaled <- all_rs_test_scaled %>%
  mutate(ORS = CRS + ERS + PRS) %>%
  mutate(ORS = as.numeric(scale(ORS)))

write.table(all_rs_test_scaled, paste0(out_dir, "/risk_scores_scaled_testing.txt"),
  quote = FALSE, row.names = FALSE
)

# Save the unscaled risk scores
all_rs_test <- Reduce(
  function(x, y) merge(x, y, by = "epr_number"),
  list(crs_test, ers_test, prs_test)
)
write.table(all_rs_test, paste0(out_dir, "/risk_scores_testing.txt"),
  quote = FALSE, row.names = FALSE
)

testing_data <- epr.bcbb.map.conv %>%
  select(epr_number, race) %>%
  inner_join(testing_data, by = "epr_number") %>%
  inner_join(all_rs_test_scaled, by = "epr_number") %>%
  mutate(race = recode_factor(race, !!!race_labels)) %>%
  mutate(gender = recode_factor(gender, !!!gender_labels))

saveRDS(testing_data, paste0(out_dir, "/testing_data.rds"))

# Plot the unscaled risk scores for training & test data
pdf(paste0(out_dir, "/risk_scores_unscaled.pdf"))
par(mfrow = c(2, 3))
sapply(
  colnames(all_rs_train %>% select(-epr_number)),
  function(x) {
    plot(density(all_rs_train[, x]), main = paste("Training data - ", x))
  }
)

sapply(
  colnames(all_rs_test %>% select(-epr_number)),
  function(x) {
    plot(density(all_rs_test[, x]), main = paste("Testing data -", x))
  }
)
dev.off()

#
# Risk score association and prediction
#

# Strata to do testing for:
strata <- list()
strata[["all"]] <- "all"
strata[["race"]] <- c("Black or African American", "White")
strata[["gender"]] <- c("Female", "Male")
print("strata to test risk scores:")
print(strata)

# Test for each level of each stratum type
for (i in seq_along(strata)) {
  strat_by <- names(strata)[i]
  strat_levels <- unname(strata[[i]])
  cat("Stratifying by:", strat_by, "\n")

  for (s_lvl in strat_levels) {
    cat("Stratum level: ", s_lvl, "\n")

    strata_dir <- paste0(out_dir, "/", make.names(s_lvl))
    cat("strata_dir = ", strata_dir, "\n")
    if (!dir.exists(strata_dir)) dir.create(strata_dir)

    strata_train_data <- training_data
    strata_test_data <- testing_data

    # For each stratum level except when its the entire dataset, subset the data,
    # recompute PCs and create the dataframe with new PCs.
    if (s_lvl != "all") {
      # Do this for both the training data and the testing data
      for (ts in c("train", "test")) {
        strata_data <- NULL
        cat(ts, "dataset", "\n")
        if (ts == "train") {
          strata_data <- strata_train_data
        } else {
          strata_data <- strata_test_data
        }

        strata_data <- strata_data %>%
          # subset the data for this stratum level
          filter(!!as.symbol(strat_by) == s_lvl) %>%
          # remove the PCs since they will be recomputed for this strata level
          select(-all_of(pc_vars))

        print("before PCA") # TODO remove
        cat("data dimensions:", dim(strata_data), "\n")
        print("Y:")
        print(table(strata_data$Y, exclude = NULL))
        print(colnames(strata_data))

        strata_pc_dir <- paste0(strata_dir, "/pca_", ts)
        if (!dir.exists(strata_pc_dir)) dir.create(strata_pc_dir)

        write.table(
          strata_data %>% select(epr_number, broad_wgs_sample_id_CHILDQ),
          paste0(strata_pc_dir, "/strata_ppts"),
          row.names = FALSE, col.names = FALSE, quote = FALSE
        )

        # Compute PCs for this strata subset
        compute_pcs(
          bfile = plink_bfileset,
          pca_dir = strata_pc_dir,
          indivs_keep = paste0(strata_pc_dir, "/strata_ppts")
        )
        strata_pcs <- read.table(
          paste0(strata_pc_dir, "/plink.eigenvec"),
          header = TRUE, quote = ""
        )
        if (ts == "train") {
          strata_train_data <- inner_join(strata_data, strata_pcs, by = "epr_number")
        } else {
          strata_test_data <- inner_join(strata_data, strata_pcs, by = "epr_number")
        }
      }
    }

    print("after PCA") # TODO remove
    cat("strata training data dimensions:", dim(strata_train_data), "\n")
    print("Y:")
    print(table(strata_train_data$Y, exclude = NULL))

    cat("strata testing data dimensions:", dim(strata_test_data), "\n")
    print("Y:")
    print(table(strata_test_data$Y, exclude = NULL))

    # Risk score association, prediction and plotting
    test_rs_assoc(
      train_data = strata_train_data,
      test_data = strata_test_data,
      rs_outdir = strata_dir
    )
    plot_rs_results(
      train_data = strata_train_data,
      test_data = strata_test_data,
      rs_outdir = strata_dir
    )

    # Risk score + E association
    test_rs_e(
      test_data = strata_test_data,
      rs_term = "CRS",
      rs_outdir = strata_dir
    )
    test_rs_e(
      test_data = strata_test_data,
      rs_term = "PRS",
      rs_outdir = strata_dir
    )

    # Risk score association and prediction with prediabetes
    test_rs_pd(
      train_data = strata_train_data,
      test_data = strata_test_data,
      rs_outdir = strata_dir
    )
  }
}


### cleanup ###

dev.off()
