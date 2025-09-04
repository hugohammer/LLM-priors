# --- 0. SCRIPT SETUP ---
cat("--- Script Setup ---\n")
start_time <- Sys.time()

# Load necessary libraries
# install.packages(c("INLA", "RColorBrewer", "correctR")) # Run once if not installed
if (!requireNamespace("INLA", quietly = TRUE)) stop("Package 'INLA' is not installed.")
library(INLA)
if (!requireNamespace("RColorBrewer", quietly = TRUE)) warning("Package 'RColorBrewer' not installed.")
if (!requireNamespace("correctR", quietly = TRUE)) {
  stop("Package 'correctR' is not installed. Please run: install.packages('correctR')")
}
library(correctR)

# --- Helper Functions ---
format_metric_with_pval <- function(point_est, p_val) {
  if(is.null(point_est) || is.na(point_est)) return("N/A")
  p_val_str <- if(is.null(p_val) || is.na(p_val)) "" else if(p_val < 0.001) "(p<0.001)" else sprintf("(p=%.3f)", p_val)
  sprintf("%.3f %s", point_est, p_val_str)
}

create_k_folds <- function(data_size, k = 5) {
  if (k <= 1 || k > data_size) stop("Number of folds k must be > 1 and <= data_size.")
  shuffled_indices <- sample(1:data_size); fold_indices <- cut(1:data_size, breaks=k, labels=F)
  folds <- list(); for (i in 1:k) folds[[i]] <- shuffled_indices[which(fold_indices == i)]; return(folds)
}
generate_colors <- function(n) {if(n==0)return(character(0));if(requireNamespace("RColorBrewer",quietly=T)&&n>0)RColorBrewer::brewer.pal(max(3,n),"Set2")[1:n] else if (n>0&&n<=8)(1:n) else if (n>0) rainbow(n) else character(0)}


# --- 1. DATA LOADING AND FORMULA DEFINITION ---
cat("\n--- 1. Data Loading and Formula Definition ---\n")
data_file_path <- "../Data/Concrete_Data.csv"
if (!file.exists(data_file_path)) stop(paste("Data file not found at:", data_file_path))
concrete_data <- read.table(data_file_path, header = TRUE, sep = ",")
if (!"CCS" %in% colnames(concrete_data)) stop("Response variable 'CCS' not found.")
model_formula <- CCS ~ Cement + Blast_Furnace_Slag + Fly_Ash + Water + Superplasticizer + Coarse_Aggregate + Fine_Aggregate + Age
cat("Concrete data loaded. N =", nrow(concrete_data), "\n")
coeff_names_to_plot <- c("Cement", "Blast_Furnace_Slag", "Fly_Ash", "Water", 
                         "Superplasticizer", "Coarse_Aggregate", "Fine_Aggregate", "Age")

# --- 2. DEFINE ALL PRIOR SETS FOR CONCRETE EXAMPLE ---
cat("\n--- 2. Define All Prior Sets for Concrete Example ---\n")
default_prior_intercept_mean <- 0
default_prior_intercept_prec <- 1 / (5^2)
prior_sets_concrete <- list()

prior_sets_concrete$Gemini_Mod <- list(
  means = list( `(Intercept)` = default_prior_intercept_mean, Cement = 0.10, Blast_Furnace_Slag = 0.06, Fly_Ash = 0.05, Water = -0.15, Superplasticizer = 0.25, Coarse_Aggregate = 0.02, Fine_Aggregate = 0.02, Age = 0.30),
  precs = list( `(Intercept)` = default_prior_intercept_prec, Cement = 400, Blast_Furnace_Slag = 625, Fly_Ash = 625, Water = 1/(0.07^2), Superplasticizer = 1/(0.15^2), Coarse_Aggregate = 2500, Fine_Aggregate = 2500, Age = 1/(0.15^2))
)
prior_sets_concrete$Gemini_Weak <- list(
  means = list( `(Intercept)` = default_prior_intercept_mean, Cement = 0, Blast_Furnace_Slag = 0, Fly_Ash = 0, Water = 0, Superplasticizer = 0, Coarse_Aggregate = 0, Fine_Aggregate = 0, Age = 0),
  precs = list( `(Intercept)` = default_prior_intercept_prec, Cement = 25, Blast_Furnace_Slag = 25, Fly_Ash = 25, Water = 25, Superplasticizer = 4, Coarse_Aggregate = 100, Fine_Aggregate = 100, Age = 1)
)
prior_sets_concrete$ChatGPT_Mod <- list(
  means = list( `(Intercept)` = default_prior_intercept_mean, Cement = 0.03, Blast_Furnace_Slag = 0.02, Fly_Ash = 0.01, Water = -0.10, Superplasticizer = 0.50, Coarse_Aggregate = 0.002, Fine_Aggregate = -0.001, Age = 0.03),
  precs = list( `(Intercept)` = default_prior_intercept_prec, Cement = 1/(0.02^2), Blast_Furnace_Slag = 1/(0.02^2), Fly_Ash = 1/(0.015^2), Water = 1/(0.05^2), Superplasticizer = 1/(0.30^2), Coarse_Aggregate = 1/(0.002^2), Fine_Aggregate = 1/(0.002^2), Age = 1/(0.01^2))
)
prior_sets_concrete$ChatGPT_Weak <- list(
  means = list( `(Intercept)` = default_prior_intercept_mean, Cement = 0, Blast_Furnace_Slag = 0, Fly_Ash = 0, Water = 0, Superplasticizer = 0, Coarse_Aggregate = 0, Fine_Aggregate = 0, Age = 0),
  precs = list( `(Intercept)` = default_prior_intercept_prec, Cement = 1/(1^2), Blast_Furnace_Slag = 1, Fly_Ash = 1, Water = 1, Superplasticizer = 1, Coarse_Aggregate = 1, Fine_Aggregate = 1, Age = 1)
)
prior_sets_concrete$Claude_Mod <- list(
  means = list(
    `(Intercept)` = default_prior_intercept_mean, Cement = 0.10, Blast_Furnace_Slag = 0.08, Fly_Ash = 0.05,
    Water = -0.30, Superplasticizer = 2.0, Coarse_Aggregate = 0.02,
    Fine_Aggregate = 0.02, Age = 0.30
  ),
  precs = list(
    `(Intercept)` = default_prior_intercept_prec, Cement = 400, Blast_Furnace_Slag = 625, Fly_Ash = 625,
    Water = 100, Superplasticizer = 1, Coarse_Aggregate = 2500,
    Fine_Aggregate = 2500, Age = 100
  )
)
prior_sets_concrete$Claude_Weak <- list(
  means = list(
    `(Intercept)` = default_prior_intercept_mean, Cement = 0.05, Blast_Furnace_Slag = 0.03, Fly_Ash = 0.02,
    Water = -0.10, Superplasticizer = 1.0, Coarse_Aggregate = 0.01,
    Fine_Aggregate = 0.01, Age = 0.20
  ),
  precs = list(
    `(Intercept)` = default_prior_intercept_prec, Cement = 25, Blast_Furnace_Slag = 25, Fly_Ash = 25,
    Water = 11.11, Superplasticizer = 0.111, Coarse_Aggregate = 100,
    Fine_Aggregate = 100, Age = 4
  )
)
cat(length(prior_sets_concrete), "prior sets defined.\n")
get_sds_from_precs_concrete <- function(precs_list) sapply(precs_list, function(p) 1/sqrt(p))

# --- 3. FIT FULL MODELS ---
cat("\n--- 3. Fitting Full Models (for parameter plots) ---\n")
full_lm_model_concrete <- lm(model_formula, data = concrete_data)
coeffs_lm_concrete <- summary(full_lm_model_concrete)$coefficients
cat("Frequentist lm model (full data) fitted.\n")

inla_models_concrete_list_fullfit <- list()
cat("Fitting INLA models (full data, all priors) for posterior distributions...\n")
for (prior_name in names(prior_sets_concrete)) {
  current_priors_inla_format <- list(mean=prior_sets_concrete[[prior_name]]$means, prec=prior_sets_concrete[[prior_name]]$precs)
  tryCatch({
    inla_models_concrete_list_fullfit[[prior_name]] <- INLA::inla(model_formula, family="gaussian", data=concrete_data, control.fixed=current_priors_inla_format, control.compute=list(config=TRUE))
  }, error=function(e) message(paste("Error INLA full fit", prior_name, ":", e$message)))
  if(is.null(inla_models_concrete_list_fullfit[[prior_name]])) warning(paste("Full INLA fit failed for", prior_name))
}
cat("Full INLA models fitted for posterior plots.\n")

# --- 4. 5-FOLD CV ANALYSIS (FULL DATASET) ---
cat("\n--- 4. 5-Fold CV Analysis (Full Dataset) ---\n")
k_folds_cv <- 5
folds_full_data <- create_k_folds(nrow(concrete_data), k = k_folds_cv)
performance_results_full <- list()
y_obs_full <- concrete_data$CCS

# 4.1 Frequentist LM: 5-fold CV
cat("Performing 5-fold CV for LM (full dataset)...\n")
preds_lm_kfold_full <- numeric(nrow(concrete_data)); sigmas_lm_kfold_full <- numeric(nrow(concrete_data))
for (k_f in 1:k_folds_cv) {
  test_indices <- folds_full_data[[k_f]]; train_indices <- setdiff(1:nrow(concrete_data), test_indices)
  lm_kfold <- lm(model_formula, data = concrete_data[train_indices,])
  preds_lm_kfold_full[test_indices] <- predict(lm_kfold, newdata = concrete_data[test_indices,])
  sigmas_lm_kfold_full[test_indices] <- summary(lm_kfold)$sigma
}
log_s_lm_ind <- -dnorm(y_obs_full, preds_lm_kfold_full, sigmas_lm_kfold_full, log=T); log_s_lm <- mean(log_s_lm_ind, na.rm=T)
rmse_lm <- sqrt(mean((y_obs_full - preds_lm_kfold_full)^2, na.rm=T)); mae_lm <- mean(abs(y_obs_full - preds_lm_kfold_full), na.rm=T)
performance_results_full[["LM_Full_5CV"]] <- list(ModelType="LM", PriorSet="N/A", LogS=log_s_lm, RMSE=rmse_lm, MAE=mae_lm, p_logS=NA, p_rmse=NA, p_mae=NA)
cat("LM 5-fold CV complete.\n")

# 4.2 INLA Models: 5-fold CV
cat("Performing 5-fold CV for INLA priors...\n")
for (prior_name in names(prior_sets_concrete)) {
  cat(paste("  Processing INLA 5-fold CV for prior set:", prior_name, "\n"))
  preds_inla_kfold <- numeric(nrow(concrete_data)); sigmas_inla_kfold <- numeric(nrow(concrete_data))
  priors_current_inla_kfold <- list(mean=prior_sets_concrete[[prior_name]]$means, prec=prior_sets_concrete[[prior_name]]$precs)
  pb_inla_kcv <- txtProgressBar(min=0, max=k_folds_cv, style=3)
  for (k_f in 1:k_folds_cv) {
    setTxtProgressBar(pb_inla_kcv, k_f)
    test_idx <- folds_full_data[[k_f]]; train_idx <- setdiff(1:nrow(concrete_data), test_idx)
    data_kcv_inla_train <- concrete_data[train_idx,]; data_kcv_inla_test_na <- concrete_data[test_idx,]; data_kcv_inla_test_na$CCS <- NA
    data_kcv_inla_comb <- rbind(data_kcv_inla_train, data_kcv_inla_test_na)
    inla_kcv_model <- NULL
    tryCatch({inla_kcv_model <- INLA::inla(model_formula, family="gaussian", data=data_kcv_inla_comb, control.fixed=priors_current_inla_kfold, control.predictor=list(compute=T), control.compute=list(config=T),verbose=F)
    }, error = function(e) message(paste("Error INLA kCV fold",k_f, prior_name, ":", e$message)))
    if(!is.null(inla_kcv_model) && !is.null(inla_kcv_model$summary.fitted.values)){
      pred_idx_kcv <- (nrow(data_kcv_inla_train)+1):nrow(data_kcv_inla_comb)
      preds_inla_kfold[test_idx] <- inla_kcv_model$summary.fitted.values[pred_idx_kcv, "mean"]
      if(!is.null(inla_kcv_model$summary.hyperpar)){
        prec_eps_sum <- inla_kcv_model$summary.hyperpar[grep("Precision for the Gaussian observations",rownames(inla_kcv_model$summary.hyperpar)),]
        if(nrow(prec_eps_sum)>0) sigmas_inla_kfold[test_idx]<-1/sqrt(prec_eps_sum[1,"mean"]) else sigmas_inla_kfold[test_idx]<-NA
      } else {sigmas_inla_kfold[test_idx]<-NA}
    } else { preds_inla_kfold[test_idx]<-NA; sigmas_inla_kfold[test_idx]<-NA }
  }
  close(pb_inla_kcv)
  
  log_s_inla_ind <- -dnorm(y_obs_full, preds_inla_kfold, sigmas_inla_kfold, log=T); log_s_inla <- mean(log_s_inla_ind, na.rm=T)
  rmse_inla <- sqrt(mean((y_obs_full - preds_inla_kfold)^2, na.rm=T)); mae_inla <- mean(abs(y_obs_full - preds_inla_kfold), na.rm=T)
  
  log_scores_inla <- -dnorm(y_obs_full, preds_inla_kfold, sigmas_inla_kfold, log=T)
  log_scores_lm <- -dnorm(y_obs_full, preds_lm_kfold_full, sigmas_lm_kfold_full, log=T)
  se_scores_inla <- (y_obs_full - preds_inla_kfold)^2
  se_scores_lm <- (y_obs_full - preds_lm_kfold_full)^2
  ae_scores_inla <- abs(y_obs_full - preds_inla_kfold)
  ae_scores_lm <- abs(y_obs_full - preds_lm_kfold_full)

  valid_log_indices <- !is.na(log_scores_inla) & !is.na(log_scores_lm)
  p_logS <- NA
  if(sum(valid_log_indices) > 1) {
      p_logS <- kfold_ttest(x = log_scores_inla[valid_log_indices], y = log_scores_lm[valid_log_indices], n = nrow(concrete_data), k = k_folds_cv, tailed = "one", greater = "y")$p.value
  }

  valid_se_indices <- !is.na(se_scores_inla) & !is.na(se_scores_lm)
  p_rmse <- NA
  if(sum(valid_se_indices) > 1) {
      p_rmse <- kfold_ttest(x = se_scores_inla[valid_se_indices], y = se_scores_lm[valid_se_indices], n = nrow(concrete_data), k = k_folds_cv, tailed = "one", greater = "y")$p.value
  }

  valid_ae_indices <- !is.na(ae_scores_inla) & !is.na(ae_scores_lm)
  p_mae <- NA
  if(sum(valid_ae_indices) > 1) {
      p_mae <- kfold_ttest(x = ae_scores_inla[valid_ae_indices], y = ae_scores_lm[valid_ae_indices], n = nrow(concrete_data), k = k_folds_cv, tailed = "one", greater = "y")$p.value
  }
  
  performance_results_full[[paste0(prior_name,"_5CV")]] <- list(ModelType="INLA", PriorSet=prior_name, LogS=log_s_inla, RMSE=rmse_inla, MAE=mae_inla, p_logS=p_logS, p_rmse=p_rmse, p_mae=p_mae)
}
cat("INLA 5-fold CV for all priors and t-test p-value calculations complete.\n")


# --- 5. REVISED PLOTTING SECTION (FULL DATASET, CONCRETE EXAMPLE) ---
cat("\n--- 5. Generating Requested Plots (Concrete Example) ---\n")
prior_categories_concrete <- list(
    Moderate = grep("_Mod$", names(prior_sets_concrete), value = TRUE),
    Weak = grep("_Weak$", names(prior_sets_concrete), value = TRUE)
)

# --- Define Manual X-Limits for PLOTS 1 & 2 (Combined Priors) ---
xlims_moderate <- list(
    Cement             = c(-0.05, 0.25),   Blast_Furnace_Slag = c(-0.1, 0.25),
    Fly_Ash            = c(-0.1, 0.20),     Water              = c(-0.7, 0.2),
    Superplasticizer   = c(-0.5, 4.0),     Coarse_Aggregate   = c(-0.025, 0.08),
    Fine_Aggregate     = c(-0.025, 0.075), Age                = c(-0.1, 0.5)
)
xlims_weak <- list(
    Cement             = c(-1.5, 1.5),     Blast_Furnace_Slag = c(-1.5, 1.5),
    Fly_Ash            = c(-1.5, 1.5),     Water              = c(-1.5, 1.5),
    Superplasticizer   = c(-3, 3),         Coarse_Aggregate   = c(-0.5, 0.5),
    Fine_Aggregate     = c(-0.5, 0.5),     Age                = c(-2, 2)
)
ylims_weak <- list(
    Cement             = c(0, 5),   Blast_Furnace_Slag = c(0, 5),
    Fly_Ash            = c(0, 5),   Water              = c(0, 5),
    Superplasticizer   = c(0, 2.5), Coarse_Aggregate   = c(0, 10),
    Fine_Aggregate     = c(0, 10),  Age                = c(0, 2)
)

# --- Define Manual X-Limits for PLOT 3 (Detailed Prior/Posterior/Likelihood) ---
manual_xlims_detailed <- list()
manual_xlims_detailed$Gemini_Mod <- list(
    Cement             = c(0.05, 0.20),    Blast_Furnace_Slag = c(0, 0.15),
    Fly_Ash            = c(0, 0.15),      Water              = c(-0.3, 0),
    Superplasticizer   = c(-0.2, 0.8),    Coarse_Aggregate   = c(-0.02, 0.06),
    Fine_Aggregate     = c(-0.02, 0.06),  Age                = c(0.08, 0.15)
)

# --- Plot 1 & 2: MLE vs. Moderate/Weak Priors ---
color_palette_concrete <- c(
  "Claude_Mod" = "#E41A1C",  # Red for Claude
  "Claude_Weak" = "#E41A1C",
  "Gemini_Mod" = "#4DAF4A",  # Green for Gemini
  "Gemini_Weak" = "#4DAF4A",
  "ChatGPT_Mod" = "#377EB8", # Blue for ChatGPT
  "ChatGPT_Weak" = "#377EB8"
)

linetype_palette_concrete <- c(
  "Claude_Mod" = 2,   # Dashed line for Claude
  "Claude_Weak" = 2,
  "Gemini_Mod" = 3,   # Dotted line for Gemini
  "Gemini_Weak" = 3,
  "ChatGPT_Mod" = 4,  # Dot-dashed line for ChatGPT
  "ChatGPT_Weak" = 4
)


# --- Plot 1 & 2: MLE vs. Moderate/Weak Priors ---
cat("Generating combined prior comparison plots...\n")
for (category_name in names(prior_categories_concrete)) {
  prior_names_in_cat <- prior_categories_concrete[[category_name]]
  if(length(prior_names_in_cat) == 0) next

  active_xlims <- if (category_name == "Moderate") xlims_moderate else xlims_weak
  active_ylims <- if (category_name == "Weak") ylims_weak else NULL

  combined_plot_fn_concrete <- paste0("concrete_dist_lm_vs_", tolower(category_name), "_priors.png")
  png(combined_plot_fn_concrete, width = 11, height = 16, units = "in", res = 300)
  op_comb_concrete <- par(mfrow=c(4, 2), mar=c(5, 5, 2, 2), oma=c(0, 0, 1, 0))
  
  for (coeff_name in coeff_names_to_plot) {
    lm_est <- coeffs_lm_concrete[coeff_name,"Estimate"]; lm_se <- coeffs_lm_concrete[coeff_name,"Std. Error"]
    plot_xlim_c <- active_xlims[[coeff_name]]
    if (is.null(plot_xlim_c)) { plot_xlim_c <- c(-1, 1) }

    if (!is.null(active_ylims)) {
        plot_ylim_c <- active_ylims[[coeff_name]]
        if(is.null(plot_ylim_c)) plot_ylim_c <- c(0,1)
    } else {
        y_maxes_calc <- c()
        x_seq_for_ylim <- seq(plot_xlim_c[1], plot_xlim_c[2], length.out = 200)
        if(!is.na(lm_est) && !is.na(lm_se) && lm_se>0) y_maxes_calc <- c(y_maxes_calc, max(dnorm(x_seq_for_ylim, lm_est, lm_se),na.rm=T))
        for(p_name in prior_names_in_cat) {
            pr_m <- prior_sets_concrete[[p_name]]$means[[coeff_name]]; pr_s <- 1/sqrt(prior_sets_concrete[[p_name]]$precs[[coeff_name]])
            if(!is.na(pr_m) && !is.na(pr_s) && pr_s > 0) y_maxes_calc <- c(y_maxes_calc, max(dnorm(x_seq_for_ylim, pr_m, pr_s),na.rm=T))
        }
        plot_ylim_c <- if(length(y_maxes_calc)>0 && any(is.finite(y_maxes_calc))) c(0,max(y_maxes_calc[is.finite(y_maxes_calc)],na.rm=T)*1.2) else c(0,1)
    }

    plot(0,type="n",xlim=plot_xlim_c,ylim=plot_ylim_c,xlab=coeff_name,ylab="Density", cex.lab=1.6, cex.axis=1.4)
    grid(col="lightgray",lty="dotted"); legend_items_c <- list()
    
    if(!is.na(lm_est) && !is.na(lm_se) && lm_se>0){
        curve(dnorm(x,lm_est,lm_se),add=T,col="black",lwd=2,lty=1)
        legend_items_c[["MLE"]] <- list(col="black",lty=1,lwd=2)
    }
    for(p_name in prior_names_in_cat){
      pr_m <- prior_sets_concrete[[p_name]]$means[[coeff_name]]
      pr_s <- 1/sqrt(prior_sets_concrete[[p_name]]$precs[[coeff_name]])
      current_color <- color_palette_concrete[p_name]
      current_linetype <- linetype_palette_concrete[p_name]

      if(!is.na(pr_m) && !is.na(pr_s) && pr_s > 0){
          curve(dnorm(x,pr_m,pr_s),add=T,col=current_color,lwd=2,lty=current_linetype)
          legend_items_c[[p_name]]<-list(col=current_color,lty=current_linetype,lwd=2)
      }
    }
    if(length(legend_items_c)>0) legend("topright",names(legend_items_c),col=sapply(legend_items_c,`[[`,"col"),lty=sapply(legend_items_c,`[[`,"lty"),lwd=sapply(legend_items_c,`[[`,"lwd"),cex=1.2,bty="n")
  }
  dev.off(); par(op_comb_concrete)
  cat(paste("Combined prior plot for",category_name,"category saved to:",combined_plot_fn_concrete,"\n"))
}

# --- Plot 3: Detailed Prior vs. Likelihood vs. Posterior for EACH prior set ---
cat("Generating detailed comparison plots for each prior set...\n")
for (prior_name in names(prior_sets_concrete)) {
  detailed_plot_filename_concrete <- paste0("concrete_dist_detailed_comparison_", prior_name, ".png")
  png(detailed_plot_filename_concrete, width = 11, height = 16, units = "in", res = 300)
  op_detail <- par(mfrow = c(4, 2), mar = c(5.1, 5.1, 2.1, 2.1), oma = c(0, 0, 1, 0))
  current_inla_model <- inla_models_concrete_list_fullfit[[prior_name]]
  current_priors_means <- prior_sets_concrete[[prior_name]]$means
  current_priors_sds <- get_sds_from_precs_concrete(prior_sets_concrete[[prior_name]]$precs)
  if(is.null(current_inla_model)){
    plot(0,type="n",axes=F,xlab="",ylab="",main=paste("INLA model for",prior_name,"missing."));dev.off();par(op_detail);next
  }

  for (coeff_name in coeff_names_to_plot) {
    prior_m <- current_priors_means[[coeff_name]]; prior_s <- current_priors_sds[[coeff_name]]
    inla_marg <- current_inla_model$marginals.fixed[[coeff_name]]
    lm_est <- coeffs_lm_concrete[coeff_name,"Estimate"]; lm_se <- coeffs_lm_concrete[coeff_name,"Std. Error"]
    
    plot_xlim <- NULL # Initialize

    if (prior_name %in% names(manual_xlims_detailed) && coeff_name %in% names(manual_xlims_detailed[[prior_name]])) {
        plot_xlim <- manual_xlims_detailed[[prior_name]][[coeff_name]]
    }

    # If no manual limit was found, fall back to the dynamic calculation
    if (is.null(plot_xlim)) {
        x_ranges <- list()
        if(!is.null(inla_marg)) x_ranges[[length(x_ranges)+1]] <- range(inla_marg[,1])
        if(!is.na(lm_est) && !is.na(lm_se) && lm_se > 0) x_ranges[[length(x_ranges)+1]] <- c(lm_est - 4*lm_se, lm_est + 4*lm_se)
        # Dynamic calculation for weak priors excludes the prior's range
        if(!prior_name %in% prior_categories_concrete$Weak) {
          if(!is.na(prior_m) && !is.na(prior_s) && prior_s > 0) x_ranges[[length(x_ranges)+1]] <- c(prior_m-4*prior_s, prior_m+4*prior_s)
        }
        plot_xlim <- if(length(x_ranges)>0) range(unlist(x_ranges),na.rm=T) else c(-1,1)
    }
    
    # Y-limit calculation based on visible range
    y_maxes <- c()
    x_seq_for_ylim <- seq(plot_xlim[1], plot_xlim[2], length.out = 200)
    if(!is.na(prior_m) && !is.na(prior_s) && prior_s > 0) y_maxes <- c(y_maxes, max(dnorm(x_seq_for_ylim, prior_m, prior_s),na.rm=T))
    if(!is.null(inla_marg)){inla_func<-approxfun(inla_marg[,1],inla_marg[,2],rule=2); y_maxes<-c(y_maxes,max(inla_func(x_seq_for_ylim),na.rm=T))}
    if(!is.na(lm_est) && !is.na(lm_se) && lm_se > 0) y_maxes <- c(y_maxes, max(dnorm(x_seq_for_ylim, lm_est, lm_se),na.rm=T))
    plot_ylim <- if(length(y_maxes)>0 && any(is.finite(y_maxes))) c(0,max(y_maxes[is.finite(y_maxes)],na.rm=T)*1.2) else c(0,1)

    # Create Plot
    plot(0,type="n",xlim=plot_xlim,ylim=plot_ylim,xlab=coeff_name,ylab="Density",cex.lab=1.6,cex.axis=1.4)
    grid(col="lightgray",lty="dotted"); legend_items <- list()
    if(!is.na(prior_m) && !is.na(prior_s) && prior_s > 0){curve(dnorm(x,prior_m,prior_s),add=T,col="darkolivegreen4",lwd=2,lty=3); legend_items$Prior<-list(col="darkolivegreen4",lty=3,lwd=2)}
    if(!is.null(inla_marg)){lines(inla_marg,col="deepskyblue3",lwd=2); legend_items$Posterior<-list(col="deepskyblue3",lty=1,lwd=2)}
    if(!is.na(lm_est) && !is.na(lm_se) && lm_se>0){curve(dnorm(x,lm_est,lm_se),add=T,col="black",lwd=2,lty=1); legend_items[["MLE"]] <- list(col="black",lty=1,lwd=2)}
    if(length(legend_items)>0) legend("topright",names(legend_items),col=sapply(legend_items,`[[`,"col"),lty=sapply(legend_items,`[[`,"lty"),lwd=sapply(legend_items,`[[`,"lwd"),cex=1.2,bty="n")
  }
  dev.off(); par(op_detail)
}
cat("Detailed comparison plots for each prior set saved.\n")


# --- 7. COMPILE AND PRINT RESULTS TABLE ---
cat("\n--- 7. Summary Table of Predictive Performance ---\n")
results_df_full_concrete <- do.call(rbind, lapply(names(performance_results_full), function(name) {
  x <- performance_results_full[[name]]
  data.frame(
    ModelName = name, ModelType = x$ModelType, PriorSet = x$PriorSet,
    MeanNegLogScore = format_metric_with_pval(x$LogS, x$p_logS),
    RMSE = format_metric_with_pval(x$RMSE, x$p_rmse),
    MAE = format_metric_with_pval(x$MAE, x$p_mae),
    stringsAsFactors = FALSE
  )
}))
cat("\n--- 5-Fold CV Performance (Full Dataset) ---\n")
print(results_df_full_concrete)
results_csv_filename <- "kfold_cv_concrete_performance_summary.csv"
tryCatch({
  write.csv(results_df_full_concrete, results_csv_filename, row.names = FALSE)
  cat(paste("\nConcrete performance summary saved to:", results_csv_filename, "\n"))
}, error = function(e){
  cat(paste("\nError saving concrete summary to CSV:", e$message, "\n"))
})


# --- 8. GENERATE COMBINED KL DIVERGENCE SUMMARY TABLE (Concrete Example) ---
cat("\n--- 8. Generating Combined KL Divergence Summary Table ---\n")
kl_divergence_norm <- function(mu1, sd1, mu2, sd2) {
  if(is.na(sd1) || is.na(sd2) || sd1 <= 0 || sd2 <= 0) return(NA)
  log(sd2 / sd1) + (sd1^2 + (mu1 - mu2)^2) / (2 * sd2^2) - 0.5
}
coeffs_lm_concrete <- summary(full_lm_model_concrete)$coefficients
lm_means <- coeffs_lm_concrete[, "Estimate"]
lm_sds   <- coeffs_lm_concrete[, "Std. Error"]
kl_divergence_table_concrete <- data.frame(matrix(
    NA,
    nrow = length(names(prior_sets_concrete)),
    ncol = length(coeff_names_to_plot),
    dimnames = list(names(prior_sets_concrete), coeff_names_to_plot)
))

cat("Calculating KL Divergence for each prior and predictor variable...\n")
for (prior_name in names(prior_sets_concrete)) {
  for (coeff in coeff_names_to_plot) {
    prior_mean <- prior_sets_concrete[[prior_name]]$means[[coeff]]
    prior_sd   <- 1/sqrt(prior_sets_concrete[[prior_name]]$precs[[coeff]])
    kl_divergence_table_concrete[prior_name, coeff] <- kl_divergence_norm(
      mu1 = lm_means[coeff], sd1 = lm_sds[coeff],
      mu2 = prior_mean,      sd2 = prior_sd
    )
  }
}

cat("Calculating summary metrics (Average KL and Average Rank)...\n")
kl_ranks_table_concrete <- apply(kl_divergence_table_concrete, 2, rank, ties.method = "min")
kl_divergence_table_concrete$`Average KL Divergence` <- rowMeans(kl_divergence_table_concrete, na.rm = TRUE)
kl_divergence_table_concrete$`Average Rank` <- rowMeans(kl_ranks_table_concrete, na.rm = TRUE)
kl_divergence_table_concrete <- kl_divergence_table_concrete[order(kl_divergence_table_concrete$`Average Rank`), ]
kl_divergence_table_concrete <- cbind(`Prior Set` = rownames(kl_divergence_table_concrete), kl_divergence_table_concrete)
rownames(kl_divergence_table_concrete) <- NULL

cat("\n--- Combined KL Divergence Table (Concrete): KL(Likelihood || Prior) ---\n")
formatted_table_concrete <- kl_divergence_table_concrete
cols_to_format_concrete <- names(formatted_table_concrete)[-1]
formatted_table_concrete[cols_to_format_concrete] <- lapply(formatted_table_concrete[cols_to_format_concrete], function(x) {
    format(round(x, 4), nsmall = 4)
})
print(formatted_table_concrete)

kl_table_concrete_filename <- "kl_divergence_summary_combined_concrete.csv"
tryCatch({
  write.csv(kl_divergence_table_concrete, kl_table_concrete_filename, row.names = FALSE)
  cat(paste("\nCombined KL Divergence summary table for concrete saved to:", kl_table_concrete_filename, "\n"))
}, error = function(e){
  cat(paste("\nError saving combined KL Divergence summary to CSV:", e$message, "\n"))
})

end_time <- Sys.time()
cat("\n--- End of Concrete Example Script ---")
cat("\nTotal execution time:", format(end_time - start_time), "\n")
