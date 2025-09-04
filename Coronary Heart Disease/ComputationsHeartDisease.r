# --- 0. SCRIPT SETUP ---
cat("--- Script Setup ---\n")
start_time <- Sys.time()

# Load necessary libraries
# install.packages(c("INLA", "pROC", "RColorBrewer", "correctR")) # Run once if not installed
if (!requireNamespace("INLA", quietly = TRUE)) {
  stop("Package 'INLA' is not installed. Please install it from https://www.r-inla.org/download-install")
}
library(INLA)
if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Package 'pROC' is not installed. AUC comparison will be skipped.")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  warning("Package 'RColorBrewer' not installed. Using default colors.")
}
if (!requireNamespace("correctR", quietly = TRUE)) {
  stop("Package 'correctR' is not installed. Please run: install.packages('correctR')")
}
library(correctR)

N_BOOT_PVAL_AUC <- 2000 # Replicates for AUC bootstrap test

# --- Helper Functions ---
format_metric_with_pval <- function(point_est, p_val) {
  if(is.null(point_est) || is.na(point_est)) return("N/A")
  p_val_str <- if(is.null(p_val) || is.na(p_val)) "" else if(p_val < 0.001) "(p<0.001)" else sprintf("(p=%.3f)", p_val)
  sprintf("%.4f %s", point_est, p_val_str)
}

generate_colors <- function(n) {if(n==0)return(character(0));if(requireNamespace("RColorBrewer",quietly=T)&&n>0)RColorBrewer::brewer.pal(max(3,n),"Set2")[1:n] else if (n>0&&n<=8)(1:n) else if (n>0) rainbow(n) else character(0)}
create_k_folds <- function(data_size, k = 5) {
  if (k <= 1 || k > data_size) stop("Number of folds k must be > 1 and <= data_size.")
  shuffled_indices <- sample(1:data_size); fold_indices <- cut(1:data_size, breaks=k, labels=F)
  folds <- list(); for (i in 1:k) folds[[i]] <- shuffled_indices[which(fold_indices == i)]; return(folds)
}


# --- 1. DATA LOADING AND FORMULA DEFINITION ---
cat("\n--- 1. Data Loading and Formula Definition ---\n")
data_file_path <- "../Data/heart_cleveland_upload.csv"
if (!file.exists(data_file_path)) stop(paste("Data file not found at:", data_file_path))
cad_data <- read.table(data_file_path, header = TRUE, sep = ",")
if (!"condition" %in% colnames(cad_data)) stop("Response variable 'condition' not found.")
formula <- condition ~ age + sex + trestbps + chol + thalach + oldpeak
cat("Data loaded and formula defined.\n")
coeff_names_to_plot <- c("age", "sex", "trestbps", "chol", "thalach", "oldpeak")

# --- 2. DEFINE ALL PRIOR SETS ---
cat("\n--- 2. Define All Prior Sets ---\n")
default_prior_intercept_mean <- 0; default_prior_intercept_prec <- 1 / (5^2)
prior_sets <- list()
prior_sets$Claude_Mod <- list(means = list(`(Intercept)`=default_prior_intercept_mean,age=0.06,sex=0.7,trestbps=0.02,chol=0.008,thalach=-0.01,oldpeak=0.5), precs = list(`(Intercept)`=default_prior_intercept_prec,age=1/(0.02^2),sex=1/(0.3^2),trestbps=1/(0.01^2),chol=1/(0.005^2),thalach=1/(0.01^2),oldpeak=1/(0.2^2)))
prior_sets$Claude_Weak <- list(means = list(`(Intercept)`=default_prior_intercept_mean,age=0.04,sex=0.5,trestbps=0.01,chol=0.005,thalach=-0.01,oldpeak=0.3), precs = list(`(Intercept)`=default_prior_intercept_prec,age=1/(0.04^2),sex=1/(0.6^2),trestbps=1/(0.02^2),chol=1/(0.01^2),thalach=1/(0.02^2),oldpeak=1/(0.4^2)))
prior_sets$Gemini_Mod <- list(means = list(`(Intercept)`=default_prior_intercept_mean,age=0.067,sex=0.69,trestbps=0.01,chol=0.005,thalach=-0.02,oldpeak=0.588), precs = list(`(Intercept)`=default_prior_intercept_prec,age=1/(0.021^2),sex=1/(0.173^2),trestbps=1/(0.005^2),chol=1/(0.0025^2),thalach=1/(0.01^2),oldpeak=1/(0.163^2)))
prior_sets$Gemini_Weak <- list(means = list(`(Intercept)`=default_prior_intercept_mean,age=0.0,sex=0.0,trestbps=0.0,chol=0.0,thalach=0.0,oldpeak=0.0), precs = list(`(Intercept)`=default_prior_intercept_prec,age=1/(0.1^2),sex=1/(1.0^2),trestbps=1/(0.05^2),chol=1/(0.02^2),thalach=1/(0.05^2),oldpeak=1/(1.0^2)))
prior_sets$ChatGPT_Mod <- list(means = list(`(Intercept)`=default_prior_intercept_mean,age=0.05,sex=0.7,trestbps=0.002,chol=0.002,thalach=-0.02,oldpeak=0.4), precs = list(`(Intercept)`=default_prior_intercept_prec,age=1/(0.02^2),sex=1/(0.3^2),trestbps=1/(0.005^2),chol=1/(0.004^2),thalach=1/(0.01^2),oldpeak=1/(0.3^2)))
prior_sets$ChatGPT_Weak <- list(means = list(`(Intercept)`=default_prior_intercept_mean,age=0,sex=0,trestbps=0,chol=0,thalach=0,oldpeak=0), precs = list(`(Intercept)`=default_prior_intercept_prec,age=1/(2.5^2),sex=1/(2.5^2),trestbps=1/(2.5^2),chol=1/(2.5^2),thalach=1/(2.5^2),oldpeak=1/(2.5^2)))
cat(length(prior_sets), "prior sets defined.\n")
get_sds_from_precs <- function(precs_list) sapply(precs_list, function(p) 1/sqrt(p))

# --- 3. FIT FULL MODELS (for plotting posteriors) ---
cat("\n--- 3. Fitting Full Models (for parameter plots) ---\n")
full_glm_model <- glm(formula, data = cad_data, family = binomial(link = "logit"))
if (!full_glm_model$converged) warning("Full GLM model did not converge.")
coeffs_glm <- summary(full_glm_model)$coefficients
cat("Frequentist GLM model (full data) fitted.\n")

inla_models_list_fullfit <- list()
cat("Fitting INLA models (full data, all priors) for posterior distributions...\n")
for (prior_name in names(prior_sets)) {
  current_priors_inla_format <- list(mean=prior_sets[[prior_name]]$means, prec=prior_sets[[prior_name]]$precs)
  tryCatch({
    inla_models_list_fullfit[[prior_name]] <- INLA::inla(formula,family="binomial",data=cad_data,control.fixed=current_priors_inla_format,control.compute=list(config=TRUE))
  }, error=function(e) message(paste("Error INLA full fit", prior_name, ":", e$message)))
  if(is.null(inla_models_list_fullfit[[prior_name]])) warning(paste("Full INLA fit failed for", prior_name))
}
cat("Full INLA models fitted for posterior plots.\n")

# --- 4. 5-FOLD CV ANALYSIS (FULL DATASET) ---
cat("\n--- 4. 5-Fold CV Analysis (Full Dataset) ---\n")
k_folds_cv <- 5
folds_full_data <- create_k_folds(nrow(cad_data), k = k_folds_cv)
performance_results_list <- list()
y_obs_full <- cad_data$condition

# 4.1 Frequentist GLM: 5-fold CV (Run this first to have baseline predictions)
cat("Performing 5-fold CV for GLM (full dataset)...\n")
phat_kfold_glm <- numeric(nrow(cad_data))
for (k_f in 1:k_folds_cv) {
  test_indices <- folds_full_data[[k_f]]; train_indices <- setdiff(1:nrow(cad_data), test_indices)
  glm_kfold <- glm(formula, data=cad_data[train_indices,], family=binomial(link="logit"))
  phat_kfold_glm[test_indices] <- predict(glm_kfold, newdata=cad_data[test_indices,], type="response")
}

brier_glm <- mean((y_obs_full - phat_kfold_glm)^2, na.rm=T)
eps <- .Machine$double.eps; phat_clamped_glm <- pmax(eps, pmin(1 - eps, phat_kfold_glm));
log_s_glm_ind <- ifelse(y_obs_full == 1, -log(phat_clamped_glm), -log(1 - phat_clamped_glm));
log_s_glm <- mean(log_s_glm_ind, na.rm=T)
auc_glm <- NA; roc_obj_glm <- NULL
if(requireNamespace("pROC", quietly=T) && length(unique(na.omit(y_obs_full))) > 1) {
    roc_obj_glm <- pROC::roc(response=y_obs_full, predictor=phat_kfold_glm, levels=c(0,1), direction="<", quiet=T)
    auc_glm <- as.numeric(pROC::auc(roc_obj_glm))
}
performance_results_list[["GLM_Frequentist"]] <- list(brier=brier_glm, log_score=log_s_glm, auc=auc_glm, p_brier=NA, p_log_score=NA, p_auc=NA, prior_set="N/A (Frequentist)", model_type="GLM")
cat("GLM 5-fold CV complete.\n")

# 4.2 INLA Models: 5-fold CV
cat("Performing 5-fold CV for INLA priors...\n")
for (prior_name in names(prior_sets)) {
  cat(paste("  Processing INLA 5-fold CV for prior set:", prior_name, "\n"))
  phat_kfold_inla <- numeric(nrow(cad_data))
  priors_current_kfold <- list(mean=prior_sets[[prior_name]]$means, prec=prior_sets[[prior_name]]$precs)
  pb_inla_kcv <- txtProgressBar(min=0, max=k_folds_cv, style=3)
  for (k_f in 1:k_folds_cv) {
    setTxtProgressBar(pb_inla_kcv, k_f)
    test_idx <- folds_full_data[[k_f]]; train_idx <- setdiff(1:nrow(cad_data), test_idx)
    data_kcv_train <- cad_data[train_idx,]; data_kcv_test_na <- cad_data[test_idx,]; data_kcv_test_na$condition <- NA
    data_kcv_comb <- rbind(data_kcv_train, data_kcv_test_na)
    inla_kcv_model <- NULL
    tryCatch({inla_kcv_model <- INLA::inla(formula, family="binomial", data=data_kcv_comb, control.fixed=priors_current_kfold, control.predictor=list(compute=T, link=1), verbose=F)
    }, error = function(e) message(paste("Error INLA kCV fold",k_f, prior_name, ":", e$message)))
    if(!is.null(inla_kcv_model) && !is.null(inla_kcv_model$summary.fitted.values)){
      pred_idx_kcv <- (nrow(data_kcv_train)+1):nrow(data_kcv_comb)
      phat_kfold_inla[test_idx] <- inla_kcv_model$summary.fitted.values[pred_idx_kcv, "mean"]
    } else { phat_kfold_inla[test_idx]<-NA }
  }
  close(pb_inla_kcv)

  brier_inla <- mean((y_obs_full - phat_kfold_inla)^2, na.rm=T)
  phat_clamped_inla <- pmax(eps, pmin(1 - eps, phat_kfold_inla));
  log_s_inla_ind <- ifelse(y_obs_full == 1, -log(phat_clamped_inla), -log(1 - phat_clamped_inla));
  log_s_inla <- mean(log_s_inla_ind, na.rm=T)

  # P-value Calculation using correctR
  brier_scores_inla <- (phat_kfold_inla - y_obs_full)^2
  brier_scores_glm <- (phat_kfold_glm - y_obs_full)^2
  phat_clamped_inla <- pmax(eps, pmin(1 - eps, phat_kfold_inla))
  phat_clamped_glm <- pmax(eps, pmin(1 - eps, phat_kfold_glm))
  log_scores_inla <- ifelse(y_obs_full == 1, -log(phat_clamped_inla), -log(1 - phat_clamped_inla))
  log_scores_glm <- ifelse(y_obs_full == 1, -log(phat_clamped_glm), -log(1 - phat_clamped_glm))

  valid_brier_indices <- !is.na(brier_scores_inla) & !is.na(brier_scores_glm)
  valid_log_indices <- !is.na(log_scores_inla) & !is.na(log_scores_glm)

  p_brier <- NA
  if (sum(valid_brier_indices) > 1) {
      p_brier <- kfold_ttest(x = brier_scores_inla[valid_brier_indices], y = brier_scores_glm[valid_brier_indices], n = nrow(cad_data), k = k_folds_cv, tailed = "one", greater = "y")$p.value
  }
  p_log_score <- NA
  if (sum(valid_log_indices) > 1) {
      p_log_score <- kfold_ttest(x = log_scores_inla[valid_log_indices], y = log_scores_glm[valid_log_indices], n = nrow(cad_data), k = k_folds_cv, tailed = "one", greater = "y")$p.value
  }

  # Paired p-value for AUC
  auc_inla <- NA; p_auc <- NA
  if(!is.null(roc_obj_glm)){
      tryCatch({
          roc_obj_inla <- pROC::roc(response=y_obs_full, predictor=phat_kfold_inla, levels=c(0,1), direction="<", quiet=T)
          auc_inla <- as.numeric(pROC::auc(roc_obj_inla))
          roc_test_res <- pROC::roc.test(roc_obj_inla, roc_obj_glm, method="bootstrap", paired=TRUE, alternative="greater", n.boot=N_BOOT_PVAL_AUC, progress="none")
          p_auc <- roc_test_res$p.value
      }, error=function(e_auc){warning(paste("AUC p-value test failed for", prior_name), call.=F)})
  }

  performance_results_list[[prior_name]] <- list(brier=brier_inla, log_score=log_s_inla, auc=auc_inla, p_brier=p_brier, p_log_score=p_log_score, p_auc=p_auc, prior_set=prior_name, model_type="INLA")
}


# --- 5. PLOTS ---
cat("\n--- 5. Generating Requested Plots (Heart Disease) ---\n")
prior_categories <- list(Moderate = grep("_Mod$", names(prior_sets), value = TRUE), Weak = grep("_Weak$", names(prior_sets), value = TRUE))

color_palette <- c(
  "Claude_Mod" = "#E41A1C",  # Red for Claude
  "Claude_Weak" = "#E41A1C",
  "Gemini_Mod" = "#4DAF4A",  # Green for Gemini
  "Gemini_Weak" = "#4DAF4A",
  "ChatGPT_Mod" = "#377EB8", # Blue for ChatGPT
  "ChatGPT_Weak" = "#377EB8"
)

linetype_palette <- c(
  "Claude_Mod" = 2,   # Dashed line for Claude
  "Claude_Weak" = 2,
  "Gemini_Mod" = 3,   # Dotted line for Gemini
  "Gemini_Weak" = 3,
  "ChatGPT_Mod" = 4,  # Dot-dashed line for ChatGPT
  "ChatGPT_Weak" = 4
)


# --- Plot 3a & 3b: GLM vs. Moderate/Weak Priors ---
cat("Generating combined prior comparison plots...\n")
for (category_name in names(prior_categories)) {
  prior_names_in_category <- prior_categories[[category_name]]
  if(length(prior_names_in_category) == 0) next

  combined_plot_filename <- paste0("heart_dist_glm_vs_", tolower(category_name), "_priors.png")
  png(combined_plot_filename, width = 11, height = 15, units = "in", res = 300)
  op_comb <- par(mfrow=c(3, 2), mar=c(5, 5, 2, 2), oma=c(0, 0, 1, 0))

  for (coeff_name in coeff_names_to_plot) {
    glm_est <- coeffs_glm[coeff_name, "Estimate"]; glm_se <- coeffs_glm[coeff_name, "Std. Error"]
    x_ranges_calc <- list(); y_maxes_calc <- c()
    if(!is.na(glm_est) && !is.na(glm_se) && glm_se > 0) x_ranges_calc[[1]]<-c(glm_est - 4*glm_se, glm_est + 4*glm_se)
    for(p_name in prior_names_in_category) {
      pr_m <- prior_sets[[p_name]]$means[[coeff_name]]; pr_s <- 1/sqrt(prior_sets[[p_name]]$precs[[coeff_name]])
      if(!is.na(pr_m) && !is.na(pr_s) && pr_s > 0) {
        if (category_name == "Weak" && p_name == "ChatGPT_Weak") {} else { x_ranges_calc[[length(x_ranges_calc)+1]] <- c(pr_m - 4*pr_s, pr_m + 4*pr_s) }
      }
    }
    plot_xlim_comb <- if(length(x_ranges_calc)>0) range(unlist(x_ranges_calc),na.rm=T) else c(-1,1)
    x_seq_for_ylim <- seq(plot_xlim_comb[1], plot_xlim_comb[2], length.out=200)
    if(!is.na(glm_est) && !is.na(glm_se) && glm_se>0) y_maxes_calc <- c(y_maxes_calc, max(dnorm(x_seq_for_ylim, glm_est, glm_se),na.rm=T))
    for(p_name in prior_names_in_category) { pr_m <- prior_sets[[p_name]]$means[[coeff_name]]; pr_s <- 1/sqrt(prior_sets[[p_name]]$precs[[coeff_name]]); if(!is.na(pr_m) && !is.na(pr_s) && pr_s > 0) y_maxes_calc <- c(y_maxes_calc, max(dnorm(x_seq_for_ylim, pr_m, pr_s),na.rm=T)) }
    plot_ylim_comb <- if(length(y_maxes_calc)>0 && any(is.finite(y_maxes_calc))) c(0,max(y_maxes_calc[is.finite(y_maxes_calc)],na.rm=T)*1.2) else c(0,1)

    plot(0, type="n", xlim=plot_xlim_comb, ylim=plot_ylim_comb, xlab=coeff_name, ylab="Density", cex.lab=1.6, cex.axis=1.4)
    grid(col="lightgray",lty="dotted"); legend_items_comb <- list()
    
    if(!is.na(glm_est) && !is.na(glm_se) && glm_se > 0){
        curve(dnorm(x,glm_est,glm_se),add=T,col="black",lwd=2,lty=1)
        legend_items_comb[["MLE"]]<-list(col="black",lty=1,lwd=2)
    }
    
    for (p_name in prior_names_in_category) {
      pr_m <- prior_sets[[p_name]]$means[[coeff_name]]; pr_s <- 1/sqrt(prior_sets[[p_name]]$precs[[coeff_name]])
      current_color <- color_palette[p_name]
      current_linetype <- linetype_palette[p_name]
      
      if(!is.na(pr_m) && !is.na(pr_s) && pr_s > 0){
          curve(dnorm(x,pr_m,pr_s),add=T,col=current_color,lwd=2,lty=current_linetype)
          legend_items_comb[[p_name]]<-list(col=current_color,lty=current_linetype,lwd=2)
      }
    }
    if(length(legend_items_comb)>0) legend("topright",names(legend_items_comb),col=sapply(legend_items_comb,`[[`,"col"),lty=sapply(legend_items_comb,`[[`,"lty"),lwd=sapply(legend_items_comb,`[[`,"lwd"),cex=1.2,bty="n")
  }
  
  dev.off(); par(op_comb)
  cat(paste("Combined prior plot for", category_name, "category saved to:", combined_plot_filename, "\n"))
}

# --- Plot 3c: Detailed Prior vs. Likelihood vs. Posterior for EACH prior set ---
cat("Generating detailed comparison plots for each prior set...\n")
for (prior_name in names(prior_sets)) {
  detailed_plot_filename <- paste0("heart_dist_detailed_comparison_", prior_name, ".png")
  png(detailed_plot_filename, width = 11, height = 15, units = "in", res = 300)
  
  op_detail <- par(mfrow = c(3, 2), mar = c(5.1, 5.1, 2.1, 2.1), oma = c(0, 0, 1, 0))
  
  current_inla_model <- inla_models_list_fullfit[[prior_name]]
  current_priors_means <- prior_sets[[prior_name]]$means
  current_priors_sds <- get_sds_from_precs(prior_sets[[prior_name]]$precs)

  if(is.null(current_inla_model)){
    plot(0, type="n", axes=FALSE, xlab="", ylab="", main=paste("INLA model for", prior_name, "missing or failed."))
    dev.off()
    par(op_detail)
    next
  }
  for (coeff_name in coeff_names_to_plot) {
    prior_m <- current_priors_means[[coeff_name]]
    prior_s <- current_priors_sds[[coeff_name]]
    inla_marg <- current_inla_model$marginals.fixed[[coeff_name]]
    glm_est <- coeffs_glm[coeff_name, "Estimate"]
    glm_se <- coeffs_glm[coeff_name, "Std. Error"]
    
    plot_xlim <- NULL 

    is_claude_prior <- prior_name == "Claude_Mod" || prior_name == "Claude_Weak"
    
    if (is_claude_prior && coeff_name == "sex") {
        plot_xlim <- c(-1.5, 3.5)
    } else if (is_claude_prior && coeff_name == "thalach") {
        plot_xlim <- c(-0.08, 0.02)
    }

    if (is.null(plot_xlim)) {
        x_ranges <- list()
        if (!is.null(inla_marg)) x_ranges[[length(x_ranges)+1]] <- range(inla_marg[,1])
        if (!is.na(glm_est) && !is.na(glm_se) && glm_se > 0) x_ranges[[length(x_ranges)+1]] <- c(glm_est-4*glm_se, glm_est+4*glm_se)
        
        weak_prior_names <- c("Claude_Weak", "Gemini_Weak", "ChatGPT_Weak")
        if (!prior_name %in% weak_prior_names) {
            if(!is.na(prior_m) && !is.na(prior_s) && prior_s > 0) x_ranges[[length(x_ranges)+1]] <- c(prior_m-4*prior_s, prior_m+4*prior_s)
        }
        plot_xlim <- if(length(x_ranges)>0) range(unlist(x_ranges), na.rm=TRUE) else c(-1,1)
    }
    
    y_maxes <- c()
    x_seq_for_ylim <- seq(plot_xlim[1], plot_xlim[2], length.out = 200)
    if(!is.na(prior_m) && !is.na(prior_s) && prior_s > 0) y_maxes <- c(y_maxes, max(dnorm(x_seq_for_ylim, prior_m, prior_s), na.rm=TRUE))
    if (!is.null(inla_marg)) { inla_func <- approxfun(inla_marg[,1], inla_marg[,2], rule=2); y_maxes <- c(y_maxes, max(inla_func(x_seq_for_ylim), na.rm=TRUE)) }
    if (!is.na(glm_est) && !is.na(glm_se) && glm_se > 0) y_maxes <- c(y_maxes, max(dnorm(x_seq_for_ylim, glm_est, glm_se), na.rm=TRUE))
    plot_ylim <- if(length(y_maxes)>0 && any(is.finite(y_maxes))) c(0, max(y_maxes[is.finite(y_maxes)], na.rm=T)*1.2) else c(0,1)

    plot(0, type="n", xlim=plot_xlim, ylim=plot_ylim, xlab=coeff_name, ylab="Density",
         cex.lab=1.6, cex.axis=1.4)
    grid(col="lightgray", lty="dotted"); legend_items <- list()
    if(!is.na(prior_m) && !is.na(prior_s) && prior_s > 0){curve(dnorm(x,prior_m,prior_s),add=TRUE,col="darkolivegreen4",lwd=2,lty=3); legend_items$Prior<-list(col="darkolivegreen4",lty=3,lwd=2)}
    if(!is.null(inla_marg)){lines(inla_marg,col="deepskyblue3",lwd=2); legend_items$Posterior<-list(col="deepskyblue3",lty=1,lwd=2)}
    if(!is.na(glm_est) && !is.na(glm_se) && glm_se>0){
        curve(dnorm(x,glm_est,glm_se),add=TRUE,col="black",lwd=2,lty=1)
        legend_items[["MLE"]]<-list(col="black",lty=1,lwd=2)
    }
    if(length(legend_items) > 0) legend("topright", names(legend_items),col=sapply(legend_items,`[[`,"col"),lty=sapply(legend_items,`[[`,"lty"),lwd=sapply(legend_items,`[[`,"lwd"),cex=1.2,bty="n")
  }
  dev.off()
  par(op_detail)
}
cat("Detailed comparison plots for each prior set saved.\n")


# --- 6. COMPILE AND PRINT RESULTS TABLE ---
cat("\n--- 6. Summary Table of 5-Fold CV Predictive Performance ---\n")
results_df_hd_5cv <- do.call(rbind, lapply(names(performance_results_list), function(name) {
  x <- performance_results_list[[name]]
  data.frame(
    ModelName = name, ModelType = x$model_type, PriorSetDescription = x$prior_set,
    BrierScore = format_metric_with_pval(x$brier, x$p_brier),
    MeanNegLogScore = format_metric_with_pval(x$log_score, x$p_log_score),
    AUC = format_metric_with_pval(x$auc, x$p_auc),
    stringsAsFactors = FALSE
  )
}))
print(results_df_hd_5cv)
results_csv_filename_hd_5cv <- "kfold_cv_heart_disease_performance_summary.csv"
tryCatch({ write.csv(results_df_hd_5cv, results_csv_filename_hd_5cv, row.names = FALSE)
           cat(paste("\nHeart Disease 5-Fold CV performance summary saved to:", results_csv_filename_hd_5cv, "\n"))
         }, error=function(e){cat(paste("\nError saving Heart Disease CV summary to CSV:",e$message,"\n"))})


# --- 10. GENERATE COMBINED KL DIVERGENCE SUMMARY TABLE ---
cat("\n--- 10. Generating Combined KL Divergence Summary Table ---\n")
kl_divergence_norm <- function(mu1, sd1, mu2, sd2) {
  if(is.na(sd1) || is.na(sd2) || sd1 <= 0 || sd2 <= 0) return(NA)
  log(sd2 / sd1) + (sd1^2 + (mu1 - mu2)^2) / (2 * sd2^2) - 0.5
}

glm_coeffs <- summary(full_glm_model)$coefficients
glm_means <- glm_coeffs[, "Estimate"]
glm_sds   <- glm_coeffs[, "Std. Error"]

kl_divergence_table <- data.frame(matrix(
    NA,
    nrow = length(names(prior_sets)),
    ncol = length(coeff_names_to_plot),
    dimnames = list(names(prior_sets), coeff_names_to_plot)
))

cat("Calculating KL Divergence for each prior and predictor variable...\n")
for (prior_name in names(prior_sets)) {
  for (coeff in coeff_names_to_plot) {
    prior_mean <- prior_sets[[prior_name]]$means[[coeff]]
    prior_sd   <- 1/sqrt(prior_sets[[prior_name]]$precs[[coeff]])
    
    kl_divergence_table[prior_name, coeff] <- kl_divergence_norm(
      mu1 = glm_means[coeff], sd1 = glm_sds[coeff],
      mu2 = prior_mean,       sd2 = prior_sd
    )
  }
}

cat("Calculating summary metrics (Average KL and Average Rank)...\n")
kl_ranks_table <- apply(kl_divergence_table, 2, rank, ties.method = "min")
kl_divergence_table$`Average KL Divergence` <- rowMeans(kl_divergence_table, na.rm = TRUE)
kl_divergence_table$`Average Rank` <- rowMeans(kl_ranks_table, na.rm = TRUE)
kl_divergence_table <- kl_divergence_table[order(kl_divergence_table$`Average KL Divergence`), ]
kl_divergence_table <- cbind(`Prior Set` = rownames(kl_divergence_table), kl_divergence_table)
rownames(kl_divergence_table) <- NULL

cat("\n--- Combined KL Divergence Table: KL(Likelihood || Prior) ---\n")
formatted_table <- kl_divergence_table
cols_to_format <- names(formatted_table)[-1]
formatted_table[cols_to_format] <- lapply(formatted_table[cols_to_format], function(x) {
    format(round(x, 4), nsmall = 4)
})
print(formatted_table)

kl_table_filename <- "kl_divergence_summary_combined.csv"
tryCatch({
  write.csv(kl_divergence_table, kl_table_filename, row.names = FALSE)
  cat(paste("\nCombined KL Divergence summary table saved to:", kl_table_filename, "\n"))
}, error = function(e){
  cat(paste("\nError saving combined KL Divergence summary to CSV:", e$message, "\n"))
})

end_time <- Sys.time()
cat("\n--- End of Script (Heart Disease Example) ---")
cat("\nTotal execution time:", format(end_time - start_time), "\n")

