# code created based on 
# https://github.com/marije-sluiskes/AccelerAge_SimulationStudy/
# and 
# https://github.com/Phil-Ber/Thesis_AccelerAge


################################################################################
#### Set-up 
################################################################################

# fall-back directory for testing and debugigng
PLOT_DIR    <- if (exists("PLOT_DIR",     inherits = TRUE)) PLOT_DIR     else "output/scenarios/bayesian/plots"
SUMMARY_DIR <- if (exists("SUMMARY_DIR",  inherits = TRUE)) SUMMARY_DIR  else "output/scenarios/bayesian/summaries"
EXP_DIR     <- if (exists("EXP_DIR",      inherits = TRUE)) EXP_DIR      else "output/scenarios/bayesian/experiments"
DIAG_DIR_HS <- if (exists("DIAG_DIR_HS",  inherits = TRUE)) DIAG_DIR_HS  else "output/scenarios/bayesian/diagnostics/hs"
DIAG_DIR_GS <- if (exists("DIAG_DIR_GS",  inherits = TRUE)) DIAG_DIR_GS  else "output/scenarios/bayesian/diagnostics/gs"

for (d in c(PLOT_DIR, SUMMARY_DIR, EXP_DIR, DIAG_DIR_HS, DIAG_DIR_GS)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

################################################################################
#### Packages & helper functions
################################################################################
library(eha)
library(brms)
library(posterior)
library(cmdstanr)
library(survival)
library(glue)
library(dplyr)
library(ggplot2)
library(mcmcse)


# ============================================================
# Generate survival data
# ============================================================

# get helper functions
source("1_dgp_helpers.R")

gen_data = function(
    n_obs,
    n_pop,
    p, 
    g = 10,
    rho_w,
    non_zero_groups,
    mu_l = -4, mu_u = 4,
    target_snr,
    gomp_a = exp(-9), 
    gomp_b = 0.085,
    ltname,
    standardise_X = T,
    X_plots = F,
    betas_plot = F,
    lt_plot = F,
    force_recalc,
    seed = 123) {
  
  
  cat("\nGenerating βs...\n")
  betas <- gen_betas(p = p,
                     g = g, 
                     non_zero_groups,
                     mu_l = mu_l, mu_u = mu_u,
                     rho_w = rho_w,
                     target_snr = target_snr,
                     gomp_a = gomp_a, 
                     gomp_b = gomp_b,
                     standardise_X = standardise_X,
                     betas_plot = betas_plot, 
                     seed = seed)
  print("Done!")
  print("----------")
  
  print("Generating X...")
  X_obj <- gen_X(n = n_pop, 
                 p = p, 
                 g = g, 
                 rho_w = rho_w,
                 standardise_X = standardise_X,
                 X_plots = X_plots,
                 seed = seed)
  print("Done!")
  print("----------")
  
  
  print("Generating Gompertz population lifetable...")
  lt <- gen_pop_lt_gompertz(n_pop = n_pop,
                            p = p,
                            non_zero_groups = non_zero_groups,
                            g = g,
                            rho_w = rho_w,
                            gomp_a = gomp_a,
                            gomp_b = gomp_b,
                            betas = betas$beta,
                            X_pop = X_obj$X,
                            force_recalc = force_recalc,
                            filename = ltname,
                            lt_plot = lt_plot,
                            seed = seed
                            )
  print("Done!")
  print("----------")
  
  
  plot = ggplot(lt, aes(x = t, y = mrl_pop)) +
    geom_line() +
    labs(title = "Gompertz Population Lifetable MRL",
         x = "Chronological age (years)", y = "Mean Residual Life") +
    theme_minimal()
  print(plot)
 
  
  print("Generating dataset from Gompertz distribution...")
  df_sim = create_dataset_gompertz(n_obs = n_obs,
                                   p = p,
                                   gomp_a = gomp_a,
                                   gomp_b = gomp_b,
                                   g = g,
                                   rho_w = rho_w,
                                   betas = betas$beta, # true betas
                                   followup = 20,
                                   standardise_X = standardise_X,
                                   X_plots = X_plots,
                                   lt = lt,
                                   seed = seed)
  print("Done!")
  print("----------")
  
  return(list(df = df_sim, true_betas = betas, lt = lt))
}



# ============================================================
# Frequentist AFT 
# ============================================================

fAFT_gomp <- function(df_sim, p, save_plots = FALSE) {
  
  true_betas = df_sim$true_betas
  lt         = df_sim$lt
  df_sim     = df_sim$df
  
  
  # train/test split 50/50
  train_idx    = sample(1:nrow(df_sim), nrow(df_sim)/2)
  df_sim_train = df_sim[train_idx,]
  df_sim_test  = df_sim[-train_idx,]
  
  
  # stop if the dataset is too small relative to the number of predictors
  if (nrow(df_sim_train) <= p + 1) return(NULL)
  
  
  # select the first M predictors (V1,...VM)
  x_names <- colnames(df_sim_train)[seq_len(p)]
  
  # include in the model only the first M predictors
  # 1=event ; 0=right-censored
  formula_fAFT <- as.formula(
    paste0("Surv(entry_age, age_end, status) ~ ",
           paste(x_names, collapse = " + ")
    ) )
  
  
  
  # -------------------
  # ---- FIT MODEL ----
  # -------------------
  
  fit_fAFT = tryCatch(
    aftreg(
      formula = formula_fAFT, 
      data    = df_sim_train, 
      dist    = "gompertz"), 
    error     = function(e) { 
      message("Model fitting error: ", e$message)
      NULL
    })
  # stop if model fit failed
  if (is.null(fit_fAFT)) return(NULL)
  
  
  print("----------")
  
  
  # extract estimated parameters 
  sigma_hat <- exp(fit_fAFT$coefficients["log(scale)"])
  cat("Scale/Sigma estimate:", sigma_hat, "| true:", 1/0.085, "\n")
  
  tau_hat <- exp(fit_fAFT$coefficients["log(shape)"])
  cat("Shape/Tau estimate:", tau_hat, "| true:", exp(-9)/0.085, "\n")
  
  # recover gomp_a, gomp_b
  gomp_a_hat <- tau_hat / sigma_hat
  gomp_b_hat <- 1 / sigma_hat
  
  # extract estimates betas with AFT model
  beta_hat <- fit_fAFT$coefficients[x_names]
  
  # get predictors from test set
  X_test <- as.matrix(df_sim_test[, x_names, drop = FALSE])
  # compute linear predictor
  linpred_hat <- drop(X_test %*% beta_hat) 
  
  
  
  # ----------------------
  # ---- ESTIMATE MRL ----
  # ----------------------
  
  pred_mrl <- numeric(nrow(df_sim_test))
  for (i in 1:nrow(df_sim_test)){
    
    # conditional survival at entry 
    s_cond <- gomp_baseline_surv(t = df_sim_test$entry_age[i] * exp(linpred_hat[i]), 
                                 gomp_a = gomp_a_hat, 
                                 gomp_b = gomp_b_hat)
    
    # integral on scaled-time axis
    t_unadj <- integrate(gomp_baseline_surv,
                         lower = df_sim_test$entry_age[i] * exp(linpred_hat[i]),
                         upper = Inf,
                         gomp_a = gomp_a_hat, gomp_b = gomp_b_hat)$value
    
    # convert back to original time scale
    t_adj <- t_unadj / exp(linpred_hat[i])
    
    # mean residual life 
    pred_mrl[i] <- (t_adj / s_cond) 
  }
  
  # uncomment to check the ranges
  # cat("pred_mrl range:", range(pred_mrl), "\n")
  # cat("lt$mrl_pop range:", range(lt$mrl_pop), "\n")
  
 
  print("----------")
  
  
  # --------------------------
  # ---- ESTIMATE BIO AGE ----
  # --------------------------
  
  # map predicted mrl to predicted bio age (via true life table)
  b_hat <- numeric(length(pred_mrl))
  for (i in 1:length(pred_mrl)) {
    b_hat[i] <- lt$t[ which.min(abs(lt$mrl_pop - pred_mrl[i])) ]
  }
  

  # ----------------------
  # ---- COMPUTE RMSE ----
  # ----------------------
  
  # TRUE RMSE: uses all values without capping
  rmse_b_all    <- sqrt(mean((b_hat    - df_sim_test$b  )^2))
  rmse_mrl_all  <- sqrt(mean((pred_mrl - df_sim_test$mrl)^2))
  rmse_beta_hat <- sqrt(mean((beta_hat - true_betas$beta)^2))
  
  # CAPPED RMSE: both predicted and true values are capped at 150 before squaring the difference.
  # When either or both sides of a pair exceed 150, each is clamped to 150,
  # to keep impossible human ages from distorting the metric while retaining all observations.
  rmse_b   <- sqrt(mean((pmin(b_hat,    150) - pmin(df_sim_test$b,   150))^2))
  rmse_mrl <- sqrt(mean((pmin(pred_mrl, 150) - pmin(df_sim_test$mrl, 150))^2))
  
  
  # save plots
  file_tag <- if (exists("tag_for_saving", inherits = FALSE) && 
                  !is.null(tag_for_saving)) tag_for_saving else "run"
  
  # mrl prediction plot
  p_mrl <- df_sim_test %>%
    ggplot(aes(x = mrl, y = pred_mrl)) +
    geom_point() +
    geom_abline(color = "red") +
    labs(
      title = "True vs Predicted MRL (Frequentist AFT Gompertz)",
      subtitle = glue("RMSE = {round(rmse_mrl, 4)} | tag = {file_tag}")
    ) +
    theme_minimal()
  
  
  # optionally, save plots
  if (save_plots == TRUE) {
    ggsave(file.path(PLOT_DIR, glue("fAFT_mrl_{file_tag}.png")),
           p_mrl, width = 7.2, height = 5.2, dpi = 180, bg = "white")
  }
  
  
  true_betas$beta_hat = beta_hat
  # betas prediction plot
  p_beta <- true_betas %>%
    mutate(index = 1:p) %>%
    ggplot(aes(x = index)) +
    geom_segment(aes(xend = index, y = beta, yend = beta_hat), color = "gray") +
    geom_point(aes(y = beta, color = "True")) +
    geom_point(aes(y = beta_hat, color = "Predicted")) +
    scale_color_manual(values = c("True" = "black", "Predicted" = "skyblue"), name = "Type") +
    labs(
      title = "True vs Predicted β (Frequentist AFT Gompertz)",
      subtitle = glue("RMSE = {round(rmse_beta_hat, 4)} | tag = {file_tag}"),
      x = "Index",
      y = "β"
    ) +
    theme_minimal()
  
  
  # optionally, save plots
  if (save_plots == TRUE) {
    ggsave(file.path(PLOT_DIR, glue("fAFT_betas_{file_tag}.png")),
           p_beta, width = 8.2, height = 4.8, dpi = 180, bg = "white")
  }
  
  
  return(list(
    test_data         = df_sim_test,
    true_betas        = true_betas,
    predicted_betas   = beta_hat,
    predicted_mrl     = pred_mrl,
    predicted_bio_age = b_hat,
    rmses = list(coef        = rmse_beta_hat,
                 mrl         = rmse_mrl,        # capped at 150
                 mrl_all     = rmse_mrl_all,    # all values, no capping
                 bio_age     = rmse_b,          # capped at 150
                 bio_age_all = rmse_b_all),     # all values, no capping
    method = "AFT Gompertz",
    distr_ests = list(shape = tau_hat, scale = sigma_hat)
  ))
}



# ============================================================
# Bayesian AFT 
# ============================================================

# get Stan functions, custom family, and priors
source("2_stan_family_priors.R")

# function to extract brms diagnostics
brms_dgs <- function(fit_bAFT) {
  # return NAs if fit is missing
  if (is.null(fit_bAFT)) {
    return(list(
      max_rhat     = NA_real_,
      min_ess_bulk = NA_real_,
      min_ess_tail = NA_real_,
      n_divergent  = NA_integer_
    ))
  }
  
  # extract posterior draws in an array format
  draws_arr <- as_draws_array(fit_bAFT)
  
  # get per-variable diagnostics
  diag_df <- summarise_draws(draws_arr, "rhat", "ess_bulk", "ess_tail")
  
  # rhat()/ess_bulk() expect a single variable
  # check documentation: https://mc-stan.org/posterior/reference/rhat.html
  
  max_rhat     <- max(diag_df$rhat,     na.rm = TRUE)
  min_ess_bulk <- min(diag_df$ess_bulk, na.rm = TRUE)
  min_ess_tail <- min(diag_df$ess_tail, na.rm = TRUE)
  
  if (!is.finite(max_rhat))     max_rhat     <- NA_real_
  if (!is.finite(min_ess_bulk)) min_ess_bulk <- NA_real_
  if (!is.finite(min_ess_tail)) min_ess_tail <- NA_real_
  
  # check divergent transitions from NUTS params
  np    <- try(nuts_params(fit_bAFT), silent = TRUE)
  n_div <- if (inherits(np, "try-error")) {
    NA_integer_
  } else {
    sum(np$Parameter == "divergent__" & np$Value == 1, na.rm = TRUE)
  }
  
  list(
    max_rhat     = max_rhat,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail,
    n_divergent  = as.integer(n_div)
  )
  
}


# function for Bayesian AFT model
bAFT_gomp = function(df_sim, 
                     p, 
                     g,
                     non_zero_groups,  # passed to priors_hs() for adaptive par_ratio
                     seed = 123, 
                     prior = c("hs", "gs"),
                     iter = 2000, 
                     chains = 4,
                     warmup = iter/2,
                     save_plots = FALSE
                     ) {
  
  set.seed(seed)

  
  # from gen_data() output
  true_betas = df_sim$true_betas
  lt         = df_sim$lt
  df_sim     = df_sim$df
  
  # train/test split 50/50
  train_idx    = sample(1:nrow(df_sim), nrow(df_sim)/2)
  df_sim_train = df_sim[train_idx,]
  df_sim_test  = df_sim[-train_idx,]
  
  
  # -------------------
  # ---- FIT MODEL ----
  # -------------------
  
  # select the first M predictors (V1,...VM)
  x_names <- colnames(df_sim_train)[seq_len(p)]
  
  linpred_formula <- as.formula(
    paste0("linpred ~ 0 + ", paste(x_names, collapse = " + ")) # instead of 1 bc multiple intercepts
  )
  
  # brms: 1 = right-censored ; 0 = event (check R documentation)
  formula_bAFT <- bf(
    age_end | cens(1 - status) + trunc(lb = entry_age) ~ 1,
    gamma ~ 1
  ) + linpred_formula 
  
  # uncomment to check if colnames(X) are correct in the formula
  # print(formula_bAFT)
  
  
  # check if prior argument is valid
  prior = match.arg(prior)
  
  # select prior specification
  priors <- switch(
    prior,
    hs = priors_hs(p, non_zero_groups),
    gs = priors_gs
  )
  
  # include additional Stan code only for Gaussian prior (to learn sigma from data)
  stanvars_used <- switch(
    prior,
    hs       = stanvars,              # use existing Stan definitions 
    gs       = make_stanvars_shared(  # generate additional Stan code to estimate variance 
    formula  = formula_bAFT,
    family   = gompertz_aft, 
    data     = df_sim_train,
    dpar     = "linpred",
    base_stanvars = stanvars, 
    include_intercept = FALSE 
    )
  )
  
  print("Fitting Bayesian Regression Model...")
  
  
  # fit Bayesian model
  fit_bAFT = brm(
    formula   = formula_bAFT,
    data      = df_sim_train,
    family    = gompertz_aft, # custom family 
    prior     = priors,
    iter      = iter, 
    chains    = chains,
    warmup    = warmup,
    control   = list(adapt_delta = 0.95, max_treedepth = 12),
    seed      = seed,
    backend   = "cmdstanr",
    #threads  = threading(2),
    save_pars = save_pars(all = FALSE),
    init      = 0, 
    stanvars  = stanvars_used
  )
  
  print("Done!")
  print("----------")
  
  
  # brms diagnostics
  dgs <- brms_dgs(fit_bAFT)
  print(dgs)
 
  # extract posterior samples 
  post <- as.data.frame(as_draws_df(fit_bAFT))
  
  # uncomment to check if names are correct
  # print(grep("linpred", colnames(post), value = TRUE))
  
  # # extract estimates 
  # mu_hat = exp(mean(post$b_Intercept)) # mu intercept
  # gamma_hat = exp(mean(post$b_gamma_Intercept)) # gamma intercept
  
  # get beta estimates
  beta_cols <- paste0("b_linpred_", x_names)
  brms_beta_hat <- colMeans(post[, beta_cols, drop = FALSE])
  
  # get predictors
  X_test <- as.matrix(df_sim_test[, x_names, drop = FALSE])
  
  # compute linear predictor 
  brms_linpred_hat <- drop(X_test %*% brms_beta_hat)
  
  # # create predicted a and b for Gomp AFT
  # brms_gomp_a_hat <- mu_hat
  # brms_gomp_b_hat <- gamma_hat
  
  # # create predicted sigma, tau 
  # brms_sigma_hat = 1/brms_gomp_b_hat
  # cat("Scale/Sigma estimate", (brms_sigma_hat), "\n")
  # brms_tau_hat = brms_gomp_a_hat / brms_gomp_b_hat
  # cat("Shape/Tau estimate:", brms_tau_hat, "\n")
  
  
  # --- (tau, sigma) parametrisation ---
  brms_tau_hat   <- exp(mean(post$b_Intercept))        
  brms_sigma_hat <- exp(mean(post$b_gamma_Intercept))

  # recover gomp_a, gomp_b
  brms_gomp_b_hat <- 1 / brms_sigma_hat
  brms_gomp_a_hat <- brms_tau_hat / brms_sigma_hat

  cat("Shape/Tau estimate:", brms_tau_hat, "| true:", exp(-9)/0.085, "\n")
  cat("Scale/Sigma estimate:", brms_sigma_hat, "| true:", 1/0.085, "\n")

  
  
  # ----------------------
  # ---- ESTIMATE MRL ----
  # ----------------------
  
  brms_pred_mrl <- numeric(nrow(df_sim_test))
  for (i in 1:nrow(df_sim_test)){
    
    # conditional survival at entry 
    s_cond <- gomp_baseline_surv(t = df_sim_test$entry_age[i] * exp(brms_linpred_hat[i]), 
                                 gomp_a = brms_gomp_a_hat, 
                                 gomp_b = brms_gomp_b_hat)
    
    # integral on scaled-time axis
    t_unadj <- integrate(gomp_baseline_surv,
                         lower = df_sim_test$entry_age[i] * exp(brms_linpred_hat[i]),
                         upper = Inf,
                         gomp_a = brms_gomp_a_hat, gomp_b = brms_gomp_b_hat)$value
    
    # convert back to original time scale (Jacobian)
    t_adj <- t_unadj / exp(brms_linpred_hat[i])
    
    # mean residual life 
    brms_pred_mrl[i] <- (t_adj / s_cond) 
  }
  
  
  
  # --------------------------
  # ---- ESTIMATE BIO AGE ----
  # --------------------------
  
  # map predicted mrl to bio age (via pop lt)
  brms_b_hat <- numeric(length(brms_pred_mrl))
  for (i in 1:length(brms_pred_mrl)) {
    brms_b_hat[i] <- lt$t[ which.min(abs(lt$mrl_pop - brms_pred_mrl[i])) ]
  }
  
  
  # ----------------------
  # ---- COMPUTE RMSE ----
  # ----------------------
  
  # TRUE RMSE: uses all values without capping
  rmse_brms_b_hat_all   <- sqrt(mean((brms_b_hat    - df_sim_test$b  )^2))
  rmse_brms_mrl_all     <- sqrt(mean((brms_pred_mrl - df_sim_test$mrl)^2))
  rmse_brms_beta_hat     <- sqrt(mean((brms_beta_hat - true_betas$beta)^2))
  
  # CAPPED RMSE: both predicted and true values are capped at 150 before squaring the difference.
  # When either or both sides of a pair exceed 150, each is clamped to 150,
  # to keep impossible human ages from distorting the metric while retaining all observations.
  rmse_brms_b_hat <- sqrt(mean((pmin(brms_b_hat,    150) - pmin(df_sim_test$b,   150))^2))
  rmse_brms_mrl   <- sqrt(mean((pmin(brms_pred_mrl, 150) - pmin(df_sim_test$mrl, 150))^2))
  
  
  # ---------------------------------------
  # ---- SIMULATION-IN-SIMULATION MCSE ----
  # ---------------------------------------
  
  beta_draws    <- as.matrix(post[, beta_cols, drop = FALSE])
  beta_mcse     <- apply(beta_draws, 2, function(x) mcse(x)$se)
  tau_draws     <- exp(post$b_Intercept)
  sigma_draws   <- exp(post$b_gamma_Intercept)
  tau_mcse      <- mcse(tau_draws)$se
  sigma_mcse    <- mcse(sigma_draws)$se
  max_beta_mcse <- max(beta_mcse, na.rm = TRUE)
  
  cat("Max beta MCSE:", round(max_beta_mcse, 6),
      "| tau MCSE:",   round(tau_mcse, 6),
      "| sigma MCSE:", round(sigma_mcse, 6), "\n")
  
  
  # save plots
  file_tag <- if (exists("tag_for_saving", inherits = FALSE) && 
                  !is.null(tag_for_saving)) tag_for_saving else "run"
  
  
  # mrl prediction plot
  p_mrl <- df_sim_test %>%
    ggplot(aes(x = mrl, y = brms_pred_mrl)) +
    geom_point() +
    geom_abline(color = "red") +
    labs(
      title = "True vs Predicted MRL (Bayesian brms AFT Gompertz)",
      subtitle = glue("RMSE = {round(rmse_brms_mrl, 4)} | prior = {prior} | tag = {file_tag}")
    ) +
    theme_minimal()
  
  
  # optionally, save plots
  if (save_plots == TRUE) {
    ggsave(file.path(PLOT_DIR, glue("bAFT_mrl_{prior}_{file_tag}.png")),
           p_mrl, width = 7.2, height = 5.2, dpi = 180, bg = "white")
  }
  
  
  
  # beta comparison plot
  true_betas$brms_beta_hat = brms_beta_hat
  
  p_beta <- true_betas %>%
    mutate(index = 1:p) %>%
    ggplot(aes(x = index)) +
    geom_segment(aes(xend = index, y = beta, yend = brms_beta_hat), color = "gray") +
    geom_point(aes(y = beta, color = "True")) +
    geom_point(aes(y = brms_beta_hat, color = "Predicted")) +
    scale_color_manual(values = c("True" = "black", "Predicted" = "skyblue"), name = "Type") +
    labs(
      title = "True vs Predicted β (Bayesian brms AFT Gompertz)",
      subtitle = glue("RMSE = {round(rmse_brms_beta_hat, 4)} | prior = {prior} | tag = {file_tag}"),
      x = "Index",
      y = "β"
    ) +
    theme_minimal()
  
  
  # optionally, save plots
  if (save_plots) {
    ggsave(file.path(PLOT_DIR, glue("bAFT_betas_{prior}_{file_tag}.png")),
           p_beta, width = 8.2, height = 4.8, dpi = 180, bg = "white")
  }
  
  
  return(list(
    test_data         = df_sim_test,
    true_betas        = true_betas,
    predicted_betas   = brms_beta_hat,
    predicted_mrl     = brms_pred_mrl,
    predicted_bio_age = brms_b_hat,
    rmses = list(beta_hat     = rmse_brms_beta_hat,
                 mrl_hat      = rmse_brms_mrl,          # capped at 150
                 mrl_hat_all  = rmse_brms_mrl_all,      # all values, no capping
                 bio_age      = rmse_brms_b_hat,        # capped at 150
                 bio_age_all  = rmse_brms_b_hat_all),   # all values, no capping
    method      = "brms Gompertz",
    distr_ests  = list(shape = brms_tau_hat, scale = brms_sigma_hat),
    fit_bAFT    = fit_bAFT,
    diagnostics = dgs,
    within_rep_mcse = list(
      beta     = beta_mcse,
      max_beta = max_beta_mcse,
      tau      = tau_mcse, 
      sigma    = sigma_mcse
    )    
  ))
  
}


# ============================================================
# Prior Predictive Check
# ============================================================

run_prior_predictive_check <- function(
    df_sim,
    p,
    g,
    non_zero_groups, 
    prior     = c("hs", "gs"),
    iter      = 2000,
    chains    = 4,
    warmup    = 1000,
    ndraws    = 100,
    xlim_vals = NULL,
    seed      = 123,
    save_plot = TRUE
) {
  
  prior <- match.arg(prior)
  set.seed(seed)       # ensure identical train split across both prior checks
  
  
  # ---- unpack data ----
  df_full <- df_sim$df
  
  train_idx <- sample(seq_len(nrow(df_full)), nrow(df_full) / 2, )
  df_train <- df_full[train_idx, ]
  
  x_names <- colnames(df_train)[seq_len(p)]
  
  
  # ---- formula (identical to bAFT_gomp) ----
  linpred_formula <- as.formula(
    paste0("linpred ~ 0 + ", paste(x_names, collapse = " + "))
  )
  formula_bAFT <- bf(
    age_end | cens(1 - status) + trunc(lb = entry_age) ~ 1,
    gamma ~ 1
  ) + linpred_formula
  
  # ---- priors & stanvars (identical to bAFT_gomp) ----
  priors <- switch(prior, 
                   hs = priors_hs(p, non_zero_groups), 
                   gs = priors_gs)
  
  # brms-level normal(0, 1) for b_linpred so brms has something to sample from.
  # (beta | sigma_beta ~ N(0, sigma_beta), sigma_beta ~ half-N(0,1)).
  if (prior == "gs") {
    priors <- c(priors, prior(normal(0, 1), class = "b", dpar = "linpred"))
  }
  
  stanvars_used <- switch(
    prior,
    hs = stanvars,
    gs = make_stanvars_shared(
      formula           = formula_bAFT,
      family            = gompertz_aft,
      data              = df_train,
      dpar              = "linpred",
      base_stanvars     = stanvars,
      include_intercept = FALSE
    )
  )
  
  cat("\n[Prior Predictive Check] Fitting model with sample_prior = 'only'...\n")
  
  
  # ---- fit prior-only model ----
  fit_prior <- brm(
    formula  = formula_bAFT,
    data     = df_train,
    family   = gompertz_aft,
    prior    = priors,
    iter     = iter,
    chains   = chains,
    warmup   = warmup,
    control  = list(adapt_delta = 0.95, max_treedepth = 12),
    seed     = seed,
    backend  = "cmdstanr",
    save_pars = save_pars(all = FALSE), 
    init     = 0,
    stanvars = stanvars_used,
    sample_prior = "only"
  )
  
  cat("[Prior Predictive Check] Done. Generating pp_check plot...\n")
  
  
  # ---- build pp_check plot ----
  p_ppc <- pp_check(fit_prior, ndraws = ndraws)
  
  if (!is.null(xlim_vals)) {
    p_ppc <- p_ppc + xlim(xlim_vals[1], xlim_vals[2])
  }
  
  p_ppc <- p_ppc +
    labs(
      title    = glue("Prior Predictive Check  |  prior = {prior}"),
      subtitle = glue("ndraws = {ndraws}  |  p = {p}  |  g = {g}")
    ) +
    theme_minimal()
  
  
  # ---- optionally, save ----
  if (save_plot) {
    file_tag  <- format(Sys.time(), "%Y%m%d_%H%M%S")
    diag_dir  <- if (prior == "hs") DIAG_DIR_HS else DIAG_DIR_GS
    save_path <- file.path(diag_dir, glue("ppc_prior_{prior}_{file_tag}.png"))
    ggsave(save_path, p_ppc, width = 8, height = 5, dpi = 200, bg = "white")
    cat("[Prior Predictive Check] Plot saved to:", save_path, "\n")
  }
  
  print(p_ppc)
  
  invisible(list(plot = p_ppc, fit = fit_prior))
}


# ============================================================
# Monte Carlo experiments
# ============================================================

grid_mc_exp = function(
    n_grid = c(250, 500),
    p_grid = c(20, 100),
    g_grid = 10,
    #beta_prior = c("hs", "gs"),
    non_zero_grid = c(0.25, 1.0),
    M = 200,
    n_pop = 1e5, # Life table population n
    rho_w = c(0.2, 0.7),
    target_snr = 2,
    seed_base = 123,
    experiment_name = NA) {
  
  
  print("Running MC experiment...")
  
  # time stamp to avoid overwriting runs on the same day
  datetime_start <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  
  # if experiment_name is NA, fall back to timestamp
  exp_name <- if (is.na(experiment_name)) datetime_start else experiment_name
  
  # create experiment folder (run-specific)
  exp_dir <- file.path(EXP_DIR, exp_name)
  dir.create(exp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # create parameter grid
  param_grid = expand.grid(
    n_obs = n_grid,
    p     = p_grid,
    g     = g_grid,
    rho_w = rho_w,
    non_zero = non_zero_grid,
    stringsAsFactors = FALSE # prior as character
  )
  
  # filter invalid combinations before saving (p must be divisible by g)
  param_grid <- subset(param_grid, p %% g == 0)
  
  # save the grid
  write.csv(
    param_grid,
    file = file.path(SUMMARY_DIR, glue("Grid_{exp_name}.csv")),
    row.names = FALSE
  )
  
  cat("\n",
      "============================================================\n",
      "Running ", nrow(param_grid),
      " parameter combinations with ", M, " replications each.\n",
      "Total experiments: ", nrow(param_grid) * M, "\n",
      "============================================================\n\n",
      sep = ""
  )
  
  # save to csv and remove if already exists
  summary_csv <- file.path(SUMMARY_DIR, glue::glue("MC_RMSE_summary_{exp_name}.csv"))
  if (file.exists(summary_csv)) file.remove(summary_csv)
  
  # iterate through the grid
  for (i in 1:nrow(param_grid)) {
    
    # skip scenario if already successfully saved
    rds_file <- file.path(exp_dir, glue("Experiment_{exp_name}_giter_{i}.rds"))
    
    if (file.exists(rds_file)) {
      cat(glue("Skipping grid item {i} - already exists.\n"))
      next
    }
    
    # --- actual loop ---
    # parameters to compare
    params = param_grid[i,]
    
    cat("\n",
        "<<<<<<<< Parameter Set ", i, "/", nrow(param_grid),
        " >>>>>>>>>>\n", sep = "")
    
    cat(glue("n = {params$n_obs} | p = {params$p} | g = {params$g} | non-zero = {params$non_zero}\n\n"))
    
    
    # --------------------------------------------------------------------------
    # Scenario-level seed
    # betas and pop lt are fixed per scenario, shared across all M replications
    # --------------------------------------------------------------------------
    
    seed_scenario <- seed_base + (i - 1)
    
    cat(glue("Generating scenario-level betas and life table (seed = {seed_scenario})...\n\n"))
    
    scenario_betas <- gen_betas(
      p                = params$p,
      g                = params$g,
      non_zero_groups  = params$non_zero,
      mu_l             = -4, mu_u = 4,
      rho_w            = params$rho_w,
      target_snr       = target_snr,
      gomp_a           = exp(-9),
      gomp_b           = 0.085,
      standardise_X    = TRUE,
      betas_plot       = FALSE,
      seed             = seed_scenario
    )
    
    scenario_X_obj <- gen_X(
      n             = n_pop,
      p             = params$p,
      g             = params$g,
      rho_w         = params$rho_w,
      standardise_X = TRUE,
      X_plots       = FALSE,
      seed          = seed_scenario
    )
    
    scenario_lt <- gen_pop_lt_gompertz(
      n_pop           = n_pop,
      p               = params$p,
      non_zero_groups = params$non_zero,
      g               = params$g,
      rho_w           = params$rho_w,
      gomp_a          = exp(-9),
      gomp_b          = 0.085,
      betas           = scenario_betas$beta,
      X_pop           = scenario_X_obj$X,
      force_recalc    = TRUE,
      filename        = glue("mc{i}_lt"), # one lt per scenario (not per replication)
      lt_plot         = FALSE,
      seed            = seed_scenario
    )
    
    cat("Scenario-level betas and life table ready.\n")
    print("----------")
    
    
    # ------------------------------------------------------------------
    # Run M Monte Carlo replications
    # ------------------------------------------------------------------
    
    # Each worker runs a full replication (fAFT + bAFT_hs + bAFT_gs) 
    # with 2 chains in parallel internally.
    
    mc_results <- future_lapply(1:M, function(m) {
      # clear the console otherwise Rstudio client freezes
      if (m %% 10 == 0) {
        cat(m, '\n')
        flush.console()
      }
      cat(
        glue::glue(
          "<<< Running grid item {i}/{nrow(param_grid)}, repetition {m}/{M} >>>\n"
        )
      )
      
      # Replication-level seed - only the observed dataset varies across replications
      # i = index of parameter grid (simulation scenarios)
      # M = number of replications per scenario (MC repetitions)
      # m = replication index (individual simulated dataset)
      seed <- seed_base + (i - 1) * M + m
      
      
      # Generate dataset only (betas and lt are fixed at scenario level),
      # one data set shared by fAFT, bAFT_hs, bAFT_gs.
      cat("\nGenerating dataset from Gompertz distribution...\n")
      
      df <- create_dataset_gompertz(
        n_obs         = params$n_obs * 2, # x2 for 50/50 train/test split
        p             = params$p,
        gomp_a        = exp(-9),
        gomp_b        = 0.085,
        g             = params$g,
        rho_w         = params$rho_w,
        betas         = scenario_betas$beta,
        followup      = 20,
        standardise_X = TRUE,
        X_plots       = FALSE,
        lt            = scenario_lt,
        seed          = seed
      )
      
      
      gres <- list(
        df          = df,
        true_betas  = scenario_betas,
        lt          = scenario_lt
      )
      
      # uncomment to check if lt is crushing
      # cat("lt rows:", nrow(gres$lt), "\n")
      # cat("lt$t length:", length(gres$lt$t), "\n")
      # cat("lt$mrl_pop length:", length(gres$lt$mrl_pop), "\n")
      
      results = list()
      
      
      # ---------------------------------
      # --- frequentist AFT model fit ---
      # ---------------------------------
      print("Running Frequentist AFT...")
      
      FA_start = Sys.time()
      
      set.seed(seed)
      
      results$fAFT = withCallingHandlers(
        fAFT_gomp(df_sim = gres, p = params$p),
        
        error = function(e) {
          message("AFT failed: ", conditionMessage(e))
          traceback()
        }
      )
      FA_stop = Sys.time()
      FA_time = difftime(FA_stop, FA_start, units = "secs")
      
      if (is.null(results$fAFT)) {
        print("Frequentist method failed!")
      } else {
        print("Frequentist method successful")
      }
      
    
      # -------------------------------------------------
      # --- Bayesian AFT model fit // Horseshoe prior ---
      # -------------------------------------------------
      print("Running bAFT (horseshoe)...")
      BHS_start = Sys.time()
      set.seed(seed)
      results$bAFT_hs = withCallingHandlers(
        bAFT_gomp(df_sim = gres, 
                  p = params$p, 
                  g = params$g,
                  non_zero_groups = params$non_zero, 
                  prior = "hs", 
                  seed = seed),
        error = function(e) { message("bAFT_hs failed: ", conditionMessage(e)); traceback() }
      )
      
      BHS_stop = Sys.time()
      BHS_time = difftime(BHS_stop, BHS_start, units = "secs")
      
      cat(if (is.null(results$bAFT_hs)) "bAFT_hs failed!\n" else "bAFT_hs successful\n") 
      
      
      # ------------------------------------------------
      # --- Bayesian AFT model fit // Gaussian prior ---
      # ------------------------------------------------
      print("Running bAFT (Gaussian)...")
      BGS_start = Sys.time()
      set.seed(seed)
      results$bAFT_gs = withCallingHandlers(
        bAFT_gomp(df_sim = gres, 
                  p = params$p, 
                  g = params$g,
                  non_zero_groups = params$non_zero, 
                  prior = "gs", 
                  seed = seed),
        error = function(e) { message("bAFT_gs failed: ", conditionMessage(e)); traceback() }
      )
      BGS_stop = Sys.time()
      BGS_time = difftime(BGS_stop, BGS_start, units = "secs")

      cat(if (is.null(results$bAFT_gs)) "bAFT_gs failed!\n" else "bAFT_gs successful\n")
    
      
      # extract diagnostics for both Bayesian fits
      extract_diag <- function(res) {
        if (!is.null(res) && !is.null(res$diagnostics)) return(res$diagnostics)
        list(max_rhat = NA_real_, min_ess_bulk = NA_real_,
             min_ess_tail = NA_real_, n_divergent = NA_integer_)
      }
      
      # log progress
      cat(glue("rep {m}/{M} done\n"), 
          file = file.path(exp_dir, "progress.log"), 
          append = TRUE)
      
      # save
      return(list(
        seed    = seed,
        results = results,
        timing  = list(
          frequentist_seconds = as.numeric(FA_time),
          bAFT_hs_seconds     = as.numeric(BHS_time),
          bAFT_gs_seconds     = as.numeric(BGS_time)
        ),
        bayes_diagnostics = list(
          hs = extract_diag(results$bAFT_hs),
          gs = extract_diag(results$bAFT_gs)
        )
      ))
      
    }, future.seed=TRUE)
    
    
    # ------------------------------------------------------------------
    # Collect per-replication RMSEs into vectors
    # ------------------------------------------------------------------
    fAFT_rmse_bio_age     <- rep(NA, M)      # capped at 150
    fAFT_rmse_bio_age_all <- rep(NA, M)      # all values
    fAFT_rmse_mrl         <- rep(NA, M)      # capped at 150
    fAFT_rmse_mrl_all     <- rep(NA, M)      # all values
    fAFT_rmse_beta_hat    <- rep(NA, M)
    fAFT_tau_hat          <- rep(NA, M)
    fAFT_sigma_hat        <- rep(NA, M)
    
    bAFT_hs_rmse_bio_age      <- rep(NA, M)  # capped at 150
    bAFT_hs_rmse_bio_age_all  <- rep(NA, M)  # all values
    bAFT_hs_rmse_mrl          <- rep(NA, M)  # capped at 150
    bAFT_hs_rmse_mrl_all      <- rep(NA, M)  # all values
    bAFT_hs_rmse_beta_hat     <- rep(NA, M)
    bAFT_hs_tau_hat           <- rep(NA, M)
    bAFT_hs_sigma_hat         <- rep(NA, M)
    
    bAFT_gs_rmse_bio_age      <- rep(NA, M)  # capped at 150
    bAFT_gs_rmse_bio_age_all  <- rep(NA, M)  # all values
    bAFT_gs_rmse_mrl          <- rep(NA, M)  # capped at 150
    bAFT_gs_rmse_mrl_all      <- rep(NA, M)  # all values
    bAFT_gs_rmse_beta_hat     <- rep(NA, M)
    bAFT_gs_tau_hat           <- rep(NA, M)
    bAFT_gs_sigma_hat         <- rep(NA, M)
    
    # within-rep MCSE vectors (Bayesian methods only)
    bAFT_hs_mcse_max_beta <- rep(NA, M)
    bAFT_hs_mcse_tau      <- rep(NA, M)
    bAFT_hs_mcse_sigma    <- rep(NA, M)
    bAFT_gs_mcse_max_beta <- rep(NA, M)
    bAFT_gs_mcse_tau      <- rep(NA, M)
    bAFT_gs_mcse_sigma    <- rep(NA, M)
    
    
    for (m in 1:M) {
      if (!is.null(mc_results[[m]]$results$fAFT)) {
        fAFT_rmse_bio_age[m]     <- mc_results[[m]]$results$fAFT$rmses$bio_age
        fAFT_rmse_bio_age_all[m] <- mc_results[[m]]$results$fAFT$rmses$bio_age_all
        fAFT_rmse_mrl[m]         <- mc_results[[m]]$results$fAFT$rmses$mrl
        fAFT_rmse_mrl_all[m]     <- mc_results[[m]]$results$fAFT$rmses$mrl_all
        fAFT_rmse_beta_hat[m]    <- mc_results[[m]]$results$fAFT$rmses$coef
        fAFT_tau_hat[m]          <- mc_results[[m]]$results$fAFT$distr_ests$shape
        fAFT_sigma_hat[m]        <- mc_results[[m]]$results$fAFT$distr_ests$scale
      }
      if (!is.null(mc_results[[m]]$results$bAFT_hs)) {
        bAFT_hs_rmse_bio_age[m]     <- mc_results[[m]]$results$bAFT_hs$rmses$bio_age
        bAFT_hs_rmse_bio_age_all[m] <- mc_results[[m]]$results$bAFT_hs$rmses$bio_age_all
        bAFT_hs_rmse_mrl[m]         <- mc_results[[m]]$results$bAFT_hs$rmses$mrl_hat
        bAFT_hs_rmse_mrl_all[m]     <- mc_results[[m]]$results$bAFT_hs$rmses$mrl_hat_all
        bAFT_hs_rmse_beta_hat[m]    <- mc_results[[m]]$results$bAFT_hs$rmses$beta_hat
        bAFT_hs_tau_hat[m]          <- mc_results[[m]]$results$bAFT_hs$distr_ests$shape
        bAFT_hs_sigma_hat[m]        <- mc_results[[m]]$results$bAFT_hs$distr_ests$scale
        bAFT_hs_mcse_max_beta[m]    <- mc_results[[m]]$results$bAFT_hs$within_rep_mcse$max_beta
        bAFT_hs_mcse_tau[m]         <- mc_results[[m]]$results$bAFT_hs$within_rep_mcse$tau
        bAFT_hs_mcse_sigma[m]       <- mc_results[[m]]$results$bAFT_hs$within_rep_mcse$sigma
      }
      if (!is.null(mc_results[[m]]$results$bAFT_gs)) {
        bAFT_gs_rmse_bio_age[m]      <- mc_results[[m]]$results$bAFT_gs$rmses$bio_age
        bAFT_gs_rmse_bio_age_all[m]  <- mc_results[[m]]$results$bAFT_gs$rmses$bio_age_all
        bAFT_gs_rmse_mrl[m]          <- mc_results[[m]]$results$bAFT_gs$rmses$mrl_hat
        bAFT_gs_rmse_mrl_all[m]      <- mc_results[[m]]$results$bAFT_gs$rmses$mrl_hat_all
        bAFT_gs_rmse_beta_hat[m]     <- mc_results[[m]]$results$bAFT_gs$rmses$beta_hat
        bAFT_gs_tau_hat[m]           <- mc_results[[m]]$results$bAFT_gs$distr_ests$shape
        bAFT_gs_sigma_hat[m]         <- mc_results[[m]]$results$bAFT_gs$distr_ests$scale
        bAFT_gs_mcse_max_beta[m]     <- mc_results[[m]]$results$bAFT_gs$within_rep_mcse$max_beta
        bAFT_gs_mcse_tau[m]          <- mc_results[[m]]$results$bAFT_gs$within_rep_mcse$tau
        bAFT_gs_mcse_sigma[m]        <- mc_results[[m]]$results$bAFT_gs$within_rep_mcse$sigma
      }
    }
    
    # true parameter values (fixed across all replications)
    true_tau   <- exp(-9) / 0.085
    true_sigma <- 1 / 0.085
    
    
    # ------------------------------------------------------------------
    # Compute MC mean and MC SE for each method and metric
    # ------------------------------------------------------------------
    n_fAFT_ok  <- sum(!is.na(fAFT_rmse_bio_age))
    n_hs_ok    <- sum(!is.na(bAFT_hs_rmse_bio_age))
    n_gs_ok    <- sum(!is.na(bAFT_gs_rmse_bio_age))
    
    mc_se_safe <- function(x, n) ifelse(n > 1, sd(x, na.rm = TRUE) / sqrt(n), NA)
    
    mc_summary <- data.frame(
      n_obs    = params$n_obs,
      p        = params$p,
      g        = params$g,
      rho_w    = params$rho_w,
      non_zero = params$non_zero,
      method   = rep(c("fAFT", "bAFT_hs", "bAFT_gs"), each = 5),
      metric   = rep(c("bio_age", "mrl", "beta_hat", "tau", "sigma"), times = 3),
      true_value = c(
        rep(NA, 3), true_tau, true_sigma,   # fAFT
        rep(NA, 3), true_tau, true_sigma,   # bAFT_hs
        rep(NA, 3), true_tau, true_sigma    # bAFT_gs
      ),
      # ---- capped RMSE (values capped at 150) ----
      mc_mean = c(
        mean(fAFT_rmse_bio_age,     na.rm = TRUE),
        mean(fAFT_rmse_mrl,         na.rm = TRUE),
        mean(fAFT_rmse_beta_hat,    na.rm = TRUE),
        mean(fAFT_tau_hat,          na.rm = TRUE),
        mean(fAFT_sigma_hat,        na.rm = TRUE),
        mean(bAFT_hs_rmse_bio_age,  na.rm = TRUE),
        mean(bAFT_hs_rmse_mrl,      na.rm = TRUE),
        mean(bAFT_hs_rmse_beta_hat, na.rm = TRUE),
        mean(bAFT_hs_tau_hat,       na.rm = TRUE),
        mean(bAFT_hs_sigma_hat,     na.rm = TRUE),
        mean(bAFT_gs_rmse_bio_age,  na.rm = TRUE),
        mean(bAFT_gs_rmse_mrl,      na.rm = TRUE),
        mean(bAFT_gs_rmse_beta_hat, na.rm = TRUE),
        mean(bAFT_gs_tau_hat,       na.rm = TRUE),
        mean(bAFT_gs_sigma_hat,     na.rm = TRUE)
      ),
      mc_se = c(
        mc_se_safe(fAFT_rmse_bio_age,     n_fAFT_ok),
        mc_se_safe(fAFT_rmse_mrl,         n_fAFT_ok),
        mc_se_safe(fAFT_rmse_beta_hat,    n_fAFT_ok),
        mc_se_safe(fAFT_tau_hat,          n_fAFT_ok),
        mc_se_safe(fAFT_sigma_hat,        n_fAFT_ok),
        mc_se_safe(bAFT_hs_rmse_bio_age,  n_hs_ok),
        mc_se_safe(bAFT_hs_rmse_mrl,      n_hs_ok),
        mc_se_safe(bAFT_hs_rmse_beta_hat, n_hs_ok),
        mc_se_safe(bAFT_hs_tau_hat,       n_hs_ok),
        mc_se_safe(bAFT_hs_sigma_hat,     n_hs_ok),
        mc_se_safe(bAFT_gs_rmse_bio_age,  n_gs_ok),
        mc_se_safe(bAFT_gs_rmse_mrl,      n_gs_ok),
        mc_se_safe(bAFT_gs_rmse_beta_hat, n_gs_ok),
        mc_se_safe(bAFT_gs_tau_hat,       n_gs_ok),
        mc_se_safe(bAFT_gs_sigma_hat,     n_gs_ok)
      ),
      # ---- true RMSE (all values, no capping; NA for tau/sigma rows) ----
      mc_mean_all = c(
        mean(fAFT_rmse_bio_age_all,     na.rm = TRUE),
        mean(fAFT_rmse_mrl_all,         na.rm = TRUE),
        NA,   # beta_hat 
        NA,   # tau
        NA,   # sigma
        mean(bAFT_hs_rmse_bio_age_all,  na.rm = TRUE),
        mean(bAFT_hs_rmse_mrl_all,      na.rm = TRUE),
        NA, NA, NA,
        mean(bAFT_gs_rmse_bio_age_all,  na.rm = TRUE),
        mean(bAFT_gs_rmse_mrl_all,      na.rm = TRUE),
        NA, NA, NA
      ),
      mc_se_all = c(
        mc_se_safe(fAFT_rmse_bio_age_all,     n_fAFT_ok),
        mc_se_safe(fAFT_rmse_mrl_all,         n_fAFT_ok),
        NA, NA, NA,
        mc_se_safe(bAFT_hs_rmse_bio_age_all,  n_hs_ok),
        mc_se_safe(bAFT_hs_rmse_mrl_all,      n_hs_ok),
        NA, NA, NA,
        mc_se_safe(bAFT_gs_rmse_bio_age_all,  n_gs_ok),
        mc_se_safe(bAFT_gs_rmse_mrl_all,      n_gs_ok),
        NA, NA, NA
      ),
      n_valid  = c(
        rep(n_fAFT_ok, 5),
        rep(n_hs_ok,   5),
        rep(n_gs_ok,   5)
      ),
      mean_within_mcse = c(
        rep(NA, 5), # no within-rep MCSE for fAFT
        mean(bAFT_hs_mcse_max_beta, na.rm = TRUE),
        mean(bAFT_hs_mcse_max_beta, na.rm = TRUE),
        mean(bAFT_hs_mcse_max_beta, na.rm = TRUE),
        mean(bAFT_hs_mcse_tau,      na.rm = TRUE),
        mean(bAFT_hs_mcse_sigma,    na.rm = TRUE),
        mean(bAFT_gs_mcse_max_beta, na.rm = TRUE),
        mean(bAFT_gs_mcse_max_beta, na.rm = TRUE),
        mean(bAFT_gs_mcse_max_beta, na.rm = TRUE),
        mean(bAFT_gs_mcse_tau,      na.rm = TRUE),
        mean(bAFT_gs_mcse_sigma,    na.rm = TRUE)
      )
    )
    
    cat("\n--- MC RMSE summary (scenario", i, ") ---\n")
    print(mc_summary)
    
    # save to csv 
    if (!file.exists(summary_csv)) {
      write.csv(mc_summary, file = summary_csv, row.names = FALSE)
    } else {
      write.table(mc_summary, file = summary_csv, sep = ",",
                  row.names = FALSE, col.names = FALSE, append = TRUE)
    }
    cat("Saved to:", summary_csv, "\n")
    
    
    # ------------------------------------------------------------------
    # Plot: true RMSE vs. capped RMSE (bio_age and mrl) per method
    # ------------------------------------------------------------------
    rmse_compare <- dplyr::bind_rows(
      data.frame(
        method      = "fAFT",
        metric      = "bio_age",
        rmse_capped = fAFT_rmse_bio_age,
        rmse_all    = fAFT_rmse_bio_age_all
      ),
      data.frame(
        method      = "fAFT",
        metric      = "mrl",
        rmse_capped = fAFT_rmse_mrl,
        rmse_all    = fAFT_rmse_mrl_all
      ),
      data.frame(
        method      = "bAFT_hs",
        metric      = "bio_age",
        rmse_capped = bAFT_hs_rmse_bio_age,
        rmse_all    = bAFT_hs_rmse_bio_age_all
      ),
      data.frame(
        method      = "bAFT_hs",
        metric      = "mrl",
        rmse_capped = bAFT_hs_rmse_mrl,
        rmse_all    = bAFT_hs_rmse_mrl_all
      ),
      data.frame(
        method      = "bAFT_gs",
        metric      = "bio_age",
        rmse_capped = bAFT_gs_rmse_bio_age,
        rmse_all    = bAFT_gs_rmse_bio_age_all
      ),
      data.frame(
        method      = "bAFT_gs",
        metric      = "mrl",
        rmse_capped = bAFT_gs_rmse_mrl,
        rmse_all    = bAFT_gs_rmse_mrl_all
      )
    )
    
    p_rmse_compare <- ggplot(
      rmse_compare,
      aes(x = rmse_all, y = rmse_capped, colour = method)
    ) +
      geom_abline(linetype = "dashed", colour = "grey50", linewidth = 0.5) +
      geom_point(alpha = 0.5, size = 1.5) +
      facet_wrap(~ metric, scales = "free",
                 labeller = labeller(metric = c(bio_age = "Bio Age", mrl = "MRL"))) +
      scale_colour_manual(
        values = c(fAFT = "#E69F00", bAFT_hs = "#56B4E9", bAFT_gs = "#009E73")
      ) +
      labs(
        title    = glue("True vs. Capped (≤150) RMSE  |  scenario {i}"),
        subtitle = glue("n = {params$n_obs} | p = {params$p} | rho_w = {params$rho_w} | non_zero = {params$non_zero}"),
        x        = "True RMSE  (all values)",
        y        = "Capped RMSE  (truncated at 150)",
        colour   = "Method"
      ) +
      theme_minimal(base_size = 11)
    
    print(p_rmse_compare)
    
    ggsave(
      file.path(PLOT_DIR, glue("RMSE_all_vs_capped_scenario_{i}.png")),
      p_rmse_compare, width = 9, height = 5, dpi = 180, bg = "white"
    )
    
    res_i <- list(
      params      = params,
      mc_results  = mc_results,
      scenario_lt = scenario_lt
    )
    
    base::saveRDS(
      res_i,
      file = file.path(exp_dir, glue::glue("Experiment_{exp_name}_giter_{i}.rds"))
    )
    
    # free memory before next scenario
    rm(mc_results, res_i)
    gc()
  }
  
  
  return("Finished experiments!")
}



