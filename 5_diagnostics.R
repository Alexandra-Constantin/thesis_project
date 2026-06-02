# ------------------------------------------------------------------------------
# Set-up
# ------------------------------------------------------------------------------
DIAG_DIR_HS <- if (exists("DIAG_DIR_HS", inherits = TRUE)) DIAG_DIR_HS else "output/scenarios/bayesian/diagnostics/hs"
DIAG_DIR_GS <- if (exists("DIAG_DIR_GS", inherits = TRUE)) DIAG_DIR_GS else "output/scenarios/bayesian/diagnostics/gs"
EXP_DIR     <- if (exists("EXP_DIR",     inherits = TRUE)) EXP_DIR     else "output/scenarios/bayesian/experiments"
SUMMARY_DIR <- if (exists("SUMMARY_DIR", inherits = TRUE)) SUMMARY_DIR else "output/scenarios/bayesian/summaries"

for (d in c(DIAG_DIR_HS, DIAG_DIR_GS, EXP_DIR, SUMMARY_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ==============================================================================
# MCMC Diagnostics for Gompertz AFT brms model
# Based on "Diagnosing Biased Inference with Divergences" (Betancourt, 2017) 
# ==============================================================================

library(brms)
library(posterior)
library(ggplot2)
library(dplyr)
library(tidyr)
library(bayesplot)

# ------------------------------------------------------------------------------
# Colours
# ------------------------------------------------------------------------------
c_light           <- "#DCBCBC"
c_light_highlight <- "#C79999"
c_mid             <- "#B97C7C"
c_mid_highlight   <- "#A25050"
c_dark            <- "#8F2727"
c_dark_highlight  <- "#7C0000"


# ==============================================================================
# Master Function    (used for any brms fit)
# ==============================================================================

run_mcmc_diagnostics <- function(
    fit,                                    # brmsfit object (from bAFT_gomp result$fit_bAFT)
    funnel_x_param  = "b_gamma_Intercept",  # x-axis param for funnel plot; NULL = first beta
    n_beta_trace    = 5,                    # how many beta chains to show in trace
    prior    = c("hs", "gs"),
    save_plots      = FALSE                 # set TRUE to write PNGs to plot_dir
) {
  
  prior    <- match.arg(prior)
  plot_dir <- if (prior == "hs") DIAG_DIR_HS else DIAG_DIR_GS
  
  if (save_plots) dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --------------------------------------------------------------------------
  # 1. PRINT SUMMARY DIAGNOSTICS
  # --------------------------------------------------------------------------
  cat("\n==============================================================\n")
  cat("  MCMC DIAGNOSTIC SUMMARY\n")
  cat("==============================================================\n\n")
  
  # extract posterior draws in an array format
  draws_arr <- as_draws_array(fit)
  
  # get per-variable diagnostics
  diag_df <- summarise_draws(draws_arr, "rhat", "ess_bulk", "ess_tail")
  
  # create named vectors
  rhat_vals     <- setNames(diag_df$rhat,     diag_df$variable)
  ess_bulk_vals <- setNames(diag_df$ess_bulk, diag_df$variable)
  ess_tail_vals <- setNames(diag_df$ess_tail, diag_df$variable)
  
  cat("Max Rhat          :", round(max(rhat_vals, na.rm = TRUE), 4), "\n")
  cat("Min ESS bulk      :", round(min(ess_bulk_vals, na.rm = TRUE), 0), "\n")
  cat("Min ESS tail      :", round(min(ess_tail_vals, na.rm = TRUE), 0), "\n")
  
  # --- NUTS diagnostics ---
  np <- nuts_params(fit)
  # count divergent transitions
  n_div  <- sum(np$Parameter == "divergent__" & np$Value == 1, na.rm = TRUE)
  # total post-warmup draws checked for divergences
  n_iter <- sum(np$Parameter == "divergent__", na.rm = TRUE)
  cat("Divergent transitions:", n_div, "/", n_iter,
      sprintf("(%.2f%%)\n", 100 * n_div / n_iter))
  
  # max treedepth
  max_td <- max(np$Value[np$Parameter == "treedepth__"], na.rm = TRUE)
  cat("Max tree depth recorded:", max_td, "\n\n")
  
  # worst Rhat
  bad_rhat <- sort(rhat_vals[rhat_vals > 1.01], decreasing = TRUE)
  if (length(bad_rhat) > 0) {
    cat("Parameters with Rhat > 1.01:\n")
    print(head(bad_rhat, 10))
  } else {
    cat("All Rhat values < 1.01.\n")
  }
  
  # worst ESS
  bad_ess <- sort(ess_bulk_vals[ess_bulk_vals < 400], decreasing = FALSE)
  if (length(bad_ess) > 0) {
    cat("\nParameters with ESS_bulk < 400:\n")
    print(head(bad_ess, 10))
  } else {
    cat("All ESS_bulk values > 400.\n")
  }
  
  
  
  # --------------------------------------------------------------------------
  # 2. EXTRACT POSTERIOR DRAWS 
  # --------------------------------------------------------------------------
  post <- as.data.frame(as_draws_df(fit))
  
  # identify beta columns
  beta_cols <- grep("^b_linpred_", colnames(post), value = TRUE)
  
  # divergence indicator per draw
  div_vec <- np$Value[np$Parameter == "divergent__"]
  
  # number of post-warmup samples in each chain
  n_draws_per_chain <- nrow(as_draws_array(fit)[,,1])
  # total number of chains
  n_chains <- dim(as_draws_array(fit))[2]
  
  # iteration number within chain
  post$iter <- rep(seq_len(n_draws_per_chain), times = n_chains)
  # chain id
  post$chain <- rep(seq_len(n_chains), each = n_draws_per_chain)
  # divergence flag
  post$divergent <- as.integer(div_vec)
  
  
  # --------------------------------------------------------------------------
  # 3. TRACE PLOTS (log of distributional params + one beta)
  # --------------------------------------------------------------------------
  cat("Generating trace plots...\n")
  
  # Gompertz distributional parameters (on log scale) to trace
  dist_params <- list(
    list(col = "b_Intercept",       label = "log(tau) [shape]",  log = FALSE),
    list(col = "b_gamma_Intercept", label = "log(sigma) [scale]", log = FALSE)
  )
  
  # loop over these parameters
  for (dp in dist_params) {
    # if it is absent from the posterior data frame, skip it
    if (!dp$col %in% colnames(post)) next
    # otherwise extract its sampled values (optionally, log transform it)
    vals <- if (dp$log) log(post[[dp$col]]) else post[[dp$col]]
    # create df for plotting
    df_tr <- data.frame(iter = post$iter, chain = factor(post$chain), val = vals)
    
    # create trace plot
    p_tr <- ggplot(df_tr, aes(x = iter, y = val, colour = chain)) +
      geom_point(size = 0.4, alpha = 0.5) +
      scale_colour_manual(values = c(c_dark, c_mid, c_light, c_dark_highlight)) +
      labs(title = paste("Trace plot -", dp$label),
           x = "Iteration", y = dp$label, colour = "Chain") +
      theme_minimal(base_size = 12)
    
    print(p_tr)
    
    if (save_plots) {
      ggsave(file.path(plot_dir, paste0("trace_", dp$col, ".png")),
             p_tr, width = 8, height = 4, dpi = 180, bg = "white")
    }
  }
  
  # trace a selection of beta coefficients
  if (length(beta_cols) > 0) {
    # step = (end-start)/nr of intervals
    idx_show <- round(seq(1, length(beta_cols), length.out = min(n_beta_trace, length(beta_cols))))
    betas_show <- beta_cols[idx_show]
    
    for (bc in betas_show) {
      df_tr <- data.frame(iter = post$iter, chain = factor(post$chain), val = post[[bc]])
      p_tr <- ggplot(df_tr, aes(x = iter, y = val, colour = chain)) +
        geom_point(size = 0.4, alpha = 0.5) +
        scale_colour_manual(values = c(c_dark, c_mid, c_light, c_dark_highlight)) +
        labs(title = paste("Trace plot -", bc),
             x = "Iteration", y = bc, colour = "Chain") +
        theme_minimal(base_size = 12)
      print(p_tr)
      
      if (save_plots) {
        ggsave(file.path(plot_dir, paste0("trace_", bc, ".png")),
               p_tr, width = 8, height = 4, dpi = 180, bg = "white")
      }
    }
  }
  
  
  # --------------------------------------------------------------------------
  # 4. RUNNING MEAN PLOT (bias check)
  # --------------------------------------------------------------------------
  cat("Generating running mean plots...\n")
  
  # distribution parameters
  params_for_bias <- c("b_Intercept", "b_gamma_Intercept")
  # select and append one regression coefficient
  if (length(beta_cols) > 0) params_for_bias <- c(params_for_bias, beta_cols[1])
  
  # compute (global) running mean
  for (pc in params_for_bias) {
    if (!pc %in% colnames(post)) next
    
    # get all posterior draws for each parameter
    vals <- post[[pc]]
    # cumulative mean after each draw
    run_mean <- cumsum(vals) / seq_along(vals)
    # overall posterior mean across all draws
    true_val <- mean(vals)    # use posterior mean as reference
    
    df_rm <- data.frame(iter = seq_along(vals), running_mean = run_mean)
    
    # plot
    p_rm <- ggplot(df_rm, aes(x = iter, y = running_mean)) +
      geom_point(size = 0.5, colour = c_dark, alpha = 0.7) +
      geom_hline(yintercept = true_val, colour = "grey50", linetype = "dashed", linewidth = 1.2) +
      labs(
        title    = paste("Running mean -", pc),
        subtitle = "Dashed line = posterior mean",
        x = "Draw", y = paste("Running mean of", pc)
      ) +
      theme_minimal(base_size = 12)
    
    print(p_rm)
    
    if (save_plots) {
      ggsave(file.path(plot_dir, paste0("runmean_", pc, ".png")),
             p_rm, width = 8, height = 4, dpi = 180, bg = "white")
    }
  }
  
  
  # --------------------------------------------------------------------------
  # 5. FUNNEL / SCATTER PLOT  (divergences highlighted in green)
  # --------------------------------------------------------------------------
  cat("Generating funnel / divergence scatter plot...\n")
  
  # y-axis: log(tau)
  if (!"b_Intercept" %in% colnames(post)) {
    cat("Skipping funnel plot: b_Intercept not found.\n")
  } else {
    y_vals <- post[["b_Intercept"]]   # already on log scale
    
    # x-axis: first beta
    x_col <- if (!is.null(funnel_x_param) && funnel_x_param %in% colnames(post)) {
      funnel_x_param
    } else if (length(beta_cols) > 0) {
      beta_cols[1]
    } else {
      NULL
    }
    
    if (!is.null(x_col)) {
      x_vals <- post[[x_col]]
      
      # build plotting data
      df_fun <- data.frame(
        x_val     = x_vals,
        y_val     = y_vals,
        divergent = factor(post$divergent, levels = c(0, 1),
                           labels = c("Non-divergent", "Divergent"))
      )
      
      # plot divergent vs non-divergent draws
      p_fun <- ggplot() +
        geom_point(data = df_fun[df_fun$divergent == "Non-divergent", ],
                   aes(x = x_val, y = y_val),
                   colour = c_dark, size = 0.7, alpha = 0.4) +
        geom_point(data = df_fun[df_fun$divergent == "Divergent", ],
                   aes(x = x_val, y = y_val),
                   colour = "green3", size = 1.5, alpha = 0.9) +
        labs(
          title    = "Divergence scatter plot",
          subtitle = "Green = divergent transitions.",
          x = if (x_col == "b_gamma_Intercept") "log(sigma)" else x_col,
          y        = "log(tau)  [b_Intercept]"
        ) +
        theme_minimal(base_size = 12)
      
      print(p_fun)
      
      if (save_plots) {
        ggsave(file.path(plot_dir, "funnel_divergences.png"),
               p_fun, width = 7, height = 6, dpi = 180, bg = "white")
      }
    }
  }
  
  
  # --------------------------------------------------------------------------
  # 6. HORSESHOE PATHOLOGY
  #    plot log(tau) vs log(lambda) to see if divergences present
  # --------------------------------------------------------------------------
  
  # find the column names that hold the global and local hs parameters.
  all_col_names <- colnames(post)
  hs_tau_cols   <- all_col_names[grepl("hs_global", all_col_names)]  # global tau columns
  hs_lam_cols   <- all_col_names[grepl("hs_local",  all_col_names)]  # local lambda columns
  
  # check if both types of column actually exist
  hs_cols_found <- length(hs_tau_cols) > 0 && length(hs_lam_cols) > 0
  
  if (hs_cols_found) {
    cat("Generating horseshoe hyperparameter funnel...\n")
    
    # extract the first global-tau column and the first local-lambda column (raw)
    global_tau_raw <- post[[ hs_tau_cols[1] ]]
    local_lam_raw  <- post[[ hs_lam_cols[1] ]]
    
    # log-transform so small values near zero are visible
    log_global_tau <- log(abs(global_tau_raw))
    log_local_lam  <- log(abs(local_lam_raw))
    
    # build data frame
    df_hs <- data.frame(
      log_tau   = log_global_tau,
      log_lam   = log_local_lam,
      divergent = post$divergent
    )
    
    # split into divergent and non-divergent rows
    df_hs_ok  <- df_hs[df_hs$divergent == 0, ]
    df_hs_div <- df_hs[df_hs$divergent == 1, ]
    
    # plot
    p_hs <- ggplot() +
      geom_point(data = df_hs_ok,  aes(x = log_lam, y = log_tau),
                 colour = c_dark,  size = 0.7, alpha = 0.4) +
      geom_point(data = df_hs_div, aes(x = log_lam, y = log_tau),
                 colour = "green3", size = 1.5) +
      labs(title    = "Horseshoe hyperparameter funnel",
           subtitle = "Green = divergent transitions.",
           x = "log(local shrinkage  lambda_1)",
           y = "log(global scale  tau)") +
      theme_minimal()
    
    print(p_hs)
    
    if (save_plots) {
      ggsave(file.path(plot_dir, "hs_funnel.png"),
             p_hs, width = 7, height = 6, dpi = 180, bg = "white")
    }
  }
  
  
  # --------------------------------------------------------------------------
  # 7. ENERGY / BAYESPLOT DIAGNOSTICS
  # --------------------------------------------------------------------------
  cat("Generating energy and pair plots...\n")
  
  # NUTS energy plot
  p_energy <- mcmc_nuts_energy(np) +
    labs(title = "NUTS Energy diagnostic",
         subtitle = "E = Marginal energy ; E\u0305 = Energy Transition") +
    theme_minimal()
  print(p_energy)
  
  if (save_plots) {
    ggsave(file.path(plot_dir, "nuts_energy.png"),
           p_energy, width = 7, height = 4, dpi = 180, bg = "white")
  }
  
  # pairs plot for distributional parameters + first beta
  pair_params <- intersect(c("b_Intercept", "b_gamma_Intercept", beta_cols[1]),
                           colnames(post))
  if (length(pair_params) >= 2) {
    p_pairs <- mcmc_pairs(
      draws_arr,
      pars       = pair_params,
      np         = np,
      np_style   = pairs_style_np(div_color = "green3", div_size = 2)
    )
    print(p_pairs)
    
    if (save_plots) {
      ggsave(file.path(plot_dir, "pairs_plot.png"),
             p_pairs, width = 8, height = 8, dpi = 180, bg = "white")
    }
  }
  
  
  # --------------------------------------------------------------------------
  # 8. ESS AND RHAT BAR PLOTS
  # --------------------------------------------------------------------------
  cat("Generating Rhat and ESS plots...\n")
  
  p_rhat <- mcmc_rhat(rhat_vals) +
    yaxis_text(hjust=1) +  # to see the names of the problematic parameters only
    labs(title = "Rhat for all parameters") +
    theme_minimal()
  print(p_rhat)
  
  p_ess <- mcmc_neff(ess_bulk_vals / length(div_vec)) +
    labs(title = "Effective sample size ratio (ESS_bulk / N)") +
    theme_minimal()
  print(p_ess)
  
  if (save_plots) {
    ggsave(file.path(plot_dir, "rhat.png"), p_rhat, width = 8, height = 5, dpi = 180, bg = "white")
    ggsave(file.path(plot_dir, "ess.png"),  p_ess,  width = 8, height = 5, dpi = 180, bg = "white")
  }
  
  
  invisible(list(n_divergent = n_div, max_rhat = max(rhat_vals, na.rm = TRUE)))
}



# ==============================================================================
# Experiment loaders
# ==============================================================================

# main folder
BASE_OUT_DIR <- if (exists("BASE_OUT_DIR", inherits = TRUE)) BASE_OUT_DIR else "output/MCexp"

# ------------------------------------------------------------------------------
# List available experiments
# ------------------------------------------------------------------------------
list_experiments <- function(exp_dir = EXP_DIR) {
  # get folders
  dirs <- list.dirs(exp_dir, recursive = FALSE, full.names = FALSE)
  if (length(dirs) == 0) {
    cat("No experiments found in", exp_dir, "\n")
  } else {
    cat("Available experiments:\n", paste(" -", dirs, collapse = "\n"), "\n")
  }
  invisible(dirs)
}


# ------------------------------------------------------------------------------
# List grid iterations (giter) files
# If no grid CSV exists it just lists the file names instead.
# ------------------------------------------------------------------------------
list_giters <- function(exp_name, 
                        exp_dir = EXP_DIR, 
                        summary_dir = SUMMARY_DIR) {
  
  # find all RDS files for this experiment
  pattern   <- paste0("Experiment_", exp_name, "_giter_*.rds")
  rds_files <- Sys.glob(file.path(exp_dir, exp_name, pattern))
  
  if (length(rds_files) == 0) {
    cat("No RDS files found for experiment:", exp_name, "\n")
    return(invisible(NULL))
  }
  
  # sort by giter index
  giter_idx <- as.integer(sub(".*_giter_([0-9]+)\\.rds$", "\\1", rds_files))
  rds_files <- rds_files[order(giter_idx)]
  giter_idx <- sort(giter_idx)
  
  # load the matching grid csv
  grid_file <- file.path(summary_dir, paste0("Grid_", exp_name, ".csv"))
  
  # attach grid 
  if (file.exists(grid_file)) {
    
    # load the grid and add the giter numbers
    grid       <- read.csv(grid_file)
    grid$giter <- giter_idx
    
    # select columns of interest
    table_to_show <- grid[ , c("giter", "n_obs", "p", "non_zero", "rho_w")]
    print(table_to_show)
    
  } else {
    # print file names
    cat("RDS files found:\n")
    cat(paste(" -", basename(rds_files), collapse = "\n"), "\n")
  }
  
  invisible(rds_files)
}


# ------------------------------------------------------------------------------
# Load a single brmsfit from a saved RDS
# ------------------------------------------------------------------------------
load_fit <- function(exp_name, giter, m = 1, prior = c("hs", "gs"), exp_dir = EXP_DIR) {
  
  prior <- match.arg(prior)
  
  path <- file.path(exp_dir, exp_name,
                    paste0("Experiment_", exp_name, "_giter_", giter, ".rds"))
  
  if (!file.exists(path)) stop("File not found:\n  ", path)
  
  res <- readRDS(path)
  
  fit_key <- paste0("bAFT_", prior)
  fit <- res$mc_results[[m]]$results[[fit_key]]$fit_bAFT
  
  if (is.null(fit)) stop(
    "fit_bAFT is NULL at giter=", giter, " m=", m, " prior=", prior,
    " - model may have failed for this replication."
  )
  
  cat("Loaded:", exp_name, "| giter =", giter, "| m =", m, "| prior =", prior, "\n")
  cat("p :", res$params$p,
      "| n :", res$params$n_obs,
      "| non_zero :", res$params$non_zero, "\n")
  
  fit
}


# ------------------------------------------------------------------------------
# Load both hs and gs fits from a single giter
# ------------------------------------------------------------------------------
load_debug_fits <- function(exp_name = "debug_run", giter = 1, m = 1, exp_dir = EXP_DIR) {
  
  list_giters(exp_name, exp_dir)
  
  fit_hs <- load_fit(exp_name, giter = giter, m = m, prior = "hs", exp_dir = exp_dir)
  fit_gs <- load_fit(exp_name, giter = giter, m = m, prior = "gs", exp_dir = exp_dir)
  
  list(hs = fit_hs, gs = fit_gs)
}

