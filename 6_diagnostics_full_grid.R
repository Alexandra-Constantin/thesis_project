source("5_diagnostics.R")

EXP_NAME <- "full_grid_run_01"

# ==============================================================================
# (1) Quick divergence scan across all giters and replications
# ==============================================================================

list_experiments()
list_giters(EXP_NAME)

# scan all 16 scenarios
rows <- vector("list", 16 * 100)
idx  <- 1L

for (g in 1:16) {
  res <- readRDS(file.path(EXP_DIR, EXP_NAME,
                           paste0("Experiment_", EXP_NAME, "_giter_", g, ".rds")))
  for (m in 1:100) {
    hs_div  <- res$mc_results[[m]]$bayes_diagnostics$hs$n_divergent
    gs_div  <- res$mc_results[[m]]$bayes_diagnostics$gs$n_divergent
    hs_rhat <- res$mc_results[[m]]$bayes_diagnostics$hs$max_rhat
    gs_rhat <- res$mc_results[[m]]$bayes_diagnostics$gs$max_rhat
    
    rows[[idx]] <- data.frame(giter   = g,
                              m       = m,
                              hs_div  = hs_div,
                              gs_div  = gs_div,
                              hs_rhat = hs_rhat,
                              gs_rhat = gs_rhat)
    idx <- idx + 1L
  }
}

diag_df <- do.call(rbind, rows)
diag_df_flagged <- diag_df[diag_df$hs_div > 0 | diag_df$gs_div > 0 |
                             (!is.na(diag_df$hs_rhat) & diag_df$hs_rhat > 1.01) |
                             (!is.na(diag_df$gs_rhat) & diag_df$gs_rhat > 1.01), ]

write.csv(diag_df,         file = "mcmc_diagnostics_all.csv",     row.names = FALSE)
write.csv(diag_df_flagged, file = "mcmc_diagnostics_flagged.csv", row.names = FALSE)


diag_df <- read.csv("mcmc_diagnostics_all.csv")

range(diag_df$hs_div)
library(dplyr)
diag_df %>%
  filter(hs_rhat > 1.01) %>%
  arrange(desc(hs_div)) %>%
  select(giter, m, hs_div, hs_rhat)


# ==============================================================================
# (2) Full diagnostics on a specific replication
# Change giter and m to whichever case the scan flagged
# ==============================================================================

# generic case
fits <- load_debug_fits(EXP_NAME, giter = 1, m = 1)

run_mcmc_diagnostics(fits$hs, prior = "hs", save_plots = TRUE)
run_mcmc_diagnostics(fits$gs, prior = "gs", save_plots = TRUE)

run_posterior_predictive_check(fits$hs, prior = "hs", save_plot = TRUE)
run_posterior_predictive_check(fits$gs, prior = "gs", save_plot = TRUE)

compare_ppc(fits$hs, fits$gs, save_plot = TRUE)


# the worst case: hs divs=584
fits2 <- load_debug_fits(EXP_NAME, giter = 2, m = 38)

run_mcmc_diagnostics(fits2$hs, prior = "hs", save_plots = TRUE)
run_posterior_predictive_check(fits2$hs, prior = "hs", save_plot = TRUE)

# second worst case: hs divs=266
fits3 <- load_debug_fits(EXP_NAME, giter = 1, m = 97)

run_mcmc_diagnostics(fits3$hs, prior = "hs", save_plots = TRUE)
run_posterior_predictive_check(fits3$hs, prior = "hs", save_plot = TRUE)

# third worst case: hs divs=106
fits4 <- load_debug_fits(EXP_NAME, giter = 9, m = 6)

run_mcmc_diagnostics(fits4$hs, prior = "hs", save_plots = TRUE)
run_posterior_predictive_check(fits4$hs, prior = "hs", save_plot = TRUE)



# ---------------------------
# COUNT BULK AND TAIL ESS
# ---------------------------
x <- readRDS("output/MCexp/experiments/full_grid_run_01/Experiment_full_grid_run_01_giter_1.rds")
x$mc_results[[1]]$results$bAFT_hs$diagnostics
extract_convergence_summary <- function(exp_dir, exp_name, n_scenarios = 16) {
  
  all_rows <- list()
  
  for (i in seq_len(n_scenarios)) {
    
    rds_path <- file.path(exp_dir, exp_name,
                          glue::glue("Experiment_{exp_name}_giter_{i}.rds"))
    if (!file.exists(rds_path)) { message("Missing: ", rds_path); next }
    
    message("Processing scenario ", i)
    x      <- readRDS(rds_path)
    params <- x$params
    
    for (m in seq_along(x$mc_results)) {
      rep_res <- x$mc_results[[m]]$results
      
      for (method in c("bAFT_hs", "bAFT_gs")) {
        dg <- rep_res[[method]]$diagnostics
        if (is.null(dg)) next
        
        all_rows <- append(all_rows, list(tibble(
          scenario      = i,
          rep           = m,
          method        = method,
          n_obs         = params$n_obs,
          p             = params$p,
          rho_w         = params$rho_w,
          non_zero      = params$non_zero,
          max_rhat      = dg$max_rhat,
          min_bulk_ess  = dg$min_ess_bulk,
          min_tail_ess  = dg$min_ess_tail,
          n_divergent   = dg$n_divergent
        )))
      }
    }
    
    rm(x); gc()
  }
  
  bind_rows(all_rows)
}

conv_df <- extract_convergence_summary(
  exp_dir     = "output/MCexp/experiments",
  exp_name    = "full_grid_run_01",
  n_scenarios = 16
)

saveRDS(conv_df, "output/MCexp/experiments/full_grid_run_01/convergence_summary.rds")

conv_df <- readRDS("output/MCexp/experiments/full_grid_run_01/convergence_summary.rds")

# count number of reps per scenario/method
conv_df %>% count(method, n_obs, p, rho_w, non_zero) %>%
  print(n = Inf)

# count number of warnings
conv_df %>%
  group_by(method, n_obs, p, rho_w, non_zero) %>%
  summarise(
    n_reps        = n(),
    pct_rhat_warn = mean(max_rhat > 1.01) * 100,
    pct_bulk_warn = mean(min_bulk_ess < 400) * 100,  # 100 × 4 chains
    pct_tail_warn = mean(min_tail_ess < 400) * 100,
    pct_diverg    = mean(n_divergent > 0) * 100,
    med_bulk      = median(min_bulk_ess),
    med_tail      = median(min_tail_ess),
    .groups = "drop"
  ) %>%
  print(n = Inf)


conv_df %>%
  filter(method == "bAFT_hs") %>%
  group_by(n_obs, p, rho_w, non_zero) %>%
  summarise(
    med_diverg  = median(n_divergent),
    p90_diverg  = quantile(n_divergent, 0.9),
    max_diverg  = max(n_divergent),
    .groups = "drop"
  ) %>%
  print(n = Inf)


conv_df %>%
  filter(method == "bAFT_hs", n_divergent > 100) %>%
  dplyr::select(scenario, rep, n_obs, p, rho_w, non_zero,
                n_divergent, max_rhat, min_bulk_ess, min_tail_ess)
