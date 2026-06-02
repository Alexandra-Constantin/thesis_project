# ==============================================================================
# Packages
# ==============================================================================
library(eha)
library(brms)
library(cmdstanr) # faster Stan
library(survival)
library(glue)
library(dplyr)
library(ggplot2)
library(bayesplot)


# ==============================================================================
# Set-up
# ==============================================================================

# clean environment
rm(list=ls())

# avoid plot windows during MC
op <- options(device = function(...) pdf(NULL))
on.exit(options(op), add = TRUE)

# uncomment after running the MC experiments to plot
# dev.off()
# options(device = "RStudioGD")

# use multiple cores (sequentially)
#options(mc.cores = parallel::detectCores())

# use multiple cores (simultaneously)
library(future.apply)
n_workers <- parallel::detectCores() %/% 2
plan(multisession, workers = n_workers)
options(mc.cores = 2)   # 2 chains per brm() call
options(future.globals.maxSize = Inf) 

# ----

# set-up output directory
BASE_OUT_DIR <- "output/MCexp"
SIM_DIR      <- file.path(BASE_OUT_DIR, "simulated")
PLOT_DIR     <- file.path(BASE_OUT_DIR, "plots")
EXP_DIR      <- file.path(BASE_OUT_DIR, "experiments")
SUMMARY_DIR  <- file.path(BASE_OUT_DIR, "summaries")
DIAG_DIR_HS  <- file.path(BASE_OUT_DIR, "diagnostics", "hs")
DIAG_DIR_GS  <- file.path(BASE_OUT_DIR, "diagnostics", "gs")

for (d in c(SIM_DIR, PLOT_DIR, EXP_DIR, SUMMARY_DIR, DIAG_DIR_HS, DIAG_DIR_GS)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# source helper functions
source("3_fit_models_parallel.R")


# ==============================================================================
# Prior Predictive Check
# run before the MC grid to validate priors
# ==============================================================================

# --- generate a small representative dataset for the check ---
ppc_data <- gen_data(
  n_obs          = 250,      # same n as debug grid
  n_pop          = 1e4,      # smaller pop just for a quick check
  p              = 20,
  g              = 10,
  rho_w          = 0.2,
  non_zero_groups = 0.5,
  target_snr     = 2,
  gomp_a         = exp(-9),
  gomp_b         = 0.085,
  ltname         = "ppc_lt",
  standardise_X  = TRUE,
  force_recalc   = TRUE,
  seed           = 123
)


# check amount of censoring vs event times
table(ppc_data$df$status)
prop.table(table(ppc_data$df$status))


# --- run the check for the horseshoe prior ---
ppc_hs <- run_prior_predictive_check(
  df_sim    = ppc_data,
  p         = 20,
  g         = 10,
  non_zero_groups = 0.5,   # match ppc_data generation above
  prior     = "hs",
  ndraws    = 100,
  xlim_vals = c(0, 150),
  seed      = 123
)


# --- run the check for the Gaussian prior ---
ppc_gs <- run_prior_predictive_check(
  df_sim    = ppc_data,
  p         = 20,
  g         = 10,
  non_zero_groups = 0.5,   # match ppc_data generation above
  prior     = "gs",
  ndraws    = 100,
  xlim_vals = c(0, 150),
  seed      = 123
)


# re-enable suppressed plot window so you can see both plots side by side
dev.off()
dev.new()
options(device = "RStudioGD")
print(ppc_hs$plot)
print(ppc_gs$plot)


# ==============================================================================
# Debug run (sanity check before the full grid)
# ==============================================================================

# get time for the entire procedure
t0 <- proc.time()

# debugging grid (not relevant for the experiment)
grid_mc_exp(
  n_grid = 250,         
  p_grid = 20,   
  g_grid = 10,
  rho_w = 0.2,           
  non_zero_grid = 0.5,  
  #beta_prior = c("hs", "gs"),
  M = 10,                 
  n_pop = 1e5,
  target_snr = 2,
  seed_base = 123,
  experiment_name = "debug_run"
)

# get time of the whole procedure
elapsed_mins <- (proc.time() - t0)[["elapsed"]] / 60
cat("\nTotal experiment time:", round(elapsed_mins, 2), "minutes\n")


# ---- check relationship predicted mrl and predicted bio age ----
res_debug <- readRDS("output/MCexp/experiments/debug_run/Experiment_debug_run_giter_1.rds")

# frequentist
plot(res_debug$mc_results[[1]]$results$fAFT$predicted_bio_age,
     res_debug$mc_results[[1]]$results$fAFT$predicted_mrl,
     xlab = "Predicted Bio Age", ylab = "Predicted MRL",
     main = "fAFT: pred MRL vs pred Bio Age")

# Bayesian Horseshoe
plot(res_debug$mc_results[[1]]$results$bAFT_hs$predicted_bio_age,
     res_debug$mc_results[[1]]$results$bAFT_hs$predicted_mrl,
     xlab = "Predicted Bio Age", ylab = "Predicted MRL",
     main = "bAFT(hs): pred MRL vs pred Bio Age")

# Bayesian Gaussian
plot(res_debug$mc_results[[1]]$results$bAFT_gs$predicted_bio_age,
     res_debug$mc_results[[1]]$results$bAFT_gs$predicted_mrl,
     xlab = "Predicted Bio Age", ylab = "Predicted MRL",
     main = "bAFT(gs): pred MRL vs pred Bio Age")

# check life table ranges
cat("range of lt$t:", range(res_debug$scenario_lt$t), "\n")
cat("MRL range in life table:", range(res_debug$scenario_lt$mrl_pop), "\n")


# ---- inspect capped vs all RMSE ----
summary_debug <- read.csv(file.path(SUMMARY_DIR, "MC_RMSE_summary_debug_run.csv"))

# show only bio_age and mrl rows where both columns are present
rmse_check <- subset(summary_debug, metric %in% c("bio_age", "mrl"),
                     select = c(method, metric, mc_mean, mc_mean_all, mc_se, mc_se_all))
print(rmse_check)



# ==============================================================================
# Full grid  (n, p, sparsity, correlation) x M=100 replications
# ==============================================================================
# restart R session before running full grid just in case

t0 <- proc.time()

grid_mc_exp(
  n_grid          = c(250, 500),
  p_grid          = c(20, 100),
  g_grid          = 10,
  rho_w           = c(0.2, 0.7),
  non_zero_grid   = c(0.25, 1.0),
  M               = 100,
  n_pop           = 1e5,
  target_snr      = 2,
  seed_base       = 123,
  experiment_name = "full_grid_run_01"
)

elapsed_mins <- (proc.time() - t0)[["elapsed"]] / 60
cat("\nFull grid time:", round(elapsed_mins, 2), "minutes\n")


# ==============================================================================
# Ultra-high-dimensional experiment 
# n = 500, p = 2000, g = 100 ; M = 1
# ==============================================================================

t0 <- proc.time()

grid_mc_exp(
  n_grid          = 500,
  p_grid          = 2000,
  g_grid          = 100,
  rho_w           = c(0.2, 0.7),
  non_zero_grid   = c(0.25, 1.0),
  M               = 1,
  n_pop           = 1e5,
  target_snr      = 2,
  seed_base       = 14,
  experiment_name = "stress_test_01"
)

elapsed_mins <- (proc.time() - t0)[["elapsed"]] / 60
cat("\nStress test time:", round(elapsed_mins, 2), "minutes\n")


# ==============================================================================
# Load and inspect results
# ==============================================================================

# read the MC RMSE summary CSVs
summary_debug    <- read.csv(file.path(SUMMARY_DIR, "MC_RMSE_summary_debug_run.csv"))
summary_full     <- read.csv(file.path(SUMMARY_DIR, "MC_RMSE_summary_full_grid_run_01.csv"))
summary_stress   <- read.csv(file.path(SUMMARY_DIR, "MC_RMSE_summary_stress_test_01.csv"))

print(summary_full)



