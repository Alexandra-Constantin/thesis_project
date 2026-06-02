# ==============================================================================
# Packages
# ============================================================================== 
source("1_dgp_helpers.R")
library(tidyverse)
library(ggh4x)
library(ggdist)
library(glue)

# ==============================================================================
# Set-up
# ==============================================================================
FIGURES_DIR <- "output/MCexp/figures"
exp_dir     <- "output/MCexp/experiments/full_grid_run_01"
exp_name    <- "full_grid_run_01"
csv_path    <- "output/MCexp/summaries/MC_RMSE_summary_full_grid_run_01.csv"
param_grid  <- read.csv("output/MCexp/summaries/Grid_full_grid_run_01.csv")
methods     <- c("fAFT", "bAFT_hs", "bAFT_gs")

# colours for the plots
group_colors <- c(
  "1"  = "#E41A1C", "2"  = "#377EB8", "3"  = "#4DAF4A",
  "4"  = "#984EA3", "5"  = "#FF7F00", "6"  = "#A65628",
  "7"  = "#F781BF", "8"  = "#666666", "9"  = "#1B9E77",
  "10" = "#E7298A"
)
p_colors <- c("20" = "#E41A1C", "100" = "#377EB8")




# ==============================================================================
# 1. Distribution of true beta values
# ==============================================================================

# --- generate data ---
scenarios <- expand.grid(
  p     = c(20, 100),
  g_nz  = c(0.25, 1),
  rho_w = c(0.2, 0.7),
  stringsAsFactors = FALSE
)

plot_df <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  sc <- scenarios[i, ]
  df <- gen_betas(p = sc$p, g = 10, non_zero_groups = sc$g_nz,
                  rho_w = sc$rho_w, target_snr = 2,
                  betas_plot = FALSE, seed = 123)
  df$p     <- sc$p
  df$g_nz  <- sc$g_nz
  df$rho_w <- sc$rho_w
  df
}))

plot_df <- plot_df %>%
  mutate(
    group_f      = factor(group, levels = 1:10),
    p            = factor(p),
    sparsity_lab = factor(
      ifelse(g_nz == 1, "Dense~(g[nz]==1)", "Sparse~(g[nz]==0.25)"),
      levels = c("Sparse~(g[nz]==0.25)", "Dense~(g[nz]==1)")
    )
  )

# --- plot function ---
make_true_beta_plot <- function(rho_val) {
  plot_df %>%
    filter(rho_w == rho_val) %>%
    ggplot(aes(x = group_f, y = beta, colour = group_f, fill = group_f)) +
    geom_boxplot(alpha = 0.35, outlier.shape = NA, width = 0.55, linewidth = 0.45) +
    geom_jitter(width = 0.18, size = 1.3, alpha = 0.75, shape = 16) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
    facet_nested(
      sparsity_lab ~ p,
      scales   = "free_y",
      switch   = "y",
      labeller = labeller(
        sparsity_lab = label_parsed,
        p            = function(x) paste0("p == ", x)
      )
    ) +
    scale_colour_manual(name = "Group", values = group_colors) +
    scale_fill_manual(  name = "Group", values = group_colors) +
    labs(
      title    = expression("Distribution of true " * bold(beta) * "-values"),
      subtitle = bquote(g == 10 ~ "|" ~ rho[b] == 0 ~ "," ~ t[snr] == 2 ~ "," ~ rho[w] == .(rho_val)),
      x = "Group",
      y = expression(beta)
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.placement    = "outside",
      strip.background   = element_rect(fill = "black"),
      strip.text         = element_text(colour = "white", face = "bold", size = 10),
      strip.text.y.left  = element_text(angle = 90),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position    = "right",
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(size = 10, colour = "grey30")
    )
}

# --- save ---
ggsave(file.path(PLOT_DIR, "true_betas_rho02.png"),
       plot = make_true_beta_plot(0.2), width = 8.27, height = 6.5, dpi = 300)

ggsave(file.path(PLOT_DIR, "true_betas_rho07.png"),
       plot = make_true_beta_plot(0.7), width = 8.27, height = 6.5, dpi = 300)




# ==============================================================================
# 2. Correlation heatmaps of X
# ==============================================================================

# --- generate data ---
scenarios_cor <- expand.grid(
  p     = c(20, 100),
  rho_w = c(0.2, 0.7),
  stringsAsFactors = FALSE
)

cor_df <- do.call(rbind, lapply(seq_len(nrow(scenarios_cor)), function(i) {
  sc      <- scenarios_cor[i, ]
  X       <- gen_X(n = 1000, p = sc$p, g = 10, rho_w = sc$rho_w,
                   standardise_X = TRUE, X_plots = FALSE, seed = 123)$X
  cor_mat <- cor(X)
  idx     <- seq_len(sc$p)
  data.frame(
    Var1  = rep((idx - 1) / (sc$p - 1), times = sc$p),
    Var2  = rep((idx - 1) / (sc$p - 1), each  = sc$p),
    value = as.vector(cor_mat),
    p     = sc$p,
    rho_w = sc$rho_w
  )
})) %>%
  mutate(p_lab = factor(paste0("p==", p), levels = paste0("p==", sort(unique(p)))))

# --- plot function ---
make_cor_plot <- function(rho_val) {
  cor_df %>%
    filter(rho_w == rho_val) %>%
    ggplot(aes(x = Var1, y = Var2, fill = value)) +
    geom_raster() +
    coord_fixed() +
    facet_wrap(~ p_lab, labeller = label_parsed) +
    scale_fill_viridis_c(name = "Correlation", limits = c(-1, 1), option = "viridis") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title    = "Correlation Structure of the Design Matrix X",
      subtitle = bquote(g == 10 ~ "|" ~ rho[b] == 0 ~ "," ~ n[obs] == 1000 ~ "," ~ rho[w] == .(rho_val)),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "black"),
      strip.text       = element_text(colour = "white", face = "bold", size = 10),
      panel.grid       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      legend.position  = "right",
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 10, colour = "grey30")
    )
}

# --- save ---
ggsave(file.path(PLOT_DIR, "X_corr_heatmaps_rho02.png"),
       plot = make_cor_plot(0.2), width = 8.27, height = 4.5, dpi = 300)

ggsave(file.path(PLOT_DIR, "X_corr_heatmaps_rho07.png"),
       plot = make_cor_plot(0.7), width = 8.27, height = 4.5, dpi = 300)




# ==============================================================================
# 3. RMSE (capped) for bio_age and mrl
# ==============================================================================

# --- load and prep ---
df_rmse <- read.csv(csv_path) %>%
  filter(metric %in% c("bio_age", "mrl")) %>%
  mutate(
    method = recode(method,
                    fAFT    = "Frequentist",
                    bAFT_hs = "Horseshoe",
                    bAFT_gs = "Gaussian"),
    method = factor(method, levels = c("Horseshoe", "Gaussian", "Frequentist")),
    metric = recode(metric,
                    bio_age = "Biological Age",
                    mrl     = "Mean Residual Life"),
    metric = factor(metric, levels = c("Biological Age", "Mean Residual Life")),
    p      = factor(p),
    n_obs  = factor(n_obs),
    sparsity_lab = factor(
      ifelse(non_zero == 1, "Dense~(g[nz]==1)", "Sparse~(g[nz]==0.25)"),
      levels = c("Sparse~(g[nz]==0.25)", "Dense~(g[nz]==1)")
    )
  )

# --- plot function ---
make_rmse_plot <- function(metric_name, rho_val) {
  df_rmse %>%
    filter(metric == metric_name, rho_w == rho_val) %>%
    ggplot(aes(x = n_obs, y = mc_mean, colour = p, group = p)) +
    geom_line(position = position_dodge(0.3), linewidth = 0.5) +
    geom_point(size = 2.5, position = position_dodge(0.3)) +
    geom_errorbar(
      aes(ymin = mc_mean - mc_se, ymax = mc_mean + mc_se),
      width = 0.08, linewidth = 1.2, position = position_dodge(0.3)
    ) +
    facet_nested(
      sparsity_lab ~ method,
      scales   = "free_y",
      switch   = "y",
      labeller = labeller(sparsity_lab = label_parsed)
    ) +
    scale_colour_manual(values = p_colors, name = "p") +
    scale_x_discrete(name = expression(n[obs])) +
    labs(
      title    = "Mean RMSE across 100 Monte Carlo Replications",
      subtitle = bquote(.(metric_name) ~ "|" ~ rho[w] == .(rho_val)),
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.placement   = "outside",
      strip.background  = element_rect(fill = "black"),
      strip.text        = element_text(colour = "white", face = "bold", size = 10),
      strip.text.y.left = element_text(angle = 90),
      panel.grid.minor  = element_blank(),
      legend.position   = "right",
      plot.title        = element_text(face = "bold", size = 12),
      plot.subtitle     = element_text(size = 10, colour = "grey30")
    )
}

# --- save ---
ggsave(file.path(FIGURES_DIR, "rmse_capped_bio_age_rho02.png"),
       plot = make_rmse_plot("Biological Age",    0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "rmse_capped_bio_age_rho07.png"),
       plot = make_rmse_plot("Biological Age",    0.7), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "rmse_capped_mrl_rho02.png"),
       plot = make_rmse_plot("Mean Residual Life", 0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "rmse_capped_mrl_rho07.png"),
       plot = make_rmse_plot("Mean Residual Life", 0.7), width = 10, height = 6.5, dpi = 300)




# ==============================================================================
# 4. Predicted vs true betas
# ==============================================================================

# --- load data ---
beta_df <- map_dfr(seq_len(nrow(param_grid)), function(i) {
  rds_file <- file.path(exp_dir, glue("Experiment_{exp_name}_giter_{i}.rds"))
  if (!file.exists(rds_file)) return(NULL)

  res    <- readRDS(rds_file)
  params <- param_grid[i, ]
  M      <- length(res$mc_results)

  map_dfr(methods, function(meth) {
    rep_data <- map_dfr(seq_len(M), function(m) {
      r <- res$mc_results[[m]]$results[[meth]]
      if (is.null(r)) return(NULL)
      tibble(
        rep       = m,
        j         = seq_along(r$predicted_betas),
        true_beta = r$true_betas$beta,
        pred_beta = as.numeric(r$predicted_betas),
        group     = r$true_betas$group
      )
    })
    if (nrow(rep_data) == 0) return(NULL)

    rep_data %>%
      group_by(j, group) %>%
      summarise(
        true_beta = mean(true_beta, na.rm = TRUE),
        pred_beta = mean(pred_beta, na.rm = TRUE),
        .groups   = "drop"
      ) %>%
      mutate(
        method   = meth,
        n_obs    = params$n_obs,
        p        = params$p,
        rho_w    = params$rho_w,
        non_zero = params$non_zero
      )
  })
})

beta_df <- beta_df %>%
  mutate(
    method = recode(method,
                    fAFT    = "Frequentist",
                    bAFT_hs = "Horseshoe",
                    bAFT_gs = "Gaussian"),
    method    = factor(method, levels = c("Horseshoe", "Gaussian", "Frequentist")),
    group     = factor(as.character(group)),
    n_obs_lab = factor(paste0("n[obs]==", n_obs),
                       levels = paste0("n[obs]==", sort(unique(n_obs))))
  )

# --- plot function ---
make_beta_plot <- function(nz, rho_val) {
  sparsity_label <- if (nz == 1) "Dense" else "Sparse"

  beta_df %>%
    filter(non_zero == nz, rho_w == rho_val) %>%
    ggplot(aes(x = true_beta, y = pred_beta, colour = group)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_point(size = 1.8, alpha = 0.8) +
    facet_nested(
      n_obs_lab ~ method,
      scales   = "free",
      switch   = "y",
      labeller = labeller(n_obs_lab = label_parsed)
    ) +
    scale_colour_manual(values = group_colors, name = "Group") +
    labs(
      title    = "Estimated vs True Coefficients across 100 Monte Carlo Replications",
      subtitle = bquote(g[nz] == .(nz) ~ "(" * .(sparsity_label) * ")" ~ "|" ~ rho[w] == .(rho_val)),
      x = expression(beta[j] ~ "(true)"),
      y = expression(hat(beta)[j] ~ "(estimated)")
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.placement   = "outside",
      strip.background  = element_rect(fill = "black"),
      strip.text        = element_text(colour = "white", face = "bold", size = 9),
      strip.text.y.left = element_text(angle = 90),
      panel.grid.minor  = element_blank(),
      legend.position   = "right",
      plot.title        = element_text(face = "bold", size = 12),
      plot.subtitle     = element_text(size = 10, colour = "grey30")
    )
}

# --- save ---
ggsave(file.path(FIGURES_DIR, "beta_scatter_sparse_rho02.png"),
       plot = make_beta_plot(0.25, 0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "beta_scatter_sparse_rho07.png"),
       plot = make_beta_plot(0.25, 0.7), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "beta_scatter_dense_rho02.png"),
       plot = make_beta_plot(1, 0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "beta_scatter_dense_rho07.png"),
       plot = make_beta_plot(1, 0.7), width = 10, height = 6.5, dpi = 300)




# ==============================================================================
# 5. Gompertz baseline parameters (tau, sigma)
# ==============================================================================

TRUE_TAU   <- 0.00145188
TRUE_SIGMA <- 11.76470588

# --- load data ---
distr_df <- map_dfr(seq_len(nrow(param_grid)), function(i) {
  rds_file <- file.path(exp_dir, glue("Experiment_{exp_name}_giter_{i}.rds"))
  if (!file.exists(rds_file)) return(NULL)

  res    <- readRDS(rds_file)
  params <- param_grid[i, ]
  M      <- length(res$mc_results)

  out <- map_dfr(methods, function(meth) {
    map_dfr(seq_len(M), function(m) {
      r <- res$mc_results[[m]]$results[[meth]]
      if (is.null(r) || is.null(r$distr_ests)) return(NULL)
      tibble(
        rep      = m,
        tau_hat  = as.numeric(r$distr_ests$shape),
        sig_hat  = as.numeric(r$distr_ests$scale),
        method   = meth,
        n_obs    = params$n_obs,
        p        = params$p,
        rho_w    = params$rho_w,
        non_zero = params$non_zero
      )
    })
  })
  rm(res); gc(); out
})

distr_long <- distr_df %>%
  mutate(
    method = recode(method,
                    fAFT    = "Frequentist",
                    bAFT_hs = "Horseshoe",
                    bAFT_gs = "Gaussian"),
    method  = factor(method, levels = c("Horseshoe", "Gaussian", "Frequentist")),
    p       = factor(p),
    n_obs   = factor(n_obs),
    sparsity_lab = factor(
      ifelse(non_zero == 1, "Dense~(g[nz]==1)", "Sparse~(g[nz]==0.25)"),
      levels = c("Sparse~(g[nz]==0.25)", "Dense~(g[nz]==1)")
    )
  ) %>%
  pivot_longer(cols = c(tau_hat, sig_hat), names_to = "param", values_to = "estimate") %>%
  mutate(
    param = recode(param, tau_hat = "Shape (tau)", sig_hat = "Scale (sigma)"),
    param = factor(param, levels = c("Shape (tau)", "Scale (sigma)"))
  )

# --- plot function ---
make_distr_plot <- function(param_name, rho_val) {
  true_val <- if (param_name == "Shape (tau)") TRUE_TAU else TRUE_SIGMA

  distr_long %>%
    filter(param == param_name, rho_w == rho_val) %>%
    ggplot(aes(x = n_obs, y = estimate, colour = p, fill = p)) +
    stat_pointinterval(
      aes(group = p),
      position       = position_dodge(0.6),
      point_size     = 2,
      linewidth      = 0.8,
      .width         = c(0.50, 0.95),
      point_interval = median_qi,
      alpha          = 0.7
    ) +
    geom_hline(yintercept = true_val,
               linetype = "dashed", colour = "grey30", linewidth = 0.6) +
    facet_nested(
      sparsity_lab ~ method,
      scales   = "free_y",
      switch   = "y",
      labeller = labeller(sparsity_lab = label_parsed)
    ) +
    scale_colour_manual(values = p_colors, name = "p") +
    scale_fill_manual(  values = p_colors, name = "p") +
    scale_x_discrete(name = expression(n[obs])) +
    labs(
      title    = "Recovery of Gompertz Baseline Parameters across 100 MC Replications",
      subtitle = bquote(.(param_name) ~ "|" ~ rho[w] == .(rho_val)),
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.placement   = "outside",
      strip.background  = element_rect(fill = "black"),
      strip.text        = element_text(colour = "white", face = "bold", size = 9),
      strip.text.y.left = element_text(angle = 90),
      panel.grid.minor  = element_blank(),
      legend.position   = "right",
      plot.title        = element_text(face = "bold", size = 12),
      plot.subtitle     = element_text(size = 10, colour = "grey30")
    )
}

# --- save ---
ggsave(file.path(FIGURES_DIR, "distr_params_tau_rho02.png"),
       plot = make_distr_plot("Shape (tau)",   0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "distr_params_tau_rho07.png"),
       plot = make_distr_plot("Shape (tau)",   0.7), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "distr_params_sigma_rho02.png"),
       plot = make_distr_plot("Scale (sigma)", 0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "distr_params_sigma_rho07.png"),
       plot = make_distr_plot("Scale (sigma)", 0.7), width = 10, height = 6.5, dpi = 300)




# ==============================================================================
# 6. Computational cost (in minutes)
# ==============================================================================

# ------------------------
# --- individual plots ---
# ------------------------

# --- load and prep ---
timing <- read.csv("output/MCexp/summaries/all_runs_timing.csv") %>%
  filter(experiment == "full_grid_run_01", !is.na(n_obs), !is.na(non_zero)) %>%
  dplyr::select(n_obs, p, rho_w, non_zero, n_reps,
                Horseshoe = bAFT_hs_mins,
                Gaussian  = bAFT_gs_mins,
                Frequentist = fAFT_mins) %>%
  pivot_longer(
    cols      = c(Horseshoe, Gaussian, Frequentist),
    names_to  = "method",
    values_to = "mean_mins"
  ) %>%
  mutate(
    mean_mins    = mean_mins / n_reps,
    method       = factor(method, levels = c("Horseshoe", "Gaussian", "Frequentist")),
    p            = factor(p),
    n_obs        = factor(n_obs),
    sparsity_lab = factor(
      ifelse(non_zero == 1, "Dense~(g[nz]==1)", "Sparse~(g[nz]==0.25)"),
      levels = c("Sparse~(g[nz]==0.25)", "Dense~(g[nz]==1)")
    )
  )

# --- plot function ---
make_timing_plot <- function(rho_val) {
  timing %>%
    filter(rho_w == rho_val) %>%
    ggplot(aes(x = n_obs, y = mean_mins, fill = p)) +
    geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
    facet_nested(
      sparsity_lab ~ method,
      labeller = labeller(sparsity_lab = label_parsed)
    ) +
    scale_fill_manual(values = p_colors, name = "p") +
    scale_x_discrete(name = expression(n[obs])) +
    scale_y_continuous(
      name   = "Mean fit time per replication (minutes)",
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title    = "Computational Cost per Replication across 100 Monte Carlo Replications",
      subtitle = bquote(rho[w] == .(rho_val))
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.placement    = "outside",
      strip.background   = element_rect(fill = "black"),
      strip.text         = element_text(colour = "white", face = "bold", size = 10),
      strip.text.y.left  = element_text(angle = 90),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position    = "right",
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(size = 10, colour = "grey30")
    )
}

# --- save ---
ggsave(file.path(FIGURES_DIR, "timing_rho02.png"),
       plot = make_timing_plot(0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "timing_rho07.png"),
       plot = make_timing_plot(0.7), width = 10, height = 6.5, dpi = 300)


# ------------------------
# ----- wide format -----
# ------------------------

timing <- read.csv("output/MCexp/summaries/all_runs_timing.csv") %>%
  filter(experiment == "full_grid_run_01", !is.na(n_obs), !is.na(non_zero)) %>%
  dplyr::select(n_obs, p, rho_w, non_zero, n_reps,
                Horseshoe   = bAFT_hs_mins,
                Gaussian    = bAFT_gs_mins,
                Frequentist = fAFT_mins) %>%
  pivot_longer(
    cols      = c(Horseshoe, Gaussian, Frequentist),
    names_to  = "method",
    values_to = "mean_mins"
  ) %>%
  mutate(
    mean_mins    = mean_mins / n_reps,
    method       = factor(method, levels = c("Horseshoe", "Gaussian", "Frequentist")),
    p            = factor(p),
    n_obs        = factor(n_obs),
    sparsity_lab = factor(
      ifelse(non_zero == 1, "Dense~(g[nz]==1)", "Sparse~(g[nz]==0.25)"),
      levels = c("Sparse~(g[nz]==0.25)", "Dense~(g[nz]==1)")
    ),
    rho_w_lab = factor(             # <- add this
      paste0("rho[w]==", rho_w),
      levels = paste0("rho[w]==", sort(unique(rho_w)))
    )
  )

make_timing_plot_combined <- function() {
  timing %>%
    ggplot(aes(x = n_obs, y = mean_mins, fill = p)) +
    geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
    facet_nested(
      sparsity_lab ~ rho_w_lab + method,
      scales   = "free_y",
      switch   = "y",
      labeller = label_parsed
    ) +
    scale_fill_manual(values = p_colors, name = "p") +
    scale_x_discrete(name = expression(n[obs])) +
    scale_y_continuous(
      name   = "Mean fit time per replication (minutes)",
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(title = "Computational Cost per Replication across 100 Monte Carlo Replications") +
    theme_bw(base_size = 11) +
    theme(
      strip.placement    = "outside",
      strip.background   = element_rect(fill = "black"),
      strip.text         = element_text(colour = "white", face = "bold", size = 10),
      strip.text.y.left  = element_text(angle = 90),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position    = "right",
      plot.title         = element_text(face = "bold", size = 12)
    )
}

# --- save ---
ggsave(file.path(FIGURES_DIR, "timing_wide.png"),
       plot = make_timing_plot_combined(),
       width = 20, height = 6.5, dpi = 300)



# ------------------------
# ----- long format -----
# ------------------------

make_timing_plot_long <- function() {
  timing %>%
    ggplot(aes(x = n_obs, y = mean_mins, fill = p)) +
    geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
    facet_nested(
      rho_w_lab + sparsity_lab ~ method,
      scales   = "free_y",
      switch   = "y",
      labeller = label_parsed
    ) +
    scale_fill_manual(values = p_colors, name = "p") +
    scale_x_discrete(name = expression(n[obs])) +
    scale_y_continuous(
      name   = "Mean fit time per replication (minutes)",
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(title = "Computational Cost per Replication across 100 Monte Carlo Replications") +
    theme_bw(base_size = 11) +
    theme(
      strip.placement    = "outside",
      strip.background   = element_rect(fill = "black"),
      strip.text         = element_text(colour = "white", face = "bold", size = 10),
      strip.text.y.left  = element_text(angle = 90),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position    = "right",
      plot.title         = element_text(face = "bold", size = 12)
    )
}

# --- save ---
ggsave(file.path(FIGURES_DIR, "timing_long.png"),
       plot = make_timing_plot_long(),
       width = 10, height = 13, dpi = 300)




# ==============================================================================
# 7. Capped vs uncapped RMSE
# ==============================================================================

# --- load and prep ---
df_comp <- read.csv(csv_path) %>%
  filter(metric %in% c("bio_age", "mrl")) %>%
  mutate(
    method = recode(method,
                    fAFT    = "Frequentist",
                    bAFT_hs = "Horseshoe",
                    bAFT_gs = "Gaussian"),
    method  = factor(method, levels = c("Horseshoe", "Gaussian", "Frequentist")),
    metric  = recode(metric,
                     bio_age = "Biological Age",
                     mrl     = "Mean Residual Life"),
    metric  = factor(metric, levels = c("Biological Age", "Mean Residual Life")),
    p       = factor(p),
    n_obs   = factor(n_obs),
    sparsity_lab = factor(
      ifelse(non_zero == 1, "Dense~(g[nz]==1)", "Sparse~(g[nz]==0.25)"),
      levels = c("Sparse~(g[nz]==0.25)", "Dense~(g[nz]==1)")
    ),
    rmse_diff = mc_mean_all - mc_mean
  )

df_long_comp <- df_comp %>%
  pivot_longer(cols = c(mc_mean, mc_mean_all), names_to = "version", values_to = "rmse") %>%
  mutate(
    version = recode(version,
                     mc_mean     = "Capped (\u2264150)",
                     mc_mean_all = "Uncapped"),
    version = factor(version, levels = c("Capped (\u2264150)", "Uncapped"))
  )

version_colors <- c("Capped (\u2264150)" = "#2166AC", "Uncapped" = "#D73027")

# --- plot functions ---
make_comp_plot <- function(metric_name, rho_val) {
  df_long_comp %>%
    filter(metric == metric_name, rho_w == rho_val) %>%
    ggplot(aes(x = n_obs, y = rmse, colour = version,
               group = interaction(p, version))) +
    geom_line(aes(linetype = p), position = position_dodge(0.4), linewidth = 0.6) +
    geom_point(aes(shape = p), size = 2.2, position = position_dodge(0.4)) +
    facet_nested(
      sparsity_lab ~ method,
      scales   = "free_y",
      switch   = "y",
      labeller = labeller(sparsity_lab = label_parsed)
    ) +
    scale_colour_manual(values = version_colors, name = "RMSE version") +
    scale_linetype_manual(values = c("20" = "solid", "100" = "dashed"), name = "p") +
    scale_shape_manual(  values = c("20" = 16,      "100" = 17),       name = "p") +
    scale_x_discrete(name = expression(n[obs])) +
    labs(
      title    = "Capped vs Uncapped RMSE: identifying extreme estimation errors",
      subtitle = bquote(.(metric_name) ~ "|" ~ rho[w] == .(rho_val)),
      y        = "RMSE"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.placement   = "outside",
      strip.background  = element_rect(fill = "black"),
      strip.text        = element_text(colour = "white", face = "bold", size = 10),
      strip.text.y.left = element_text(angle = 90),
      panel.grid.minor  = element_blank(),
      legend.position   = "right",
      plot.title        = element_text(face = "bold", size = 12),
      plot.subtitle     = element_text(size = 10, colour = "grey30")
    )
}

make_diff_plot <- function(metric_name, rho_val) {
  df_comp %>%
    filter(metric == metric_name, rho_w == rho_val) %>%
    ggplot(aes(x = n_obs, y = rmse_diff, fill = p)) +
    geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
    facet_nested(
      sparsity_lab ~ method,
      scales   = "free_y",
      switch   = "y",
      labeller = labeller(sparsity_lab = label_parsed)
    ) +
    scale_fill_manual(values = c("20" = "#E41A1C", "100" = "#377EB8"), name = "p") +
    scale_x_discrete(name = expression(n[obs])) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title    = "Impact of capping: difference between uncapped and capped RMSE",
      subtitle = bquote(.(metric_name) ~ "|" ~ rho[w] == .(rho_val)),
      y        = "Uncapped \u2212 Capped RMSE"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.placement    = "outside",
      strip.background   = element_rect(fill = "black"),
      strip.text         = element_text(colour = "white", face = "bold", size = 10),
      strip.text.y.left  = element_text(angle = 90),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position    = "right",
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(size = 10, colour = "grey30")
    )
}

# --- save ---
ggsave(file.path(FIGURES_DIR, "rmse_capped_vs_all_bio_age_rho02.png"),
       plot = make_comp_plot("Biological Age",    0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "rmse_capped_vs_all_bio_age_rho07.png"),
       plot = make_comp_plot("Biological Age",    0.7), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "rmse_capped_vs_all_mrl_rho02.png"),
       plot = make_comp_plot("Mean Residual Life", 0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "rmse_capped_vs_all_mrl_rho07.png"),
       plot = make_comp_plot("Mean Residual Life", 0.7), width = 10, height = 6.5, dpi = 300)

ggsave(file.path(FIGURES_DIR, "rmse_diff_bio_age_rho02.png"),
       plot = make_diff_plot("Biological Age",    0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "rmse_diff_bio_age_rho07.png"),
       plot = make_diff_plot("Biological Age",    0.7), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "rmse_diff_mrl_rho02.png"),
       plot = make_diff_plot("Mean Residual Life", 0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "rmse_diff_mrl_rho07.png"),
       plot = make_diff_plot("Mean Residual Life", 0.7), width = 10, height = 6.5, dpi = 300)




# ==============================================================================
# 8. Predicted MRL vs Predicted Biological Age
# ==============================================================================

# --- load data ---
scatter_df <- map_dfr(seq_len(nrow(param_grid)), function(i) {
  rds_file <- file.path(exp_dir, glue("Experiment_{exp_name}_giter_{i}.rds"))
  if (!file.exists(rds_file)) return(NULL)

  res    <- readRDS(rds_file)
  params <- param_grid[i, ]
  M      <- length(res$mc_results)

  out <- map_dfr(methods, function(meth) {
    map_dfr(seq_len(M), function(m) {
      r <- res$mc_results[[m]]$results[[meth]]
      if (is.null(r)) return(NULL)
      td <- r$test_data
      tibble(
        rep           = m,
        pred_bio_age  = as.numeric(r$predicted_bio_age),
        pred_mrl      = as.numeric(r$predicted_mrl),
        true_bio_age  = as.numeric(td$b),
        true_mrl      = as.numeric(td$mrl),
        bio_age_error = as.numeric(r$predicted_bio_age) - as.numeric(td$b),
        mrl_error     = as.numeric(r$predicted_mrl)     - as.numeric(td$mrl),
        method        = meth,
        n_obs         = params$n_obs,
        p             = params$p,
        rho_w         = params$rho_w,
        non_zero      = params$non_zero
      )
    })
  })
  rm(res); gc(); return(out)
})

scatter_df <- scatter_df %>%
  mutate(
    method = recode(method,
                    fAFT    = "Frequentist",
                    bAFT_hs = "Horseshoe",
                    bAFT_gs = "Gaussian"),
    method  = factor(method, levels = c("Horseshoe", "Gaussian", "Frequentist")),
    p       = factor(p),
    n_obs   = factor(n_obs),
    sparsity_lab = factor(
      ifelse(non_zero == 1, "Dense~(g[nz]==1)", "Sparse~(g[nz]==0.25)"),
      levels = c("Sparse~(g[nz]==0.25)", "Dense~(g[nz]==1)")
    ),
    bio_age_error_cap = pmax(pmin(bio_age_error, 150), -150),
    mrl_error_cap     = pmax(pmin(mrl_error,     150), -150)
  )

# --- plot function ---
make_scatter_plot <- function(nobs_val, p_val, rho_val) {
  scatter_df %>%
    filter(n_obs == nobs_val, p == p_val, rho_w == rho_val, rep <= 20) %>%
    ggplot(aes(x = pred_bio_age, y = pred_mrl, colour = bio_age_error_cap)) +
    geom_point(size = 1, alpha = 0.4) +
    facet_nested(
      sparsity_lab ~ method,
      scales   = "free",
      switch   = "y",
      labeller = labeller(sparsity_lab = label_parsed)
    ) +
    scale_colour_gradient2(
      low      = "#2166AC",
      mid      = "grey90",
      high     = "#D73027",
      midpoint = 0,
      name     = "Bio Age Error\n(pred \u2212 true)"
    ) +
    labs(
      title    = "Predicted MRL vs Predicted Biological Age",
      subtitle = bquote(n[obs] == .(nobs_val) ~ "," ~ p == .(p_val) ~ "," ~ rho[w] == .(rho_val)),
      x = expression(hat(b) ~ "(predicted bio age)"),
      y = expression(hat(MRL) ~ "(predicted MRL)")
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.placement   = "outside",
      strip.background  = element_rect(fill = "black"),
      strip.text        = element_text(colour = "white", face = "bold", size = 9),
      strip.text.y.left = element_text(angle = 90),
      panel.grid.minor  = element_blank(),
      legend.position   = "right",
      plot.title        = element_text(face = "bold", size = 12),
      plot.subtitle     = element_text(size = 10, colour = "grey30")
    )
}

# --- save ---
ggsave(file.path(FIGURES_DIR, "scatter_mrl_bioage_n250_p20_rho02.png"),
       plot = make_scatter_plot(250, 20,  0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "scatter_mrl_bioage_n250_p20_rho07.png"),
       plot = make_scatter_plot(250, 20,  0.7), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "scatter_mrl_bioage_n250_p100_rho02.png"),
       plot = make_scatter_plot(250, 100, 0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "scatter_mrl_bioage_n250_p100_rho07.png"),
       plot = make_scatter_plot(250, 100, 0.7), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "scatter_mrl_bioage_n500_p20_rho02.png"),
       plot = make_scatter_plot(500, 20,  0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "scatter_mrl_bioage_n500_p20_rho07.png"),
       plot = make_scatter_plot(500, 20,  0.7), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "scatter_mrl_bioage_n500_p100_rho02.png"),
       plot = make_scatter_plot(500, 100, 0.2), width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "scatter_mrl_bioage_n500_p100_rho07.png"),
       plot = make_scatter_plot(500, 100, 0.7), width = 10, height = 6.5, dpi = 300)

