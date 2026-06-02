# code created based on 
# https://github.com/marije-sluiskes/AccelerAge_SimulationStudy/
# and 
# https://github.com/Phil-Ber/Thesis_AccelerAge

# ==============================================================================
# Set-up
# ==============================================================================

# fall-back directory for testing and debugigng
SIM_DIR  <- if (exists("SIM_DIR",  inherits = TRUE)) SIM_DIR  else "output/scenarios/simulated"
PLOT_DIR <- if (exists("PLOT_DIR", inherits = TRUE)) PLOT_DIR else "output/scenarios/plots"
dir.create(SIM_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# Packages
# ==============================================================================
library(MASS) 
library(Matrix)
library(eha)
library(survival)
library(scam)
library(glue)
library(dplyr)
library(ggplot2)


# ==============================================================================
# Gompertz AFT - helper functions
# ==============================================================================

# Gompertz baseline survival probability 
gomp_baseline_surv <- function(t, gomp_a, gomp_b){
  exp( (gomp_a / gomp_b) * (1 - exp(gomp_b * t)))
}

# AFT survival at baseline 
S0 <- function(t, tau, sigma) {
  exp(-tau * (exp(t / sigma) - 1))
}

# AFT survival S(t | linpred) = S0(t * exp(linpred)) ; linpred = beta'X
S_gomp_aft <- function(t, linpred, tau, sigma) {
  exp(-tau * (exp((t * exp(linpred)) / sigma) - 1))
}

# draw T from Gompertz AFT (age-at-event from birth) via inverse-CDF using
# S(t | linpred) = exp(-tau * (exp( (t * exp(linpred))/sigma)-1) )
rgompertz_aft <- function(n, sigma, tau, linpred) {
  u <- runif(n)
  sigma * log(1 - log(u) / tau) / exp(linpred)
}


# baseline AFT Gompertz draw with linpred = 0
rgompertz0_aft <- function(n, gomp_a = exp(-9), gomp_b = 0.085) {
  
  sigma <- 1 / gomp_b           # scale 
  tau   <- gomp_a / gomp_b      # shape 
  
  u <- runif(n)
  T0 = sigma * log(1 - log(u) / tau)
  
  return(T0)
}

# epsilon variance = Var(log(T0))
epsilon_var_gompertz <- function(n_sim, gomp_a = exp(-9), gomp_b = 0.085) {
  T0 <- rgompertz0_aft(n_sim, gomp_a = gomp_a, gomp_b = gomp_b)
  var(log(T0))
}





# ==============================================================================
# X generation
# (g blocks, within-block corr = rho_w, between-block = 0)
# ==============================================================================

gen_X <- function(n, 
                  p, 
                  g = 10, 
                  rho_w,
                  standardise_X = T,
                  X_plots = F,
                  seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  

  # create variance-covariance matrix, Sigma (within-block = rho_w, between = 0)
  make_block_cov <- function(p, g, rho_w) {
    if (p %% g != 0) stop("p must be divisible by g")
    
    # block size
    m <- p / g
    # initialise
    Sigma <- matrix(0, nrow = p, ncol = p)
    # loop over all groups
    for (grp in 1:g) {
      idx <- ((grp - 1) * m + 1):(grp * m)
      Sigma[idx, idx] <- rho_w
      diag(Sigma)[idx] <- 1
    }
    Sigma
  }
  Sigma <- make_block_cov(p = p, g = g, rho_w = rho_w)
  
  # generate matrix of predictors from multivariate normal 
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  # set colnames
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  
  # optionally, standardise X
  if (isTRUE(standardise_X)) X <- scale(X)
  
  
  # group membership vector 
  p_g <- p / g
  group_membership <- rep(1:g, each = p_g)
  
  
  # plots
  if (isTRUE(X_plots)) {
    
    # visualise correlation structure
    heatmap(Sigma, Rowv = NA, Colv = NA, scale = "none",
            main = "Sigma structure (block corr)")
    
    # visualise all rows and columns
    heatmap(X, Rowv = NA, Colv = NA, scale = "none",
            main = "X matrix")
    
    X_long <- data.frame(
      Value = as.vector(X),
      Group = rep(rep(1:g, each = p_g), times = n),
      Variable = rep(1:p, each = n)
    )
    
    # boxplot of all variables by group
    boxplot(Value ~ Group, data = X_long,
            main = "Distribution of X values by group",
            xlab = "Group", ylab = "Value")
    abline(h = 0, lty = 2, col = "gray50")
    
    heatmap(cor(X), Rowv = NA, Colv = NA, scale = "none",
            main = "Empirical cor(X)")
    
    # heathmap 
    cor(X) %>% heatmap(Rowv = NA, Colv = NA, 
                       scale = "none", main = "X Matrix Correlations")
  }
  
  
  list(
    X = X,
    Sigma = Sigma,
    group_membership = group_membership
  )
  
}



# ==============================================================================
# Beta generation
# (g blocks, within-block corr = rho_w, between-block = 0)
# ==============================================================================

gen_betas <- function(p, 
                      g = 10, 
                      non_zero_groups,
                      mu_l = -4, mu_u = 4,
                      rho_w,
                      target_snr = 2,
                      gomp_a = exp(-9), 
                      gomp_b = 0.085,
                      standardise_X = T,
                      betas_plot = T, 
                      seed = 123) { 
  
  if (!is.null(seed)) set.seed(seed)
  if (p %% g != 0) stop("p must be divisible by g")
  
  p_g <- p / g
  group_membership <- rep(1:g, each = p_g)
  
  # randomly select which groups will have non-zero coefficients
  active_groups <- sample(1:g, size = round(non_zero_groups * g))
  
  # sample group means
  mu_groups <- runif(g, min = mu_l, max = mu_u)
  
  # mean vector mu_full for betas (only active groups get a group mean)
  mu_full <- numeric(p)
  for (grp in active_groups) {
    idx <- ((grp - 1) * p_g + 1):(grp * p_g)
    mu_full[idx] <- mu_groups[grp]
  }
  
  # create variance-covariance matrix, Sigma (within-block = rho_w, between = 0)
  make_block_cov <- function(p, g, rho_w) {
    if (p %% g != 0) stop("p must be divisible by g")
    
    # block size
    m <- p / g
    # initialise
    Sigma <- matrix(0, nrow = p, ncol = p)
    # loop over all groups
    for (grp in 1:g) {
      idx <- ((grp - 1) * m + 1):(grp * m)
      Sigma[idx, idx] <- rho_w
      diag(Sigma)[idx] <- 1
    }
    Sigma
  }
  Sigma_beta <- make_block_cov(p = p, g = g, rho_w = rho_w)
  
  # sample betas ~ N(mu_full, Sigma_beta)
  betas <- as.numeric(mvrnorm(n = 1, mu = mu_full, Sigma = Sigma_beta))
  
  # set inactive groups exactly to 0
  inactive <- setdiff(1:g, active_groups)
  for (grp in inactive) {
    idx <- ((grp - 1) * p_g + 1):(grp * p_g)
    betas[idx] <- 0
  }
  
  print(glue("Mean of betas pre-scaling: { round(mean(betas), 3) }"))
  print(glue("{c('Lower', 'Upper')} range of beta values pre-scaling: { round(range(betas), 3) }"))
  
  # Scale betas 
  #betas = betas * beta_scale
  
  # initialise 
  epsilon_var  <- NA_real_
  scale_factor <- NA_real_
  achieved_snr <- NA_real_
  
  # SNR scaling
  if (!is.null(target_snr) && target_snr != 0) {
    n_sim <- 1e5
    
    # baseline noise variance: Var(log(T0))
    baseline_T  <- rgompertz0_aft(n_sim, gomp_a = gomp_a, gomp_b = gomp_b)
    epsilon_var <- var(log(baseline_T))
    
    # simulate X_scale under the same block structure (within=rho, between=0)
    X_scale <- gen_X(n = n_sim, p = p, g = g, rho_w = rho_w,
                     seed = seed, standardise_X = standardise_X)
    
    
    # centre active coefficients so E[beta] = 0
    active_idx <- which(group_membership %in% active_groups)
    betas[active_idx] <- betas[active_idx] - mean(betas[active_idx])
    
    
    # compute linear predictor
    linpred_raw <- as.numeric(X_scale$X %*% betas)
    
    # set target SNR and rescale
    target_var_linpred  <- target_snr * epsilon_var
    current_var_linpred <- var(linpred_raw)
    
    # compute scaling factor, sf
    scale_factor <- as.numeric(sqrt(target_var_linpred / current_var_linpred))
    
    # rescale betas
    betas <- betas * scale_factor
    
    # checked achieved SNR
    linpred_final <- as.numeric(X_scale$X %*% betas)
    achieved_snr  <- var(linpred_final) / epsilon_var
    
    print(glue("Target SNR: {target_snr}"))
    print(glue("Epsilon variance: {round(epsilon_var, 4)}"))
    print(glue("Achieved SNR: {round(achieved_snr, 2)}"))
  }
  
  # output  
  beta_df <- data.frame(beta = betas, group = group_membership)
  
  
  # plot betas
  if (isTRUE(betas_plot)) {

    heatmap(Sigma_beta, Rowv = NA, Colv = NA, scale = "none", main = "Beta covariance (Sigma_beta)")
    
    ((betas) %*% t(betas)) %>% 
      heatmap(Rowv = NA, Colv = NA, scale = "none", main = "Beta Gram Matrix")
    
    hist(betas[betas != 0], breaks = 50, main = "Non-zero betas", xlab = "beta")
    abline(v = 0, lty = 2)
    
    beta_by_group <- split(beta_df$beta, beta_df$group)
    boxplot(beta_by_group,
            main = "Beta distribution by group",
            xlab = "Group", ylab = "Beta",
            col = rainbow(g),
            border = "black",
            outline = TRUE)
    # add a horizontal line at y=0 for reference
    abline(h = 0, lty = 2, col = "gray50")
  }
  
  return(beta_df)
  
}


# ==============================================================================
# Population life table (Gompertz AFT)
# ==============================================================================

gen_pop_lt_gompertz = function(
    n_pop = 1e5,
    p,
    g = 10,
    non_zero_groups,
    rho_w,
    gomp_a = exp(-9),
    gomp_b = 0.085,
    betas,
    X_pop,
    force_recalc = TRUE,
    filename = NULL,
    lt_plot = F,
    seed = 123) {
  
  # file path
  path <- file.path(SIM_DIR, glue("pop_lifetable_{filename}.rds"))
  
  # if file does not exist or force_recalc=T, recompute pop lt
  if (!file.exists(path) || force_recalc) {
    
    set.seed(seed)
    
    sigma <- 1 / gomp_b
    tau   <- gomp_a / gomp_b
    
    
    # computed linear predictor (numeric)
    linpred <- drop(X_pop %*% betas)
    
    # diagnostics
    stopifnot(length(linpred) == n_pop)
    stopifnot(is.numeric(linpred))
    
    
    # simulate population event times
    t = numeric(n_pop) 
    for (i in 1:n_pop){
      t[i] <- rgompertz_aft(1, sigma = sigma, tau = tau, linpred = linpred[i])
    }
    
    
    lifetable_pop <- as.data.frame(X_pop)
    lifetable_pop$t <- t
    lifetable_pop <- lifetable_pop[order(lifetable_pop$t), ]
    
    # mean residual life
    mrl_pop <- numeric(n_pop)
    for (j in 1:n_pop) {
      mrl_pop[j] <- mean(lifetable_pop$t[j:n_pop]) - lifetable_pop$t[j]
    }
    lifetable_pop$mrl_pop <- mrl_pop
    
    # diagnostics
    cat("\nSummary: MRL Population\n")
    print(summary(mrl_pop))
    cat("\nSummary: Survival Times (t)\n")
    print(summary(t))
    
    # # smooth - gives negative values for lt$mrl_pop
    # fit4 <- scam(mrl_pop ~ s(t, bs = "mpd"), data = lifetable_pop) 
    # xx <- seq(0, max(lifetable_pop$t), by = 0.1) 
    # lt <- as.data.frame(cbind(t = xx, mrl_pop = predict(fit4, data.frame(t=xx))))
    
    
    # smooth - fix negative mrl values
    eps <- 1e-6 # tiny constant to avoid log(0)
    fit4 <- scam(log(mrl_pop + eps) ~ s(t, bs = "mpd"),
                 data = lifetable_pop)
    xx <- seq(0, max(lifetable_pop$t), by = 0.1)
    lt <- data.frame(
      t = xx,
      # transform back
      mrl_pop = exp(predict(fit4, data.frame(t = xx))) - eps
    )
    
    # save life table 
    saveRDS(lt, path)
 
    } else {
    # if life table exists, read it
    lt = readRDS(path)
  }
  
  # plot when smoothing
  if (isTRUE(lt_plot)) {
    print(
      ggplot(lt, aes(x = t, y = mrl_pop)) +
        geom_line() +
        labs(title = "Gompertz Population Lifetable MRL",
             x = "Time (t)", y = "Mean Residual Life") +
        theme_minimal()
    )
  }
  
  
  return(lt)
}


# ==============================================================================
# Create data sets from Gompertz AFT 
# (new individuals, left truncation + right censoring)
# ==============================================================================

create_dataset_gompertz = function(
    n_obs,
    p,
    gomp_a = exp(-9),
    gomp_b = 0.085,
    g = 10,
    rho_w,
    betas,
    followup = 20,
    standardise_X = T,
    X_plots = F,
    lt,
    seed = 123) {
  
  sigma = 1/gomp_b
  tau = gomp_a / gomp_b
  
  # sample 5x as many to ensure enough observations, because for some T < C => not observed
  n_gen <- 5 * n_obs
  
  X_obj = gen_X(n = n_gen, 
                p = p, 
                g = g, 
                rho_w = rho_w,
                standardise_X = standardise_X,
                X_plots = X_plots,
                seed = seed)
  
  # predictors
  X = as.matrix(X_obj$X)
  
  x_names <- colnames(X)[seq_len(p)]
  
  # ages at entry
  entry_age <- runif(n_gen, 20, 80)       
  
  # compute linear predictor (numeric)
  linpred <- drop(X %*% betas)
  
  
  # diagnostics
  stopifnot(length(linpred) == nrow(X))       
  stopifnot(is.numeric(linpred))
  
  
  # get age of death
  age_death_t <- numeric(n_gen)
  for (i in 1:n_gen){
    age_death_t[i] <- rgompertz_aft(1, sigma = sigma, tau = tau, linpred = linpred[i])
  }
  
  
  # remove observations that are left-truncated
  ok <- which(entry_age < age_death_t)  # get all eligible
  stopifnot(length(ok) >= n_obs)        # safety check if not enough n_test
  indx_obs <- sample(ok, n_obs)
  
  # create data frame with valid observations only
  df_sim <- as.data.frame(
    cbind(X, age_death_t, entry_age, linpred)[indx_obs,]
  )
  
  cat("Range of linpred values:", range(linpred), "\n")
  if(any(is.infinite(exp(linpred)))) {
    warning("Some exp(linpred) values are Inf")
  }
  
  
  # get mean residual life (mrl) 
  df_sim$mrl <- numeric(nrow(df_sim)) 
  for (i in 1:nrow(df_sim)){
    
    # conditional survival at entry: denominator S(c | X)
    s_cond <- gomp_baseline_surv(df_sim$entry_age[i] * exp(df_sim$linpred[i]), 
                                 gomp_a = gomp_a, gomp_b = gomp_b)
    
    # integral on scaled-time axis
    t_unadj <- integrate(gomp_baseline_surv,
                            lower = df_sim$entry_age[i] * exp(df_sim$linpred[i]),
                            upper = Inf, 
                            gomp_a = gomp_a, gomp_b = gomp_b)$value
    
    # convert back to original time scale 
    t_adj <- t_unadj / exp(df_sim$linpred[i])
    
    # mean residual life
    df_sim$mrl[i] <- (t_adj / s_cond)
  }
  
  
  # get biological age (via population life table)
  bio_age <- numeric(nrow(df_sim))
  for (i in 1:nrow(df_sim)) {
    bio_age[i] <- lt$t[ which.min(abs(lt$mrl_pop - df_sim$mrl[i])) ]
  }
  df_sim$b <- bio_age
  
  
  # add censoring
  df_sim$yrs_rem <- df_sim$age_death_t - df_sim$entry_age     # true, remaining years
  wh <- which(df_sim$yrs_rem > followup)                      # TRUE=censored
  
  # censoring for frequentist AFT (aftreg::eha)
  df_sim$status     <- as.integer(df_sim$yrs_rem <= followup) # 1=event ; 0=censored
  # censoring for Bayesian AFT (brms cens())
  df_sim$cens_ind <- 1L - df_sim$status                       # 0=none/event, 1=right-censored
  
  
  # df_sim$status <- 1                                        # event
  # df_sim$status[wh] <- 0                                    # no event
  df_sim$follow_up_time <- df_sim$yrs_rem
  df_sim$follow_up_time[wh] <- followup
  df_sim$age_end <- df_sim$entry_age + df_sim$follow_up_time  # age at the end of observation
  
  # save everything 
  file_name <- paste0(
    "gompertz_p", p,
    "_n", n_obs,
    "_seed", seed,
    ".rds"
  )
  
  saveRDS(df_sim, file.path(SIM_DIR, file_name))
  
  return(df_sim)
}


  