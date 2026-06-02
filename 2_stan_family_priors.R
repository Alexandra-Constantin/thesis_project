# code created based on 
# https://metrumrg.com/wp-content/uploads/2024/12/2024-Yoder-cure-rate-brms.pdf


################################################################################
# This file defines:
#  - Stan custom functions (stan_funs)
#  - brms stanvars object (stanvars)
#  - custom_family for Gompertz AFT (gompertz_aft)
#  - prior sets for Horseshoe (priors_hs) and Gaussian (priors_gs)
#  - Stan additions for a learned global Gaussian scale (stanvars_shared)
################################################################################

library(brms)

# ==============================================================================
# Stan custom functions for Gompertz AFT
# (log-pdf, log-CDF, log-CCDF, RNG)
# ==============================================================================

stan_funs <- "

  // ============================================================
  // BASELINE GOMPERTZ (tau, sigma)
  // ============================================================

  // log density
  real gompertz_base_lpdf(real t, real tau, real sigma) {
  return log(tau) - log(sigma) + (t / sigma) + tau - tau * exp(t / sigma);
  }

  // log CDF
  real gompertz_base_lcdf(real t, real tau, real sigma) {
    // log F(t) = log(1 - exp(log S(t))) = log1m_exp(log S(t))
    real logS = tau * (1 - exp(t / sigma));
    return log1m_exp(logS);
  }

  // log survival = log CCDF
  real gompertz_base_lccdf(real t, real tau, real sigma) {
    return tau * (1 - exp(t / sigma));
  }

  // RNG via inverse-CDF
  real gompertz_base_rng(real tau, real sigma) {
    real u = uniform_rng(0, 1);
    real t = sigma * log(1 - (1 / tau) * log1m(u));
    if (is_nan(t))
      t = positive_infinity();
    return t;
  }


  // ============================================================
  // GOMPERTZ AFT MODEL (tau, sigma, linpred)
  // S(t|linpred) = S0(t * exp(linpred))
  // ============================================================

  real gompertz_aft_lpdf(real t, real tau, real sigma, real linpred) {
    if (t <= 0) return negative_infinity();
    real tt = t * exp(linpred);
    // baseline density under a transformation of variables
    return gompertz_base_lpdf(tt | tau, sigma) + linpred;
  }

  real gompertz_aft_lcdf(real t, real tau, real sigma, real linpred) {
    if (t <= 0) return negative_infinity();
    real tt = t * exp(linpred);
    return gompertz_base_lcdf(tt | tau, sigma);
  }

  real gompertz_aft_lccdf(real t, real tau, real sigma, real linpred) {
    if (t <= 0) return 0;
    real tt = t * exp(linpred);
    return gompertz_base_lccdf(tt | tau, sigma);
  }

  real gompertz_aft_rng(real tau, real sigma, real linpred) {
    real t0 = gompertz_base_rng(tau, sigma);
    return t0 / exp(linpred);
  }

"


stanvars <- stanvar(scode = stan_funs, block = "functions")


# ==============================================================================
# R-level RNG functions ( for pp_check() )
# ==============================================================================

# Stan's log1m(u) = log(1 - u) = R's log1p(-u)
gompertz_base_rng <- function(tau, sigma) {
  u <- runif(length(tau))
  t <- sigma * log(1 - (1 / tau) * log1p(-u))
  t[is.nan(t) | !is.finite(t)] <- Inf
  t
}


# AFT scaling: T = T0 / exp(linpred)
gompertz_aft_rng <- function(tau, sigma, linpred) {
  t0 <- gompertz_base_rng(tau, sigma)
  t0 / exp(linpred)
}


# ==============================================================================
# Define custom family in brms for Gompertz AFT
# ==============================================================================

gompertz_aft <- custom_family(
  name  = "gompertz_aft",
  # mu = tau (shape), gamma = sigma (scale)
  dpars = c("mu", "gamma", "linpred"), # Error: All families must have a 'mu' parameter.
  links = c(mu = "log", gamma = "log", linpred = "identity"),
  #lb    = c(0, 0, NA), 
  type  = "real",

  # log-likelihood contribution for one observation i
  log_lik = function(i, prep) {

    tau     <- get_dpar(prep, "mu", i = i)    # local variable renamed (brms doesnt know tau)
    sigma   <- get_dpar(prep, "gamma", i = i)
    linpred <- get_dpar(prep, "linpred", i = i)

    t    <- prep$data$Y[i]
    cens <- prep$data$cens[i]
    lb   <- prep$data$lb[i]   # entry_age, set by trunc(lb = entry_age)
    
    if (is.na(cens)) stop("cens is NA at i = ", i)
    
    log_trunc <- gompertz_aft_lccdf(lb, tau, sigma, linpred)  # log S(entry_age)

    if (cens == 0) {
      x <- gompertz_aft_lpdf(t, tau, sigma, linpred) - log_trunc
    } else if (cens == 1) {
      x <- gompertz_aft_lccdf(t, tau, sigma, linpred) - log_trunc
    } else {
      stop("cens must be 0 or 1. Found: ", cens, " at i = ", i)
    }

    return(x)
  },

  # generate a random draw of the response for observation i from the model
  posterior_predict = function(i, prep, ...) {

    tau     <- get_dpar(prep, "mu", i = i)
    sigma   <- get_dpar(prep, "gamma", i = i)
    linpred <- get_dpar(prep, "linpred", i = i)
    lb      <- prep$data$lb[i]

    return(gompertz_aft_rng(tau, sigma, linpred))
  }
)

# ==============================================================================
# Priors
# ==============================================================================

# (1) (Regularised) Horseshoe prior (Piironen & Vehtari, 2017)

priors_hs <- function(p, non_zero_groups) {
  
  # builds the hs prior adaptaviley: 
  # scale_global is set via par_ratio - as advised by Piironen & Vehtari (2017),
  # where par_ratio = p0 / (p - p0) and scale_global = par_ratio/ sqrt(N)
  # with p0 = expected number of non-zero coefficients
  # and p-p0 =  the expected number of zero coefficients
  
  cat(glue::glue("[hs prior] par_ratio = {non_zero_groups} (p = {p})\n"))
  
  c(
    # Gompertz parameters - informative priors
    # --- shape / tau ---
    prior(normal(-6.535, 0.6),  class = "Intercept", dpar = "mu"),
    
    # --- scale / sigma ---
    prior(normal(2.465,  0.03), class = "Intercept", dpar = "gamma"),
    
    # --- linpred ---
    set_prior(
      paste0("horseshoe(df = 3, df_global = 1, par_ratio = ", non_zero_groups,
             ", scale_slab = 1, df_slab = 4)"),
      class = "b", dpar = "linpred"
    )
  )
}



# (2) Gaussian prior (Pavlou et al., 2015; Van Wieringen, 2023)

priors_gs <- c(
 
  # Gompertz parameters - informative priors
  # --- shape / tau ---
  prior(normal(-6.535, 0.6), class = "Intercept", dpar = "mu"), 
  
  # --- scale / sigma ---
  prior(normal(2.465, 0.03), class = "Intercept", dpar = "gamma") 
)


# ---------------------------------------
# data-driven sigma_beta^2 ~ half-N(0, s)
# ---------------------------------------

stan_shared_sigma <- "
  real<lower=0> sigma_beta;
"

stan_shared_sigma_prior <- "
  sigma_beta ~ normal(0, 1);  // half-normal due to lower=0 
"

stan_override_betas <- "
  // Replace the default prior on AFT slopes with shared-sigma
  b_linpred ~ normal(0, sigma_beta);
"

stanvars_shared <- c(
  stanvars,
  stanvar(scode = stan_shared_sigma,       block = "parameters"),
  stanvar(scode = stan_shared_sigma_prior, block = "model"),
  stanvar(scode = stan_override_betas,     block = "model")
)


# helper function to create a shared-sigma stanvars object with a different slope-vector name
make_stanvars_shared <- function(formula, data,
                                 family,  
                                 dpar = "linpred",
                                 include_intercept = FALSE,
                                 hyperprior = "sigma_beta ~ normal(0, 1);",
                                 base_stanvars = stanvars) {
  
  stopifnot(is.character(dpar), length(dpar) == 1)
  stopifnot(is.logical(include_intercept), length(include_intercept) == 1)
  
  pr <- brms::get_prior(formula = formula, data = data, family = family) 
  
  # get population-level (class == "b") coefficients for the chosen dpar
  coefs <- pr$coef[pr$class == "b" & pr$dpar == dpar]
  coefs <- unique(coefs)
  coefs <- coefs[!is.na(coefs) & nzchar(coefs)]
  
  if (!include_intercept) {
    coefs <- setdiff(coefs, "Intercept")
  }
  
  if (length(coefs) == 0) {
    stop("No coefficients found for dpar = '", dpar, "'. Check your formula/dpar name.")
  }
  
  
  # hierarchical shrinkage with shared, data-estimated sigma_beta
  override <- paste0(
    "\n  // shared-sigma prior on all slopes of ", dpar, "\n  ",
    "target += normal_lpdf(b_", dpar, " | 0, sigma_beta);\n"
  )
  
  
  c(
    base_stanvars,
    stanvar(scode = stan_shared_sigma, block = "parameters"),
    stanvar(scode = paste0("\n  ", hyperprior, "\n"), block = "model"),
    stanvar(scode = override, block = "model")
  )
}
