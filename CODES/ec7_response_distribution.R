################################################################
## EC7 -- Response distribution and model diagnostics
## Editorial comment 7 (Schlogl): overdispersion diagnostics,
## posterior predictive checks, zero inflation, alternative count
## models, sensitivity to influential high-count catchments, and
## the decisive test: do the spatial random effects shrink once the
## likelihood accommodates extra-Poisson variability?
##
## Setup replicates CODES/INLA.R exactly (Queen contiguity, same
## standardised predictors, same Leroux Cmatrix, log(area) offset).
## Run from the PAPER_BHGLM root.
################################################################

suppressMessages({
  library(sf); library(spdep); library(INLA); library(dplyr); library(Matrix)
})
set.seed(20260922)
INLA::inla.setOption(num.threads = "4:1")

OUT <- "RESULTS/EC7"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf(...), "\n", sep = "")

## ---------------------------------------------------------------
## 1. Data, exactly as in INLA.R
## ---------------------------------------------------------------
aoi <- st_read("DATA/df_catchments_kmeans.gpkg", quiet = TRUE)
aoi <- aoi %>% mutate_at(c('elev_mean', 'slope_mean', 'RainfallDaysmean'),
                         ~ (scale(.) %>% as.vector))
aoi <- aoi %>% mutate(cuenca_num = case_when(cuenca == 'Atrato' ~ 1,
                                             cuenca == 'Cauca' ~ 2,
                                             cuenca == 'Magdalena' ~ 3))
N <- nrow(aoi)

aoi.nb    <- poly2nb(aoi)
aoi.mat   <- as(nb2mat(aoi.nb, style = "B"), "Matrix")
aoi.listw <- nb2listw(aoi.nb)          # style "W" (row-standardised), as in INLA.R

d <- as.data.frame(aoi)
## Graph index must address rows of aoi.mat. INLA.R used the file's `id`;
## verify it is 1..N in row order, otherwise use the row index.
if (!identical(as.integer(d$id), 1:N)) {
  log_line("NOTE: file `id` is not 1..N in row order; using row index for the graph.")
  d$id <- 1:N
}
y_obs <- d$lands_rec

## Leroux structure matrix, exactly as in INLA.R
Lap     <- Diagonal(nrow(aoi.mat), apply(aoi.mat, 1, sum)) - aoi.mat
Cmatrix <- Diagonal(N, 1) - Lap

## ---------------------------------------------------------------
## 2. Marginal overdispersion diagnostics (EC7a)
## ---------------------------------------------------------------
od <- data.frame(
  n_catchments  = N,
  total_count   = sum(y_obs),
  mean_count    = mean(y_obs),
  var_count     = var(y_obs),
  var_mean_ratio= var(y_obs) / mean(y_obs),
  median_count  = median(y_obs),
  max_count     = max(y_obs),
  n_zero_obs    = sum(y_obs == 0),
  pct_zero_obs  = 100 * mean(y_obs == 0),
  ## zeros expected under a Poisson with the observed marginal mean
  n_zero_pois_marginal = N * dpois(0, mean(y_obs))
)
write.csv(od, file.path(OUT, "ec7a_overdispersion.csv"), row.names = FALSE)
log_line("[EC7a] var/mean = %.1f | zeros obs = %d (%.1f%%) | zeros expected under Poisson = %.2g",
         od$var_mean_ratio, od$n_zero_obs, od$pct_zero_obs, od$n_zero_pois_marginal)

## ---------------------------------------------------------------
## 3. Model fitting helper
## ---------------------------------------------------------------
CC <- list(dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = FALSE)

fit_model <- function(form, family, data = d, ctrl.family = list()) {
  inla(form, family = family, offset = log(data$area), data = data,
       control.predictor = list(compute = TRUE),
       control.compute   = CC,
       control.family    = ctrl.family,
       control.inla      = list(strategy = "adaptive"))
}

f_M1 <- lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean
f_M2 <- update(f_M1, . ~ . + f(cuenca_num, model = "iid"))
f_M3 <- update(f_M2, . ~ . + f(id, model = "besag",    graph = aoi.mat))
f_M4 <- update(f_M2, . ~ . + f(id, model = "bym",      graph = aoi.mat))
f_M5 <- update(f_M2, . ~ . + f(id, model = "generic1", Cmatrix = Cmatrix))

## Pearson residuals and residual Moran's I, defined once for all models.
## Residual: (y - fitted)/sqrt(fitted); weights: row-standardised Queen.
resid_moran <- function(fit) {
  mu <- fit$summary.fitted.values$mean[1:N]
  r  <- (y_obs - mu) / sqrt(mu)
  mt <- moran.test(r, listw = aoi.listw, randomisation = TRUE, alternative = "two.sided")
  list(resid = r,
       MI    = unname(mt$estimate[1]),
       p     = mt$p.value,
       max_abs_resid = max(abs(r)))
}

summarise_fit <- function(fit, label, family) {
  rm <- resid_moran(fit)
  data.frame(model = label, family = family,
             DIC   = fit$dic$dic,
             WAIC  = fit$waic$waic,
             p_eff = fit$dic$p.eff,
             MI_res = rm$MI, MI_p = rm$p,
             max_abs_pearson = rm$max_abs_resid,
             mlik = fit$mlik[1, 1])
}

## ---------------------------------------------------------------
## 4. Baseline Poisson progression M1-M5 (reproduction; also settles EC9)
## ---------------------------------------------------------------
log_line("\n[4] Fitting Poisson progression M1-M5 ...")
pois <- list(
  M1 = fit_model(f_M1, "poisson"),
  M2 = fit_model(f_M2, "poisson"),
  M3 = fit_model(f_M3, "poisson"),
  M4 = fit_model(f_M4, "poisson"),
  M5 = fit_model(f_M5, "poisson")
)
tab_pois <- do.call(rbind, Map(function(f, l) summarise_fit(f, l, "poisson"),
                               pois, names(pois)))
write.csv(tab_pois, file.path(OUT, "ec9_poisson_progression.csv"), row.names = FALSE)
print(tab_pois, row.names = FALSE)

## ---------------------------------------------------------------
## 5. Alternative count models (EC7b)
## ---------------------------------------------------------------
log_line("\n[5] Fitting negative-binomial, zero-inflated and hurdle variants ...")
alt <- list()
alt[["M1_nb"]]     <- fit_model(f_M1, "nbinomial")
alt[["M5_nb"]]     <- fit_model(f_M5, "nbinomial")
alt[["M4_nb"]]     <- fit_model(f_M4, "nbinomial")
## zero-inflated type 1: p*1{y=0} + (1-p)*f(y)
alt[["M1_zip"]]    <- fit_model(f_M1, "zeroinflatedpoisson1")
alt[["M5_zip"]]    <- fit_model(f_M5, "zeroinflatedpoisson1")
alt[["M1_zinb"]]   <- fit_model(f_M1, "zeroinflatednbinomial1")
alt[["M5_zinb"]]   <- fit_model(f_M5, "zeroinflatednbinomial1")
## type 0 = hurdle (zero-truncated count for y>0)
alt[["M1_hurdle"]] <- fit_model(f_M1, "zeroinflatedpoisson0")
alt[["M5_hurdle"]] <- fit_model(f_M5, "zeroinflatednbinomial0")

fam_of <- c(M1_nb = "nbinomial", M5_nb = "nbinomial", M4_nb = "nbinomial",
            M1_zip = "ZIP", M5_zip = "ZIP", M1_zinb = "ZINB", M5_zinb = "ZINB",
            M1_hurdle = "hurdle-P", M5_hurdle = "hurdle-NB")
tab_alt <- do.call(rbind, Map(function(f, l) summarise_fit(f, l, fam_of[[l]]),
                              alt, names(alt)))
tab_all <- rbind(tab_pois, tab_alt)
write.csv(tab_all, file.path(OUT, "ec7b_model_comparison.csv"), row.names = FALSE)
print(tab_alt, row.names = FALSE)

## Dispersion parameter and zero-inflation probability where applicable
hyp_report <- function(fit, label) {
  h <- fit$summary.hyperpar
  if (is.null(h) || nrow(h) == 0) return(NULL)
  data.frame(model = label, parameter = rownames(h),
             mean = h[, "mean"], sd = h[, "sd"],
             q025 = h[, "0.025quant"], q975 = h[, "0.975quant"],
             row.names = NULL)
}
hyp <- do.call(rbind, c(Map(hyp_report, pois, names(pois)),
                        Map(hyp_report, alt,  names(alt))))
write.csv(hyp, file.path(OUT, "ec7b_hyperparameters.csv"), row.names = FALSE)

## ---------------------------------------------------------------
## 6. Posterior predictive checks (EC7c)
## ---------------------------------------------------------------
log_line("\n[6] Posterior predictive checks ...")
NSAMP <- 500

ppc_stats <- function(fit, label, family, nsamp = NSAMP) {
  s <- inla.posterior.sample(nsamp, fit, verbose = FALSE)
  pred_rows <- grep("^Predictor", rownames(s[[1]]$latent))[1:N]
  reps <- matrix(NA_real_, nrow = nsamp, ncol = 4)
  for (k in seq_len(nsamp)) {
    eta <- s[[k]]$latent[pred_rows]         # includes the log(area) offset
    mu  <- exp(eta)
    hp  <- s[[k]]$hyperpar
    size_nm <- grep("size for", names(hp), value = TRUE)
    zprob_nm <- grep("zero-probability|zero probability", names(hp), value = TRUE)
    yrep <- if (length(size_nm)) rnbinom(N, size = hp[[size_nm[1]]], mu = mu) else rpois(N, mu)
    if (length(zprob_nm)) {
      p0 <- hp[[zprob_nm[1]]]
      yrep[runif(N) < p0] <- 0
    }
    reps[k, ] <- c(sum(yrep), sum(yrep == 0), max(yrep), var(yrep))
  }
  colnames(reps) <- c("total", "n_zero", "max", "variance")
  obs <- c(total = sum(y_obs), n_zero = sum(y_obs == 0),
           max = max(y_obs), variance = var(y_obs))
  data.frame(
    model = label, family = family, statistic = names(obs), observed = as.numeric(obs),
    ppc_mean   = apply(reps, 2, mean),
    ppc_q025   = apply(reps, 2, quantile, 0.025),
    ppc_q975   = apply(reps, 2, quantile, 0.975),
    ## Bayesian p-value: P(T(y_rep) >= T(y_obs)); ~0 or ~1 indicates misfit
    bayes_p    = sapply(seq_along(obs), function(j) mean(reps[, j] >= obs[j])),
    row.names = NULL
  )
}

ppc_targets <- list(M1_poisson = list(pois$M1, "poisson"),
                    M5_poisson = list(pois$M5, "poisson"),
                    M4_poisson = list(pois$M4, "poisson"),
                    M1_nb      = list(alt$M1_nb, "nbinomial"),
                    M5_nb      = list(alt$M5_nb, "nbinomial"),
                    M5_zinb    = list(alt$M5_zinb, "ZINB"))
ppc_tab <- do.call(rbind, Map(function(t, l) ppc_stats(t[[1]], l, t[[2]]),
                              ppc_targets, names(ppc_targets)))
write.csv(ppc_tab, file.path(OUT, "ec7c_posterior_predictive_checks.csv"), row.names = FALSE)
print(ppc_tab, row.names = FALSE)

## ---------------------------------------------------------------
## 7. Sensitivity to influential high-count catchments (EC7d)
## ---------------------------------------------------------------
log_line("\n[7] Leave-out-five sensitivity ...")
top5 <- order(y_obs, decreasing = TRUE)[1:5]
log_line("    excluded counts: %s", paste(y_obs[top5], collapse = ", "))
d_lo <- d
d_lo$lands_rec[top5] <- NA        # drop from likelihood, keep in the graph
fit_lo <- function(form, family) {
  inla(form, family = family, offset = log(d_lo$area), data = d_lo,
       control.predictor = list(compute = TRUE), control.compute = CC,
       control.inla = list(strategy = "adaptive"))
}
lo <- list(M5_poisson_lo5 = fit_lo(f_M5, "poisson"),
           M5_nb_lo5      = fit_lo(f_M5, "nbinomial"))

fixed_of <- function(fit, label) {
  fx <- fit$summary.fixed
  data.frame(model = label, term = rownames(fx), mean = fx[, "mean"], sd = fx[, "sd"],
             q025 = fx[, "0.025quant"], q975 = fx[, "0.975quant"], row.names = NULL)
}
fixed_tab <- do.call(rbind, c(
  Map(fixed_of, pois, names(pois)),
  Map(fixed_of, alt,  names(alt)),
  Map(fixed_of, lo,   names(lo))))
write.csv(fixed_tab, file.path(OUT, "ec7_fixed_effects_all.csv"), row.names = FALSE)

sens <- do.call(rbind, Map(function(f, l) summarise_fit(f, l, "leave-out-5"),
                           lo, names(lo)))
write.csv(sens, file.path(OUT, "ec7d_sensitivity_top5.csv"), row.names = FALSE)
print(sens, row.names = FALSE)

## ---------------------------------------------------------------
## 8. THE DECISIVE TEST (EC7e)
##    Does the structured spatial component shrink once the
##    likelihood accommodates extra-Poisson variability?
##    generic1: hyperpars are Precision and Beta (= rho, the mixing
##    parameter of the Leroux specification).
## ---------------------------------------------------------------
log_line("\n[8] Decisive test: Leroux latent field under Poisson vs negative binomial ...")
leroux_field <- function(fit, label) {
  h <- fit$summary.hyperpar
  prec_row <- grep("Precision for id", rownames(h))
  beta_row <- grep("Beta for id",      rownames(h))
  data.frame(
    model = label,
    precision_id = if (length(prec_row)) h[prec_row[1], "mean"] else NA,
    variance_id  = if (length(prec_row)) 1 / h[prec_row[1], "mean"] else NA,
    rho          = if (length(beta_row)) h[beta_row[1], "mean"] else NA,
    rho_q025     = if (length(beta_row)) h[beta_row[1], "0.025quant"] else NA,
    rho_q975     = if (length(beta_row)) h[beta_row[1], "0.975quant"] else NA,
    sd_latent_field = sd(fit$summary.random$id$mean),
    DIC = fit$dic$dic, WAIC = fit$waic$waic
  )
}
decisive <- rbind(leroux_field(pois$M5, "M5 Leroux, Poisson"),
                  leroux_field(alt$M5_nb, "M5 Leroux, negative binomial"),
                  leroux_field(alt$M5_zinb, "M5 Leroux, zero-inflated NB"))
write.csv(decisive, file.path(OUT, "ec7e_decisive_test.csv"), row.names = FALSE)
print(decisive, row.names = FALSE)

saveRDS(list(pois = lapply(pois, function(f) f$summary.hyperpar),
             alt  = lapply(alt,  function(f) f$summary.hyperpar)),
        file.path(OUT, "ec7_hyperpar_raw.rds"))

log_line("\nDone. Outputs written to %s", OUT)
