################################################################
## EC4 (part 2) -- Sensitivity of the model to the aggregation
## function of the slope predictor.
## The published catchment-mean slope is replaced, one at a time,
## by alternatives computed from the native 12.5 m DEM
## (CODES/ec4_zonal_12m.R). M1, M4 and M5 (Poisson) are refitted
## with everything else identical to CODES/INLA.R and
## CODES/ec7_response_distribution.R (Queen contiguity, standardised
## predictors, Leroux Cmatrix, log(area) offset, same Moran's I).
## Reported: fixed effects, DIC/WAIC, residual Moran's I, BYM
## variance partition, Leroux rho and latent-field variance, and the
## correlation of the M5 latent field with the published one.
## Run from the PAPER_BHGLM root.
################################################################

suppressMessages({
  library(sf); library(spdep); library(INLA); library(dplyr); library(Matrix)
})
set.seed(20260929)
INLA::inla.setOption(num.threads = "4:1")

OUT <- "RESULTS/EC4"
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

## ---------------------------------------------------------------
## 1. Data, as in INLA.R, plus the 12.5 m terrain alternatives
## ---------------------------------------------------------------
aoi <- st_read("DATA/df_catchments_kmeans.gpkg", quiet = TRUE)
ter <- read.csv(file.path(OUT, "ec4_catchment_terrain_12m.csv"))
stopifnot(all(as.numeric(ter$id) == as.numeric(aoi$id)))
alt_vars <- c("slope_mean12", "slope_p90", "frac_gt20", "frac_gt25", "frac_gt30", "slope_sd")
aoi[alt_vars] <- ter[alt_vars]

aoi <- aoi %>% mutate_at(c("elev_mean", "slope_mean", "RainfallDaysmean", alt_vars),
                         ~ (scale(.) %>% as.vector))
aoi <- aoi %>% mutate(cuenca_num = case_when(cuenca == "Atrato" ~ 1,
                                             cuenca == "Cauca" ~ 2,
                                             cuenca == "Magdalena" ~ 3))
N <- nrow(aoi)
aoi.nb    <- poly2nb(aoi)
aoi.mat   <- as(nb2mat(aoi.nb, style = "B"), "Matrix")
aoi.listw <- nb2listw(aoi.nb)
d <- as.data.frame(aoi)
if (!identical(as.integer(d$id), 1:N)) d$id <- 1:N
y_obs <- d$lands_rec

Lap     <- Diagonal(nrow(aoi.mat), apply(aoi.mat, 1, sum)) - aoi.mat
Cmatrix <- Diagonal(N, 1) - Lap

## ---------------------------------------------------------------
## 2. Helpers (identical definitions to ec7_response_distribution.R)
## ---------------------------------------------------------------
CC <- list(dic = TRUE, waic = TRUE, return.marginals.predictor = FALSE)
fit_model <- function(form, threads = "4:1") {
  fit <- try(inla(form, family = "poisson", offset = log(d$area), data = d,
                  control.predictor = list(compute = TRUE), control.compute = CC,
                  control.inla = list(strategy = "adaptive"), num.threads = threads))
  if (inherits(fit, "try-error") && threads != "1:1") {
    log_line("  inla crashed; retrying single-threaded")
    fit <- fit_model(form, threads = "1:1")
  }
  fit
}
resid_moran <- function(fit) {
  mu <- fit$summary.fitted.values$mean[1:N]
  r  <- (y_obs - mu) / sqrt(mu)
  mt <- moran.test(r, listw = aoi.listw, randomisation = TRUE, alternative = "two.sided")
  c(MI = unname(mt$estimate[1]), MI_p = mt$p.value)
}
hp <- function(fit, pat, what = "mean") {
  h <- fit$summary.hyperpar; j <- grep(pat, rownames(h))
  if (length(j)) h[j[1], what] else NA
}

forms_for <- function(slope_terms) {
  fx <- paste("lands_rec ~ 1 + RainfallDaysmean + elev_mean +",
              paste(slope_terms, collapse = " + "))
  re <- "+ f(cuenca_num, model = 'iid')"
  list(
    M1 = as.formula(fx),
    M4 = as.formula(paste(fx, re, "+ f(id, model = 'bym', graph = aoi.mat)")),
    M5 = as.formula(paste(fx, re, "+ f(id, model = 'generic1', Cmatrix = Cmatrix)")))
}

variants <- list(
  "mean (published)"     = "slope_mean",
  "mean (12.5 m)"        = "slope_mean12",
  "p90"                  = "slope_p90",
  "area > 20 deg"        = "frac_gt20",
  "area > 25 deg"        = "frac_gt25",
  "area > 30 deg"        = "frac_gt30",
  "mean + within-SD"     = c("slope_mean", "slope_sd"))

## ---------------------------------------------------------------
## 3. Fit all variants
## ---------------------------------------------------------------
summ <- list(); fixed <- list(); fields <- list()
for (v in names(variants)) {
  fs <- forms_for(variants[[v]])
  for (m in names(fs)) {
    log_line("fitting %s | %s", m, v)
    fit <- fit_model(fs[[m]])
    mi  <- resid_moran(fit)
    summ[[length(summ) + 1]] <- data.frame(
      variant = v, model = m,
      DIC = fit$dic$dic, WAIC = fit$waic$waic,
      MI_res = mi[["MI"]], MI_p = mi[["MI_p"]],
      ## BYM: iid and spatial precisions; Leroux: precision and rho
      prec_iid     = if (m == "M4") hp(fit, "iid component") else NA,
      prec_spatial = if (m == "M4") hp(fit, "spatial component") else NA,
      prec_id      = if (m == "M5") hp(fit, "Precision for id") else NA,
      rho          = if (m == "M5") hp(fit, "Beta for id") else NA,
      rho_q025     = if (m == "M5") hp(fit, "Beta for id", "0.025quant") else NA,
      rho_q975     = if (m == "M5") hp(fit, "Beta for id", "0.975quant") else NA,
      sd_latent_field = if (m != "M1") sd(fit$summary.random$id$mean[1:N]) else NA)
    fx <- fit$summary.fixed
    fixed[[length(fixed) + 1]] <- data.frame(
      variant = v, model = m, term = rownames(fx), mean = fx[, "mean"],
      q025 = fx[, "0.025quant"], q975 = fx[, "0.975quant"], row.names = NULL)
    if (m == "M5") fields[[v]] <- fit$summary.random$id$mean[1:N]
  }
}
summ  <- do.call(rbind, summ);  fixed <- do.call(rbind, fixed)
summ$var_iid     <- 1 / summ$prec_iid
summ$var_spatial <- 1 / summ$prec_spatial
summ$var_id      <- 1 / summ$prec_id
write.csv(summ,  file.path(OUT, "ec4_sensitivity_summary.csv"), row.names = FALSE)
write.csv(fixed, file.path(OUT, "ec4_sensitivity_fixed.csv"),   row.names = FALSE)

## Similarity of the M5 latent field to the published specification
F <- do.call(cbind, fields)
fc <- data.frame(variant = colnames(F),
                 cor_with_published = as.numeric(cor(F)[, "mean (published)"]))
write.csv(fc, file.path(OUT, "ec4_latent_field_similarity.csv"), row.names = FALSE)

options(width = 200)
print(summ[, c("variant", "model", "DIC", "WAIC", "MI_res", "var_iid", "var_spatial",
               "var_id", "rho", "sd_latent_field")], row.names = FALSE, digits = 4)
print(fixed[fixed$term != "(Intercept)", ], row.names = FALSE, digits = 3)
print(fc, row.names = FALSE, digits = 4)
log_line("done")
