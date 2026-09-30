################################################################
## Second editor round, additional analyses on the full data set
##   - non-spatial overdispersed baseline (Poisson + catchment iid)
##   - BYM2 with PC priors (phi = share of marginal latent variance
##     that is spatially structured) and PC-prior sensitivity of M5
##   - basin level: iid (default), fixed effect, or removed
##   - lithology and land cover as factors in M4/M5 (Reviewer 1)
##   - non-linear (rw2) slope and elevation in M4/M5 (Reviewer 1)
##   - negative control for the EC4 latent-field comparison
##   - two-sided Monte Carlo Moran's I, pD, p_WAIC, CPO log score
##   - M1 linearity tests corrected for overdispersion
##   - PCA of the candidate predictors (Appendix A2)
## Setup identical to CODES/INLA.R and CODES/ec7_response_distribution.R.
## Run from the PAPER_BHGLM root. Outputs: RESULTS/LB/
################################################################

suppressMessages({
  library(sf); library(spdep); library(INLA); library(dplyr); library(Matrix)
})
set.seed(20260930)
INLA::inla.setOption(num.threads = "4:1")
OUT <- "RESULTS/LB"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

## ---------------------------------------------------------------
## Data
## ---------------------------------------------------------------
aoi <- st_read("DATA/df_catchments_kmeans.gpkg", quiet = TRUE)
aoi <- aoi %>% mutate_at(c("elev_mean", "slope_mean", "RainfallDaysmean"), ~ (scale(.) %>% as.vector))
aoi <- aoi %>% mutate(cuenca_num = case_when(cuenca == "Atrato" ~ 1, cuenca == "Cauca" ~ 2, cuenca == "Magdalena" ~ 3))
N <- nrow(aoi)
aoi.nb    <- poly2nb(aoi)
aoi.mat   <- as(nb2mat(aoi.nb, style = "B"), "Matrix")
aoi.listw <- nb2listw(aoi.nb)
comp <- n.comp.nb(aoi.nb)
write.csv(data.frame(component = as.integer(names(table(comp$comp.id))),
                     n_catchments = as.integer(table(comp$comp.id))),
          file.path(OUT, "graph_components.csv"), row.names = FALSE)
d <- as.data.frame(aoi)
if (!identical(as.integer(d$id), 1:N)) d$id <- 1:N
d$geo <- factor(d$geomedian, levels = c("granitic", "sediment", "volcanic"))
d$lc  <- factor(d$landcovermedian, levels = c("forest", "grass"))
d$basin_f <- factor(d$cuenca)
d$slope_g <- inla.group(d$slope_mean, n = 20)
d$elev_g  <- inla.group(d$elev_mean, n = 20)
set.seed(1); d$slope_perm <- sample(d$slope_mean)
y_obs <- d$lands_rec
Cmatrix <- Diagonal(N, 1) - (Diagonal(N, apply(aoi.mat, 1, sum)) - aoi.mat)

## ---------------------------------------------------------------
## Model specifications
## ---------------------------------------------------------------
pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
fx   <- "lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean"
fxr2 <- "lands_rec ~ 1 + RainfallDaysmean + f(elev_g, model = 'rw2', scale.model = TRUE) + f(slope_g, model = 'rw2', scale.model = TRUE)"
B    <- " + f(cuenca_num, model = 'iid')"
ICAR <- " + f(id, model = 'besag', graph = aoi.mat)"
BYM  <- " + f(id, model = 'bym', graph = aoi.mat)"
LER  <- " + f(id, model = 'generic1', Cmatrix = Cmatrix)"
LERpc<- " + f(id, model = 'generic1', Cmatrix = Cmatrix, hyper = pc_prec)"
BYM2 <- " + f(id, model = 'bym2', graph = aoi.mat, scale.model = TRUE, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)), phi = list(prior = 'pc', param = c(0.5, 2/3))))"
IID  <- " + f(id, model = 'iid')"
specs <- list(
  M1 = fx, M2 = paste0(fx, B), M2iid = paste0(fx, B, IID),
  M3 = paste0(fx, B, ICAR), M4 = paste0(fx, B, BYM), M5 = paste0(fx, B, LER),
  BYM2 = paste0(fx, B, BYM2), M5_pcprior = paste0(fx, B, LERpc),
  M4_nobasin = paste0(fx, BYM), M5_nobasin = paste0(fx, LER),
  M4_basinfixed = paste0(fx, " + basin_f", BYM), M5_basinfixed = paste0(fx, " + basin_f", LER),
  M4_lith_lc = paste0(fx, " + geo + lc", B, BYM), M5_lith_lc = paste0(fx, " + geo + lc", B, LER),
  M4_rw2 = paste0(fxr2, B, BYM), M5_rw2 = paste0(fxr2, B, LER),
  M5_noslope = paste0("lands_rec ~ 1 + RainfallDaysmean + elev_mean", B, LER),
  M5_permslope = paste0("lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_perm", B, LER))

CC <- list(dic = TRUE, waic = TRUE, cpo = TRUE, config = FALSE, return.marginals.predictor = FALSE)
fit_model <- function(form, threads = "4:1") {
  f <- try(inla(as.formula(form), family = "poisson", offset = log(d$area), data = d,
                control.predictor = list(compute = TRUE), control.compute = CC,
                control.inla = list(strategy = "adaptive"), num.threads = threads))
  if (inherits(f, "try-error") && threads != "1:1") f <- fit_model(form, "1:1")
  f
}
hp <- function(fit, pat, what = "mean") {
  h <- fit$summary.hyperpar; j <- grep(pat, rownames(h)); if (length(j)) h[j[1], what] else NA
}

## ---------------------------------------------------------------
## Fit and summarise
## ---------------------------------------------------------------
fits <- list(); summ <- list(); fixed <- list(); hyper <- list()
set.seed(20260930)
for (k in names(specs)) {
  log_line("fitting %s", k)
  f <- fit_model(specs[[k]]); fits[[k]] <- f
  mu <- f$summary.fitted.values$mean[1:N]
  r  <- (y_obs - mu) / sqrt(mu)
  mc <- moran.mc(r, aoi.listw, nsim = 999, alternative = "two.sided")
  lat <- if (!is.null(f$summary.random$id)) f$summary.random$id$mean[1:N] else rep(NA, N)
  summ[[k]] <- data.frame(
    model = k, DIC = f$dic$dic, pD = f$dic$p.eff, WAIC = f$waic$waic, pWAIC = f$waic$p.eff,
    mlik = f$mlik[1, 1], logscore_CPO = -mean(log(f$cpo$cpo), na.rm = TRUE),
    cpo_failures = sum(f$cpo$failure > 0, na.rm = TRUE),
    MI = unname(mc$statistic), MI_p_two_sided = mc$p.value, max_abs_pearson = max(abs(r)),
    prec_id = hp(f, "Precision for id"), rho_leroux = hp(f, "Beta for id"),
    rho_q025 = hp(f, "Beta for id", "0.025quant"), rho_q975 = hp(f, "Beta for id", "0.975quant"),
    phi_bym2 = hp(f, "Phi for id"), phi_q025 = hp(f, "Phi for id", "0.025quant"), phi_q975 = hp(f, "Phi for id", "0.975quant"),
    prec_basin = hp(f, "cuenca_num"), prec_basin_sd = hp(f, "cuenca_num", "sd"),
    sd_latent_field = sd(lat))
  fx_s <- f$summary.fixed
  fixed[[k]] <- data.frame(model = k, term = rownames(fx_s), mean = fx_s[, "mean"], sd = fx_s[, "sd"],
                           q025 = fx_s[, "0.025quant"], q975 = fx_s[, "0.975quant"], row.names = NULL)
  h <- f$summary.hyperpar
  if (!is.null(h) && nrow(h) > 0) hyper[[k]] <- data.frame(model = k, parameter = rownames(h), mean = h[, "mean"], sd = h[, "sd"],
                                        q025 = h[, "0.025quant"], q975 = h[, "0.975quant"], row.names = NULL)
}
summ  <- do.call(rbind, summ)
lat_of <- function(k) fits[[k]]$summary.random$id$mean[1:N]
summ$r_field_vs_M5 <- sapply(summ$model, function(k) if (grepl("^M5", k)) cor(lat_of(k), lat_of("M5")) else NA)
summ$r_field_vs_M4 <- sapply(summ$model, function(k) if (grepl("^M4", k)) cor(lat_of(k), lat_of("M4")) else NA)
write.csv(summ, file.path(OUT, "lb_model_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, fixed), file.path(OUT, "lb_fixed_effects.csv"), row.names = FALSE)
write.csv(do.call(rbind, hyper), file.path(OUT, "lb_hyperparameters.csv"), row.names = FALSE)

## Latent field of M5 by lithology and land-cover class
f5 <- lat_of("M5")
lm_geo <- lm(f5 ~ d$geo); lm_lc <- lm(f5 ~ d$lc)
lf <- rbind(
  data.frame(variable = "lithology", class = levels(d$geo), n = as.integer(table(d$geo)),
             mean_field = as.numeric(tapply(f5, d$geo, mean)), R2 = summary(lm_geo)$r.squared),
  data.frame(variable = "land cover", class = levels(d$lc), n = as.integer(table(d$lc)),
             mean_field = as.numeric(tapply(f5, d$lc, mean)), R2 = summary(lm_lc)$r.squared))
write.csv(lf, file.path(OUT, "lb_latent_field_by_class.csv"), row.names = FALSE)

## Non-linear effects in M4/M5 (for plotting / reporting)
for (k in c("M4_rw2", "M5_rw2")) for (v in c("slope_g", "elev_g")) {
  r <- fits[[k]]$summary.random[[v]]
  write.csv(data.frame(model = k, term = v, x = r$ID, mean = r$mean, q025 = r$`0.025quant`, q975 = r$`0.975quant`),
            file.path(OUT, sprintf("lb_rw2_%s_%s.csv", k, v)), row.names = FALSE)
}

## ---------------------------------------------------------------
## M1 linearity, Poisson and quasi-Poisson (overdispersion-corrected)
## ---------------------------------------------------------------
g0 <- glm(lands_rec ~ RainfallDaysmean + elev_mean + slope_mean + offset(log(area)), family = poisson, data = d)
phi <- sum(residuals(g0, "pearson")^2) / g0$df.residual
lin <- do.call(rbind, lapply(c("RainfallDaysmean", "elev_mean", "slope_mean"), function(v) {
  g1 <- update(g0, as.formula(paste(". ~ . + I(", v, "^2)")))
  q0 <- update(g0, family = quasipoisson); q1 <- update(g1, family = quasipoisson)
  a  <- anova(q0, q1, test = "F")
  data.frame(term = v, dAIC_poisson = AIC(g0) - AIC(g1),
             LR = deviance(g0) - deviance(g1), p_LR_poisson = pchisq(deviance(g0) - deviance(g1), 1, lower.tail = FALSE),
             quad_coef = unname(coef(g1)[length(coef(g1))]),
             dispersion = phi, F_quasi = a$F[2], p_F_quasi = a$`Pr(>F)`[2])
}))
write.csv(lin, file.path(OUT, "lb_linearity_M1.csv"), row.names = FALSE)

## ---------------------------------------------------------------
## PCA of the candidate predictors (as in the original screening)
## ---------------------------------------------------------------
cand <- read.csv("DATA/catchments_with_new_predictors.csv")
pv <- c(area = "Area", hypso_inte = "Hypsometric integral", Densidad = "Lineament density",
        rainfallAnnual_mean = "Mean annual rainfall", elev_mean = "Mean elevation (H)",
        slope_mean = "Mean slope (S)", rel_mean = "Mean relief", RainfallDaysmean = "Rainfall days (R_d)")
e <- eigen(cor(cand[, names(pv)]), symmetric = TRUE)   # identical to PCA on standardized variables
V <- e$vectors; for (j in 1:ncol(V)) if (V[which.max(abs(V[, j])), j] < 0) V[, j] <- -V[, j]
write.csv(data.frame(component = paste0("PC", seq_along(e$values)), eigenvalue = e$values,
                     explained = e$values / sum(e$values), cumulative = cumsum(e$values) / sum(e$values)),
          file.path(OUT, "pca_variance.csv"), row.names = FALSE)
ld <- data.frame(variable = unname(pv), V[, 1:3]); names(ld)[2:4] <- c("PC1", "PC2", "PC3")
write.csv(ld, file.path(OUT, "pca_loadings.csv"), row.names = FALSE)
sc <- scale(cand[, names(pv)]) %*% V[, 1:2]
write.csv(data.frame(id = cand$id, PC1 = sc[, 1], PC2 = sc[, 2], lands = cand$lands_rec, basin = cand$cuenca),
          file.path(OUT, "pca_scores.csv"), row.names = FALSE)

options(width = 220)
print(summ[, c("model", "DIC", "pD", "WAIC", "pWAIC", "mlik", "logscore_CPO", "cpo_failures", "MI",
               "MI_p_two_sided", "rho_leroux", "phi_bym2", "prec_basin", "sd_latent_field", "r_field_vs_M5", "r_field_vs_M4")],
      digits = 4, row.names = FALSE)
print(lf); print(lin, digits = 4)
log_line("done")
