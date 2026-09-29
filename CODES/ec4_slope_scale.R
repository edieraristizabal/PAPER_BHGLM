################################################################
## EC4 (part 3) -- Does aggregation to catchments mask the
## slope-scale control on landslide occurrence?
## Editorial comment 4 / Reviewer 1: landslides occur at hillslope
## scale, predictors are catchment aggregates. Three analyses, all on
## the native 12.5 m ALOS-PALSAR DEM (Sect. 3) and the point inventory
## (DATA/inventarioV1.shp; reproduces lands_rec in all 526 catchments):
##
##  (A) How much slope information does aggregation discard?
##      Pixel-level variance of slope split into between- and
##      within-catchment components.
##  (B) Where do landslides sit within their own catchment?
##      Percentile of the crown-pixel slope in the slope distribution
##      of its catchment (uniform on [0,1] if the within-catchment
##      slope distribution were irrelevant).
##  (C) Process-consistent aggregation. The slope-scale response f(s)
##      is estimated from WITHIN-catchment contrasts only:
##        n_cb ~ Poisson(A_cb * exp(alpha_c + f(s_b)))
##      with a fixed effect alpha_c per catchment, so f uses how
##      landslides are distributed over slope classes inside each
##      catchment and not how many each catchment has. It is then
##      integrated over every pixel of each catchment,
##        S_c = log( sum_b A_cb exp(f(s_b)) / A_c ),
##      which is the catchment-level expectation implied by a purely
##      slope-scale control (the correct change of support; by Jensen
##      it differs from f(mean slope)). M1, M4 and M5 are refitted
##      with S_c replacing, or added to, the catchment-mean slope, and
##      with S_c entered as a fixed offset.
## Run from the PAPER_BHGLM root.
################################################################

suppressMessages({
  library(sf); library(terra); library(spdep); library(INLA)
  library(dplyr); library(Matrix); library(splines)
})
sf_use_s2(FALSE)
set.seed(20260929)
terraOptions(progress = 0, memfrac = 0.5)
INLA::inla.setOption(num.threads = "4:1")

DEM_PATH <- "C:/Users/edier/Documents/INVESTIGACION/PAPERS/PUBLICADOS/PAPER_SAR/Data/DEM_filled_12m.tif"
OUT <- "RESULTS/EC4"
TMP <- file.path(tempdir(), "ec4")
dir.create(TMP, recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")
tabulate_w <- function(id, w, n) {
  r <- rowsum(w, id, reorder = FALSE)
  out <- numeric(n); out[as.integer(rownames(r))] <- r[, 1]; out
}

aoi <- st_read("DATA/df_catchments_kmeans.gpkg", quiet = TRUE)
aoi$rid <- seq_len(nrow(aoi))
N <- nrow(aoi)
pts <- st_transform(st_read("DATA/inventarioV1.shp", quiet = TRUE), st_crs(aoi))

## ---------------------------------------------------------------
## 1. Slope raster, catchment raster, per-catchment slope histograms
## ---------------------------------------------------------------
dem <- crop(rast(DEM_PATH), ext(vect(aoi)) + 250, snap = "out",
            filename = file.path(TMP, "dem.tif"), overwrite = TRUE)
slp <- terrain(dem, v = "slope", neighbors = 8, unit = "degrees",
               filename = file.path(TMP, "slope.tif"), overwrite = TRUE)
idr <- rasterize(vect(aoi), dem, field = "rid", datatype = "INT2U",
                 filename = file.path(TMP, "id.tif"), overwrite = TRUE)
log_line("rasters ready")

BW <- 0.05; NB <- ceiling(90 / BW)
hist <- numeric(N * NB); ssum <- numeric(N); sss <- numeric(N); sn <- numeric(N)
readStart(slp); readStart(idr)
nr <- nrow(slp)
for (r0 in seq(1L, nr, by = 1000L)) {
  nrw <- min(1000L, nr - r0 + 1L)
  s  <- readValues(slp, r0, nrw); id <- readValues(idr, r0, nrw)
  ok <- !is.na(id) & !is.na(s); if (!any(ok)) next
  s <- s[ok]; id <- as.integer(id[ok])
  hist <- hist + tabulate((id - 1L) * NB + pmin(as.integer(s %/% BW) + 1L, NB), nbins = N * NB)
  ssum <- ssum + tabulate_w(id, s, N); sss <- sss + tabulate_w(id, s^2, N)
  sn   <- sn + tabulate(id, nbins = N)
}
readStop(slp); readStop(idr)
H <- matrix(hist, nrow = N, ncol = NB, byrow = TRUE)   # cell counts, 0.05 deg bins
log_line("histograms done")

## ---------------------------------------------------------------
## (A) Between- vs within-catchment variance of pixel slope
## ---------------------------------------------------------------
m_c  <- ssum / sn
v_c  <- sss / sn - m_c^2
m_all <- sum(ssum) / sum(sn)
var_between <- sum(sn * (m_c - m_all)^2) / sum(sn)
var_within  <- sum(sn * v_c) / sum(sn)
A <- data.frame(
  n_pixels = sum(sn), mean_slope = m_all,
  var_total = var_between + var_within,
  var_between = var_between, var_within = var_within,
  pct_within = 100 * var_within / (var_between + var_within),
  sd_between_catchment_means = sd(m_c),
  median_within_sd = median(sqrt(v_c)))
write.csv(A, file.path(OUT, "ec4A_variance_decomposition.csv"), row.names = FALSE)
log_line("[A] within-catchment share of pixel slope variance = %.1f%%", A$pct_within)

## ---------------------------------------------------------------
## (B) Slope at landslide crowns relative to their catchment
## ---------------------------------------------------------------
pv  <- vect(pts)
p_s <- terra::extract(slp, pv)[, 2]
p_c <- terra::extract(idr, pv)[, 2]
keep <- !is.na(p_s) & !is.na(p_c)
p_s <- p_s[keep]; p_c <- as.integer(p_c[keep])
stopifnot(sum(keep) > 10000)
cdf <- t(apply(H, 1, function(h) cumsum(h) / sum(h)))
p_bin <- pmin(as.integer(p_s %/% BW) + 1L, NB)
pct <- cdf[cbind(p_c, p_bin)]                    # percentile within own catchment
B <- data.frame(
  n_points = length(pct),
  median_percentile = median(pct),
  mean_percentile = mean(pct),
  frac_above_catchment_median = mean(pct > 0.5),
  frac_above_catchment_p75 = mean(pct > 0.75),
  frac_above_catchment_p90 = mean(pct > 0.90),
  median_slope_at_crowns = median(p_s),
  median_diff_from_catchment_mean = median(p_s - m_c[p_c]),
  ks_uniform_D = unname(suppressWarnings(ks.test(pct, "punif"))$statistic))
write.csv(B, file.path(OUT, "ec4B_crown_percentiles.csv"), row.names = FALSE)
write.csv(data.frame(catchment = p_c, slope = p_s, percentile = pct),
          file.path(OUT, "ec4B_crown_points.csv"), row.names = FALSE)
log_line("[B] median within-catchment percentile of crown slope = %.2f", B$median_percentile)

## ---------------------------------------------------------------
## (C1) Within-catchment slope response f(s), 1-degree classes
## ---------------------------------------------------------------
K <- 1 / BW                                      # 20 fine bins per degree
NB1 <- 70                                        # 0-70 deg; steeper cells pooled
H1 <- sapply(seq_len(NB1), function(b) {
  j <- ((b - 1) * K + 1):(if (b < NB1) b * K else NB)
  rowSums(H[, j, drop = FALSE])
})
p_b1 <- pmin(as.integer(p_s %/% 1) + 1L, NB1)
L1 <- matrix(tabulate((p_c - 1L) * NB1 + p_b1, nbins = N * NB1),
             nrow = N, ncol = NB1, byrow = TRUE)
mid1 <- seq_len(NB1) - 0.5
long <- data.frame(c = rep(seq_len(N), times = NB1),
                   s = rep(mid1, each = N),
                   area = as.vector(H1) * 12.5^2 / 1e6,       # km2
                   n = as.vector(L1))
long <- long[long$area > 0, ]
## only catchments with >= 1 landslide carry within-catchment information
long_fit <- long[long$c %in% which(tabulate(p_c, N) > 0), ]
g <- glm(n ~ factor(c) + ns(s, df = 5), family = poisson,
         offset = log(area), data = long_fit)
disp <- sum(residuals(g, "pearson")^2) / g$df.residual
log_line("[C1] within-catchment glm fitted; Pearson dispersion = %.2f", disp)

## f(s) on a fine grid, centred at the study-area mean slope, with SE
grid <- data.frame(s = seq(0.5, 69.5, by = 0.5))
Xg   <- predict(ns(long_fit$s, df = 5), grid$s)
bi   <- grep("ns\\(s", names(coef(g)))
fg   <- as.vector(Xg %*% coef(g)[bi])
seg  <- sqrt(pmax(rowSums((Xg %*% vcov(g)[bi, bi]) * Xg), 0)) * sqrt(max(disp, 1))
ref  <- approx(grid$s, fg, xout = m_all)$y
## pooled (marginal) frequency ratio for comparison: ignores catchments
fr <- tapply(long$n, cut(long$s, seq(0, 70, 1)), sum) /
      tapply(long$area, cut(long$s, seq(0, 70, 1)), sum)
fr <- fr / (sum(long$n) / sum(long$area))
Cf <- data.frame(slope = grid$s, f = fg - ref, lo = fg - ref - 1.96 * seg,
                 hi = fg - ref + 1.96 * seg,
                 pooled_log_FR = log(approx(mid1, as.numeric(fr), xout = grid$s)$y))
write.csv(Cf, file.path(OUT, "ec4C_slope_response.csv"), row.names = FALSE)

## ---------------------------------------------------------------
## (C2) Process-consistent catchment covariate S_c
## ---------------------------------------------------------------
f1  <- as.vector(predict(ns(long_fit$s, df = 5), mid1) %*% coef(g)[bi]) - ref
S_c <- log(as.vector(H1 %*% exp(f1)) / rowSums(H1))
f_at_mean <- approx(mid1, f1, xout = m_c, rule = 2)$y
C2 <- data.frame(id = aoi$id, slope_mean = aoi$slope_mean, S_c = S_c,
                 f_at_mean_slope = f_at_mean, jensen_gap = S_c - f_at_mean)
write.csv(C2, file.path(OUT, "ec4C_process_covariate.csv"), row.names = FALSE)
log_line("[C2] cor(S_c, slope_mean) = %.3f; median Jensen gap = %.3f",
         cor(S_c, aoi$slope_mean), median(C2$jensen_gap))

## ---------------------------------------------------------------
## (C3) Refit M1, M4, M5 (setup identical to INLA.R / EC7)
## ---------------------------------------------------------------
aoi$S_c_raw <- S_c
aoi$S_c     <- S_c
aoi <- aoi %>% mutate_at(c("elev_mean", "slope_mean", "RainfallDaysmean", "S_c"),
                         ~ (scale(.) %>% as.vector))
aoi <- aoi %>% mutate(cuenca_num = case_when(cuenca == "Atrato" ~ 1,
                                             cuenca == "Cauca" ~ 2,
                                             cuenca == "Magdalena" ~ 3))
aoi.nb <- poly2nb(aoi); aoi.mat <- as(nb2mat(aoi.nb, style = "B"), "Matrix")
aoi.listw <- nb2listw(aoi.nb)
d <- as.data.frame(aoi); if (!identical(as.integer(d$id), 1:N)) d$id <- 1:N
y_obs <- d$lands_rec
Lap <- Diagonal(N, apply(aoi.mat, 1, sum)) - aoi.mat
Cmatrix <- Diagonal(N, 1) - Lap

CC <- list(dic = TRUE, waic = TRUE, return.marginals.predictor = FALSE)
fit_model <- function(form, off, threads = "4:1") {
  fit <- try(inla(form, family = "poisson", offset = off, data = d,
                  control.predictor = list(compute = TRUE), control.compute = CC,
                  control.inla = list(strategy = "adaptive"), num.threads = threads))
  if (inherits(fit, "try-error") && threads != "1:1") fit <- fit_model(form, off, "1:1")
  fit
}
resid_moran <- function(fit) {
  mu <- fit$summary.fitted.values$mean[1:N]; r <- (y_obs - mu) / sqrt(mu)
  unname(moran.test(r, aoi.listw, randomisation = TRUE, alternative = "two.sided")$estimate[1])
}
hp <- function(fit, pat, what = "mean") {
  h <- fit$summary.hyperpar; j <- grep(pat, rownames(h)); if (length(j)) h[j[1], what] else NA
}
forms_for <- function(terms) {
  fx <- paste("lands_rec ~ 1 +", paste(terms, collapse = " + "))
  re <- "+ f(cuenca_num, model = 'iid')"
  list(M1 = as.formula(fx),
       M4 = as.formula(paste(fx, re, "+ f(id, model = 'bym', graph = aoi.mat)")),
       M5 = as.formula(paste(fx, re, "+ f(id, model = 'generic1', Cmatrix = Cmatrix)")))
}
base <- c("RainfallDaysmean", "elev_mean")
variants <- list(
  "mean slope (published)"      = list(terms = c(base, "slope_mean"), off = log(d$area)),
  "S_c (slope-scale integrated)"= list(terms = c(base, "S_c"),        off = log(d$area)),
  "mean slope + S_c"            = list(terms = c(base, "slope_mean", "S_c"), off = log(d$area)),
  "S_c as fixed offset"         = list(terms = base, off = log(d$area) + d$S_c_raw))

summ <- list(); fixed <- list(); fields <- list()
for (v in names(variants)) {
  fs <- forms_for(variants[[v]]$terms)
  for (m in names(fs)) {
    log_line("fitting %s | %s", m, v)
    fit <- fit_model(fs[[m]], variants[[v]]$off)
    summ[[length(summ) + 1]] <- data.frame(
      variant = v, model = m, DIC = fit$dic$dic, WAIC = fit$waic$waic,
      MI_res = resid_moran(fit),
      var_iid = if (m == "M4") 1 / hp(fit, "iid component") else NA,
      var_spatial = if (m == "M4") 1 / hp(fit, "spatial component") else NA,
      var_id = if (m == "M5") 1 / hp(fit, "Precision for id") else NA,
      rho = if (m == "M5") hp(fit, "Beta for id") else NA,
      sd_latent_field = if (m != "M1") sd(fit$summary.random$id$mean[1:N]) else NA)
    fx <- fit$summary.fixed
    fixed[[length(fixed) + 1]] <- data.frame(variant = v, model = m, term = rownames(fx),
      mean = fx[, "mean"], q025 = fx[, "0.025quant"], q975 = fx[, "0.975quant"], row.names = NULL)
    if (m == "M5") fields[[v]] <- fit$summary.random$id$mean[1:N]
  }
}
summ <- do.call(rbind, summ); fixed <- do.call(rbind, fixed)
F <- do.call(cbind, fields)
summ$cor_M5_field_with_published <- NA
summ$cor_M5_field_with_published[summ$model == "M5"] <- cor(F)[, 1]
write.csv(summ,  file.path(OUT, "ec4C_refit_summary.csv"), row.names = FALSE)
write.csv(fixed, file.path(OUT, "ec4C_refit_fixed.csv"),   row.names = FALSE)

options(width = 200)
print(A, digits = 4); print(B, digits = 3)
print(Cf[Cf$slope %in% c(2.5, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60), ], digits = 3, row.names = FALSE)
print(summary(C2[, -1]))
print(summ, digits = 4, row.names = FALSE)
print(fixed[fixed$term != "(Intercept)", ], digits = 3, row.names = FALSE)
log_line("done")
