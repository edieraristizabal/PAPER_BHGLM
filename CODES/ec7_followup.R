################################################################
## EC7 follow-up
## (a) Leave-out-five sensitivity, with Moran's I and residual
##     diagnostics computed on the 521 RETAINED catchments only.
##     (In ec7_response_distribution.R the held-out units were still
##     entering the residual vector, which contaminated the statistic.)
## (b) Fig. S5: posterior predictive checks.
## Run from the PAPER_BHGLM root.
################################################################

suppressMessages({
  library(sf); library(spdep); library(INLA); library(dplyr); library(Matrix)
})
set.seed(20260922)
INLA::inla.setOption(num.threads = "4:1")
OUT <- "RESULTS/EC7"

aoi <- st_read("DATA/df_catchments_kmeans.gpkg", quiet = TRUE)
aoi <- aoi %>% mutate_at(c('elev_mean', 'slope_mean', 'RainfallDaysmean'),
                         ~ (scale(.) %>% as.vector))
aoi <- aoi %>% mutate(cuenca_num = case_when(cuenca == 'Atrato' ~ 1,
                                             cuenca == 'Cauca' ~ 2,
                                             cuenca == 'Magdalena' ~ 3))
N <- nrow(aoi)
aoi.nb  <- poly2nb(aoi)
aoi.mat <- as(nb2mat(aoi.nb, style = "B"), "Matrix")
d <- as.data.frame(aoi); d$id <- 1:N
y_obs <- d$lands_rec
Lap     <- Diagonal(nrow(aoi.mat), apply(aoi.mat, 1, sum)) - aoi.mat
Cmatrix <- Diagonal(N, 1) - Lap
CC <- list(dic = TRUE, waic = TRUE, config = TRUE)

f_M5 <- lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean +
  f(cuenca_num, model = "iid") + f(id, model = "generic1", Cmatrix = Cmatrix)

## ---------------------------------------------------------------
## (a) Leave-out-five, diagnostics on retained units only
## ---------------------------------------------------------------
top5 <- order(y_obs, decreasing = TRUE)[1:5]
keep <- setdiff(seq_len(N), top5)

## Subset the neighbour structure to the retained catchments so that
## Moran's I is computed on a graph that matches the residual vector.
nb_keep <- subset(aoi.nb, seq_len(N) %in% keep)
lw_keep <- nb2listw(nb_keep, zero.policy = TRUE)

d_lo <- d; d_lo$lands_rec[top5] <- NA   # drop from likelihood, keep in the graph

fit_one <- function(data, family) {
  inla(f_M5, family = family, offset = log(data$area), data = data,
       control.predictor = list(compute = TRUE), control.compute = CC,
       control.inla = list(strategy = "adaptive"))
}

diag_keep <- function(fit, label, family) {
  mu <- fit$summary.fitted.values$mean[1:N]
  r  <- (y_obs - mu) / sqrt(mu)
  mt <- moran.test(r[keep], listw = lw_keep, alternative = "two.sided",
                   zero.policy = TRUE)
  h  <- fit$summary.hyperpar
  pr <- grep("Precision for id", rownames(h)); be <- grep("Beta for id", rownames(h))
  fx <- fit$summary.fixed
  data.frame(model = label, family = family,
             DIC = fit$dic$dic, WAIC = fit$waic$waic,
             MI_retained = unname(mt$estimate[1]), MI_p = mt$p.value,
             max_abs_pearson_retained = max(abs(r[keep])),
             variance_id = if (length(pr)) 1 / h[pr[1], "mean"] else NA,
             rho         = if (length(be)) h[be[1], "mean"] else NA,
             b_rain = fx["RainfallDaysmean", "mean"],
             b_rain_q025 = fx["RainfallDaysmean", "0.025quant"],
             b_rain_q975 = fx["RainfallDaysmean", "0.975quant"],
             b_elev = fx["elev_mean", "mean"],
             b_slope = fx["slope_mean", "mean"],
             ## out-of-sample prediction for the five excluded catchments
             pred_top5_mean = if (label == "M5 full, Poisson") NA else
               paste(round(mu[top5], 1), collapse = "; "),
             row.names = NULL)
}

fits <- list("M5 full, Poisson"       = fit_one(d,    "poisson"),
             "M5 full, neg. binomial" = fit_one(d,    "nbinomial"),
             "M5 leave-out-5, Poisson"       = fit_one(d_lo, "poisson"),
             "M5 leave-out-5, neg. binomial" = fit_one(d_lo, "nbinomial"))
fams <- c("poisson", "nbinomial", "poisson", "nbinomial")
sens <- do.call(rbind, Map(function(f, l, fm) diag_keep(f, l, fm),
                           fits, names(fits), fams))
write.csv(sens, file.path(OUT, "ec7d_sensitivity_top5_corrected.csv"), row.names = FALSE)
cat("\n[EC7d] Leave-out-five, diagnostics on the 521 retained catchments\n")
print(sens[, c("model", "DIC", "MI_retained", "variance_id", "rho",
               "b_rain", "b_rain_q025", "b_rain_q975", "b_elev", "b_slope")],
      row.names = FALSE)
cat("\nObserved counts of the five excluded catchments:",
    paste(y_obs[top5], collapse = ", "), "\n")

## ---------------------------------------------------------------
## (b) Fig. S5 -- posterior predictive checks
## ---------------------------------------------------------------
ppc <- read.csv(file.path(OUT, "ec7c_posterior_predictive_checks.csv"))
ppc$label <- factor(ppc$model,
  levels = c("M1_poisson", "M4_poisson", "M5_poisson", "M1_nb", "M5_nb", "M5_zinb"),
  labels = c("M1 Poisson", "M4 BYM Poisson", "M5 Leroux Poisson",
             "M1 neg. binomial", "M5 Leroux neg. bin.", "M5 Leroux ZINB"))
stat_lab <- c(total = "Total count", n_zero = "Zero-count catchments",
              max = "Maximum count", variance = "Count variance")
ppc$stat <- factor(stat_lab[ppc$statistic], levels = stat_lab)

## Perceptually uniform, colour-vision-safe (Okabe-Ito), per EC10a
col_ok   <- "#0072B2"   # blue: observed value inside the predictive interval
col_fail <- "#D55E00"   # vermillion: observed value outside it
ppc$inside <- ppc$observed >= ppc$ppc_q025 & ppc$observed <= ppc$ppc_q975

png("FIGURES/FigA4_posterior_predictive_checks.png",
    width = 2400, height = 1800, res = 260)
op <- par(mfrow = c(2, 2), mar = c(4.2, 10.5, 3.0, 1.4), cex.axis = 0.82,
          cex.main = 1.0, font.main = 1, las = 1)
for (s in levels(ppc$stat)) {
  sub <- ppc[ppc$stat == s, ]
  sub <- sub[order(sub$label, decreasing = TRUE), ]
  k <- nrow(sub)
  xr <- range(c(sub$ppc_q025, sub$ppc_q975, sub$observed))
  logx <- s %in% c("Maximum count", "Count variance")
  if (logx) xr <- pmax(xr, 1)
  plot(NA, xlim = xr, ylim = c(0.5, k + 0.5), yaxt = "n", xlab = "", ylab = "",
       main = s, log = if (logx) "x" else "")
  axis(2, at = 1:k, labels = sub$label, tick = FALSE, line = -0.4)
  abline(v = sub$observed[1], col = "grey55", lty = 2, lwd = 1.4)
  for (i in 1:k) {
    cl <- if (sub$inside[i]) col_ok else col_fail
    segments(sub$ppc_q025[i], i, sub$ppc_q975[i], i, col = cl, lwd = 4, lend = 1)
    points(sub$ppc_mean[i], i, pch = 21, bg = "white", col = cl, cex = 1.1, lwd = 2)
  }
  points(rep(sub$observed[1], 1), k + 0.42, pch = 25, bg = "grey25",
         col = "grey25", cex = 0.9, xpd = NA)
  mtext(sprintf("observed = %s", format(sub$observed[1], big.mark = ",")),
        side = 3, line = 0.05, adj = 1, cex = 0.66, col = "grey35")
}
par(op)
dev.off()
cat("\nFig. A4 written to FIGURES/FigA4_posterior_predictive_checks.png\n")
