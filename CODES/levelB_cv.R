################################################################
## Second editor round: cross-validation redone with proper
## predictive scoring.
##   - partitions: three spatially blocked k-means partitions
##     (k = 5, different seeds; one reproduces CODES/spatial_cross_validation.R)
##     and one random 5-fold partition (interpolation task)
##   - models: M1, M2, M2iid, M3, M4, M5, BYM2 (as in levelB_models.R)
##   - scores per held-out catchment, from 300 posterior samples of
##     the linear predictor: log predictive density (log score),
##     CRPS (sample based), squared/absolute error of the predictive
##     median, and the squared error of the posterior mean of lambda
##     (the metric used in the first revision)
##   - records whether a held-out catchment's graph component has no
##     observed catchment in the training set
## Run from the PAPER_BHGLM root. Outputs: RESULTS/LB/cv_*.csv
################################################################

suppressMessages({
  library(sf); library(spdep); library(INLA); library(dplyr); library(Matrix)
})
INLA::inla.setOption(num.threads = "4:1")
OUT <- "RESULTS/LB"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")
S <- 300

aoi <- st_read("DATA/df_catchments_kmeans.gpkg", quiet = TRUE)
aoi <- aoi %>% mutate_at(c("elev_mean", "slope_mean", "RainfallDaysmean"), ~ (scale(.) %>% as.vector))
aoi <- aoi %>% mutate(cuenca_num = case_when(cuenca == "Atrato" ~ 1, cuenca == "Cauca" ~ 2, cuenca == "Magdalena" ~ 3))
N <- nrow(aoi)
aoi.nb  <- poly2nb(aoi)
aoi.mat <- as(nb2mat(aoi.nb, style = "B"), "Matrix")
comp_id <- n.comp.nb(aoi.nb)$comp.id
d0 <- as.data.frame(aoi); if (!identical(as.integer(d0$id), 1:N)) d0$id <- 1:N
y_obs <- d0$lands_rec
Cmatrix <- Diagonal(N, 1) - (Diagonal(N, apply(aoi.mat, 1, sum)) - aoi.mat)

fx <- "lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean"
B  <- " + f(cuenca_num, model = 'iid')"
specs <- list(
  M1 = fx, M2 = paste0(fx, B), M2iid = paste0(fx, B, " + f(id, model = 'iid')"),
  M3 = paste0(fx, B, " + f(id, model = 'besag', graph = aoi.mat)"),
  M4 = paste0(fx, B, " + f(id, model = 'bym', graph = aoi.mat)"),
  M5 = paste0(fx, B, " + f(id, model = 'generic1', Cmatrix = Cmatrix)"),
  BYM2 = paste0(fx, B, " + f(id, model = 'bym2', graph = aoi.mat, scale.model = TRUE, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)), phi = list(prior = 'pc', param = c(0.5, 2/3))))"))

## Partitions
cent <- st_coordinates(st_centroid(st_geometry(aoi)))
parts <- list()
for (s in c(20260820, 101, 202)) { set.seed(s); parts[[paste0("spatial_seed", s)]] <- kmeans(cent, centers = 5, nstart = 25)$cluster }
set.seed(303); parts[["random"]] <- sample(rep(1:5, length.out = N))

logmeanexp <- function(x) { m <- max(x); if (!is.finite(m)) return(m); m + log(mean(exp(x - m))) }

rows <- list()
for (pn in names(parts)) {
  fold <- parts[[pn]]
  for (k in 1:5) {
    test <- which(fold == k)
    d <- d0; d$lands_rec[test] <- NA
    orphan <- !(comp_id[test] %in% comp_id[-test])      # component with no training catchment
    for (m in names(specs)) {
      log_line("%s | fold %d | %s", pn, k, m)
      fit <- try(inla(as.formula(specs[[m]]), family = "poisson", offset = log(d$area), data = d,
                      control.predictor = list(compute = TRUE, link = 1),
                      control.compute = list(config = TRUE), num.threads = "4:1"), silent = TRUE)
      if (inherits(fit, "try-error")) { log_line("  FAILED"); next }
      smp <- try(inla.posterior.sample(S, fit, selection = list(Predictor = test)), silent = TRUE)
      if (inherits(smp, "try-error")) { log_line("  sampling FAILED"); next }
      eta <- sapply(smp, function(s) s$latent[, 1])        # length(test) x S
      if (is.null(dim(eta))) eta <- matrix(eta, nrow = length(test))
      lam <- exp(eta)
      yt  <- y_obs[test]
      lpd <- sapply(seq_along(test), function(i) logmeanexp(dpois(yt[i], lam[i, ], log = TRUE)))
      yrep <- matrix(rpois(length(lam), pmin(lam, 1e7)), nrow = nrow(lam))
      crps <- sapply(seq_along(test), function(i) {
        yy <- yrep[i, ]; mean(abs(yy - yt[i])) - 0.5 * mean(abs(yy - sample(yy)))
      })
      med <- apply(yrep, 1, median)
      pm  <- fit$summary.fitted.values$mean[test]
      rows[[length(rows) + 1]] <- data.frame(
        partition = pn, fold = k, model = m, id = test, observed = yt, orphan = orphan,
        lpd = lpd, crps = crps, pred_median = med, post_mean_lambda = pm)
    }
  }
  write.csv(do.call(rbind, rows), file.path(OUT, "cv_per_catchment.csv"), row.names = FALSE)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT, "cv_per_catchment.csv"), row.names = FALSE)

## Summaries
by_fold <- res %>% group_by(partition, fold, model) %>%
  summarise(n = n(), n_orphan = sum(orphan), log_score = -mean(lpd), crps = mean(crps),
            rmse_median = sqrt(mean((pred_median - observed)^2)), mae_median = mean(abs(pred_median - observed)),
            rmse_postmean = sqrt(mean((post_mean_lambda - observed)^2)), .groups = "drop")
write.csv(by_fold, file.path(OUT, "cv_by_fold.csv"), row.names = FALSE)
res$type <- ifelse(res$partition == "random", "random", "spatial")
by_type <- res %>% group_by(type, model) %>%
  summarise(n = n(), log_score = -mean(lpd), crps = mean(crps),
            rmse_median = sqrt(mean((pred_median - observed)^2)), mae_median = mean(abs(pred_median - observed)),
            median_rmse_postmean_fold = NA_real_, .groups = "drop")
by_type$median_rmse_postmean_fold <- sapply(seq_len(nrow(by_type)), function(i) {
  sel <- by_fold$model == by_type$model[i] & ((by_fold$partition == "random") == (by_type$type[i] == "random"))
  median(by_fold$rmse_postmean[sel])
})
write.csv(by_type, file.path(OUT, "cv_summary.csv"), row.names = FALSE)
options(width = 200); print(as.data.frame(by_type), digits = 4)
print(res %>% filter(type == "spatial") %>% group_by(model, orphan) %>%
        summarise(n = n(), log_score = -mean(lpd), crps = mean(crps), .groups = "drop") %>% as.data.frame(), digits = 4)
log_line("done")
