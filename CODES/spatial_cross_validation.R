################################################################
# Spatial k-fold cross-validation for M1-M5 (RC2 response,
# comment 2.5: train/test split). Uses spatially-blocked folds
# (k-means on catchment centroids) rather than random folds, to
# avoid leakage from strong spatial autocorrelation between
# neighboring catchments (Brenning, 2005, NHESS).
#
# For each fold, the response for held-out catchments is set to
# NA (INLA computes posterior predictive there) while the full
# adjacency graph is kept intact, so CAR/BYM/Leroux models can
# still borrow strength from neighboring catchments as they
# would for genuinely unobserved units.
################################################################

library(sf)
library(spdep)
library(INLA)
library(dplyr)
library(Matrix)

set.seed(20260820)

aoi <- st_read("DATA/df_catchments_kmeans.gpkg", quiet = TRUE)
aoi <- aoi %>% mutate_at(c('elev_mean','slope_mean','RainfallDaysmean'), ~(scale(.) %>% as.vector))
aoi <- aoi %>% mutate(cuenca_num = case_when(cuenca == 'Atrato' ~ 1, cuenca == 'Cauca' ~ 2, cuenca == 'Magdalena' ~ 3))

aoi.nb <- poly2nb(aoi)
aoi.mat <- as(nb2mat(aoi.nb, style = "B"), "Matrix")
colnames(aoi.mat) <- rownames(aoi.mat)

D <- Diagonal(nrow(aoi.mat), apply(aoi.mat, 1, sum)) - aoi.mat
Cmatrix <- Diagonal(nrow(aoi), 1) - D

# --- spatial folds via k-means on centroids ---
K <- 5
cent <- st_coordinates(st_centroid(aoi))
km <- kmeans(cent, centers = K, nstart = 25)
fold <- km$cluster
cat("Fold sizes:", table(fold), "\n")

y_full <- aoi$lands_rec
df_base <- as.data.frame(aoi)

fit_model <- function(model_name, df) {
  form_fixed <- lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean
  if (model_name == "M1") {
    m <- inla(form_fixed, offset = log(area), data = df, family = "poisson",
              control.predictor = list(compute = TRUE, link = 1),
              control.compute = list(config = TRUE))
  } else if (model_name == "M2") {
    m <- inla(update(form_fixed, . ~ . + f(cuenca, model = "iid")),
              offset = log(area), data = df, family = "poisson",
              control.predictor = list(compute = TRUE, link = 1),
              control.compute = list(config = TRUE))
  } else if (model_name == "M3") {
    m <- inla(update(form_fixed, . ~ . + f(cuenca, model = "iid") + f(id, model = "besag", graph = aoi.mat)),
              offset = log(area), data = df, family = "poisson",
              control.predictor = list(compute = TRUE, link = 1),
              control.compute = list(config = TRUE))
  } else if (model_name == "M4") {
    m <- inla(update(form_fixed, . ~ . + f(cuenca, model = "iid") + f(id, model = "bym", graph = aoi.mat)),
              offset = log(area), data = df, family = "poisson",
              control.predictor = list(compute = TRUE, link = 1),
              control.compute = list(config = TRUE))
  } else if (model_name == "M5") {
    m <- inla(update(form_fixed, . ~ . + f(cuenca, model = "iid") + f(id, model = "generic1", Cmatrix = Cmatrix)),
              offset = log(area), data = df, family = "poisson",
              control.predictor = list(compute = TRUE, link = 1),
              control.compute = list(config = TRUE))
  }
  m
}

models <- c("M1", "M2", "M3", "M4", "M5")
results <- data.frame()

for (k in 1:K) {
  test_idx <- which(fold == k)
  df_k <- df_base
  df_k$lands_rec[test_idx] <- NA  # hold out: INLA predicts these

  for (mname in models) {
    cat("Fold", k, "-", mname, "...\n")
    m <- fit_model(mname, df_k)
    lambda_hat <- m$summary.fitted.values$mean[test_idx]
    y_obs <- y_full[test_idx]

    loglik <- dpois(y_obs, lambda = pmax(lambda_hat, 1e-6), log = TRUE)
    rmse <- sqrt(mean((y_obs - lambda_hat)^2))
    mae <- mean(abs(y_obs - lambda_hat))

    results <- rbind(results, data.frame(
      fold = k, model = mname, n_test = length(test_idx),
      mean_loglik = mean(loglik), rmse = rmse, mae = mae
    ))
  }
}

out_dir <- "RESULTS"
write.csv(results, file.path(out_dir, "spatial_cv_results.csv"), row.names = FALSE)

summary_tab <- results %>%
  group_by(model) %>%
  summarise(mean_loglik = mean(mean_loglik), sd_loglik = sd(mean_loglik),
            mean_rmse = mean(rmse), sd_rmse = sd(rmse),
            mean_mae = mean(mae), sd_mae = sd(mae))
write.csv(summary_tab, file.path(out_dir, "spatial_cv_summary.csv"), row.names = FALSE)

cat("\n=== Per-fold results ===\n")
print(results)
cat("\n=== Summary across folds ===\n")
print(summary_tab)
cat("\nDONE\n")
