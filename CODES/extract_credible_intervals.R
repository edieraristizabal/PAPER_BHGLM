################################################################
# Re-fit M1-M5 (same specification as INLA.R) with a local data
# path, and extract 95% credible/confidence intervals for the
# fixed effects, plus DIC/WAIC/Moran's I, to verify Table 1 and
# answer RC1 Additional Comment #3 (report CI/CrI consistently).
################################################################

library(sf)
library(spdep)
library(INLA)
library(dplyr)
library(Matrix)

aoi <- st_read("DATA/df_catchments_kmeans.gpkg", quiet = TRUE)
aoi <- aoi %>% mutate_at(c('elev_mean','slope_mean','RainfallDaysmean'), ~(scale(.) %>% as.vector))
aoi <- aoi %>% mutate(cuenca_num = case_when(cuenca == 'Atrato' ~ 1, cuenca == 'Cauca' ~ 2, cuenca == 'Magdalena' ~ 3))

aoi.nb <- poly2nb(aoi)
aoi.mat <- as(nb2mat(aoi.nb, style = "B"), "Matrix")
aoi.listw <- nb2listw(aoi.nb)
colnames(aoi.mat) <- rownames(aoi.mat)
mat <- as.matrix(aoi.mat[1:dim(aoi.mat)[1], 1:dim(aoi.mat)[1]])

summarize_model <- function(model, name, resid_pearson) {
  fx <- as.data.frame(model$summary.fixed)
  fx$model <- name
  fx$term <- rownames(fx)
  moran <- moran.mc(x = resid_pearson, listw = aoi.listw, nsim = 999, alternative = "greater")
  list(
    fixed = fx,
    hyperpar = as.data.frame(model$summary.hyperpar),
    dic = model$dic$dic,
    waic = model$waic$waic,
    moran_I = moran$statistic
  )
}

cat("Fitting M1...\n")
m1_inla <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean,
                 offset = log(area), data = as.data.frame(aoi), family = "poisson",
                 control.predictor = list(compute = TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE))
res_pearson_m1 <- (aoi$lands_rec - m1_inla$summary.fitted.values$mean) / sqrt(m1_inla$summary.fitted.values$mean)
r1 <- summarize_model(m1_inla, "M1", res_pearson_m1)

cat("Fitting M2...\n")
m2_cuenca <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + f(cuenca, model = "iid"),
                   offset = log(area), data = as.data.frame(aoi), family = "poisson",
                   control.predictor = list(compute = TRUE),
                   control.compute = list(dic = TRUE, waic = TRUE))
res_pearson_m2 <- (aoi$lands_rec - m2_cuenca$summary.fitted.values$mean) / sqrt(m2_cuenca$summary.fitted.values$mean)
r2 <- summarize_model(m2_cuenca, "M2", res_pearson_m2)

cat("Fitting M3 (ICAR)...\n")
m3_icar <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean +
                   f(cuenca, model = "iid") + f(id, model = "besag", graph = aoi.mat),
                 offset = log(area), data = as.data.frame(aoi), family = "poisson",
                 control.predictor = list(compute = TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE))
res_pearson_m3 <- (aoi$lands_rec - m3_icar$summary.fitted.values$mean) / sqrt(m3_icar$summary.fitted.values$mean)
r3 <- summarize_model(m3_icar, "M3", res_pearson_m3)

cat("Fitting M4 (BYM)...\n")
m4_bym <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean +
                  f(cuenca, model = "iid") + f(id, model = "bym", graph = aoi.mat),
                offset = log(area), data = as.data.frame(aoi), family = "poisson",
                control.compute = list(dic = TRUE, waic = TRUE),
                control.predictor = list(compute = TRUE))
res_pearson_m4 <- (aoi$lands_rec - m4_bym$summary.fitted.values$mean) / sqrt(m4_bym$summary.fitted.values$mean)
r4 <- summarize_model(m4_bym, "M4", res_pearson_m4)

cat("Fitting M5 (Leroux)...\n")
D <- Diagonal(nrow(aoi.mat), apply(aoi.mat, 1, sum)) - aoi.mat
Cmatrix <- Diagonal(nrow(aoi), 1) - D
m5_ler <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean +
                 f(cuenca, model = "iid") + f(id, model = "generic1", Cmatrix = Cmatrix),
               offset = log(area), data = as.data.frame(aoi), family = "poisson",
               control.predictor = list(compute = TRUE),
               control.compute = list(dic = TRUE, waic = TRUE))
res_pearson_m5 <- (aoi$lands_rec - m5_ler$summary.fitted.values$mean) / sqrt(m5_ler$summary.fitted.values$mean)
r5 <- summarize_model(m5_ler, "M5", res_pearson_m5)

all_results <- list(M1 = r1, M2 = r2, M3 = r3, M4 = r4, M5 = r5)

out_dir <- "RESULTS"

fixed_all <- do.call(rbind, lapply(all_results, function(x) x$fixed))
write.csv(fixed_all, file.path(out_dir, "credible_intervals_fixed.csv"), row.names = FALSE)

summary_tab <- data.frame(
  model = names(all_results),
  DIC = sapply(all_results, function(x) x$dic),
  WAIC = sapply(all_results, function(x) x$waic),
  MoranI = sapply(all_results, function(x) x$moran_I)
)
write.csv(summary_tab, file.path(out_dir, "credible_intervals_summary.csv"), row.names = FALSE)

cat("\n=== FIXED EFFECTS (mean, sd, 0.025, 0.5, 0.975) ===\n")
print(fixed_all)
cat("\n=== DIC / WAIC / Moran's I per model ===\n")
print(summary_tab)

cat("\n=== Hyperparameters ===\n")
for (nm in names(all_results)) {
  cat("\n--", nm, "--\n")
  print(all_results[[nm]]$hyperpar)
}

cat("\nDONE\n")
