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
colnames(aoi.mat) <- rownames(aoi.mat)

m2_cuenca <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + f(cuenca, model = "iid"),
                   offset = log(area), data = as.data.frame(aoi), family = "poisson",
                   control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE))

m3_icar <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean +
                   f(cuenca, model = "iid") + f(id, model = "besag", graph = aoi.mat),
                 offset = log(area), data = as.data.frame(aoi), family = "poisson",
                 control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE))

m4_bym <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean +
                  f(cuenca, model = "iid") + f(id, model = "bym", graph = aoi.mat),
                offset = log(area), data = as.data.frame(aoi), family = "poisson",
                control.compute = list(dic = TRUE, waic = TRUE), control.predictor = list(compute = TRUE))

D <- Diagonal(nrow(aoi.mat), apply(aoi.mat, 1, sum)) - aoi.mat
Cmatrix <- Diagonal(nrow(aoi), 1) - D
m5_ler <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean +
                 f(cuenca, model = "iid") + f(id, model = "generic1", Cmatrix = Cmatrix),
               offset = log(area), data = as.data.frame(aoi), family = "poisson",
               control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE))

cat("=== M2 cuenca ===\n"); print(m2_cuenca$summary.random$cuenca)
cat("=== M3 cuenca ===\n"); print(m3_icar$summary.random$cuenca)
cat("=== M4 cuenca ===\n"); print(m4_bym$summary.random$cuenca)
cat("=== M5 cuenca ===\n"); print(m5_ler$summary.random$cuenca)

out_dir <- "RESULTS"
rbind(
  cbind(model="M2", m2_cuenca$summary.random$cuenca),
  cbind(model="M3", m3_icar$summary.random$cuenca),
  cbind(model="M4", m4_bym$summary.random$cuenca),
  cbind(model="M5", m5_ler$summary.random$cuenca)
) -> all_cuenca
write.csv(all_cuenca, file.path(out_dir, "cuenca_random_effects.csv"), row.names = FALSE)
cat("DONE\n")
