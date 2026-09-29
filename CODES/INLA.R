################################################
############## PAPER BHGLM #####################
################################################

library(sf) 
library(spdep)
library(INLA)
library(sp)
library(viridis)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(ggspatial) #scale

#Data
## Run from the PAPER_BHGLM root
aoi = st_read("DATA/df_catchments_kmeans.gpkg",quiet = TRUE)
aoi <- aoi %>% mutate_at(c('elev_mean','slope_mean','RainfallDaysmean'), ~(scale(.) %>% as.vector))
aoi <- aoi %>% mutate(cuenca_num = case_when(cuenca == 'Atrato' ~ 1, cuenca == 'Cauca' ~ 2, cuenca == 'Magdalena' ~ 3))

#Spatial matrix
aoi.nb <- poly2nb(aoi) #Queen matrix
aoi.mat <- as(nb2mat(aoi.nb, style = "B"), "Matrix")
aoi.listw = nb2listw(aoi.nb)
colnames(aoi.mat) <- rownames(aoi.mat) 
mat <- as.matrix(aoi.mat[1:dim(aoi.mat)[1], 1:dim(aoi.mat)[1]])


#Standard Poisson model in INLA
m1_inla <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean,
           offset=log(area), data = as.data.frame(aoi),family = "poisson",
           control.predictor = list(compute = TRUE),
           control.compute = list(dic = TRUE, waic = TRUE))

summary(m1_inla)
aoi$m1_inla <- m1_inla$summary.fitted[1:nrow(aoi), "mean"]
res_pearson_m1 <- (aoi$lands_rec - m1_inla$summary.fitted.values$mean) / sqrt(m1_inla$summary.fitted.values$mean)
aoi$res_pearson_m1 <- res_pearson_m1
moran_m1 <- moran.mc(x = res_pearson_m1, listw = aoi.listw,nsim = 999, alternative = "greater")

#basic random intercept (cuenca)
m2_cuenca <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + 
           f(cuenca, model = "iid"),
           offset=log(area), data = as.data.frame(aoi), family = "poisson", 
           control.predictor = list(compute = TRUE),
           control.compute = list(dic = TRUE, waic = TRUE))

summary(m2_cuenca)
m2_cuenca$summary.random$cuenca
aoi$m2_cuenca <- m2_cuenca$summary.fitted[1:nrow(aoi), "mean"]
res_pearson_m2 <- (aoi$lands_rec - m2_cuenca$summary.fitted.values$mean) / sqrt(m2_cuenca$summary.fitted.values$mean)
aoi$res_pearson_m2 <- res_pearson_m2
moran_m2 <- moran.mc(x = res_pearson_m2, listw = aoi.listw,nsim = 999, alternative = "greater")

# ICAR model
m3_icar <- inla(formula = lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + 
                  f(cuenca, model = "iid") +
                  f(id, model = "besag", graph = aoi.mat),
                  offset=log(area), data = as.data.frame(aoi), family ="poisson",
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE))

summary(m3_icar)
m3_icar$summary.random$cuenca
aoi$ICAR <- m3_icar$summary.fitted[1:nrow(aoi), "mean"]
res_pearson_m3 <- (aoi$lands_rec - m3_icar$summary.fitted.values$mean) / sqrt(m3_icar$summary.fitted.values$mean)
aoi$res_pearson_m3 <- res_pearson_m3
moran_m3 <- moran.mc(x = res_pearson_m3, listw = aoi.listw,nsim = 999, alternative = "greater")

# BYM model
m4_bym<-inla(formula = lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean +
           f(cuenca, model = "iid") +
           f(id, model = "bym",graph = aoi.mat), 
           offset=log(area), data = as.data.frame(aoi), family = "poisson",
           control.compute = list(dic = TRUE, waic=T), 
           control.predictor = list(compute = TRUE))

summary(m4_bym)
m4_bym$summary.random$cuenca
aoi$BYM <- m4_bym$summary.fitted[1:nrow(aoi), "mean"]
res_pearson_m4 <- (aoi$lands_rec - m4_bym$summary.fitted.values$mean) / sqrt(m4_bym$summary.fitted.values$mean)
aoi$res_pearson_m4 <- res_pearson_m4
moran_m4 <- moran.mc(x = res_pearson_m4, listw = aoi.listw,nsim = 999, alternative = "greater")

#Leroux et al. model
m5_leroux <- Diagonal(nrow(aoi.mat), apply(aoi.mat, 1, sum)) - aoi.mat
Cmatrix <- Diagonal(nrow(aoi), 1) -  m5_leroux
max(eigen(Cmatrix)$values) #just to check =1

m5_ler = inla(formula = lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean +
                f(cuenca, model = "iid") +
                f(id, model = "generic1", Cmatrix = Cmatrix),
                offset=log(area), data = as.data.frame(aoi), family ="poisson",
                control.predictor = list(compute = TRUE),
                control.compute = list(dic = TRUE, waic = TRUE))

summary(m5_ler)
m5_ler$summary.random$cuenca
aoi$LEROUX <- m5_ler$summary.fitted[1:nrow(aoi), "mean"]
res_pearson_m5 <- (aoi$lands_rec - m5_ler$summary.fitted.values$mean) / sqrt(m5_ler$summary.fitted.values$mean)
aoi$res_pearson_m5 <- res_pearson_m5
moran_m5 <- moran.mc(x = res_pearson_m5, listw = aoi.listw,nsim = 999, alternative = "greater")

###################Results###############
## Tables: see extract_credible_intervals.R, extract_random_effects.R and spatial_cross_validation.R.
## Final figures (Figs. 2b, 4-8 and A2): see figures_main.R, which refits these models.
