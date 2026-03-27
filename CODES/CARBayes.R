###################################
####CARBayes
################################

library(spBayes)
library(maps)
library(RANN)
library(gjam)
library(CARBayes)
library(CARBayesdata)
library(mgcv)
library(sf) 
library(spdep)
library(sp)
library(viridis)
library(dplyr)

aoi = st_read("G:/My Drive/INVESTIGACION/POSDOC/Data/Vector/df_catchments_kmeans.gpkg",quiet = TRUE)
aoi2 <- aoi %>% mutate_at(c('elev_mean','slope_mean','RainfallDaysmean'), ~(scale(.) %>% as.vector))
aoi_sp <- as_Spatial(aoi)

#Spatial matrix
W.nb <- poly2nb(aoi)
W.mat <- nb2mat(W.nb, style="B")
W.list <- nb2listw(W.nb, style="B")

#Standard Poisson model in CARBayes
m1_carbayes <- S.glm(formula=lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + offset(log(area)),
                    data=aoi, family="poisson", 
                    burnin=100000, n.sample=300000, thin=100, n.chains=3, n.cores=3)


m1_carbayes
m1_carbayes$modelfit
aoi_sp$carm1 <- m1_carbayes$fitted.values
moran.mc(x=residuals(m1_carbayes), listw=W.list, nsim=1000)



#Leroux=proper
cariid <- S.CARleroux(formula=lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + offset(log(area)),
                      data=aoi, family="poisson", W=W.mat,rho=0.5,
                      burnin=100000, n.sample=300000, thin=100, n.chains=3, n.cores=3)
cariid
cariid$modelfit
plot( cariid$samples$rho, bty = 'n' )
summary(cariid$samples)
summary.beta <- summary(cariid$samples$beta, quantiles=c(0.025, 0.975))
beta.mean <- summary.beta$statistics[ ,"Mean"]
beta.ci <- summary.beta$quantiles
beta.results <- cbind(beta.mean, beta.ci)
rownames(beta.results) <- colnames(cariid$X)
round(beta.results, 5)
aoi_sp$cariid <- cariid$fitted.values
aoi_sp$cariid_res <- cariid$residuals
moran.mc(x=residuals(cariid), listw=W.list, nsim=1000)

#Leroux=ICAR
caricar <- S.CARleroux(formula=lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + offset(log(area)),
                       data=aoi, family="poisson", W=W.mat,rho=1,
                       burnin=100000, n.sample=300000, thin=100, n.chains=3, n.cores=3)
caricar
caricar$modelfit
summary(caricar$samples)
summary.beta <- summary(caricar$samples$beta, quantiles=c(0.025, 0.975))
beta.mean <- summary.beta$statistics[ ,"Mean"]
beta.ci <- summary.beta$quantiles
beta.results <- cbind(beta.mean, beta.ci)
rownames(beta.results) <- colnames(caricar$X)
round(beta.results, 5)
aoi_sp$caricar <- caricar$fitted.values
aoi_sp$caricar_res <- caricar$residuals
moran.mc(x=residuals(caricar), listw=W.list, nsim=1000)


#BYM
carbym <- S.CARbym(formula=lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + offset(log(area)),
                   data=aoi, family="poisson", W=W.mat,
                   burnin=100000, n.sample=300000, thin=100, n.chains=3, n.cores=3)
carbym
carbym$modelfit
plot( carbym$samples$rho, bty = 'n' )
summary(carbym$samples)
summary.beta <- summary(carbym$samples$beta, quantiles=c(0.025, 0.975))
beta.mean <- summary.beta$statistics[ ,"Mean"]
beta.ci <- summary.beta$quantiles
beta.results <- cbind(beta.mean, beta.ci)
rownames(beta.results) <- colnames(carbym$X)
round(beta.results, 5)
aoi$carbym <- carbym$fitted.values
aoi$carbym_res <- carbym$residuals
moran.mc(x=residuals(carbym), listw=W.list, nsim=1000)


#Leroux
carleroux <- S.CARleroux(formula=lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + offset(log(area)),
                  data=aoi, family="poisson", W=W.mat,
                  burnin=100000, n.sample=300000, thin=100, n.chains=3, n.cores=3)
carleroux
carleroux$modelfit
plot( carleroux$samples$rho, bty = 'n' )
summary(carleroux$samples)
summary.beta <- summary(carleroux$samples$beta, quantiles=c(0.025, 0.975))
beta.mean <- summary.beta$statistics[ ,"Mean"]
beta.ci <- summary.beta$quantiles
beta.results <- cbind(beta.mean, beta.ci)
rownames(beta.results) <- colnames(carleroux$X)
round(beta.results, 5)
aoi_sp$carleroux <- carleroux$fitted.values
aoi$carleroux_res <- carleroux$residuals
moran.mc(x=residuals(carleroux), listw=W.list, nsim=1000)

############################

ggplot() + geom_sf(data=aoi,aes(fill=carbym),color = "black") +
  annotation_scale(location="br",style = "ticks") +
  annotation_north_arrow(location = "tr",which_north = "true", height = unit(1, "cm"), width = unit(1, "cm"),) +
  scale_fill_gradientn(colors=c("white","red"),name = "Landslides") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,size = 1),
    axis.ticks.length=unit(-0.1, "cm"),
    axis.text.x = element_text(size = 10, margin = unit(c(t = 1, r = 0, b = 0, l = 0), "mm")),
    axis.text.y = element_text(size = 10, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
    legend.text = element_text(size=12),
    legend.title.align = 0,
    legend.position = c(0.34,0.9), 
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size=14, vjust = .8, hjust = .5))

ggplot() + geom_sf(data=aoi,aes(fill=carbym_res$response),color = "black") +
  annotation_scale(location="br",style = "ticks") +
  annotation_north_arrow(location = "tr",which_north = "true", height = unit(1, "cm"), width = unit(1, "cm"),) +
  scale_fill_gradientn(colors=c("white","orange"),name = "Residuals") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,size = 1),
    axis.ticks.length=unit(-0.1, "cm"),
    axis.text.x = element_text(size = 10, margin = unit(c(t = 1, r = 0, b = 0, l = 0), "mm")),
    axis.text.y = element_text(size = 10, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
    legend.text = element_text(size=12),
    legend.title.align = 0,
    legend.position = c(0.34,0.9), 
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size=14, vjust = .8, hjust = .5)) 

##############################

#Multilevel
carmultilevel <- S.CARmultilevel(formula=lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + offset(log(area)),
                                 data=aoi,family="poisson", ind.area=aoi$cuenca_num, W=W.mat,verbose=FALSE,
                                 burnin=100000, n.sample=300000,thin=100, n.chains=3, n.cores=3)


###############################

spplot(aoi_sp, c("lands_rec", "caricar", "carbym", "carleroux"),col.regions = rev(magma(16)))
spplot(aoi_sp, c("cariid", "caricar", "carbym", "carleroux"),col.regions = rev(magma(16)))
spplot(aoi_sp, c("carm1","cariid"),col.regions = rev(magma(16)))


