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
aoi = st_read("G:/My Drive/INVESTIGACION/POSDOC/Data/Vector/df_catchments_kmeans.gpkg",quiet = TRUE)
aoi <- aoi %>% mutate_at(c('elev_mean','slope_mean','RainfallDaysmean'), ~(scale(.) %>% as.vector))
aoi <- aoi %>% mutate(cuenca_num = case_when(cuenca == 'Atrato' ~ 1, cuenca == 'Cauca' ~ 2, cuenca == 'Magdalena' ~ 3))

#Spatial matrix
aoi.nb <- poly2nb(aoi) #Queen matrix
aoi.mat <- as(nb2mat(aoi.nb, style = "B"), "Matrix")
aoi.listw = nb2listw(aoi.nb)
colnames(aoi.mat) <- rownames(aoi.mat) 
mat <- as.matrix(aoi.mat[1:dim(aoi.mat)[1], 1:dim(aoi.mat)[1]])

#Graph
nb2INLA("cl_graph",aoi.nb)
am_adj <-paste(getwd(),"G:/My Drive/INVESTIGACION/POSDOC/Figuras3/inla.graph",sep="")
H<-inla.read.graph(filename="G:/My Drive/INVESTIGACION/POSDOC/Figuras3/inla.graph")
image(inla.graph2matrix(H), xlab = "", ylab = "")

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
aoi_sp$BYM <- m4_bym$summary.fitted.values[, "mean"]
res_pearson_m4 <- (aoi$lands_rec - m4_bym$summary.fitted.values$mean) / sqrt(m4_bym$summary.fitted.values$mean)
aoi$res_pearson_m4 <- res_pearson_m4
moran_m4 <- moran.mc(x = res_pearson_m4, listw = aoi.listw,nsim = 999, alternative = "greater")

#Leroux et al. model
m5_leroux <- Diagonal(nrow(aoi.mat), apply(aoi.mat, 1, sum)) - aoi.mat
Cmatrix <- Diagonal(nrow(aoi), 1) -  m5_leroux
max(eigen(Cmatrix)$values) #just to check =1

m5_ler = inla(formula = lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean +
                f(cuenca, elev_mean, model = "iid") +
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

ggplot() + geom_sf(data=aoi,aes(fill=lands_rec),color = "black") +
  annotation_scale(location="br",style = "ticks") +
  annotation_north_arrow(location = "tr",which_north = "true", height = unit(1, "cm"), width = unit(1, "cm"),) +
  scale_fill_gradientn(colors=c("white","orange"),name = "Landslides") +
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

################################################################
#Figure M1
cuenca_colors <- c("Atrato" = "green", "Cauca" = "red", "Magdalena" = "blue")

p1 <- ggplot() + 
  geom_sf(data = aoi, aes(fill = res_pearson_m1, color = cuenca), size = 0.7) +
  scale_color_manual(values = cuenca_colors, guide = "none") +
  scale_fill_gradient2(
    low = "yellow", mid = "white", high = "purple", midpoint = 0,
    name = "Res. Pearson"
  ) +
  annotation_scale(
    location = "br", style = "ticks",
    text_cex = 1.0  # Tamaño de fuente en escala
  ) +
  annotation_north_arrow(
    location = "tr", which_north = "true",
    height = unit(1, "cm"), width = unit(1, "cm"),
    pad_x = unit(0.1, "cm"), pad_y = unit(0.1, "cm")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.3, 0.9),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.5, "cm"),
    axis.text = element_text(size = 11),          # Coordenadas del mapa
    axis.title = element_text(size = 13),         # Si decides añadir títulos
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

p2 <- ggplot(aoi, aes(x = cuenca, y = res_pearson_m1, fill = cuenca)) +
  geom_boxplot(notch = TRUE, outlier.shape = 1) +
  scale_fill_manual(values = cuenca_colors) +
  theme_bw() +
  labs(
    y = "Res. Pearson", x = ""  # Etiqueta eje Y
  ) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(size = 12),         # Etiquetas cuenca
    axis.text.y = element_text(size = 10),         # Ticks eje Y
    axis.title.y = element_text(size = 12),  # Título eje Y
    aspect.ratio = 1.95
  )

combined_plot <- plot_grid(p1, p2, ncol = 2, rel_widths = c(1.2, 0.8))

final_plot <- ggdraw(combined_plot) +
  draw_plot_label(
    label = c("A", "B"),
    x = c(0.1, 0.7), y = c(0.92, 0.92),
    hjust = 0, vjust = 1,
    fontface = "bold", size = 14)

print(final_plot)
ggsave("G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/PAPER_BHGLM/FIGURES/m1_map.jpg",plot = final_plot,dpi = 300,width = 8,height = 6,units = "in")


#Figure M2
cuenca_colors <- c("Atrato" = "green", "Cauca" = "red", "Magdalena" = "blue")

p1 <- ggplot() + 
  geom_sf(data = aoi, aes(fill = res_pearson_m2, color = cuenca), size = 0.7) +
  scale_color_manual(values = cuenca_colors, guide = "none") +
  scale_fill_gradient2(
    low = "yellow", mid = "white", high = "purple", midpoint = 0,
    name = "Res. Pearson"
  ) +
  annotation_scale(
    location = "br", style = "ticks",
    text_cex = 1.0  # Tamaño de fuente en escala
  ) +
  annotation_north_arrow(
    location = "tr", which_north = "true",
    height = unit(1, "cm"), width = unit(1, "cm"),
    pad_x = unit(0.1, "cm"), pad_y = unit(0.1, "cm")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.3, 0.9),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.5, "cm"),
    axis.text = element_text(size = 11),          # Coordenadas del mapa
    axis.title = element_text(size = 13),         # Si decides añadir títulos
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

p2 <- ggplot(aoi, aes(x = cuenca, y = res_pearson_m2, fill = cuenca)) +
  geom_boxplot(notch = TRUE, outlier.shape = 1) +
  scale_fill_manual(values = cuenca_colors) +
  theme_bw() +
  labs(
    y = "Res. Pearson", x = ""  # Etiqueta eje Y
  ) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(size = 12),         # Etiquetas cuenca
    axis.text.y = element_text(size = 10),         # Ticks eje Y
    axis.title.y = element_text(size = 12),  # Título eje Y
    aspect.ratio = 1.95
  )

combined_plot <- plot_grid(p1, p2, ncol = 2, rel_widths = c(1.2, 0.8))

final_plot <- ggdraw(combined_plot) +
  draw_plot_label(
    label = c("A", "B"),
    x = c(0.1, 0.7), y = c(0.92, 0.92),
    hjust = 0, vjust = 1,
    fontface = "bold", size = 14
  )

print(final_plot)
ggsave("G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/PAPER_BHGLM/FIGURES/m2_map.jpg",plot = final_plot,dpi = 300,width = 8,height = 6,units = "in")
###################################################

#Figura M3

cuenca_colors <- c("Atrato" = "green", "Cauca" = "red", "Magdalena" = "blue")
p1 <- ggplot() + 
  geom_sf(data = aoi, aes(fill = ICAR, color = cuenca), size = 0.7) +
  scale_color_manual(values = cuenca_colors, guide = "none") +
  scale_fill_gradientn(
    colours = c("white", "orange", "red", "darkred"),
    name = "Landslides"
  ) +
  annotation_scale(location = "br", style = "ticks", text_cex = 1.0) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         pad_x = unit(0.1, "cm"), pad_y = unit(0.1, "cm")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.3, 0.9),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.5, "cm"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

p2 <- ggplot() + 
  geom_sf(data = aoi, aes(fill = res_pearson_m3, color = cuenca), size = 0.7) +
  scale_color_manual(values = cuenca_colors, guide = "none") +
  scale_fill_gradient2(
    low = "yellow", mid = "white", high = "purple", midpoint = 0,
    name = "Res. Pearson"
  ) +
  annotation_scale(location = "br", style = "ticks", text_cex = 1.0) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         pad_x = unit(0.1, "cm"), pad_y = unit(0.1, "cm")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.3, 0.9),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.5, "cm"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank())

combined_plot <- plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1))
final_plot <- ggdraw(combined_plot) +
  draw_plot_label(
    label = c("A", "B"),
    x = c(0.1, 0.6), y = c(0.96, 0.96),
    hjust = 0, vjust = 1,
    fontface = "bold", size = 14
  )
print(final_plot)
ggsave("G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/PAPER_BHGLM/FIGURES/m3_icar.jpg",
       plot = final_plot, dpi = 300, width = 12, height = 6, units = "in")

###################################################

#Figura M4

cuenca_colors <- c("Atrato" = "green", "Cauca" = "red", "Magdalena" = "blue")
p1 <- ggplot() + 
  geom_sf(data = aoi, aes(fill = BYM, color = cuenca), size = 0.7) +
  scale_color_manual(values = cuenca_colors, guide = "none") +
  scale_fill_gradientn(
    colours = c("white", "orange", "red", "darkred"),
    name = "Landslides"
  ) +
  annotation_scale(location = "br", style = "ticks", text_cex = 1.0) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         pad_x = unit(0.1, "cm"), pad_y = unit(0.1, "cm")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.3, 0.9),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.5, "cm"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

p2 <- ggplot() + 
  geom_sf(data = aoi, aes(fill = res_pearson_m4, color = cuenca), size = 0.7) +
  scale_color_manual(values = cuenca_colors, guide = "none") +
  scale_fill_gradient2(
    low = "yellow", mid = "white", high = "purple", midpoint = 0,
    name = "Res. Pearson"
  ) +
  annotation_scale(location = "br", style = "ticks", text_cex = 1.0) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         pad_x = unit(0.1, "cm"), pad_y = unit(0.1, "cm")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.3, 0.9),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.7, "cm"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank())

combined_plot <- plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1))
final_plot <- ggdraw(combined_plot) +
  draw_plot_label(
    label = c("A", "B"),
    x = c(0.1, 0.6), y = c(0.96, 0.96),
    hjust = 0, vjust = 1,
    fontface = "bold", size = 14
  )
print(final_plot)
ggsave("G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/PAPER_BHGLM/FIGURES/m4_bym.jpg",
       plot = final_plot, dpi = 300, width = 12, height = 6, units = "in")
###################################################

#Figura M5

cuenca_colors <- c("Atrato" = "green", "Cauca" = "red", "Magdalena" = "blue")
p1 <- ggplot() + 
  geom_sf(data = aoi, aes(fill = LEROUX, color = cuenca), size = 0.7) +
  scale_color_manual(values = cuenca_colors, guide = "none") +
  scale_fill_gradientn(
    colours = c("white", "orange", "red", "darkred"),
    name = "Landslides"
  ) +
  annotation_scale(location = "br", style = "ticks", text_cex = 1.0) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         pad_x = unit(0.1, "cm"), pad_y = unit(0.1, "cm")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.3, 0.9),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.5, "cm"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

p2 <- ggplot() + 
  geom_sf(data = aoi, aes(fill = res_pearson_m5, color = cuenca), size = 0.7) +
  scale_color_manual(values = cuenca_colors, guide = "none") +
  scale_fill_gradient2(
    low = "yellow", mid = "white", high = "purple", midpoint = 0,
    name = "Res. Pearson"
  ) +
  annotation_scale(location = "br", style = "ticks", text_cex = 1.0) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         pad_x = unit(0.1, "cm"), pad_y = unit(0.1, "cm")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.3, 0.9),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.7, "cm"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank())

combined_plot <- plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1))
final_plot <- ggdraw(combined_plot) +
  draw_plot_label(
    label = c("A", "B"),
    x = c(0.1, 0.6), y = c(0.96, 0.96),
    hjust = 0, vjust = 1,
    fontface = "bold", size = 14
  )
print(final_plot)
ggsave("G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/PAPER_BHGLM/FIGURES/m5_ler.jpg",
       plot = final_plot, dpi = 300, width = 12, height = 6, units = "in")


