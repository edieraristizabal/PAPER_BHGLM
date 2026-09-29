################################################################
## Final main-text figures (editorial comment EC10)
##   Fig02b  observed landslide counts per catchment (Fig. 2b)
##   Fig04   M1: Pearson residuals (A) and by basin (B)
##   Fig05   M2: Pearson residuals (A) and by basin (B)
##   Fig06-08 M3-M5: predicted counts (A) and Pearson residuals (B)
##   FigA2_slope_partial_residual.png (appendix, replaces the Python version)
## Changes with respect to CODES/INLA.R (models are identical):
##  - Pearson residuals: perceptually uniform diverging HCL palette
##    (colorspace "Blue-Red 3") centred on zero;
##  - counts: sequential multi-hue HCL palette ("YlOrRd"), with a common
##    scale across maps so that they are directly comparable;
##  - basin boundaries drawn as dark outlines with names instead of
##    saturated green/red/blue catchment borders; box plots use a
##    qualitative HCL palette of equal luminance ("Dark 3");
##  - panels composed with patchwork on a fixed grid so that A and B have
##    identical plot regions; all labels in English.
## Run from the PAPER_BHGLM root.
################################################################

suppressMessages({
  library(sf); library(spdep); library(INLA); library(dplyr); library(Matrix)
  library(ggplot2); library(patchwork); library(ggspatial); library(colorspace)
})
sf_use_s2(FALSE)
INLA::inla.setOption(num.threads = "4:1")
REV <- "FIGURES"

## ---------------------------------------------------------------
## Data and models, exactly as in CODES/INLA.R
## ---------------------------------------------------------------
aoi <- st_read("DATA/df_catchments_kmeans.gpkg", quiet = TRUE)
raw <- aoi
aoi <- aoi %>% mutate_at(c("elev_mean", "slope_mean", "RainfallDaysmean"), ~ (scale(.) %>% as.vector))
aoi.nb  <- poly2nb(aoi)
aoi.mat <- as(nb2mat(aoi.nb, style = "B"), "Matrix")
N <- nrow(aoi)
d <- as.data.frame(aoi); if (!identical(as.integer(d$id), 1:N)) d$id <- 1:N
Cmatrix <- Diagonal(N, 1) - (Diagonal(N, apply(aoi.mat, 1, sum)) - aoi.mat)

fit <- function(f) inla(f, family = "poisson", offset = log(d$area), data = d,
                        control.predictor = list(compute = TRUE),
                        control.compute = list(dic = TRUE, waic = TRUE))
f1 <- lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean
m <- list(
  M1 = fit(f1),
  M2 = fit(update(f1, . ~ . + f(cuenca, model = "iid"))),
  M3 = fit(update(f1, . ~ . + f(cuenca, model = "iid") + f(id, model = "besag", graph = aoi.mat))),
  M4 = fit(update(f1, . ~ . + f(cuenca, model = "iid") + f(id, model = "bym", graph = aoi.mat))),
  M5 = fit(update(f1, . ~ . + f(cuenca, model = "iid") + f(id, model = "generic1", Cmatrix = Cmatrix))))
for (k in names(m)) {
  mu <- m[[k]]$summary.fitted.values$mean[1:N]
  aoi[[paste0("pred_", k)]] <- mu
  aoi[[paste0("res_", k)]]  <- (aoi$lands_rec - mu) / sqrt(mu)
}
cat("DIC:", round(sapply(m, function(x) x$dic$dic)), "\n")

## ---------------------------------------------------------------
## Common graphical elements
## ---------------------------------------------------------------
## Basin perimeters for display only: morphological closing (+/-1.5 km) and
## removal of interior holes (valley floors not assigned to any catchment)
drop_holes <- function(g) st_multipolygon(lapply(st_cast(st_sfc(g), "POLYGON"), function(p) list(p[[1]])))
basins <- aoi %>% st_buffer(1500) %>% group_by(cuenca) %>% summarise(.groups = "drop") %>%
  st_buffer(-1500) %>% st_make_valid()
basins$geom <- st_sfc(lapply(st_geometry(basins), drop_holes), crs = st_crs(aoi))
st_geometry(basins) <- "geom"
lab_xy <- st_coordinates(st_point_on_surface(basins)); basins$lx <- lab_xy[, 1]; basins$ly <- lab_xy[, 2]
basin_pal <- setNames(qualitative_hcl(3, palette = "Dark 3"), c("Atrato", "Cauca", "Magdalena"))
count_max <- max(c(aoi$lands_rec, aoi$pred_M3, aoi$pred_M4, aoi$pred_M5))

map_theme <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 9), legend.title.position = "top",
        legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.3, "cm"),
        axis.title = element_blank(), plot.tag = element_text(face = "bold", size = 12))

base_map <- function(fill_var) {
  ggplot() +
    geom_sf(data = aoi, aes(fill = .data[[fill_var]]), colour = "grey55", linewidth = 0.08) +
    geom_sf(data = basins, fill = NA, colour = "grey10", linewidth = 0.45) +
    geom_label(data = basins, aes(lx, ly, label = cuenca), size = 2.6, fontface = "italic",
               linewidth = 0, fill = alpha("white", 0.7), label.padding = unit(0.1, "lines")) +
    annotation_scale(location = "bl", style = "ticks", text_cex = 0.7, line_width = 0.5) +
    annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_orienteering(text_size = 7),
                           height = unit(0.8, "cm"), width = unit(0.6, "cm")) +
    coord_sf(expand = FALSE) + map_theme
}
res_map <- function(k) {
  base_map(paste0("res_", k)) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", mid = 0, name = "Pearson residual")
}
count_map <- function(var, title) {
  base_map(var) +
    scale_fill_continuous_sequential(palette = "YlOrRd", limits = c(0, count_max), name = title)
}
res_box <- function(k) {
  ggplot(as.data.frame(aoi), aes(cuenca, .data[[paste0("res_", k)]], fill = cuenca)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
    geom_boxplot(notch = TRUE, outlier.shape = 1, outlier.size = 1, linewidth = 0.35, width = 0.6) +
    scale_fill_manual(values = basin_pal, guide = "none") +
    labs(x = NULL, y = "Pearson residual") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
                                     plot.tag = element_text(face = "bold", size = 12))
}
two_panel <- function(a, b, file) {
  p <- (a | b) + plot_layout(widths = c(1, 1)) + plot_annotation(tag_levels = "A")
  ggsave(file.path(REV, file), p, width = 180, height = 118, units = "mm", dpi = 400, bg = "white")
}

## ---------------------------------------------------------------
## Figures
## ---------------------------------------------------------------
ggsave(file.path(REV, "Fig02b_landslide_counts.png"), count_map("lands_rec", "Mapped landslides"),
       width = 90, height = 118, units = "mm", dpi = 400, bg = "white")
two_panel(res_map("M1"), res_box("M1"), "Fig04_M1_poisson.jpg")
two_panel(res_map("M2"), res_box("M2"), "Fig05_M2_basin_intercepts.jpg")
two_panel(count_map("pred_M3", "Predicted landslides"), res_map("M3"), "Fig06_M3_icar.jpg")
two_panel(count_map("pred_M4", "Predicted landslides"), res_map("M4"), "Fig07_M4_bym.jpg")
two_panel(count_map("pred_M5", "Predicted landslides"), res_map("M5"), "Fig08_M5_leroux.jpg")

## Appendix Fig. A2: partial-residual plot for mean slope in M1 (Poisson GLM)
g1 <- glm(lands_rec ~ RainfallDaysmean + elev_mean + slope_mean + offset(log(area)),
          family = poisson, data = d)
pr <- data.frame(x = d$slope_mean,
                 y = residuals(g1, type = "working") + coef(g1)["slope_mean"] * d$slope_mean)
pal2 <- qualitative_hcl(2, palette = "Dark 3")
pA2 <- ggplot(pr, aes(x, y)) +
  geom_point(shape = 1, size = 1, colour = "grey45") +
  geom_smooth(aes(colour = "LOWESS smooth"), method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8, span = 0.6) +
  geom_abline(aes(intercept = 0, slope = coef(g1)["slope_mean"], colour = "Linear fit"), linetype = "22", linewidth = 0.7) +
  scale_colour_manual(values = c("LOWESS smooth" = pal2[1], "Linear fit" = pal2[2]), name = NULL) +
  labs(x = "Standardized mean catchment slope", y = "Partial residual (working scale)") +
  theme_classic(base_size = 10) + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
ggsave(file.path(REV, "FigA2_slope_partial_residual.png"), pA2, width = 110, height = 90, units = "mm", dpi = 400, bg = "white")

## Colour-vision check: deuteranope and protanope simulations of the palettes
pal_check <- list(diverging = diverging_hcl(9, "Blue-Red 3"), sequential = sequential_hcl(9, "YlOrRd"),
                  qualitative = qualitative_hcl(3, "Dark 3"))
png("RESULTS/figs_colour_check.png", width = 1400, height = 900, res = 150)
swatchplot(pal_check, cvd = c("deutan", "protan"))
dev.off()
cat("figures written\n")
