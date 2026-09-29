################################################################
## EC4 (part 1) -- Alternative aggregation of terrain predictors
## Editorial comment 4 (Schlogl): demonstrate that the catchment
## support / aggregation is appropriate. Computes, from the native
## 12.5 m ALOS-PALSAR DEM used in the manuscript (Sect. 3), per
## catchment:
##   slope_mean12, elev_mean12  -> validation against the published
##                                 slope_mean / elev_mean
##   slope_p90                  -> 90th percentile of slope
##   frac_gt20/25/30            -> proportion of area steeper than
##                                 20/25/30 degrees (coarse analogue
##                                 of a potential-process-area
##                                 restriction, Steger et al.)
##   slope_sd, elev_sd          -> within-catchment variability that
##                                 the catchment mean discards
## Slope: terra::terrain, 8 neighbours (Horn), degrees.
## Read block-wise and accumulated as per-catchment histograms, so
## the ~6e8 cells are never held in memory.
## Run from the PAPER_BHGLM root.
################################################################

suppressMessages({ library(sf); library(terra) })
terraOptions(progress = 0, memfrac = 0.5)

DEM_PATH <- "C:/Users/edier/Documents/INVESTIGACION/PAPERS/PUBLICADOS/PAPER_SAR/Data/DEM_filled_12m.tif"
OUT      <- "RESULTS/EC4"
TMP      <- file.path(tempdir(), "ec4")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TMP, recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")
tabulate_w <- function(id, w, n) {           # weighted tabulate
  r <- rowsum(w, id, reorder = FALSE)
  out <- numeric(n); out[as.integer(rownames(r))] <- r[, 1]; out
}

aoi <- st_read("DATA/df_catchments_kmeans.gpkg", quiet = TRUE)
aoi$rid <- seq_len(nrow(aoi))
N <- nrow(aoi)

## 1. Crop DEM to the catchments (+250 m so edge slopes are defined)
dem  <- rast(DEM_PATH)
dem  <- crop(dem, ext(vect(aoi)) + 250, snap = "out",
             filename = file.path(TMP, "dem.tif"), overwrite = TRUE)
log_line("DEM cropped: %d x %d", nrow(dem), ncol(dem))

## 2. Slope (degrees) and catchment-id raster on the same grid
slp <- terrain(dem, v = "slope", neighbors = 8, unit = "degrees",
               filename = file.path(TMP, "slope.tif"), overwrite = TRUE)
log_line("Slope computed")
idr <- rasterize(vect(aoi), dem, field = "rid", datatype = "INT2U",
                 filename = file.path(TMP, "id.tif"), overwrite = TRUE)
log_line("Catchment raster done")

## 3. Block-wise accumulation
BW   <- 0.05                          # histogram bin width (degrees)
NB   <- ceiling(90 / BW)
hist <- numeric(N * NB)               # per-catchment slope histograms
zsum <- numeric(N); zss <- numeric(N); zn <- numeric(N)
ssum <- numeric(N); sss <- numeric(N)

readStart(dem); readStart(slp); readStart(idr)
nr <- nrow(dem); step <- 1000L
starts <- seq(1L, nr, by = step)
for (k in seq_along(starts)) {
  r0 <- starts[k]; nrw <- min(step, nr - r0 + 1L)
  z  <- readValues(dem, r0, nrw)
  s  <- readValues(slp, r0, nrw)
  id <- readValues(idr, r0, nrw)
  ok <- !is.na(id) & !is.na(s) & !is.na(z) & z > -1000
  if (!any(ok)) next
  z <- z[ok]; s <- s[ok]; id <- as.integer(id[ok])
  bin  <- pmin(as.integer(s %/% BW) + 1L, NB)
  hist <- hist + tabulate((id - 1L) * NB + bin, nbins = N * NB)
  zsum <- zsum + tabulate_w(id, z,   N)
  zss  <- zss  + tabulate_w(id, z^2, N)
  ssum <- ssum + tabulate_w(id, s,   N)
  sss  <- sss  + tabulate_w(id, s^2, N)
  zn   <- zn   + tabulate(id, nbins = N)
  log_line("rows %d-%d of %d", r0, r0 + nrw - 1L, nr)
}
readStop(dem); readStop(slp); readStop(idr)

## 4. Per-catchment statistics
H    <- matrix(hist, nrow = N, ncol = NB, byrow = TRUE)
mids <- (seq_len(NB) - 0.5) * BW
q_from_hist <- function(h, p) {
  cs <- cumsum(h) / sum(h); j <- which(cs >= p)[1]
  lo <- (j - 1) * BW; prev <- if (j > 1) cs[j - 1] else 0
  lo + BW * (p - prev) / (cs[j] - prev)       # linear within bin
}
frac_gt <- function(th) rowSums(H[, mids > th, drop = FALSE]) / rowSums(H)

res <- data.frame(
  id           = aoi$id,
  ID_CUENCA    = aoi$ID_CUENCA,
  cuenca       = aoi$cuenca,
  n_cells      = zn,
  slope_mean   = aoi$slope_mean,           # published predictor
  elev_mean    = aoi$elev_mean,
  slope_mean12 = ssum / zn,
  slope_sd     = sqrt(pmax(sss / zn - (ssum / zn)^2, 0)),
  slope_p50    = apply(H, 1, q_from_hist, p = 0.50),
  slope_p90    = apply(H, 1, q_from_hist, p = 0.90),
  frac_gt20    = frac_gt(20),
  frac_gt25    = frac_gt(25),
  frac_gt30    = frac_gt(30),
  elev_mean12  = zsum / zn,
  elev_sd      = sqrt(pmax(zss / zn - (zsum / zn)^2, 0))
)
write.csv(res, file.path(OUT, "ec4_catchment_terrain_12m.csv"), row.names = FALSE)

## 5. Validation against the published predictors
v <- data.frame(
  check = c("slope_mean12 vs slope_mean", "elev_mean12 vs elev_mean"),
  pearson_r = c(cor(res$slope_mean12, res$slope_mean), cor(res$elev_mean12, res$elev_mean)),
  mean_abs_diff = c(mean(abs(res$slope_mean12 - res$slope_mean)),
                    mean(abs(res$elev_mean12 - res$elev_mean))),
  max_abs_diff = c(max(abs(res$slope_mean12 - res$slope_mean)),
                   max(abs(res$elev_mean12 - res$elev_mean))))
write.csv(v, file.path(OUT, "ec4_validation.csv"), row.names = FALSE)
print(v)

## 6. Correlation among the alternative slope aggregations
cm <- cor(res[, c("slope_mean", "slope_p50", "slope_p90", "frac_gt20",
                  "frac_gt25", "frac_gt30", "slope_sd", "elev_mean", "elev_sd")])
write.csv(round(cm, 3), file.path(OUT, "ec4_aggregation_correlations.csv"))
print(round(cm, 3))
print(summary(res[, -(1:4)]))
log_line("done")
