################################################################
# Predictor screening for RC1 response (Major Concern 4)
# Computes candidate covariates not in the original model (TWI,
# curvature, drainage density), joins them with the catchment
# dataset, and produces the full Pearson correlation matrix used
# to check multicollinearity (|r| > 0.75) against the retained
# predictors (RainfallDaysmean, elev_mean, slope_mean).
#
# Requires the 'geospatial' conda environment (rasterio, pysheds,
# geopandas, scipy, matplotlib).
################################################################

import numpy as np
if not hasattr(np, 'in1d'):
    # pysheds 0.5 calls np.in1d, removed in numpy>=2.0
    np.in1d = np.isin

import rasterio
from rasterio.enums import Resampling
from rasterio.mask import mask as rio_mask
from scipy import ndimage
from pysheds.grid import Grid
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt

DEM_PATH = r"c:/Users/edier/Documents/INVESTIGACION/PAPERS/PUBLICADOS/PAPER_SAR/Data/DEM_filled_12m.tif"
DRAINAGE_NETWORK_PATH = r"c:/Users/edier/Documents/INVESTIGACION/PAPERS/PUBLICADOS/PAPER_SAR/Data/drainage_network_5km2.gpkg"
CATCHMENTS_PATH = r"DATA/df_catchments_kmeans.gpkg"

WORKDIR = r"DATA"
DEM_25M = WORKDIR + "/dem_25m.tif"
TWI_25M = WORKDIR + "/twi_25m.tif"
CURV_25M = WORKDIR + "/curvature_25m.tif"

RESAMPLED_RES = 25.0  # metres; coarser than the native 12.5 m ALOS-PALSAR DEM used for
                       # slope/elevation, chosen only to make flow accumulation tractable
                       # over the full three-basin extent (~1e9 cells at 12.5 m).

RETAINED = ['RainfallDaysmean', 'elev_mean', 'slope_mean']
NEW_CANDIDATES = ['TWI', 'curvature', 'drainage_density']
ORIGINAL_SCREENED = ['area', 'hypso_inte', 'Densidad', 'rel_mean', 'rainfallAnnual_mean']
R_THRESHOLD = 0.75


def resample_dem():
    with rasterio.open(DEM_PATH) as src:
        scale = RESAMPLED_RES / src.res[0]
        new_width = int(src.width / scale)
        new_height = int(src.height / scale)
        data = src.read(1, out_shape=(new_height, new_width), resampling=Resampling.average)
        new_transform = src.transform * src.transform.scale(src.width / new_width, src.height / new_height)
        profile = src.profile.copy()
        profile.update(height=new_height, width=new_width, transform=new_transform, dtype='float32')
        with rasterio.open(DEM_25M, 'w', **profile) as dst:
            dst.write(data.astype('float32'), 1)


def compute_curvature():
    with rasterio.open(DEM_25M) as src:
        z = src.read(1).astype('float64')
        profile = src.profile.copy()
        cellsize = src.res[0]
        nodata = src.nodata
    mask_arr = (z == nodata) if nodata is not None else None
    if mask_arr is not None:
        z = np.where(mask_arr, np.nan, z)
    # Discrete Laplacian (general curvature): (Z_left+Z_right+Z_up+Z_down-4*Z_center)/cellsize^2
    kernel = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]], dtype='float64')
    curv = ndimage.convolve(np.nan_to_num(z, nan=0.0), kernel, mode='nearest') / (cellsize ** 2)
    if mask_arr is not None:
        curv = np.where(mask_arr, np.nan, curv)
    out_nodata = -9999.0
    curv = np.where(np.isnan(curv), out_nodata, curv).astype('float32')
    profile.update(dtype='float32', nodata=out_nodata)
    with rasterio.open(CURV_25M, 'w', **profile) as dst:
        dst.write(curv, 1)


def compute_twi():
    grid = Grid.from_raster(DEM_25M)
    dem = grid.read_raster(DEM_25M)
    flooded = grid.fill_depressions(dem)
    inflated = grid.resolve_flats(flooded)
    fdir = grid.flowdir(inflated)
    acc = grid.accumulation(fdir)
    cellsize = grid.affine[0]
    contrib_area = (np.asarray(acc) + 1) * (cellsize ** 2)
    z = np.asarray(inflated, dtype='float64')
    dzdy, dzdx = np.gradient(z, cellsize)
    slope_rad = np.arctan(np.sqrt(dzdx ** 2 + dzdy ** 2))
    tan_slope = np.tan(slope_rad)
    tan_slope = np.where(tan_slope < 0.001, 0.001, tan_slope)
    twi = np.log(contrib_area / tan_slope).astype('float32')
    with rasterio.open(DEM_25M) as src:
        profile = src.profile.copy()
    profile.update(dtype='float32', nodata=-9999.0)
    with rasterio.open(TWI_25M, 'w', **profile) as dst:
        dst.write(twi, 1)


def zonal_mean(raster_path, geoms):
    means = []
    with rasterio.open(raster_path) as src:
        nodata = src.nodata
        for geom in geoms:
            try:
                out_image, _ = rio_mask(src, [geom], crop=True, nodata=nodata)
                arr = out_image[0]
                valid = arr[arr != nodata]
                means.append(float(np.mean(valid)) if valid.size > 0 else np.nan)
            except Exception:
                means.append(np.nan)
    return means


def build_dataset():
    cat = gpd.read_file(CATCHMENTS_PATH)
    cat['TWI'] = zonal_mean(TWI_25M, cat.geometry)
    cat['curvature'] = zonal_mean(CURV_25M, cat.geometry)

    dn = gpd.read_file(DRAINAGE_NETWORK_PATH, layer='drenaje')
    dd = dn.groupby('ID_CUENCA')['length_m'].sum().reset_index()
    cat = cat.merge(dd, on='ID_CUENCA', how='left')
    cat['length_m'] = cat['length_m'].fillna(0)
    cat['drainage_density'] = cat['length_m'] / cat['area'] * 1000  # km/km^2
    return cat


def correlation_and_figure(cat):
    cols = RETAINED + NEW_CANDIDATES + ORIGINAL_SCREENED
    corr = cat[cols].corr(method='pearson')

    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(corr.values, vmin=-1, vmax=1, cmap='RdBu_r')
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=45, ha='right')
    ax.set_yticks(range(len(cols)))
    ax.set_yticklabels(cols)
    for i in range(len(cols)):
        for j in range(len(cols)):
            ax.text(j, i, f"{corr.values[i, j]:.2f}", ha='center', va='center', fontsize=8)
    fig.colorbar(im, ax=ax, label="Pearson r")
    ax.set_title("Table S1. Correlation matrix: retained predictors, new candidates\n"
                  "(TWI, curvature, drainage density), and originally screened variables")
    fig.tight_layout()

    out_fig = r"RESULTS/TableS1_correlation_matrix.png"
    fig.savefig(out_fig, dpi=300)

    out_csv = r"RESULTS/TableS1_correlation_matrix.csv"
    corr.to_csv(out_csv)

    print("Correlations exceeding |r| >", R_THRESHOLD, "with a retained predictor:")
    for c in NEW_CANDIDATES + ORIGINAL_SCREENED:
        for r in RETAINED:
            v = corr.loc[c, r]
            if abs(v) > R_THRESHOLD:
                print(f"  {c} vs {r}: r={v:.3f}")

    print("saved", out_fig)
    print("saved", out_csv)
    return corr


if __name__ == "__main__":
    print("1/5 resampling DEM to", RESAMPLED_RES, "m ...")
    resample_dem()
    print("2/5 computing curvature ...")
    compute_curvature()
    print("3/5 computing TWI (fill depressions, flow accumulation) ...")
    compute_twi()
    print("4/5 building catchment dataset (zonal stats + drainage density) ...")
    cat = build_dataset()
    cat.drop(columns='geometry').to_csv(WORKDIR + "/catchments_with_new_predictors.csv", index=False)
    print("5/5 correlation matrix and Table S1 figure ...")
    correlation_and_figure(cat)
    print("done")
