################################################################
# Nonlinearity check for RC1 response (Major Concern 4)
# Tests whether the log-linear Poisson assumption for the fixed
# effects (in particular mean slope) is supported, by comparing
# a linear model against models with a quadratic term and a
# penalized spline (GAM) term for each predictor, and by
# producing partial-residual plots.
################################################################

import geopandas as gpd
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.gam.api import GLMGam, BSplines
import matplotlib.pyplot as plt

CATCHMENTS_PATH = r"DATA/df_catchments_kmeans.gpkg"
OUT_DIR = r"RESULTS"

cat = gpd.read_file(CATCHMENTS_PATH)
df = pd.DataFrame(cat.drop(columns='geometry'))
df['log_area'] = np.log(df['area'])

# standardize predictors as in the manuscript's M2-M5 models
for col in ['RainfallDaysmean', 'elev_mean', 'slope_mean']:
    df[col + '_z'] = (df[col] - df[col].mean()) / df[col].std()

# --- 1. Linear (baseline) Poisson GLM ---
lin_formula = "lands_rec ~ RainfallDaysmean_z + elev_mean_z + slope_mean_z"
lin_model = smf.glm(lin_formula, data=df, family=sm.families.Poisson(),
                     offset=df['log_area']).fit()

# --- 2. Quadratic terms added one predictor at a time ---
results = {"linear": lin_model.aic}
for var in ['RainfallDaysmean_z', 'elev_mean_z', 'slope_mean_z']:
    df[var + '_sq'] = df[var] ** 2
    quad_formula = lin_formula + f" + {var}_sq"
    quad_model = smf.glm(quad_formula, data=df, family=sm.families.Poisson(),
                          offset=df['log_area']).fit()
    lr_stat = 2 * (quad_model.llf - lin_model.llf)
    from scipy.stats import chi2
    p_value = chi2.sf(lr_stat, df=1)
    results[f"+{var}^2"] = quad_model.aic
    print(f"Quadratic term for {var}: coef={quad_model.params[var+'_sq']:.4f}, "
          f"p={quad_model.pvalues[var+'_sq']:.4g}, LR p={p_value:.4g}, "
          f"AIC {lin_model.aic:.1f} -> {quad_model.aic:.1f}")

# --- 3. GAM with penalized spline for slope (the variable RC1 flagged) ---
x_spline = df[['slope_mean_z']].values
bs = BSplines(x_spline, df=[6], degree=[3])
gam_model = GLMGam.from_formula(
    "lands_rec ~ RainfallDaysmean_z + elev_mean_z",
    data=df, smoother=bs, family=sm.families.Poisson(), offset=df['log_area']
).fit()
print("\nGAM (spline on slope) AIC:", gam_model.aic, "vs linear AIC:", lin_model.aic)

# --- 4. Partial residual plot for slope ---
eta_others = (lin_model.params['Intercept']
              + lin_model.params['RainfallDaysmean_z'] * df['RainfallDaysmean_z']
              + lin_model.params['elev_mean_z'] * df['elev_mean_z']
              + df['log_area'])
partial_resid = (df['lands_rec'] - np.exp(eta_others)) / np.exp(eta_others) \
                 + lin_model.params['slope_mean_z'] * df['slope_mean_z']

fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(df['slope_mean_z'], partial_resid, s=15, alpha=0.5)
order = np.argsort(df['slope_mean_z'].values)
lowess_x = df['slope_mean_z'].values[order]
from statsmodels.nonparametric.smoothers_lowess import lowess as lowess_fn
sm_line = lowess_fn(partial_resid.values[order], lowess_x, frac=0.5)
ax.plot(sm_line[:, 0], sm_line[:, 1], color='red', lw=2, label='LOWESS smooth')
ax.plot(lowess_x, lin_model.params['slope_mean_z'] * lowess_x, color='black',
        linestyle='--', label='Linear fit')
ax.set_xlabel('Standardized mean slope')
ax.set_ylabel('Partial residual (working scale)')
ax.legend()
ax.set_title('Partial-residual plot: mean slope (M1 baseline)')
fig.tight_layout()
fig_path = OUT_DIR + "/FigS2_slope_partial_residual.png"
fig.savefig(fig_path, dpi=300)
print("saved", fig_path)

print("\nAIC summary:", results)
