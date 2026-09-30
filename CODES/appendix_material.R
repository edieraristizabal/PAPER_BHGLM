################################################################
## Appendix A of the revised manuscript: figures and LaTeX tables
## generated from the archived analysis outputs, so that no number
## in the appendix is transcribed by hand.
##   Fig. A1  spatial support (descriptive)          -> new
##   Tables A1-A6 -> RESULTS/appendix_tables/tabA*.tex (inlined in Appendix A of the manuscript)
## Run from the PAPER_BHGLM root after the EC4/EC7 scripts.
################################################################

suppressMessages({ library(ggplot2); library(patchwork); library(colorspace) })

REV <- "RESULTS"
FIG <- "FIGURES"
OUT <- file.path(REV, "appendix_tables")
dir.create(OUT, showWarnings = FALSE)
fmt  <- function(x, d = 2) formatC(x, format = "f", digits = d, big.mark = "\\,")
fmtn <- function(x) formatC(round(x), format = "d", big.mark = "\\,")
write_tex <- function(lines, f) writeLines(lines, file.path(OUT, f), useBytes = TRUE)

## ---------------------------------------------------------------
## Fig. A1 and Table A1: spatial support
## ---------------------------------------------------------------
ter <- read.csv(file.path(REV, "EC4/ec4_catchment_terrain_12m.csv"))
cat <- read.csv("DATA/catchments_with_new_predictors.csv")
stopifnot(all(cat$id == ter$id))
ter$area_km2  <- cat$area / 1e6
ter$lands     <- cat$lands_rec
ter$rain_days <- cat$RainfallDaysmean

pal <- qualitative_hcl(2, palette = "Dark 3")
th  <- theme_classic(base_size = 9) + theme(plot.tag = element_text(face = "bold"))
h   <- function(df, x, xlab, bins = 30, logx = FALSE) {
  p <- ggplot(df, aes(.data[[x]])) +
    geom_histogram(bins = bins, fill = pal[1], colour = "white", linewidth = 0.15) +
    labs(x = xlab, y = "Catchments") + th
  if (logx) p <- p + scale_x_log10()
  p
}
f_a1 <- (h(ter, "area_km2", expression("Catchment area (km"^2*")"), logx = TRUE) |
         h(ter, "lands", "Mapped landslides per catchment", bins = 40)) /
        (h(ter, "slope_sd", "Within-catchment SD of slope (°)") |
         h(ter, "elev_sd", "Within-catchment SD of elevation (m)")) +
        plot_annotation(tag_levels = "a")
ggsave(file.path(FIG, "FigA1_spatial_support.png"), f_a1,
       width = 170, height = 120, units = "mm", dpi = 300)

row_q <- function(lbl, x, d, dm = d) {
  q <- quantile(x, c(0, .25, .5, .75, 1))
  v <- c(fmt(q[1:3], d), fmt(mean(x), dm), fmt(q[4:5], d))
  paste(lbl, "&", paste(v, collapse = " & "), "\\\\")
}
write_tex(c(
  "\\begin{table}[h]",
  "\\caption{Descriptive statistics of the 526 catchments. Terrain statistics computed from the 12.5~m ALOS-PALSAR DEM; within-catchment standard deviations (SD) quantify the variability that aggregation to catchment means discards.}\\label{tabA1}",
  "\\footnotesize",
  "\\begin{tabular}{@{}lrrrrrr@{}}",
  "\\toprule",
  "Variable & Min & Q1 & Median & Mean & Q3 & Max \\\\",
  "\\midrule",
  row_q("Catchment area (km$^2$)", ter$area_km2, 1),
  row_q("Mapped landslides", ter$lands, 0, dm = 1),
  row_q("Rainfall days, $R_d$", ter$rain_days, 0),
  row_q("Mean elevation, $H$ (m)", ter$elev_mean, 0),
  row_q("Within-catchment SD of elevation (m)", ter$elev_sd, 0),
  row_q("Mean slope, $S$ ($^{\\circ}$)", ter$slope_mean, 1),
  row_q("Within-catchment SD of slope ($^{\\circ}$)", ter$slope_sd, 1),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), "tabA1.tex")

## ---------------------------------------------------------------
## Table A2: correlation matrix of the candidate predictors (lower triangle)
## ---------------------------------------------------------------
cm <- as.matrix(read.csv(file.path(REV, "TableS1_correlation_matrix.csv"), row.names = 1))
lab <- c(RainfallDaysmean = "$R_d$", elev_mean = "$H$", slope_mean = "$S$", TWI = "TWI",
         curvature = "Curv.", drainage_density = "DD", area = "Area",
         hypso_inte = "HI", Densidad = "LD", rel_mean = "Rel.", rainfallAnnual_mean = "MAR")
k <- nrow(cm)
body <- sapply(seq_len(k), function(i) {
  v <- ifelse(seq_len(k) < i, fmt(cm[i, ], 2), ifelse(seq_len(k) == i, "1", ""))
  v <- ifelse(seq_len(k) < i & abs(cm[i, ]) > 0.75, paste0("\\textbf{", v, "}"), v)
  paste(lab[rownames(cm)[i]], "&", paste(v, collapse = " & "), "\\\\")
})
write_tex(c(
  "\\begin{table*}[h]",
  "\\caption{Pearson correlation matrix of the candidate catchment-level predictors. Values with $|r| > 0.75$, the multicollinearity threshold, are in bold. $R_d$: days with rainfall $>$20~mm; $H$: mean elevation; $S$: mean slope; TWI: topographic wetness index; Curv.: curvature index; DD: drainage density; HI: hypsometric integral; LD: lineament density; Rel.: mean relief; MAR: mean annual rainfall.}\\label{tabA2}",
  "\\footnotesize",
  paste0("\\begin{tabular}{@{}l", strrep("r", k), "@{}}"),
  "\\toprule",
  paste("&", paste(lab[colnames(cm)], collapse = " & "), "\\\\"),
  "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\end{table*}"), "tabA2.tex")

## ---------------------------------------------------------------
## Table A3: posterior summaries of the fixed effects
## ---------------------------------------------------------------
fx <- read.csv(file.path(REV, "EC7/ec7_fixed_effects_all.csv"))   # single archived run (EC9)
fx <- fx[fx$model %in% paste0("M", 1:5), ]
term_lab <- c("(Intercept)" = "Intercept", RainfallDaysmean = "Rainfall days",
              elev_mean = "Mean elevation", slope_mean = "Mean slope")
rows <- unlist(lapply(names(term_lab), function(t) {
  s <- fx[fx$term == t, ]; s <- s[match(paste0("M", 1:5), s$model), ]
  c(paste0(term_lab[[t]], " & ", paste(sprintf("%s (%s)", fmt(s$mean, 3), fmt(s$sd, 3)), collapse = " & "), " \\\\"),
    paste0(" & ", paste(sprintf("[%s,\\, %s]", fmt(s$q025, 3), fmt(s$q975, 3)), collapse = " & "), " \\\\"))
}))
write_tex(c(
  "\\begin{table*}[h]",
  "\\caption{Posterior mean (posterior standard deviation) and 95\\% credible interval of the fixed effects in M1--M5, on standardized predictors (mean 0, SD 1).}\\label{tabA3}",
  "\\footnotesize",
  "\\begin{tabular}{@{}lccccc@{}}",
  "\\toprule",
  "Term & M1 & M2 & M3-ICAR & M4-BYM & M5-Leroux \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}", "\\end{table*}"), "tabA3.tex")

## ---------------------------------------------------------------
## Table A4: spatially blocked cross-validation, RMSE per fold
## ---------------------------------------------------------------
cv <- read.csv(file.path(REV, "spatial_cv_results.csv"))
fr <- function(x) ifelse(x >= 1e4, sprintf("$%.1f \\times 10^{%d}$", x / 10^floor(log10(x)), floor(log10(x))), fmt(x, 1))
mods <- paste0("M", 1:5)
cv_rows <- sapply(sort(unique(cv$fold)), function(f) {
  s <- cv[cv$fold == f, ]; s <- s[match(mods, s$model), ]
  paste(f, "&", s$n_test[1], "&", paste(fr(s$rmse), collapse = " & "), "\\\\")
})
med <- sapply(mods, function(m) median(cv$rmse[cv$model == m]))
write_tex(c(
  "\\begin{table}[h]",
  "\\caption{Spatially blocked 5-fold cross-validation: root-mean-square error (RMSE) of the predicted landslide counts in the held-out catchments of each fold. Folds are geographically compact groups of catchments obtained by $k$-means clustering of catchment centroids.}\\label{tabA4}",
  "\\footnotesize",
  "\\begin{tabular}{@{}ccrrrrr@{}}",
  "\\toprule",
  "Fold & $n$ & M1 & M2 & M3-ICAR & M4-BYM & M5-Leroux \\\\",
  "\\midrule", cv_rows, "\\midrule",
  paste("Median & &", paste(fr(med), collapse = " & "), "\\\\"),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), "tabA4.tex")

## ---------------------------------------------------------------
## Table A5: alternative likelihoods (EC7)
## ---------------------------------------------------------------
mc <- read.csv(file.path(REV, "EC7/ec7b_model_comparison.csv"))
sel <- c(M1 = "M1, Poisson", M1_nb = "M1, negative binomial",
         M1_zip = "M1, zero-inflated Poisson", M1_zinb = "M1, zero-inflated neg.\\ binomial",
         M1_hurdle = "M1, hurdle Poisson",
         M4 = "M4-BYM, Poisson", M4_nb = "M4-BYM, negative binomial",
         M5 = "M5-Leroux, Poisson", M5_zip = "M5-Leroux, zero-inflated Poisson",
         M5_nb = "M5-Leroux, negative binomial", M5_zinb = "M5-Leroux, zero-inflated neg.\\ binomial",
         M5_hurdle = "M5-Leroux, hurdle neg.\\ binomial")
s5 <- mc[match(names(sel), mc$model), ]
write_tex(c(
  "\\begin{table}[h]",
  "\\caption{Model comparison under alternative count likelihoods. All specifications share the predictors, Queen adjacency and $\\log(A_i)$ offset of the main models. MI: residual Moran's I (Eq.~\\ref{eq:moran}); max $|r|$: largest absolute Pearson residual.}\\label{tabA5}",
  "\\footnotesize",
  "\\begin{tabular}{@{}lrrrr@{}}",
  "\\toprule",
  "Specification & DIC & WAIC & MI & max $|r|$ \\\\",
  "\\midrule",
  paste(sel, "&", fmtn(s5$DIC), "&", fmtn(s5$WAIC), "&", fmt(s5$MI_res, 3), "&", fmt(s5$max_abs_pearson, 1), "\\\\"),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), "tabA5.tex")

## ---------------------------------------------------------------
## Table A6: aggregation and slope-scale sensitivity (EC4)
## ---------------------------------------------------------------
sa <- read.csv(file.path(REV, "EC4/ec4_sensitivity_summary.csv"))
fa <- read.csv(file.path(REV, "EC4/ec4_sensitivity_fixed.csv"))
la <- read.csv(file.path(REV, "EC4/ec4_latent_field_similarity.csv"))
sc <- read.csv(file.path(REV, "EC4/ec4C_refit_summary.csv"))
fc <- read.csv(file.path(REV, "EC4/ec4C_refit_fixed.csv"))
slope_terms <- c("slope_mean", "slope_p90", "frac_gt20", "frac_gt25", "frac_gt30", "S_c")
row6 <- function(lbl, S, Fx, v, corr) {
  m5 <- S[S$variant == v & S$model == "M5", ]
  m4 <- S[S$variant == v & S$model == "M4", ]
  m1 <- S[S$variant == v & S$model == "M1", ]
  b  <- Fx[Fx$variant == v & Fx$model == "M5" & Fx$term %in% slope_terms, ]
  bs <- if (nrow(b) == 1) sprintf("%s [%s, %s]", fmt(b$mean), fmt(b$q025), fmt(b$q975)) else "--"
  paste(lbl, "&", fmtn(m1$DIC), "&", fmtn(m4$DIC), "&", fmtn(m5$DIC), "&", fmt(m5$MI_res, 3), "&",
        fmt(m5$var_id), "&", fmt(m5$rho, 3), "&", fmt(corr, 3), "&", bs, "\\\\")
}
cr <- function(v) la$cor_with_published[la$variant == v]
cs <- function(v) sc$cor_M5_field_with_published[sc$variant == v & sc$model == "M5"]
write_tex(c(
  "\\begin{table*}[h]",
  "\\caption{Sensitivity of the models to the aggregation of the slope predictor. Upper block: the catchment-mean slope is replaced by alternative within-catchment summaries computed from the 12.5~m DEM. Lower block: slope enters through $S_c$, the catchment integral of the within-catchment slope response (Sect.~\\ref{secA6}). $\\sigma^2_u$ and $\\rho$: variance and mixing parameter of the Leroux latent field; $r_u$: correlation of the M5 latent field with that of the published specification; $\\beta_S$: standardized slope coefficient in M5 with 95\\% credible interval.}\\label{tabA6}",
  "\\footnotesize",
  "\\begin{tabular}{@{}lrrrrrrrl@{}}",
  "\\toprule",
  " & M1 & M4 & \\multicolumn{6}{c}{M5-Leroux} \\\\",
  "\\cmidrule(lr){4-9}",
  "Slope predictor & DIC & DIC & DIC & MI & $\\sigma^2_u$ & $\\rho$ & $r_u$ & $\\beta_S$ \\\\",
  "\\midrule",
  row6("Mean slope (published)", sa, fa, "mean (published)", 1),
  row6("90th percentile", sa, fa, "p90", cr("p90")),
  row6("Area fraction $>20^{\\circ}$", sa, fa, "area > 20 deg", cr("area > 20 deg")),
  row6("Area fraction $>25^{\\circ}$", sa, fa, "area > 25 deg", cr("area > 25 deg")),
  row6("Area fraction $>30^{\\circ}$", sa, fa, "area > 30 deg", cr("area > 30 deg")),
  row6("Mean + within-catchment SD", sa, fa, "mean + within-SD", cr("mean + within-SD")),
  "\\midrule",
  row6("$S_c$ (slope-scale integrated)", sc, fc, "S_c (slope-scale integrated)", cs("S_c (slope-scale integrated)")),
  row6("$S_c$ as fixed offset", sc, fc, "S_c as fixed offset", cs("S_c as fixed offset")),
  "\\bottomrule", "\\end{tabular}", "\\end{table*}"), "tabA6.tex")

## ---------------------------------------------------------------
## Table A2b and Fig. A2: PCA of the candidate predictors
## ---------------------------------------------------------------
pv <- read.csv(file.path(REV, "LB/pca_variance.csv"))
pl <- read.csv(file.path(REV, "LB/pca_loadings.csv"))
ps <- read.csv(file.path(REV, "LB/pca_scores.csv"))
write_tex(c(
  "\\begin{table}[h]",
  "\\caption{Principal component analysis of the eight continuous candidate predictors of the original screening (standardized): loadings on the first three components and share of variance explained.}\\label{tabA2pca}",
  "\\footnotesize",
  "\\begin{tabular}{@{}lrrr@{}}",
  "\\toprule",
  "Variable & PC1 & PC2 & PC3 \\\\",
  "\\midrule",
  paste(gsub("R_d", "$R_d$", pl$variable), "&", fmt(pl$PC1), "&", fmt(pl$PC2), "&", fmt(pl$PC3), "\\\\"),
  "\\midrule",
  paste("Variance explained (\\%) &", paste(fmt(100 * pv$explained[1:3], 1), collapse = " & "), "\\\\"),
  paste("Cumulative (\\%) &", paste(fmt(100 * pv$cumulative[1:3], 1), collapse = " & "), "\\\\"),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), "tabA2pca.tex")

sc_ar <- 3.2
arr <- data.frame(v = pl$variable, x = pl$PC1 * sc_ar, y = pl$PC2 * sc_ar)
arr$lab <- sub(" \\(.*\\)", "", arr$v)
pA2 <- ggplot(ps, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(aes(colour = basin, size = lands), alpha = 0.55, shape = 16) +
  geom_segment(data = arr, aes(x = 0, y = 0, xend = x, yend = y), inherit.aes = FALSE,
               arrow = arrow(length = unit(1.8, "mm")), linewidth = 0.4, colour = "grey15") +
  ggrepel::geom_label_repel(data = arr, aes(x, y, label = lab), inherit.aes = FALSE, size = 2.5,
             linewidth = 0, fill = alpha("white", 0.8), label.padding = unit(0.1, "lines"),
             min.segment.length = 0, segment.colour = "grey40", segment.size = 0.25, box.padding = 0.35, seed = 1) +
  scale_colour_manual(values = setNames(qualitative_hcl(3, palette = "Dark 3"), c("Atrato", "Cauca", "Magdalena")), name = "Basin") +
  scale_size_area(max_size = 5, name = "Mapped landslides", breaks = c(10, 100, 300)) +
  labs(x = sprintf("PC1 (%.1f%%)", 100 * pv$explained[1]), y = sprintf("PC2 (%.1f%%)", 100 * pv$explained[2])) +
  coord_equal() + theme_bw(base_size = 9) + theme(panel.grid = element_blank())
ggsave(file.path(FIG, "FigA2_pca_biplot.png"), pA2, width = 170, height = 120, units = "mm", dpi = 400, bg = "white")

## ---------------------------------------------------------------
## Table A7: model-specification checks (second editor round)
## ---------------------------------------------------------------
lb <- read.csv(file.path(REV, "LB/lb_model_summary.csv"))
lab_m <- c(M1 = "M1", M2 = "M2", M2iid = "M2 + catchment iid (non-spatial)", M3 = "M3-ICAR", M4 = "M4-BYM", M5 = "M5-Leroux",
           BYM2 = "BYM2 (PC priors)", M5_pcprior = "M5, PC prior on precision",
           M4_nobasin = "M4, no basin effect", M5_nobasin = "M5, no basin effect",
           M4_basinfixed = "M4, basin as fixed effect", M5_basinfixed = "M5, basin as fixed effect",
           M4_lith_lc = "M4 + lithology + land cover", M5_lith_lc = "M5 + lithology + land cover",
           M4_rw2 = "M4, rw2 slope and elevation", M5_rw2 = "M5, rw2 slope and elevation",
           M5_noslope = "M5 without slope", M5_permslope = "M5, permuted slope")
lb <- lb[match(names(lab_m), lb$model), ]
fp <- function(p) ifelse(p <= 0.002, "$\\leq$0.002", fmt(p, 3))
mix <- ifelse(!is.na(lb$phi_bym2), sprintf("$\\phi$ = %s", fmt(lb$phi_bym2, 2)),
              ifelse(!is.na(lb$rho_leroux), sprintf("$\\rho$ = %s", fmt(lb$rho_leroux, 3)), "--"))
ru <- ifelse(!is.na(lb$r_field_vs_M5), fmt(lb$r_field_vs_M5, 3), ifelse(!is.na(lb$r_field_vs_M4), fmt(lb$r_field_vs_M4, 3), "--"))
sdl <- ifelse(is.na(lb$sd_latent_field), "--", fmt(lb$sd_latent_field, 2))
write_tex(c(
  "\\begin{table*}[h]",
  "\\caption{Model-specification checks. DIC is not interpretable for M1 and M2 because their effective number of parameters $p_D$ is negative. MI: residual Moran's I with two-sided Monte Carlo $p$-value (999 permutations). $\\rho$: Leroux mixing parameter; $\\phi$: BYM2 share of the marginal latent variance that is spatially structured. SD: standard deviation of the posterior mean latent field. $r_u$: correlation of the latent field with that of M5 (M5 variants) or M4 (M4 variants).}\\label{tabA7}",
  "\\footnotesize",
  "\\begin{tabular}{@{}lrrrrrrlrr@{}}",
  "\\toprule",
  "Model & DIC & $p_D$ & WAIC & $p_{\\mathrm{WAIC}}$ & log mlik & MI ($p$) & Mixing & SD & $r_u$ \\\\",
  "\\midrule",
  paste(lab_m, "&", fmtn(lb$DIC), "&", fmtn(lb$pD), "&", fmtn(lb$WAIC), "&", fmtn(lb$pWAIC), "&", fmtn(lb$mlik), "&",
        sprintf("%s (%s)", fmt(lb$MI, 3), fp(lb$MI_p_two_sided)), "&", mix, "&", sdl, "&", ru, "\\\\"),
  "\\bottomrule", "\\end{tabular}", "\\end{table*}"), "tabA7.tex")

## ---------------------------------------------------------------
## Table A4b: cross-validation with predictive scores (second round)
## ---------------------------------------------------------------
if (file.exists(file.path(REV, "LB/cv_summary.csv"))) {
  cvs <- read.csv(file.path(REV, "LB/cv_summary.csv"))
  mods <- c("M1", "M2", "M2iid", "M3", "M4", "M5", "BYM2")
  labs <- c("M1", "M2", "M2 + catchment iid", "M3-ICAR", "M4-BYM", "M5-Leroux", "BYM2")
  g <- function(tp, m, v) { x <- cvs[cvs$type == tp & cvs$model == m, v]; if (length(x)) x else NA }
  fr1 <- function(x) ifelse(is.na(x), "--", ifelse(x >= 1e4, sprintf("$%.1f \\times 10^{%d}$", x / 10^floor(log10(x)), floor(log10(x))), fmt(x, 1)))
  rowsb <- sapply(seq_along(mods), function(i) paste(labs[i],
    "&", fmt(g("spatial", mods[i], "log_score"), 2), "&", fr1(g("spatial", mods[i], "crps")), "&", fr1(g("spatial", mods[i], "rmse_median")),
    "&", fmt(g("random", mods[i], "log_score"), 2), "&", fr1(g("random", mods[i], "crps")), "&", fr1(g("random", mods[i], "rmse_median")), "\\\\"))
  write_tex(c(
    "\\begin{table*}[h]",
    "\\caption{Cross-validation with predictive scores. Spatially blocked: three $k$-means partitions of catchment centroids ($k = 5$, different seeds; 1578 held-out predictions per model); random: one random 5-fold partition (526 predictions). Scores are computed from 300 posterior samples of the linear predictor of each held-out catchment: log score (negative mean log predictive density) and continuous ranked probability score (CRPS), both negatively oriented, and the root-mean-square error of the predictive median.}\\label{tabA4b}",
    "\\footnotesize",
    "\\begin{tabular}{@{}lrrrrrr@{}}",
    "\\toprule",
    " & \\multicolumn{3}{c}{Spatially blocked} & \\multicolumn{3}{c}{Random} \\\\",
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
    "Model & Log score & CRPS & RMSE & Log score & CRPS & RMSE \\\\",
    "\\midrule", rowsb, "\\bottomrule", "\\end{tabular}", "\\end{table*}"), "tabA4b.tex")
}

## Figs. A3-A5 are written directly by figures_main.R, ec7_followup.R and ec4_figure.R
cat("appendix material written\n")
