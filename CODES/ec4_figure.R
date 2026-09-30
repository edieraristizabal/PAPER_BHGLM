################################################################
## EC4 -- Supplementary figure: slope-scale control vs catchment
## aggregation. (a) Within-catchment slope response f(s) against the
## pooled frequency ratio; (b) percentile of landslide-crown slope in
## the slope distribution of its own catchment.
## Reads outputs of CODES/ec4_slope_scale.R. Run from PAPER_BHGLM root.
################################################################

suppressMessages({ library(ggplot2); library(patchwork); library(colorspace) })

OUT <- "RESULTS/EC4"
cf  <- read.csv(file.path(OUT, "ec4C_slope_response.csv"))
pts <- read.csv(file.path(OUT, "ec4B_crown_points.csv"))
B   <- read.csv(file.path(OUT, "ec4B_crown_percentiles.csv"))

pal <- qualitative_hcl(2, palette = "Dark 3")
th  <- theme_classic(base_size = 10) +
  theme(legend.position = c(0.03, 0.97), legend.justification = c(0, 1),
        legend.background = element_blank(), legend.title = element_blank(),
        plot.tag = element_text(face = "bold"))

cf <- cf[cf$slope <= 60, ]
pa <- ggplot(cf, aes(slope)) +
  geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)), fill = pal[1], alpha = 0.25) +
  geom_line(aes(y = exp(f), colour = "Within-catchment response"),
            linewidth = 0.9) +
  geom_line(aes(y = exp(pooled_log_FR), colour = "Pooled frequency ratio"),
            linewidth = 0.6, linetype = "22") +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3) +
  scale_colour_manual(values = c(pal[2], pal[1])) +
  scale_y_log10() +
  labs(x = "Pixel slope (°)", y = "Relative landslide intensity") + th +
  theme(legend.position = c(0.98, 0.03), legend.justification = c(1, 0))

pb <- ggplot(pts, aes(percentile)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, 0.05),
                 fill = pal[1], colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = 1, colour = "grey30", linetype = "22", linewidth = 0.5) +
  annotate("text", x = 0.02, y = 1.06, hjust = 0, vjust = 0, size = 3,
           label = "Uniform: no preference") +
  annotate("text", x = 0.02, y = 1.9, hjust = 0, vjust = 1, size = 3,
           label = sprintf("n = %s crowns\nmedian percentile = %.2f\nabove catchment P90: %.1f%%",
                           format(B$n_points, big.mark = ","), B$median_percentile,
                           100 * B$frac_above_catchment_p90)) +
  coord_cartesian(ylim = c(0, 1.9)) +
  labs(x = "Percentile of crown slope within its own catchment", y = "Density") + th

fig <- (pa | pb) + plot_annotation(tag_levels = "a")
ggsave("FIGURES/FigA5_slope_scale_control.png", fig,
       width = 180, height = 80, units = "mm", dpi = 300)
