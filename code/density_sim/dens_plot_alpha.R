library(ggplot2)
load(here("res", "sims", "alpha", "100alpha_density_results.Rdata"))
dens_plot = function(alpha_dens) {
  # Ensure threshold is a factor (discrete)
  alpha_dens$threshold_f <- factor(alpha_dens$threshold, levels = sort(unique(alpha_dens$threshold)))
  
  # Generate a visually progressive color palette (e.g., viridis)
  num_thresholds <- length(levels(alpha_dens$threshold_f))
  palette_colors <- viridis::viridis(num_thresholds)
  
  # Create 10^x expression labels for the legend
  threshold_labels <- parse(text = paste0("10^", log10(as.numeric(levels(alpha_dens$threshold_f)))))
  
  ggplot(alpha_dens, aes(x = alpha, y = mean_density, 
                         group = threshold_f, 
                         color = threshold_f)) +
    geom_line(size = 0.5, alpha = 0.6) +
    geom_point(size = 2, alpha = 1) +
    scale_color_manual(
      values = palette_colors,
      labels = threshold_labels,
      guide = guide_legend(nrow = 1)
    ) +
    scale_x_continuous(
      breaks = sort(unique(alpha_dens$alpha)),
      labels = function(x) ifelse(abs(x - round(x)) < .Machine$double.eps^0.5, as.character(round(x)), as.character(x)),
      minor_breaks = NULL
    ) +
    labs(
      x     = expression(alpha),
      y     = expression(E * "[" * d(H[n]) * "]"),
      color = "Network size"
    ) +
    theme_minimal() +
    theme(
      legend.position    = "top",
      legend.direction   = "horizontal",
      legend.box         = "horizontal",
      legend.text        = element_text(margin = margin(r = 10)),
      legend.key.width   = unit(1.5, "lines")
    )
}


dens_plot(alpha_dens100)

ggsave("res/pics/density_analysis/alpha/100net_alpha_density_DP.pdf",
       width = 8,
       height = 5,
       dpi = 600,
       bg = "white")
