library(ggplot2)
library(ggrepel)
library(here)
library(tidyverse)

load(here("res", "sims", "alpha", "100alpha_density_results.Rdata"))

dens_plot <- function(alpha_dens) {
  # Ensure threshold is a factor (discrete)
  alpha_dens$threshold_f <- factor(alpha_dens$threshold, levels = sort(unique(alpha_dens$threshold)))

  # Create 10^x expression labels for the legend
  threshold_labels <- parse(text = paste0("10^", log10(as.numeric(levels(alpha_dens$threshold_f)))))

  ggplot(alpha_dens, aes(
    x = alpha,
    y = mean_density,
    group = threshold_f,
    color = threshold_f
  )) +
    geom_line(size = 0.5, alpha = 0.6) +
    geom_point(size = 2, alpha = 1) +
    scale_color_viridis_d(
      option = "B",
      begin = 0.2,
      end = 0.8,
      direction = -1,
      labels = threshold_labels,
      guide = guide_legend(nrow = 1)
    ) +
    scale_x_continuous(
      breaks = sort(unique(alpha_dens$alpha)),
      labels = function(x) ifelse(abs(x - round(x)) < .Machine$double.eps^0.5, as.character(round(x)), as.character(x)),
      minor_breaks = NULL
    ) +
    scale_y_continuous(
      breaks = c(0, 0.2, 0.4, 0.6, 0.8),
      limit = c(0, NA)
    ) +
    labs(
      x     = expression(alpha),
      y     = expression(density),
      color = "Networks size"
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

alpha_plot <- dens_plot(alpha_dens100)

ggsave(
  filename = "res/pics/density_analysis/alpha/100net_alpha_density_DP.pdf",
  plot = alpha_plot,
  width = 7,
  height = 4,
  bg = "white"
)


label_fun <- function(x) {
  # Check if the number is effectively an integer, accounting for floating point issues
  is_integer <- abs(x - round(x)) < .Machine$double.eps^0.5
  ifelse(is_integer, as.character(round(x)), as.character(x))
}

# Create the data for the labels at the end of the lines
label_data <- alpha_dens100 %>%
  filter(alpha == max(alpha)) %>%
  mutate(
    network_label = paste0("10^", log10(threshold))
  )

# Now build the plot
log_alpha_plot <- ggplot(alpha_dens100, aes(x = alpha, y = mean_density, color = factor(threshold))) +
  geom_line(aes(group = threshold), alpha = 0.4, linewidth = 0.5) +
  geom_point(size = 2.5) +

  # --- THE MODIFIED PART ---
  # Use the new custom function for the labels argument
  scale_x_log10(
    breaks = unique(alpha_dens100$alpha),
    labels = label_fun
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +

  # Use ggrepel to add non-overlapping labels
  geom_text_repel(
    data = label_data,
    aes(label = network_label),
    parse = TRUE,
    nudge_x = 0.15,
    hjust = 0,
    direction = "y",
    segment.color = NA,
    size = 4,
    inherit.aes = TRUE
  ) +

  # A good color scheme for sequential data
  scale_color_viridis_d(
    option = "plasma", end = 0.9,
    direction = -1
  ) +

  # Final labels and theme
  labs(
    x = expression(alpha),
    y = "density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "res/pics/density_analysis/alpha/100net_log_alpha_density_DP.pdf",
  plot = log_alpha_plot,
  width = 7,
  height = 4,
  bg = "white"
)
