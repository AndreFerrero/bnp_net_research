dens <- function(edges){
  # edges is alist, where every element is an edge data.frame
  
  N <- seq(50, 10000, by = 50)
  nsim <- length(edges)
  
  # Create a list to store results for the dynamic analysis.
  results <- list()
  
  for (i in 1:nsim) {
    # Get the full network from the simulation list
    full_net <- edges[[i]]
    
    # Loop over each Size
    for (n in N) {
      # Take the first 'n' edges
      sub_net <- full_net[1:n, ]
      
      # Ensure uniqueness
      sub_unique <- unique(sub_net)
      n_sub_unique <- nrow(sub_unique)
      
      # Count unique nodes in each partition from the subnetwork
      nA_sub <- length(unique(sub_unique$X_A))
      nB_sub <- length(unique(sub_unique$X_B))
      
      # Maximum possible edges for this subnetwork
      max_possible_sub <- nA_sub * nB_sub
      
      # Calculate density for this Size
      density_sub <- n_sub_unique / max_possible_sub
      
      # Save the results into a data frame
      results[[length(results) + 1]] <- data.frame(
        Simulation = i,
        Size = n,
        UniqueEdges = n_sub_unique,
        nA = nA_sub,
        nB = nB_sub,
        MaxPossible = max_possible_sub,
        Density = density_sub
      )
    }
  }
  
  # Combine dynamic results into one data frame
  res_df <- do.call(rbind, results)
  
  return(res_df)
}

# dens_PY25 <- dens(sim_edgesPY552525)
# dens_PY75 <- dens(sim_edgesPY557575)

library(here)

dens_folder <- here("res", "density_PY")

# save(dens_PY25, file = here(dens_folder, "dens_PY_0.25.Rdata"))
# save(dens_PY75, file = here(dens_folder, "dens_PY_0.75.Rdata"))

load(here(dens_folder, "dens_PY_0.25.Rdata"))
load(here(dens_folder, "dens_PY_0.5.Rdata"))
load(here(dens_folder, "dens_PY_0.75.Rdata"))

library(dplyr)
summarise_density <- function(data,
                              size_col  = "Size",
                              dens_col  = "Density",
                              ci_level  = 0.95,
                              log_x     = TRUE,
                              log_y     = TRUE) {
  data %>%
    mutate(
      x = if (log_x) log(.data[[size_col]]) else .data[[size_col]],
      y = if (log_y) log(.data[[dens_col]]) else .data[[dens_col]]
    ) %>%
    group_by(x) %>%
    summarise(
      mean    = mean(y),
      sd      = sd(y),
      n       = n(),
      se      = sd / sqrt(n),
      ci_low  = mean + qt((1 - ci_level)/2, df = n - 1) * se,
      ci_high = mean + qt(1 - (1 - ci_level)/2, df = n - 1) * se,
      .groups = "drop"
    )
}

dens_sum_PY25 <- summarise_density(dens_PY25, log_x = F,log_y = F)
dens_sum_PY5 <- summarise_density(dens_PY5, log_x = F,log_y = F)
dens_sum_PY75 <- summarise_density(dens_PY75, log_x = F,log_y = F)
loglog_dens_sum_PY25 <- summarise_density(dens_PY25, log_x = T,log_y = T)
loglog_dens_sum_PY5 <- summarise_density(dens_PY5, log_x = T,log_y = T)
loglog_dens_sum_PY75 <- summarise_density(dens_PY75, log_x = T,log_y = T)

####
dens_pics <- here("res", "pics", "density_analysis")

library(ggplot2)
# 
# ggplot(dens_PY25, aes(x = log(Size), y = log(Density),
#                       group = factor(Simulation),
#                       color = factor(Simulation))) +
#   geom_point(alpha = 0.4, size = 0.2) +
#   scale_color_viridis_d(option = "plasma", begin = 0, end = 1) +
#   labs(x = expression(log(n)),
#        y = expression(log(d(H[n]))),
#        color = "Simulation") +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "none",
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank()
#   )


# ggsave(filename = here(dens_pics, "PY_5__0.25__density.png"),
#        width = 11,
#        height = 8,
#        dpi = 600,
#        bg = "white")

# Mean and CI
plot_dens_sum <- function(data, 
                          xlab = expression(n), 
                          ylab = expression(E * "[" * d(H[n]) * "]"),
                          log_xlab = expression(log(n)),
                          log_ylab = expression(log(E * "[" * d(H[n]) * "]")),
                          loglog = FALSE,
                          add_hline0 = FALSE,
                          add_lm = FALSE,
                          title = NULL) {
  library(scales)
  
  # If loglog, use log axis labels
  if (loglog) {
    xlab <- log_xlab
    ylab <- log_ylab
  }
  
  # Base plot
  p <- ggplot(data, aes(x = x, y = mean)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
                fill = "steelblue", alpha = 0.2) +
    geom_line(color = "steelblue", size = 1) +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # Optional elements
  if (add_hline0) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  }
  
  if (add_lm) {
    p <- p + geom_smooth(method = "lm", se = FALSE, color = "darkred",
                         linetype = "dotted")
  }
  
  return(p)
}




dens_sum_plot_py25 <- plot_dens_sum(dens_sum_PY25, add_hline0 = T,
                                    title = expression(sigma == 0.25))
loglog_dens_sum_plot_py25 <- plot_dens_sum(loglog_dens_sum_PY25, loglog = T,
                                           add_lm = F)


# ggplot(dens_sum_PY25, aes(x = x, y = mean)) +
#   geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
#               fill = "steelblue", alpha = 0.2) +
#   geom_line(color = "steelblue", size = 1) +
#   labs(x = expression(n),
#        y = expression(E(d(H[n])))) +
#   theme_minimal(base_size = 14) +
#   theme(panel.grid.minor     = element_blank(),
#         panel.grid.major.x   = element_blank())
# 
# 
# 
# ggplot(loglog_dens_sum_PY25, aes(x = x, y = mean)) +
#   geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
#               fill = "steelblue", alpha = 0.2) +
#   geom_line(color = "steelblue", size = 1) +
#   labs(x = expression(log(n)),
#        y = expression(log(E(d(H[n]))))) +
#   theme_minimal(base_size = 14) +
#   theme(panel.grid.minor     = element_blank(),
#         panel.grid.major.x   = element_blank())

#####
dens_sum_plot_py5 <- plot_dens_sum(dens_sum_PY5, add_hline0 = T, ylab = NULL,
                                   title = expression(sigma == 0.5))
loglog_dens_sum_plot_py5 <- plot_dens_sum(loglog_dens_sum_PY5, loglog = T,
                                          add_lm = F, log_ylab = NULL)

#####
# ggplot(dens_PY75, aes(x = log(Size), y = log(Density),
#                       group = factor(Simulation),
#                       color = factor(Simulation))) +
#   geom_point(alpha = 0.4, size = 0.2) +
#   scale_color_viridis_d(option = "plasma", begin = 0, end = 1) +
#   labs(x = expression(log(n)),
#        y = expression(log(d(H[n]))),
#        color = "Simulation") +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "none",
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank()
#   )
# 
# ggsave(filename = here(dens_pics, "PY_5__0.75__density.png"),
#        width = 11,
#        height = 8,
#        dpi = 600,
#        bg = "white")


# Mean and CI

dens_sum_plot_py75 <- plot_dens_sum(dens_sum_PY75, add_hline0 = T,
                                           ylab = NULL,
                                    title = expression(sigma == 0.75))
loglog_dens_sum_plot_py75 <- plot_dens_sum(loglog_dens_sum_PY75, loglog = T,
                                           add_lm = F, log_ylab = NULL)

# ggplot(dens_sum_PY75, aes(x = x, y = mean)) +
#   geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
#               fill = "steelblue", alpha = 0.2) +
#   geom_line(color = "steelblue", size = 1) +
#   labs(x = expression(n),
#        y = expression(E(d(H[n])))) +
#   theme_minimal(base_size = 14) +
#   theme(panel.grid.minor     = element_blank(),
#         panel.grid.major.x   = element_blank())
# 
# ggplot(loglog_dens_sum_PY75, aes(x = x, y = mean)) +
#   geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
#               fill = "steelblue", alpha = 0.2) +
#   geom_line(color = "steelblue", size = 1) +
#   labs(x = expression(log(n)),
#        y = expression(log(E(d(H[n]))))) +
#   theme_minimal(base_size = 14) +
#   theme(panel.grid.minor     = element_blank(),
#         panel.grid.major.x   = element_blank())

dens_py <- gridExtra::grid.arrange(dens_sum_plot_py25, dens_sum_plot_py5,
                        dens_sum_plot_py75,
                        loglog_dens_sum_plot_py25, loglog_dens_sum_plot_py5,
                        loglog_dens_sum_plot_py75,
                        nrow = 2,
                        ncol = 3)

dens_py <- gridExtra::grid.arrange(dens_sum_plot_py25, dens_sum_plot_py5,
                                   dens_sum_plot_py75,
                                   ncol = 3)
ggsave(dens_py,
       filename = here(dens_pics, "py", "dens_py.pdf"),
       width = 9,
       height = 5,
       dpi = 300,
       bg = "white")

