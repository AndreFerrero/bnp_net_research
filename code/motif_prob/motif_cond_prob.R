source("code/setup.R")

set.seed(123)  # For reproducibility

# N = c(100, 500, 5000, 10000, 100000)
N = 10000

for(n in N) {
  d = max(1000, 10*n)
  df_results <- m_cond_dist(N = n, t = 4, d, alpha = 5)
  
  # Compute proportion of non-zero values per simulation
  nonzero_proportions <- df_results %>%
    group_by(Simulation) %>%
    summarize(nonzero_prop = mean(Value > 0) * 100)  # Compute % of non-zero values
  
  # Create the plot
  p <- ggplot(df_results, aes(x = factor(Value), y = ..count../d * 100)) +  
    geom_bar(stat = "count", fill = "steelblue", color = "black", alpha = 0.7) +  
    theme_bw() +
    facet_wrap(~Simulation, scales = "free_x") +  
    labs(title = paste("Distribution of M_{N+1} given X_N,", "N =", n, ", d =", d),
         x = "New motifs",
         y = "Proportion (%)") +
    scale_x_discrete() +
    geom_text(data = nonzero_proportions, 
              aes(x = Inf, y = Inf, label = paste0("Non-zero: ", round(nonzero_prop, 1), "%")), 
              vjust = 1.5, hjust = 1.5, inherit.aes = FALSE, size = 4, color = "red")  # Add text annotation
  
  # Explicitly print the plot
  print(p)
}

