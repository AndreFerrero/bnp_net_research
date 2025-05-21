# Function to process the list output and plot boxplots
plot_boxplots <- function(sim_list, m_values, Nsim, p = 1) {
  
  means = list()
  vars = list()
  
  par(mfrow = c(1,p))
  for (i in 1:length(m_values)) {
    # Extract the i-th motif across all simulations
    motif_data <- lapply(sim_list, function(mat) mat[i, ])
    motif_data <- do.call(cbind, motif_data)  # Convert list to matrix
    
    means[[i]] = colMeans(motif_data)
    vars[[i]] = apply(motif_data, 2, var)
    
    # Boxplot for each motif
    boxplot(motif_data, names = Nsim, main = paste("Motif", m_values[i]),
            xlab = "Nsim", ylab = "Frequency", col = "lightblue")
    points(means[[i]], col = "red", pch = 16, cex = 0.7)
  }
  
  par(mfrow = c(1,1))
  
  return(list(means = means,
              vars = vars))
}

sizes = seq(500, 100000, 500)

m = c(3,5,10)
sim_res = plot_boxplots(sim_to_100000$m_freq, m, sizes)

mm = 1

plot(sizes, sim_res$means[[mm]], pch = 19, main = paste("Motif", m[mm]),
     ylab = "Average frequency")

plot(sqrt(sizes), sim_res$means[[mm]], pch = 19, main =  paste("Motif", m[mm]))

plot(log(sizes), sim_res$means[[mm]], pch = 19, main =  paste("Motif", m[mm]),
     ylab = "Average frequency")

plot(log(sizes), log(sim_res$means[[mm]]), pch = 19, main =  paste("Motif", m[mm]),
     ylab = "Log average frequency")

plot(sizes, sim_res$vars[[mm]], pch = 19, main = paste("Motif", m[mm]),
     ylab = "Variability of frequency")

plot(log(sizes), log(sim_res$vars[[mm]]), pch = 19, main = paste("Motif", m[mm]),
     ylab = "Variability of frequency")

plot(5*sizes*log(sizes), sim_res$means[[mm]], pch = 19, main =  paste("Motif", m[mm]),
     ylab = "Average frequency")

# CURVES

# sqrt_scale_factor <- max(sim_res[[mm]]) / sqrt(max(sizes))
# 
# curve(sqrt_scale_factor * sqrt(x), from = min(sizes), to = max(sizes), 
#       add = TRUE, col = "red", lwd = 2)
# 
# log_scale_factor <- max(sim_res[[mm]]) / log(max(sizes))
# 
# y_shift <- min(sim_res[[mm]]) - log_scale_factor * log(min(sizes))
# 
# curve(log_scale_factor * log(x) + y_shift, from = min(sizes), to = max(sizes), 
#       add = TRUE, col = "blue", lwd = 2)
# 
# legend("bottomright", legend = c("sqrt(sizes)", "log(sizes)"), 
#        col = c("red", "blue"), 
#        lty = c(1, 1), lwd = c(2, 2))

plot(sizes, sim_to_100000$time)
