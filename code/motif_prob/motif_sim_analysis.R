source("code/setup.R")

## RESULTS ####

sim = simulation_networks_parallel(S = 1000,
                                   N = 100000, alphaA = 5, alphaB = 50,
                                   ncores = 6)
# sim$varM = apply(sim$M, 2, var)



plot_sim_motif = function(sim, N, alpha, p0 = FALSE, pg0 = FALSE,
                          counts = FALSE, log_tot = FALSE, var = FALSE, tot = FALSE){
  
  # Titles with dynamic alpha values using bquote()
  title_alpha <- bquote( ~ alpha[A] ~ "=" ~ .(alpha[1]) ~ "," ~ alpha[B] ~ "=" ~ .(alpha[2]))
  
  if(p0){
    plot(1:N, sim$avgZero[1:N], type = "l",
         ylab = "", xlab = "Network size",
         main = bquote(P(M[n]^"*" == 0) ~ "," ~ .(title_alpha)),
         ylim = c(0,1))
  }
  
  if(pg0){
    plot(1:N, 1 - sim$avgZero[1:N], type = "l",
         ylab = "", xlab = "Network size",
         main = bquote(P(M[n]^"*" > 0) ~ "," ~ .(title_alpha)),
         ylim = c(0,1))
  }
  
  if(counts){
    plot(1:N, sim$avgM[1:N], type = "l",
         ylab = "Average Counts", xlab = "Network size",
         main = bquote("E" * "[" * M[n]^"*" * "]" ~ "," ~ .(title_alpha)))
  }
  
  if(var){
    plot(1:N, sim$varM[1:N], type = "l",
         ylab = "Counts", xlab = "Network size",
         main = bquote(Var(M(N)^"*") ~ "," ~ .(title_alpha)))
  }
  
  if(tot){
    plot(1:N, sim$avgCumulativeMotifs[1:N], type = "l",
         ylab = "Counts", xlab = "Network size",
         main = bquote("E" * "[" * M^"*" * "(" * N * ")" * "]" ~ "," ~ .(title_alpha)))
  }
  
  if(log_tot){
    b <- 0.4
    plot(log(1:N), log(sim$avgCumulativeMotifs[1:N]), type = "l",
         main = bquote("E" * "[" * M^"*" * "(" * N * ")" * "]" ~ "," ~ .(title_alpha) ~ "," ~ "Slope =" ~ .(b)),
         ylab = "Log average counts",
         xlab = "Log Network size")
    abline(a = 4, b = b, col = "red")
  }
}


plot_sim_motif(sim, 100000, alpha = c(5,5), tot = T)

plot_sim_motif(sim, 100000, alpha = c(5,5), counts = T)

plot_sim_motif(sim, 1000, alpha = c(5,5), p0 = T)

plot_sim_motif(sim, 1000, alpha = c(5,5), pg0 = T)

# plot_sim_motif(sim, 100000, alpha = c(5,5), var = T)

plot_sim_motif(sim, 100000, alpha = c(5,5), log_tot = T)

sizes = 1:100000

selected_columns <- c(500, 1000, 5000, 10000, 50000, 100000)

# Extract the data for these columns
data_to_plot <- sim$M[, selected_columns, drop = FALSE]

# Create the boxplot
boxplot(data_to_plot, 
        names = selected_columns,  # Set x-axis labels to chosen n-values
        main = "Distribution of Motifs at Different Network Sizes",
        xlab = "Network Size n", 
        ylab = "Number of Motifs",
        col = "lightblue", 
        border = "darkblue")
