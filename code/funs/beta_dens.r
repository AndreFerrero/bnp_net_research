plot_beta_density <- function(alpha, beta, from = 0, to = 1, n = 1000, ...) {
  if (alpha <= 0 || beta <= 0) {
    stop("Alpha and Beta must be positive.")
  }
  
  # Sequence of x values between 0 and 1
  x <- seq(from, to, length.out = n)
  
  # Beta density values
  y <- dbeta(x, alpha, beta)
  
  # Plot
  plot(x, y, type = "l", lwd = 2, col = "blue",
       main = bquote("Beta Density with" ~ alpha == .(alpha) ~ "and" ~ beta == .(beta)),
       xlab = "x", ylab = "Density", ...)
  
  # Optionally add a grid
  grid()
}

plot_beta_density(0.5, 1)
