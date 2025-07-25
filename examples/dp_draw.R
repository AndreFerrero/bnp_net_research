x <- seq(-4, 4, by = 0.01)
cdf <- diff(c(0, pnorm(x)))

# Dirichlet Distribution Sampling Function
rdiric <- function(n, a) {
  # Generates samples from Dirichlet distribution Dir(a)
  p <- length(a)
  m <- matrix(nrow = n, ncol = p)
  for (i in 1:p) {
    m[, i] <- rgamma(n, a[i])
  }
  sumvec <- rowSums(m)
  m / sumvec
}

# Dirichlet process from definition
dp <- function(N, G0, alpha, plot.cdf = FALSE) {
  alpha_vec <- alpha * G0 # This is the vector used for sampling

  dp_sample <- VGAM::rdiric(N, alpha_vec)

  # dp.sample.cdf <- t(apply(dp.sample, 1, cumsum))

  if (plot.cdf) {
    plot(x, cumsum(G0),
      type = "l", col = "black", lwd = 2,
      ylab = "CDF", xlab = NA,
      main = bquote(alpha == .(alpha))
    ) # Use scalar alpha here

    for (i in 1:N) {
      lines(x, cumsum(dp_sample[i, ]), col = "red", lty = 2)
    }
  }
}


pdf("dp_draw_smaller.pdf", width = 5, height = 3)
par(mfrow = c(1, 3))
dp(N = 10, G0 = cdf, alpha = 1, plot.cdf = TRUE)
dp(N = 10, G0 = npdf, alpha = 10, plot.cdf = TRUE)
dp(N = 10, G0 = npdf, alpha = 100, plot.cdf = TRUE)

dev.off()
