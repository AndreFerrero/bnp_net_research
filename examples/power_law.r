library(here)

source(here("code", "funs", "py_sample.r"))

sample <- py_sample(10000, 5, 0.6)

pdf("examples/power_law.pdf", width = 3, height = 2.5)

par(mar = c(0.4, 1, 1.3, 1))

hist(sample$x,
  main = "",
  xlab = "",
  ylab = "",
  axes = FALSE,
  col = "aquamarine2"
)

dev.off()
