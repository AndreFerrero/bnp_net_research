set.seed(123)

# Parameters
alpha <- 1
N <- 30

# Simulate CRP
table_assignments <- numeric(N)
tables <- list()

for (i in 1:N) {
  if (i == 1) {
    table_assignments[i] <- 1
    tables[[1]] <- 1
  } else {
    probs <- c(sapply(tables, length), alpha)
    probs <- probs / sum(probs)
    choice <- sample(1:(length(tables) + 1), 1, prob = probs)
    if (choice > length(tables)) {
      tables[[choice]] <- i
    } else {
      tables[[choice]] <- c(tables[[choice]], i)
    }
    table_assignments[i] <- choice
  }
}

# Save to PDF
pdf("crp.pdf", width = 4, height = 3)

# Prepare plot
n_tables <- length(tables)
colors <- rainbow(n_tables)
par(mar = c(0, 0, 0, 0)) # bottom, left, top, right
plot(NA,
  xlim = c(-10, 10), ylim = c(-10, 10),
  xlab = "", ylab = "", asp = 1, axes = FALSE,
  xaxs = "i", yaxs = "i"
)

# Function to draw circular table and customer dots (without labels)
draw_table <- function(center_x, center_y, radius, customers, color) {
  symbols(center_x, center_y,
    circles = radius, inches = FALSE,
    add = TRUE, bg = adjustcolor(color, alpha.f = 0.2), fg = NA
  )

  n <- length(customers)
  if (n == 1) {
    points(center_x, center_y + radius, pch = 21, bg = color, cex = 1.5)
  } else {
    angles <- seq(0, 2 * pi, length.out = n + 1)[-1]
    for (i in 1:n) {
      x <- center_x + radius * cos(angles[i])
      y <- center_y + radius * sin(angles[i])
      points(x, y, pch = 21, bg = color, cex = 1.5)
    }
  }
}

# Arrange tables in circular layout
angle_step <- 2 * pi / n_tables
layout_radius <- 6
for (i in 1:n_tables) {
  theta <- angle_step * (i - 1)
  cx <- layout_radius * cos(theta)
  cy <- layout_radius * sin(theta)
  draw_table(cx, cy, radius = 1.5, customers = tables[[i]], color = colors[i])
}

# Close PDF device
dev.off()
