# density_only.R

# --- 0. Load dependencies ---
if (!requireNamespace("pbapply", quietly = TRUE)) 
  install.packages("pbapply")
library(pbapply)

# --- 1. Settings ---
out_dir      <- "sims"
alpha_values <- c(0.5, 1, 1.5, 2, 3, 5, 20)
thresholds   <- c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7)
numCores     <- parallel::detectCores()

# --- 2. Compute densities by alpha (loading only that alpha's networks) ---
# List all .Rdata files in the sims directory
files <- list.files(out_dir,
                    full.names = TRUE)

# Map each file to its alpha value
file_alpha <- vapply(files, function(f) {
  fn  <- basename(f)
  # extract suffix between 'dp' and '.Rdata'
  suf <- sub("^dp(.*)\.Rdata$", "\1", fn)
  # build suffixes from alpha_values
  a_strs   <- gsub("\.", "", as.character(alpha_values))
  suffixes <- paste0(a_strs, a_strs)
  idx <- match(suf, suffixes)
  if (is.na(idx)) stop("Unrecognized file suffix: ", suf)
  alpha_values[idx]
}, numeric(1))

# Parallel over alpha files
density_list <- pbmclapply(
  seq_along(files),
  mc.cores       = min(length(files), numCores),
  mc.preschedule = FALSE,
  FUN = function(i) {
    f     <- files[i]
    a_val <- file_alpha[i]
    message(sprintf("[Density] α=%.1f loading %s", a_val, basename(f)))
    env   <- new.env()
    load(f, envir = env)
    obj_name <- ls(env)[1]
    nets     <- env[[obj_name]]
    rm(env)
    
    # compute mean density for each threshold
    dens_vals <- vapply(thresholds, function(th) {
      dens_sim <- vapply(nets, function(net) {
        sb <- unique(net[1:th, , drop = FALSE])
        ne <- nrow(sb)
        nA <- length(unique(sb$X_A))
        nB <- length(unique(sb$X_B))
        ne / (nA * nB)
      }, numeric(1))
      mean(dens_sim)
    }, numeric(1))
    
    data.frame(
      alpha        = a_val,
      threshold    = thresholds,
      mean_density = dens_vals,
      row.names    = NULL
    )
  }
)

# --- 3. Combine and save density results ---
alpha_dens <- do.call(rbind, density_list)
save(alpha_dens, file = file.path(out_dir, "alpha_density_results.Rdata"))
message("=== Density results saved to ",
        file.path(out_dir, "alpha_density_results.Rdata"),
        " ===")
