library(parallel)
library(pbapply)
library(ggplot2)

# assume sims_by_alpha is your named list of sim_edgesDP… objects
# and thresholds is your vector of sizes

avg_density_by_alpha <- function(sims_by_alpha,
                                          thresholds = c(1e2, 1e3, 1e4, 1e5, 1e6),
                                          ncores = 6) {
  # 1) build a data.frame of all (alpha_name, threshold) jobs
  jobs <- expand.grid(
    alpha_name = names(sims_by_alpha),
    threshold  = thresholds,
    stringsAsFactors = FALSE
  )
  
  # 2) launch a PSOCK cluster
  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl), add = TRUE)
  
  # 3) show a text‐mode progress bar
  pboptions(type = "txt", style = 3)
  
  # 4) run all jobs in parallel with progress
  results_list <- pblapply(
    seq_len(nrow(jobs)),
    cl = cl,
    FUN = function(i) {
      a_name <- jobs$alpha_name[i]
      th     <- jobs$threshold[i]
      sims   <- sims_by_alpha[[a_name]]
      
      # compute density for each simulation at this threshold
      dens_vec <- vapply(sims, function(net) {
        sb <- unique(net[1:th, , drop=FALSE])
        ne <- nrow(sb)
        nA <- length(unique(sb$X_A))
        nB <- length(unique(sb$X_B))
        ne / (nA * nB)
      }, numeric(1))
      
      data.frame(
        alpha        = as.numeric(a_name),
        threshold    = th,
        mean_density = mean(dens_vec)
      )
    }
  )
  
  # 5) combine and return
  avg_df <- do.call(rbind, results_list)
  rownames(avg_df) <- NULL
  avg_df
}

# ---------------------
# Usage
# ---------------------
sims_by_alpha <- list(
  "0.5" = sim_edgesDP0505,
  "1"   = sim_edgesDP11,
  "2"   = sim_edgesDP22,
  "3"   = sim_edgesDP33,
  "5"   = sim_edgesDP55,
  "20"  = sim_edgesDP2020
)

thresholds <- c(1e2, 1e3, 1e4, 1e5, 1e6)

alpha_dens <- avg_density_by_alpha(sims_by_alpha, thresholds)
