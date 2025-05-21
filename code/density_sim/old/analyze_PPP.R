library(ggplot2)
library(viridis)

analyze_sim_edges <- function(sim_edges,
                              sizes = 1:10000,
                              name = NULL,
                              save_plots = FALSE,
                              width = 10,
                              height = 5,
                              dpi = 600) {
  # 1) Determine name for title/parsing
  if (is.null(name)) name <- deparse(substitute(sim_edges))
  info    <- sub(".*sim_edges", "", name)
  process <- regmatches(info, regexpr("^[A-Z]+", info))
  digits  <- sub("^[A-Z]+", "", info)
  
  # 2) Extract alpha
  if (process == "DP") {
    if (nchar(digits) <= 2) {
      # simple single-digit case
      a1    <- as.numeric(substr(digits, 1, 1))
      a2    <- as.numeric(substr(digits, 2, 2))
      alpha <- if (a1 == a2) a1 else paste0(a1, ",", a2)
    } else if (nchar(digits) == 4) {
      # head two digits are the alpha code
      head <- substr(digits, 1, 2)
      if (substr(head, 1, 1) == "0") {
        # alpha < 1
        alpha <- as.numeric(paste0("0.", substr(head, 2, 2)))
      } else {
        # alpha >= 10
        alpha <- as.numeric(head)
      }
    } else {
      stop("Cannot parse DP alpha code: ", digits)
    }
  } else {
    # PY: first two single-digit alphas (assumed equal)
    a1    <- as.numeric(substr(digits, 1, 1))
    a2    <- as.numeric(substr(digits, 2, 2))
    alpha <- if (a1 == a2) a1 else paste0(a1, ",", a2)
  }
  
  # 3) Extract sigma only for PY
  sigma <- NA
  if (process == "PY" && nchar(digits) > 2) {
    tail <- substring(digits, 3)
    if (nchar(tail) == 2) {
      s1 <- as.numeric(substr(tail, 1, 1)) / 10
      s2 <- as.numeric(substr(tail, 2, 2)) / 10
    } else if (nchar(tail) == 4) {
      s1 <- as.numeric(substr(tail, 1, 2)) / 100
      s2 <- as.numeric(substr(tail, 3, 4)) / 100
    } else {
      stop("Cannot parse PY sigma code: ", tail)
    }
    sigma <- if (!is.na(s1) && !is.na(s2) && s1 == s2) s1 else paste0(s1, ",", s2)
  }
  
  # 4) Build process label
  process_label <- if (process == "DP") {
    paste0("DP(", alpha, ")")
  } else {
    paste0("PY(", alpha, if (!is.na(sigma)) paste0(", ", sigma), ")")
  }
  
  # 5) Build dynamic results
  nsim     <- length(sim_edges)
  dyn_list <- vector("list", nsim * length(sizes))
  idx      <- 1
  for (i in seq_len(nsim)) {
    net <- sim_edges[[i]]
    for (th in sizes) {
      sb   <- net[1:th, , drop = FALSE]
      u    <- unique(sb)
      ne   <- nrow(u)
      nA   <- length(unique(u$X_A))
      nB   <- length(unique(u$X_B))
      mp   <- nA * nB
      dens <- ne / mp
      
      dyn_list[[idx]] <- data.frame(
        Simulation  = i,
        Threshold   = th,
        UniqueEdges = ne,
        nA          = nA,
        nB          = nB,
        MaxPossible = mp,
        Density     = dens
      )
      idx <- idx + 1
    }
  }
  dyn_results <- do.call(rbind, dyn_list)
  title_txt   <- paste0(process_label, ", ", nsim, " networks")
  
  # 6) Generate five plots
  p1 <- ggplot(dyn_results, aes(Threshold, Density,
                                group = factor(Simulation),
                                color = factor(Simulation))) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_viridis_d(option = "plasma", begin = 0, end = 1) +
    labs(x = "Observed Edges", y = "Density", title = title_txt) +
    theme_minimal() + theme(legend.position = "none")
  
  p2 <- ggplot(dyn_results, aes(log(Threshold), log(UniqueEdges),
                                group = factor(Simulation),
                                color = factor(Simulation))) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_viridis_d(option = "plasma", begin = 0, end = 1) +
    labs(x = "log(Observed Edges)", y = "log(Unique Edges)", title = title_txt) +
    theme_minimal() + theme(legend.position = "none")
  
  p3 <- ggplot(dyn_results, aes(log(Threshold), log(MaxPossible),
                                group = factor(Simulation),
                                color = factor(Simulation))) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_viridis_d(option = "plasma", begin = 0, end = 1) +
    labs(x = "log(Observed Edges)", y = "log(nA * nB)", title = title_txt) +
    theme_minimal() + theme(legend.position = "none")
  
  p4 <- ggplot(dyn_results, aes(log(MaxPossible), log(UniqueEdges),
                                group = factor(Simulation),
                                color = factor(Simulation))) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_viridis_d(option = "plasma", begin = 0, end = 1) +
    labs(x = "log(Max Possible)", y = "log(Unique Edges)", title = title_txt) +
    theme_minimal() + theme(legend.position = "none")
  
  p5 <- ggplot(dyn_results, aes(log(nA + nB), log(UniqueEdges),
                                group = factor(Simulation),
                                color = factor(Simulation))) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_viridis_d(option = "plasma", begin = 0, end = 1) +
    labs(x = "log(nA + nB)", y = "log(Unique Edges)", title = title_txt) +
    theme_minimal() + theme(legend.position = "none")
  
  # 7) Optionally save
  if (save_plots) {
    dir.create("sims/pics/density_analysis", recursive = TRUE, showWarnings = FALSE)
    fname_base <- gsub("[ (),]", "_", process_label)
    
    ggsave(
      filename = file.path("sims/pics/density_analysis", paste0(fname_base, "_density.png")),
      plot     = p1, width = width, height = height, dpi = dpi, bg = "white"
    )
    ggsave(
      filename = file.path("sims/pics/density_analysis", paste0(fname_base, "_log_unique_edges.png")),
      plot     = p2, width = width, height = height, dpi = dpi, bg = "white"
    )
    ggsave(
      filename = file.path("sims/pics/density_analysis", paste0(fname_base, "_log_max_possible.png")),
      plot     = p3, width = width, height = height, dpi = dpi, bg = "white"
    )
    ggsave(
      filename = file.path("sims/pics/density_analysis", paste0(fname_base, "_log_unique_vs_max.png")),
      plot     = p4, width = width, height = height, dpi = dpi, bg = "white"
    )
    ggsave(
      filename = file.path("sims/pics/density_analysis", paste0(fname_base, "_log_unique_vs_nAplusnB.png")),
      plot     = p5, width = width, height = height, dpi = dpi, bg = "white"
    )
  }
  
  # 8) Print plots & return data invisibly
  print(p1); print(p2); print(p3); print(p4); print(p5)
  invisible(dyn_results)
}

analyze_sim_edges(sim_edgesDP0505[1:5],
                  name = "sim_edgesDP0505",
                  save = T)

analyze_sim_edges(sim_edgesDP55[1:5],
                  name = "sim_edgesDP55",
                  save = T)

analyze_sim_edges(sim_edgesDP2020[1:5],
                  name = "sim_edgesDP2020",
                  save = T)

analyze_sim_edges(sim_edgesPY552525[1:5],
                  name = "sim_edgesPY552525",
                  save = T)

analyze_sim_edges(sim_edgesPY5566[1:5],
                  name = "sim_edgesPY5566",
                  save = T)

analyze_sim_edges(sim_edgesPY557575[1:5],
                  name = "sim_edgesPY557575",
                  save_plots = T)
