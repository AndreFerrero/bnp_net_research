# Provided m_xi function (not used in simulation functions below, but retained for reference)
m_xi <- function(xi, xA, xB) {
  xiA <- xi[1]
  xiB <- xi[2]
  
  ms <- 0
  
  # Only proceed if xiB is in xB$x
  if (xiB %in% xB$x) {
    if (!any(xA$x == xiA & xB$x == xiB)) {
      ms <- xB$counts_x[xB$idx_list_x[xiB]]
    }
  }
  
  return(ms)
}