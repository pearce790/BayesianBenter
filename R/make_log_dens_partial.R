# Log collapsed conditionals for the Metropolis-Hastings updates of alpha_r under partial rankings.
# Internal; not exported.

#' Collapsed log-conditional factory for alpha_r (partial rankings)
#'
#' Identical to `make_log_dens()` but restricts to unique classes that observe stage `r`
#' (`pi_unique[, J + 2]` holds r_max). If no class observes stage `r`, the
#' conditional collapses to the Beta(c, d) prior.
#' @noRd
make_log_dens_partial <- function(pi_unique, loglambda, r, J, c, d) {
  ## pi_unique layout: cols 1..J = completed ranking; col J+1 = count;
  ## col J+2 = r_max. Only voters with r_max >= r contribute at stage r.
  active <- pi_unique[, J + 2L] >= r
  if (!any(active)) {
    return(function(a) (c - 1) * log(a) + (d - 1) * log1p(-a))
  }
  pu        <- pi_unique[active, , drop = FALSE]
  risk_r    <- pu[, r:J, drop = FALSE]
  logL_r    <- matrix(loglambda[risk_r], nrow = nrow(risk_r))
  w         <- pu[, J + 1L]
  sum_w_win <- sum(w * logL_r[, 1])
  
  function(a) {
    (c - 1) * log(a) + (d - 1) * log1p(-a) +
      a * sum_w_win -
      sum(w * matrixStats::rowLogSumExps(a * logL_r))
  }
}
