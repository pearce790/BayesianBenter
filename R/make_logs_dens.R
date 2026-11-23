# Log collapsed conditionals for the Metropolis-Hastings updates of alpha_r under complete rankings.
# Internal; not exported.

#' Collapsed log-conditional factory for alpha_r (complete rankings)
#'
#' Returns a closure of a single argument `a` giving the (unnormalized)
#' log full conditional of `alpha_r` at stage `r`, using the unique rankings
#' and their multiplicities in `pi_unique` (cols 1..J) and counts (col J+1).
#' @noRd
make_log_dens <- function(pi_unique, loglambda, r, J, c, d) {
  risk_r    <- pi_unique[, r:J, drop = FALSE]                  # U x m
  logL_r    <- matrix(loglambda[risk_r], nrow = nrow(risk_r))  # U x m
  w         <- pi_unique[, J + 1L]                             # U
  sum_w_win <- sum(w * logL_r[, 1])                            # scalar

  function(a) {
    (c - 1) * log(a) + (d - 1) * log1p(-a) +
      a * sum_w_win -
      sum(w * matrixStats::rowLogSumExps(a * logL_r))
  }
}
