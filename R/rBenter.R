#' Simulate complete rankings from the Benter model
#'
#' Draws `n` independent complete rankings of `J = length(lambda)` items from
#' the Benter model with worths `lambda` and dampening exponents `alpha`. At
#' each stage `r`, the next item is chosen from those still available with
#' probability proportional to `lambda^alpha[r]`.
#'
#' @param n Number of rankings to simulate.
#' @param lambda Positive vector of item worths, length `J`. Need not sum to
#'   one (only relative magnitudes matter).
#' @param alpha Vector of dampening exponents, length `J - 1`, with
#'   `alpha[1] == 1`.
#'
#' @return An `n x J` integer matrix; row `i` is a complete ranking, with
#'   `pi[i, r]` the item placed `r`th.
#'
#' @seealso [gibbsBenter()] to fit the model.
#'
#' @examples
#' set.seed(1)
#' lambda <- c(5, 3, 2, 1.5, 1) / sum(c(5, 3, 2, 1.5, 1))
#' alpha  <- c(1, 0.8, 0.5, 0.3)
#' rBenter(5, lambda, alpha)
#'
#' @export
rBenter <- function(n, lambda, alpha) {
  J <- length(lambda)
  stopifnot(length(alpha) == J - 1L, alpha[1] == 1)
  t(replicate(n, {
    avail <- seq_len(J); out <- integer(J)
    for (r in seq_len(J - 1L)) {
      w      <- sample.int(length(avail), 1L,
                           prob = lambda[avail]^alpha[r])
      out[r] <- avail[w]; avail <- avail[-w]
    }
    out[J] <- avail
    out
  }))
}
