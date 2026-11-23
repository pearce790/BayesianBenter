#' Bayesian estimation of the Benter ranking model
#'
#' Fits the Benter (1994) model for ranking data by Gibbs sampling. The Benter
#' model generalizes the Plackett-Luce model by introducing stage-specific
#' *dampening* exponents \eqn{\alpha_r \in (0, 1]} that down-weight the
#' information in later ranking stages:
#' \deqn{p(\pi_i \mid \lambda, \alpha) = \prod_{r=1}^{J-1}
#'   \frac{\lambda_{w_{ir}}^{\alpha_r}}{\sum_{j \in R_{ir}} \lambda_j^{\alpha_r}},}
#' where \eqn{\lambda_j > 0} are item worths, \eqn{w_{ir}} is the item placed
#' \eqn{r}th by voter \eqn{i}, and \eqn{R_{ir}} is the set of items still at risk
#' at stage \eqn{r}. Inference uses a positive-stable frailty data augmentation
#' that restores conjugacy, so the worths are drawn exactly from Gamma full
#' conditionals while each \eqn{\alpha_r} is updated by Metropolis-Hastings on
#' its collapsed conditional.
#'
#' `gibbsBenter()` is the user-facing entry point. It validates and
#' standardizes `pi`, then dispatches to the appropriate internal sampler for
#' complete or partial rankings. A length-(J-1) ranking pins down its final
#' item and is completed automatically (this is deduction, not imputation); a
#' row is treated as genuinely partial only if it leaves two or more items
#' unplaced.
#'
#' For partial rankings, only the dampening parameters up to the deepest
#' observed rank are identifiable, so `alpha_r` for stages beyond the longest
#' observed ranking are neither updated nor returned (reported as `NA`). The
#' item worths \eqn{\lambda} for all `J` items are always estimated.
#'
#' @param pi A ranking data set: an integer matrix (or data frame, or a single
#'   vector interpreted as one ranking) with one ranking per row. Entry
#'   `pi[i, r]` is the item placed `r`th by voter `i`; ranked entries are
#'   left-aligned and any unranked tail positions are `NA`. All non-missing
#'   entries must be integers in `1:J` with no repeats within a row.
#' @param J Integer (`>= 3`). The total number of items being ranked.
#' @param a,b Shape and rate of the Gamma(`a`, `b`) prior on the item worths
#'   \eqn{\lambda_j}. On the identified (sum-to-one) scale, `a` alone controls
#'   the induced symmetric Dirichlet prior on the normalized worths; `b` only
#'   sets the (marginalized-out) overall scale.
#' @param c,d Shape parameters of the Beta(`c`, `d`) prior on each dampening
#'   exponent \eqn{\alpha_r}, `r`=2,...,J-1. (`alpha_1` is fixed at 1 for
#'   identifiability.)
#' @param alpha_min Lower floor for the Metropolis-Hastings proposals of
#'   \eqn{\alpha_r}. Proposals at or below `alpha_min` are rejected. A small
#'   positive value (e.g. `0.05`) keeps `alpha_r` inside the range where the
#'   tilted-stable sampler is numerically stable; the default `0` imposes no
#'   floor beyond the open boundary.
#' @param iters Number of Gibbs iterations (including burn-in).
#' @param mh_iters Number of Metropolis-Hastings sub-steps per \eqn{\alpha_r}
#'   per sweep.
#' @param burnin Fraction of `iters` discarded as burn-in. The adaptive tuning
#'   of the proposal scales is active during burn-in and frozen afterwards.
#' @param verbose Logical. If `TRUE`, show a text progress bar.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`lambda`}{A `(iters - burn) x J` matrix of posterior draws of the
#'       item worths, normalized to sum to one in each row.}
#'     \item{`alpha`}{A `(iters - burn) x (J-1)` matrix of posterior draws of
#'       the dampening exponents. Column 1 is fixed at 1; for partial data,
#'       columns beyond the deepest observed rank are `NA`.}
#'     \item{`parameters`}{A list recording the prior and sampler settings and
#'       the Metropolis-Hastings acceptance rates.}
#'     \item{`time_secs`}{Wall-clock run time in seconds.}
#'   }
#'
#' @references
#' Benter, W. (1994). Computer-based horse race handicapping and wagering
#' systems: A report. In *Efficiency of Racetrack Betting Markets*.
#'
#' @seealso [rBenter()] to simulate rankings from the model.
#'
#' @examples
#' set.seed(1)
#' lambda <- c(5, 3, 2, 1.5, 1)
#' lambda <- lambda / sum(lambda)
#' alpha  <- c(1, 0.9, 0.8, 0.7)
#' pi <- rBenter(200, lambda, alpha)
#'
#' \donttest{
#' fit <- gibbsBenter(pi, J = 5, a = 1, b = 1, c = 1, d = 1,
#'                    iters = 400, burnin = 0.5)
#' colMeans(fit$lambda)   # posterior mean worths (vs. truth `lambda`)
#' colMeans(fit$alpha)    # posterior mean dampening (vs. truth `alpha`)
#' }
#'
#' @export
gibbsBenter <- function(pi, J, a, b, c, d, alpha_min = 0,
                        iters = 1000, mh_iters = 5, burnin = 0.5,
                        verbose = FALSE){

  # Primary dispatcher: validate/standardize pi, then route to the complete or
  # partial sampler. Does by construction whatever needs no inference (matrix
  # coercion, padding to J columns, filling the determined last place of a
  # length-(J-1) ranking); errors clearly on anything that cannot be assumed.
  #
  # Guarantees handed to the samplers:
  #   * pi is an integer matrix with exactly J columns
  #   * every entry is an integer in 1:J, or NA
  #   * within each row, ranked entries are left-aligned and NAs (if any) trail
  #   * no item repeats within a row
  #   * every row has >= 1 ranked entry
  #   * any ranking that pins down all J items (length J-1 or J) is completed,
  #     so a row is left with NAs only if it is genuinely partial

  ## --- J is a single integer >= 3 ------------------------------------------
  if(length(J) != 1L || is.na(J) || J < 3 || J %% 1 != 0)
    stop("J must be a single integer >= 3.")
  J <- as.integer(J)

  ## --- coerce pi to a matrix (accept data frame, or a bare vector = 1 row) --
  if(is.data.frame(pi)){pi <- as.matrix(pi)}
  if(!is.matrix(pi)){
    if(is.atomic(pi) && is.null(dim(pi))) pi <- matrix(pi, nrow = 1L)
    else stop("pi must be a matrix, data frame, or atomic vector.")
  }

  ## --- coerce entries to numeric (catches character/factor columns) --------
  if(!is.numeric(pi)){
    raw_na <- is.na(pi)
    num    <- suppressWarnings(as.numeric(pi))
    if(any(is.na(num) & !raw_na))
      stop("pi contains non-numeric entries; entries must be integers in 1:J or NA.")
    pi <- matrix(num, nrow = nrow(pi))
  }

  n <- nrow(pi)
  if(n < 1L) stop("pi has no rows.")

  ## --- column count cannot exceed J ----------------------------------------
  if(ncol(pi) > J)
    stop(sprintf("pi has %d columns but J = %d; a ranking of J items needs at most J columns.",
                 ncol(pi), J))

  ## --- entries are integers in 1:J (NA allowed) ----------------------------
  vals <- pi[!is.na(pi)]
  if(length(vals) == 0L) stop("pi has no ranking entries (all values are NA).")
  if(any(vals %% 1 != 0) || any(vals < 1) || any(vals > J))
    stop("All non-missing entries of pi must be integers in 1:J.")
  storage.mode(pi) <- "integer"

  ## --- per-row structure: >=1 entry, left-aligned (NAs trail), no repeats ---
  m <- ncol(pi)
  r_len <- integer(n)
  for(i in seq_len(n)){
    nna <- !is.na(pi[i, ])
    k   <- sum(nna)
    if(k < 1L)
      stop(sprintf("Row %d has no ranking entries; every ranking needs at least one item.", i))
    if(!all(nna[seq_len(k)]))
      stop(sprintf("Row %d has an interior NA; ranked items must be left-aligned with NA padding only at the end.", i))
    if(anyDuplicated(pi[i, seq_len(k)]))
      stop(sprintf("Row %d repeats an item; items within a ranking must be distinct.", i))
    r_len[i] <- k
  }

  ## --- pad to J columns (no inference) -------------------------------------
  if(m < J) pi <- cbind(pi, matrix(NA_integer_, n, J - m))

  ## --- fill the determined last place of any length-(J-1) ranking ----------
  ## (the final item is logically pinned down; this is not inference)
  fill <- which(r_len == (J - 1L))
  for(i in fill) pi[i, J] <- setdiff(1:J, pi[i, 1:(J - 1L)])

  ## --- route -----------------------------------------------------------------
  if(anyNA(pi)){
    if(verbose) message("Partial rankings detected (", sum(is.na(pi[,J])),
                        " of ", n, " rankings are incomplete).")
    result <- gibbsBenter_partial(pi, J, a, b, c, d, alpha_min, iters, mh_iters, burnin, verbose)
  } else {
    if(verbose) message("Complete rankings detected.")
    result <- gibbsBenter_complete(pi, J, a, b, c, d, alpha_min, iters, mh_iters, burnin, verbose)
  }
  result
}
