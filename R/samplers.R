# Internal Gibbs samplers dispatched by gibbsBenter(). Not exported.
# See gibbsBenter() for the model, arguments, and return value.

#' Complete-ranking Gibbs sampler
#' @noRd
gibbsBenter_complete <- function(pi, J, a, b, c, d, alpha_min, iters, mh_iters, burnin, verbose){

  t0 <- Sys.time()
  ## Pre-Calculation
  n <- nrow(pi)
  prop_sd <- rep(0.05, J - 1)
  n_j <- sapply(1:J, function(j){sum(pi[, 1:(J - 1)] == j)})
  pi_unique <- suppressMessages(as.matrix(
    as.data.frame(pi) %>% group_by_all() %>% summarize(count = n())))
  unique_partials <- lapply(2:(J - 1), function(r){
    sorted_pi <- t(apply(pi[, r:J], 1, sort))
    tmp <- unique(sorted_pi)
    lapply(1:nrow(tmp), function(ind){
      list(sorted_partial_pi = tmp[ind, ],
           which_voters = which(apply(t(sorted_pi) == tmp[ind, ], 2, all)))
    })
  })
  burn <- floor(iters * burnin)

  ## Step 1 Pre-Calculation
  rank_mat <- matrix(NA_integer_, n, J)
  for (r in 1:J) rank_mat[cbind(1:n, pi[, r])] <- r
  mask <- array(FALSE, dim = c(n, J - 1L, J))
  for (r in 1:(J - 1L)) mask[, r, ] <- rank_mat >= r

  ## Initializers
  S <- array(1, c(n, J - 1, J))
  V <- matrix(1, nrow = n, ncol = J - 1)
  alpha <- c(1, rep(.9, J - 2))

  sum_sv <- colSums(matrix(S * as.vector(V) * mask, nrow = n * (J - 1L)))
  lambda <- rgamma(J, a + n_j, b + sum_sv)

  ## Storage
  steps <- accepts <- rep(0, J - 1)
  lambda_out <- matrix(NA, nrow = iters - burn, ncol = J)
  alpha_out <- matrix(NA, nrow = iters - burn, ncol = J - 1)

  if(verbose) pb <- utils::txtProgressBar(min = 0, max = iters, style = 3)

  ## Gibbs Sweeps
  for(iter in 1:iters){

    ## Step 1(a)
    V[, 1] <- rexp(n, sum(lambda))

    ## Step 1(b)
    loglambda <- log(lambda)

    ### 1(b)(i)
    if(iter %% 50 == 0 & iter <= burn){
      acc <- accepts[-1] / steps[-1]
      prop_sd[2:(J - 1)] <- prop_sd[2:(J - 1)] * ifelse(acc < 0.35, 0.8,
                                                        ifelse(acc > 0.50, 1.25, 1))
      prop_sd[prop_sd < 0.001] <- 0.001
      steps <- accepts <- rep(0, J - 1)
    }
    par_alpha <- lapply(2:(J - 1L), function(r) {
      log_dens <- make_log_dens(pi_unique, loglambda, r, J, c, d)
      a        <- alpha[r]
      log_den  <- log_dens(a)
      accept   <- 0L

      for (t in 1:mh_iters) {
        a_prop <- a + rnorm(1, sd = prop_sd[r])
        if (a_prop <= alpha_min || a_prop >= 1) next
        log_num <- log_dens(a_prop)
        diff <- log_num - log_den
        if (!is.nan(diff) & log(runif(1)) < diff) {
          a       <- a_prop
          log_den <- log_num
          accept  <- accept + 1L
        }
      }
      c(mh_iters, accept, a)
    })
    for(r in 2:(J - 1)){
      steps[r] <- steps[r] + par_alpha[[r - 1]][1]
      accepts[r] <- accepts[r] + par_alpha[[r - 1]][2]
      alpha[r] <- par_alpha[[r - 1]][3]
    }

    ### 1(b)(ii)
    for(r in 2:(J - 1)){
      for(ind in 1:length(unique_partials[[r - 1]])){
        R_ir <- unique_partials[[r - 1]][[ind]]$sorted_partial_pi
        which_voters <- unique_partials[[r - 1]][[ind]]$which_voters
        V[which_voters, r] <- rweibull(length(which_voters), alpha[r],
                                       exp(-logSumExp(alpha[r] * loglambda[R_ir]) / alpha[r]))
      }
    }

    ### 1(b)(iii)
    for (r in 2:(J - 1L)) {
      al       <- alpha[r]
      risk     <- pi[, r:J, drop = FALSE]                        # n x (J - r + 1)
      lam_risk <- matrix(lambda[risk], nrow = n)                 # same shape
      h        <- V[, r] * lam_risk                              # V[,r] broadcast

      ## One retstable call for every (i, j in R_ir) at stage r
      S_ir <- matrix(retstable(alpha = al, V0 = as.vector(h^al)), nrow = n) / h

      ## Winner (column 1 of risk = pi[, r]) picks up the size-biasing Gamma term
      S_ir[, 1] <- S_ir[, 1] + rgamma(n, shape = 1 - al, rate = h[, 1])

      ## Scatter back into the 3D array S[i, r, j]
      S[cbind(rep(seq_len(n), J - r + 1L), r, as.vector(risk))] <- as.vector(S_ir)
    }

    ## Step 2
    sum_sv <- colSums(matrix(S * as.vector(V) * mask, nrow = n * (J - 1L)))
    lambda <- rgamma(J, a + n_j, b + sum_sv)

    ## Saving Samples
    if(iter > burn){
      lambda_out[iter - burn, ] <- lambda / sum(lambda)
      alpha_out[iter - burn, ] <- alpha
    }

    if(verbose) utils::setTxtProgressBar(pb, iter)
  }

  if(verbose) close(pb)
  t1 <- Sys.time()

  ## Returns
  list(lambda = lambda_out,
       alpha = alpha_out,
       parameters = list(a = a, b = b, c = c, d = d,
                         alpha_min = alpha_min, iters = iters, mh_iters = mh_iters,
                         burnin = burnin, prop_sd = prop_sd,
                         accept_prob = accepts[-1] / steps[-1]),
       time_secs = as.numeric(difftime(t1, t0, units = "secs")))
}

#' Partial-ranking Gibbs sampler
#'
#' Mirrors the structure of the complete sampler, with partial-ranking
#' bookkeeping (completed risk sets, active-voter masks) and depth truncation:
#' only `alpha_2, ..., alpha_{R_obs}` are estimated/kept, where
#' `R_obs = min(max(r_max), J - 1)` is the deepest observed rank. `alpha_r` for
#' `r` > `R_obs` are returned as `NA`.
#' @noRd
gibbsBenter_partial <- function(pi, J, a, b, c, d, alpha_min, iters, mh_iters, burnin, verbose){

  t0 <- Sys.time()
  ## Pre-Calculation
  n <- nrow(pi)
  prop_sd <- rep(0.05, J - 1)

  ## --- Partial-ranking bookkeeping ------------------------------------------
  ## r_max[i] = number of items voter i ranks (assumes NAs sit at the tail)
  r_max <- rowSums(!is.na(pi))
  stopifnot(all(r_max >= 1))

  ## Highest observed rank. Only stages 2, ..., R_obs have winners, so only
  ## those alpha_r are identified / estimated / kept. Active voters at stage r
  ## are {i : r_max >= r} (monotone in r), so the stages 2:R_obs have no gaps.
  R_obs        <- min(max(r_max), J - 1L)
  alpha_stages <- if (R_obs >= 2L) 2:R_obs else integer(0)

  ## pi_completed: append unranked items in setdiff(1:J, observed) order so
  ## each row is a permutation of 1:J. For any active voter (r_i >= r), the
  ## slice pi_completed[i, r:J] is still exactly the risk set R_ir.
  pi_completed <- pi
  for (i in 1:n) {
    if (r_max[i] < J) {
      unranked <- setdiff(1:J, pi[i, 1:r_max[i]])
      pi_completed[i, (r_max[i] + 1):J] <- unranked
    }
  }
  storage.mode(pi_completed) <- "integer"

  ## Win counts: position r of voter i is a "win" iff r <= min(r_i, J-1).
  r_seq  <- matrix(1:J, n, J, byrow = TRUE)
  win_pos <- r_seq <= pmin(r_max, J - 1L)                       # n x J logical
  n_j    <- tabulate(pi_completed[win_pos], nbins = J)

  ## pi_unique: group voters by their observed ranking (NAs preserved), then
  ## complete each unique class and store r_max (col J+2) for filtering.
  pi_df        <- as.data.frame(pi)
  pi_unique_df <- suppressMessages(
    pi_df %>% group_by_all() %>% summarize(count = n(), .groups = "drop"))
  pi_unique_partial <- as.matrix(pi_unique_df[, 1:J])
  pi_unique_counts  <- pi_unique_df$count
  pi_unique_rmax    <- rowSums(!is.na(pi_unique_partial))
  pi_unique_comp    <- pi_unique_partial
  for (i in 1:nrow(pi_unique_comp)) {
    r_i <- pi_unique_rmax[i]
    if (r_i < J) {
      unranked <- setdiff(1:J, pi_unique_comp[i, 1:r_i])
      pi_unique_comp[i, (r_i + 1):J] <- unranked
    }
  }
  storage.mode(pi_unique_comp) <- "integer"
  pi_unique <- cbind(pi_unique_comp, pi_unique_counts, pi_unique_rmax)
  colnames(pi_unique) <- NULL

  ## unique_partials: built only for the identified stages 2, ..., R_obs,
  ## restricted to voters active at each stage r (r_max >= r). Element [[r-1]]
  ## holds stage r. (Empty when R_obs < 2, i.e. only top-1 data.)
  unique_partials <- lapply(alpha_stages, function(r) {
    active_idx <- which(r_max >= r)
    if (length(active_idx) == 0) return(list())
    sorted_pi <- t(apply(pi_completed[active_idx, r:J, drop = FALSE], 1, sort))
    tmp <- unique(sorted_pi)
    lapply(1:nrow(tmp), function(ind) {
      in_class <- which(apply(t(sorted_pi) == tmp[ind, ], 2, all))
      list(sorted_partial_pi = tmp[ind, ],
           which_voters = active_idx[in_class])
    })
  })
  burn <- floor(iters * burnin)

  ## Step 1 Pre-Calculation ---------------------------------------------------
  rank_mat <- matrix(NA_integer_, n, J)
  for (r in 1:J) rank_mat[cbind(1:n, pi_completed[, r])] <- r

  ## Mask combines "at risk" with "active at stage r". Stages r > R_obs have no
  ## active voters, so they contribute nothing to the lambda-rate sum.
  active_ir <- outer(r_max, 1:(J - 1L), FUN = ">=")             # n x (J-1)
  mask <- array(FALSE, dim = c(n, J - 1L, J))
  for (r in 1:(J - 1L)) mask[, r, ] <- (rank_mat >= r) & active_ir[, r]

  ## Initializers
  S <- array(1, c(n, J - 1, J))
  V <- matrix(1, nrow = n, ncol = J - 1)
  alpha <- c(1, rep(.9, J - 2))

  ## Seed lambda so Step 1(a) has something to draw V[,1] from
  sum_sv <- colSums(matrix(S * as.vector(V) * mask, nrow = n * (J - 1L)))
  lambda <- rgamma(J, a + n_j, b + sum_sv)

  ## Storage
  steps <- accepts <- rep(0, J - 1)
  lambda_out <- matrix(NA, nrow = iters - burn, ncol = J)
  alpha_out  <- matrix(NA, nrow = iters - burn, ncol = J - 1)

  if(verbose) pb <- utils::txtProgressBar(min = 0, max = iters, style = 3)

  ## Gibbs Sweeps
  for (iter in 1:iters) {

    ## Step 1(a) -- every voter has r_max >= 1, so all participate at stage 1
    V[, 1] <- rexp(n, sum(lambda))

    ## Step 1(b)
    loglambda <- log(lambda)

    ### 1(b)(i) -- alpha update over identified stages 2, ..., R_obs only
    if (iter %% 50 == 0 & iter <= burn & length(alpha_stages) > 0) {
      acc <- accepts[alpha_stages] / steps[alpha_stages]
      prop_sd[alpha_stages] <- prop_sd[alpha_stages] *
        ifelse(acc < 0.35, 0.8, ifelse(acc > 0.50, 1.25, 1))
      prop_sd[prop_sd < 0.001] <- 0.001
      steps <- accepts <- rep(0, J - 1)
    }

    par_alpha <- lapply(alpha_stages, function(r) {
      log_dens <- make_log_dens_partial(pi_unique, loglambda, r, J, c, d)
      a        <- alpha[r]
      log_den  <- log_dens(a)
      accept   <- 0L

      for (t in 1:mh_iters) {
        a_prop <- a + rnorm(1, sd = prop_sd[r])
        if (a_prop <= alpha_min || a_prop >= 1) next
        log_num <- log_dens(a_prop)
        diff <- log_num - log_den
        if (!is.nan(diff) & log(runif(1)) < diff) {
          a       <- a_prop
          log_den <- log_num
          accept  <- accept + 1L
        }
      }
      c(mh_iters, accept, a)
    })
    for (r in alpha_stages) {
      steps[r]   <- steps[r] + par_alpha[[r - 1]][1]
      accepts[r] <- accepts[r] + par_alpha[[r - 1]][2]
      alpha[r]   <- par_alpha[[r - 1]][3]
    }

    ### 1(b)(ii) -- V draws for active voters, identified stages only
    for (r in alpha_stages) {
      up_r <- unique_partials[[r - 1]]
      if (length(up_r) == 0) next
      for (ind in seq_along(up_r)) {
        R_ir         <- up_r[[ind]]$sorted_partial_pi
        which_voters <- up_r[[ind]]$which_voters
        V[which_voters, r] <- rweibull(
          length(which_voters), alpha[r],
          exp(-logSumExp(alpha[r] * loglambda[R_ir]) / alpha[r]))
      }
    }

    ### 1(b)(iii) -- S draws for active voters, identified stages only
    for (r in alpha_stages) {
      al         <- alpha[r]
      active_idx <- which(r_max >= r)
      n_active   <- length(active_idx)
      if (n_active == 0) next

      risk     <- pi_completed[active_idx, r:J, drop = FALSE]   # n_a x (J-r+1)
      lam_risk <- matrix(lambda[risk], nrow = n_active)
      h        <- V[active_idx, r] * lam_risk                   # V broadcast

      S_ir <- matrix(retstable(alpha = al, V0 = as.vector(h^al)),
                     nrow = n_active) / h
      S_ir[, 1] <- S_ir[, 1] +
        rgamma(n_active, shape = 1 - al, rate = h[, 1])

      S[cbind(rep(active_idx, J - r + 1L), r, as.vector(risk))] <-
        as.vector(S_ir)
    }

    ## Step 2
    sum_sv <- colSums(matrix(S * as.vector(V) * mask, nrow = n * (J - 1L)))
    lambda <- rgamma(J, a + n_j, b + sum_sv)

    ## Saving Samples (normalized at save time)
    if (iter > burn) {
      lambda_out[iter - burn, ] <- lambda / sum(lambda)
      alpha_out[iter - burn, ]  <- alpha
    }

    if(verbose) utils::setTxtProgressBar(pb, iter)
  }

  if(verbose) close(pb)

  ## Drop the unidentified dampening parameters: alpha_r, r > R_obs, had no
  ## data and were never updated, so return them as NA.
  if (R_obs < (J - 1L)) alpha_out[, (R_obs + 1L):(J - 1L)] <- NA_real_

  t1 <- Sys.time()

  ## Per-stage acceptance (NA where no MH steps were taken, i.e. r > R_obs)
  accept_prob <- ifelse(steps[-1] > 0, accepts[-1] / steps[-1], NA_real_)

  ## Returns
  list(lambda = lambda_out,
       alpha  = alpha_out,
       parameters = list(a = a, b = b, c = c, d = d,
                         alpha_min = alpha_min, R_obs = R_obs,
                         iters = iters, mh_iters = mh_iters, burnin = burnin,
                         prop_sd = prop_sd, accept_prob = accept_prob),
       time_secs = as.numeric(difftime(t1, t0, units = "secs")))
}
