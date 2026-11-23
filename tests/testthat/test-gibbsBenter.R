test_that("rBenter returns valid permutations", {
  set.seed(1)
  lambda <- c(4, 2, 1, 1) / 8
  alpha  <- c(1, 0.7, 0.4)
  pi <- rBenter(50, lambda, alpha)

  expect_equal(dim(pi), c(50, 4))
  expect_true(all(apply(pi, 1, function(z) identical(sort(z), 1:4))))
})

test_that("gibbsBenter returns the expected structure on complete data", {
  set.seed(1)
  lambda <- c(4, 2, 1, 1) / 8
  alpha  <- c(1, 0.9, 0.8)
  pi <- rBenter(300, lambda, alpha)

  fit <- gibbsBenter(pi, J = 4, a = 1, b = 1, c = 1, d = 1,
                     iters = 120, burnin = 0.5)

  expect_named(fit, c("lambda", "alpha", "parameters", "time_secs"))
  expect_equal(dim(fit$lambda), c(60, 4))
  expect_equal(dim(fit$alpha),  c(60, 3))
  # worths are normalized to sum to one
  expect_equal(unname(rowSums(fit$lambda)), rep(1, 60), tolerance = 1e-8)
  # alpha_1 is fixed at 1
  expect_true(all(fit$alpha[, 1] == 1))
})

test_that("partial rankings drop unidentified dampening parameters", {
  set.seed(2)
  lambda <- c(4, 2, 1, 1, 1) / 9
  alpha  <- c(1, 0.7, 0.5, 0.3)
  pi <- rBenter(500, lambda, alpha)
  pi[, 4:5] <- NA                     # top-3 rankings: deepest observed rank 3

  fit <- gibbsBenter(pi, J = 5, a = 1, b = 1, c = 1, d = 1,
                     iters = 120, burnin = 0.5)

  expect_equal(fit$parameters$R_obs, 3)
  # alpha_4 (column 4) is unidentified -> NA; alpha_1..3 are finite
  expect_true(all(is.na(fit$alpha[, 4])))
  expect_true(all(is.finite(fit$alpha[, 1:3])))
  # all J worths are still estimated
  expect_equal(dim(fit$lambda), c(60, 5))
})

test_that("input validation rejects malformed rankings", {
  expect_error(gibbsBenter(matrix(c(1, 2, 2), 1), J = 3, a = 1, b = 1, c = 1, d = 1),
               "repeats")
  expect_error(gibbsBenter(matrix(1:3, 1), J = 2, a = 1, b = 1, c = 1, d = 1),
               "J must be")
})
