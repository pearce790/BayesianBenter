# Bayesian estimation of the Benter ranking model

Fits the Benter (1994) model for ranking data by Gibbs sampling. The
Benter model generalizes the Plackett-Luce model by introducing
stage-specific *dampening* exponents \\\alpha_r \in (0, 1\]\\ that
down-weight the information in later ranking stages: \$\$p(\pi_i \mid
\lambda, \alpha) = \prod\_{r=1}^{J-1}
\frac{\lambda\_{w\_{ir}}^{\alpha_r}}{\sum\_{j \in R\_{ir}}
\lambda_j^{\alpha_r}},\$\$ where \\\lambda_j \> 0\\ are item worths,
\\w\_{ir}\\ is the item placed \\r\\th by voter \\i\\, and \\R\_{ir}\\
is the set of items still at risk at stage \\r\\. Inference uses a
positive-stable frailty data augmentation that restores conjugacy, so
the worths are drawn exactly from Gamma full conditionals while each
\\\alpha_r\\ is updated by Metropolis-Hastings on its collapsed
conditional.

## Usage

``` r
gibbsBenter(
  pi,
  J,
  a,
  b,
  c,
  d,
  alpha_min = 0,
  iters = 1000,
  mh_iters = 5,
  burnin = 0.5,
  verbose = FALSE
)
```

## Arguments

- pi:

  A ranking data set: an integer matrix (or data frame, or a single
  vector interpreted as one ranking) with one ranking per row. Entry
  `pi[i, r]` is the item placed `r`th by voter `i`; ranked entries are
  left-aligned and any unranked tail positions are `NA`. All non-missing
  entries must be integers in `1:J` with no repeats within a row.

- J:

  Integer (`>= 3`). The total number of items being ranked.

- a, b:

  Shape and rate of the Gamma(`a`, `b`) prior on the item worths
  \\\lambda_j\\. On the identified (sum-to-one) scale, `a` alone
  controls the induced symmetric Dirichlet prior on the normalized
  worths; `b` only sets the (marginalized-out) overall scale.

- c, d:

  Shape parameters of the Beta(`c`, `d`) prior on each dampening
  exponent \\\alpha_r\\, `r`=2,...,J-1. (`alpha_1` is fixed at 1 for
  identifiability.)

- alpha_min:

  Lower floor for the Metropolis-Hastings proposals of \\\alpha_r\\.
  Proposals at or below `alpha_min` are rejected. A small positive value
  (e.g. `0.05`) keeps `alpha_r` inside the range where the tilted-stable
  sampler is numerically stable; the default `0` imposes no floor beyond
  the open boundary.

- iters:

  Number of Gibbs iterations (including burn-in).

- mh_iters:

  Number of Metropolis-Hastings sub-steps per \\\alpha_r\\ per sweep.

- burnin:

  Fraction of `iters` discarded as burn-in. The adaptive tuning of the
  proposal scales is active during burn-in and frozen afterwards.

- verbose:

  Logical. If `TRUE`, show a text progress bar.

## Value

A list with components:

- `lambda`:

  A `(iters - burn) x J` matrix of posterior draws of the item worths,
  normalized to sum to one in each row.

- `alpha`:

  A `(iters - burn) x (J-1)` matrix of posterior draws of the dampening
  exponents. Column 1 is fixed at 1; for partial data, columns beyond
  the deepest observed rank are `NA`.

- `parameters`:

  A list recording the prior and sampler settings and the
  Metropolis-Hastings acceptance rates.

- `time_secs`:

  Wall-clock run time in seconds.

## Details

`gibbsBenter()` is the user-facing entry point. It validates and
standardizes `pi`, then dispatches to the appropriate internal sampler
for complete or partial rankings. A length-(J-1) ranking pins down its
final item and is completed automatically (this is deduction, not
imputation); a row is treated as genuinely partial only if it leaves two
or more items unplaced.

For partial rankings, only the dampening parameters up to the deepest
observed rank are identifiable, so `alpha_r` for stages beyond the
longest observed ranking are neither updated nor returned (reported as
`NA`). The item worths \\\lambda\\ for all `J` items are always
estimated.

## References

Benter, W. (1994). Computer-based horse race handicapping and wagering
systems: A report. In *Efficiency of Racetrack Betting Markets*.

## See also

[`rBenter()`](https://pearce790.github.io/BayesianBenter/reference/rBenter.md)
to simulate rankings from the model.

## Examples

``` r
set.seed(1)
lambda <- c(5, 3, 2, 1.5, 1)
lambda <- lambda / sum(lambda)
alpha  <- c(1, 0.9, 0.8, 0.7)
pi <- rBenter(200, lambda, alpha)

# \donttest{
fit <- gibbsBenter(pi, J = 5, a = 1, b = 1, c = 1, d = 1,
                   iters = 400, burnin = 0.5)
colMeans(fit$lambda)   # posterior mean worths (vs. truth `lambda`)
#> [1] 0.41990139 0.20843699 0.18133789 0.10639968 0.08392406
colMeans(fit$alpha)    # posterior mean dampening (vs. truth `alpha`)
#> [1] 1.0000000 0.8617744 0.7392397 0.6400031
# }
```
