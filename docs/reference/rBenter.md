# Simulate complete rankings from the Benter model

Draws `n` independent complete rankings of `J = length(lambda)` items
from the Benter model with worths `lambda` and dampening exponents
`alpha`. At each stage `r`, the next item is chosen from those still
available with probability proportional to `lambda^alpha[r]`.

## Usage

``` r
rBenter(n, lambda, alpha)
```

## Arguments

- n:

  Number of rankings to simulate.

- lambda:

  Positive vector of item worths, length `J`. Need not sum to one (only
  relative magnitudes matter).

- alpha:

  Vector of dampening exponents, length `J - 1`, with `alpha[1] == 1`.

## Value

An `n x J` integer matrix; row `i` is a complete ranking, with
`pi[i, r]` the item placed `r`th.

## See also

[`gibbsBenter()`](https://pearce790.github.io/BayesianBenter/reference/gibbsBenter.md)
to fit the model.

## Examples

``` r
set.seed(1)
lambda <- c(5, 3, 2, 1.5, 1) / sum(c(5, 3, 2, 1.5, 1))
alpha  <- c(1, 0.8, 0.5, 0.3)
rBenter(5, lambda, alpha)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    3    4    5    2
#> [2,]    1    5    4    3    2
#> [3,]    2    1    3    4    5
#> [4,]    3    1    5    2    4
#> [5,]    3    5    1    4    2
```
