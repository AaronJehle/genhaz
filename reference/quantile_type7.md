# Armadillo-based reimplementation of `stats::quantile()` for `type=7`.

Allows for NAs.

## Usage

``` r
quantile_type7(x, probs)
```

## Arguments

- x:

  a vec of values for which to find the quantiles

- probs:

  a vec of probilities that specify which quantiles

## Value

a vec of quantiles

## Examples

``` r
set.seed(12345)
y <- rnorm(1e5)
# Should be near zero for all quantiles
quantile(y) - quantile_type7(y, seq(0, 1, length.out = 5))
#>      [,1]
#> [1,]    0
#> [2,]    0
#> [3,]    0
#> [4,]    0
#> [5,]    0
quantile(y, 0.975) - quantile_type7(y, 0.975)
#>      [,1]
#> [1,]    0
# NA probs are allowed (unlike base quantile()):
quantile_type7(y, probs = c(0.1, 0.5, 1, 2, 5, 10, 50, NA) / 100)
#>              [,1]
#> [1,] -3.105140881
#> [2,] -2.588471928
#> [3,] -2.332039437
#> [4,] -2.056505903
#> [5,] -1.642808669
#> [6,] -1.282869942
#> [7,]  0.003099443
#> [8,]          NaN
```
