# Pointwise confidence bands for hazard and survival functions

Computes pointwise 95 % (or other level) confidence intervals for the
cumulative hazard H(t), survival S(t) = exp(-H(t)), and hazard h(t), at
a given covariate pattern, using the delta method on the log scale.

## Usage

``` r
CI(fit, t, covariate, alpha = 0.05)
```

## Arguments

- fit:

  A fitted model object from
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

- t:

  Numeric vector of evaluation time points.

- covariate:

  Numeric vector of covariate values (one per covariate in the model).

- alpha:

  Significance level. Default 0.05 (95 % CI).

## Value

A data frame with columns `time`, `H`, `lower_H`, `upper_H`, `S`,
`lower_S`, `upper_S`, `h`, `lower_h`, `upper_h`.

## Examples

``` r
if (FALSE) { # \dontrun{
t_grid <- seq(0.01, 8, length.out = 200)
ci     <- CI(fit, t_grid, covariate = 0)
plot(ci$time, ci$h, type = "l")
lines(ci$time, ci$lower_h, lty = 2)
lines(ci$time, ci$upper_h, lty = 2)
} # }
```
