# Wald confidence interval for a single parameter

Wald confidence interval for a single parameter

## Usage

``` r
waldCI(fit, param, alpha = 0.05)
```

## Arguments

- fit:

  A fitted model object from
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

- param:

  Character; name of the parameter (must be in `fit$parnames`).

- alpha:

  Significance level. Default 0.05 (95 % CI).

## Value

Named numeric vector with `lower` and `upper` bounds.

## Examples

``` r
if (FALSE) { # \dontrun{
waldCI(fit, "beta2_X")
} # }
```
