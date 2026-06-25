# Wald confidence interval for the difference of two parameters

Uses the delta method to compute a confidence interval for
`param1 - param2` accounting for their covariance.

## Usage

``` r
waldCI_minus(fit, param1, param2, alpha = 0.05)
```

## Arguments

- fit:

  A fitted model object from
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

- param1:

  Character; name of the first parameter.

- param2:

  Character; name of the second parameter.

- alpha:

  Significance level. Default 0.05 (95 % CI).

## Value

Named numeric vector with `lower` and `upper` bounds.

## Examples

``` r
if (FALSE) { # \dontrun{
waldCI_minus(fit, "beta1_X", "beta2_X")
} # }
```
