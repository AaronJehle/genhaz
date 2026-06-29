# Wald confidence intervals for model parameters

Computes Wald confidence intervals for the parameters of a fitted GH
model, following the
[`stats::confint()`](https://rdrr.io/r/stats/confint.html) convention.
With `diff = TRUE` it instead returns a delta-method interval for the
difference of two parameters, accounting for their covariance.

## Usage

``` r
# S3 method for class 'genhaz_fit'
confint(object, parm, level = 0.95, diff = FALSE, ...)
```

## Arguments

- object:

  A `"genhaz_fit"` object from
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

- parm:

  Character vector of parameter names (must be in `object$parnames`). If
  missing, all parameters are used. When `diff = TRUE`, exactly two
  names must be supplied.

- level:

  Confidence level. Default `0.95`.

- diff:

  Logical. If `TRUE`, return a single interval for `parm[1] - parm[2]`
  via the delta method. Default `FALSE`.

- ...:

  Currently unused.

## Value

A matrix with one row per parameter (or a single row for the difference)
and two columns giving the lower and upper bounds.

## Examples

``` r
if (FALSE) { # \dontrun{
confint(fit)                                       # all parameters
confint(fit, "beta2_X")                            # one parameter
confint(fit, c("beta1_X", "beta2_X"), diff = TRUE) # difference
} # }
```
