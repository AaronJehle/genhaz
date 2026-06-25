# Print a fitted GH model

Displays a compact summary: model metadata and a covariate coefficient
table (spline parameters are suppressed as they are not directly
interpretable).

## Usage

``` r
# S3 method for class 'genhaz_fit'
print(x, digits = 4, alpha = 0.05, ...)
```

## Arguments

- x:

  A `"genhaz_fit"` object from
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

- digits:

  Number of significant digits for numeric output. Default 4.

- alpha:

  Significance level for Wald confidence intervals. Default 0.05.

- ...:

  Currently unused.

## Value

Invisibly returns `x`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_genhaz(Surv(dat$time, dat$event), ~X, data = dat,
                  model_type = "GH", profile = TRUE, n_knots = 6)
print(fit)
} # }
```
