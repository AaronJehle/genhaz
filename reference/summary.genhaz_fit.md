# Summarise a fitted GH model

Returns a `"summary.genhaz_fit"` object containing model metadata and a
covariate coefficient table with Wald confidence intervals and
exponentiated estimates.

## Usage

``` r
# S3 method for class 'genhaz_fit'
summary(object, alpha = 0.05, ...)
```

## Arguments

- object:

  A `"genhaz_fit"` object from
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

- alpha:

  Significance level for confidence intervals. Default 0.05.

- ...:

  Currently unused.

## Value

An object of class `"summary.genhaz_fit"` (printed by
[`print.summary.genhaz_fit()`](https://aaronjehle.github.io/genhaz/reference/print.summary.genhaz_fit.md)).

## Examples

``` r
if (FALSE) { # \dontrun{
summary(fit)
} # }
```
