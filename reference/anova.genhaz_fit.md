# Likelihood ratio tests between nested GH models

Compares two or more nested `"genhaz_fit"` models with chi-squared
likelihood-ratio tests, in the style of
[`stats::anova()`](https://rdrr.io/r/stats/anova.html) for `lm`/`glm`
objects. Models are compared sequentially in the order supplied, so list
the most restricted model first.

## Usage

``` r
# S3 method for class 'genhaz_fit'
anova(object, ...)
```

## Arguments

- object:

  A `"genhaz_fit"` object from
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md)
  (the first, usually most restricted, model).

- ...:

  One or more further `"genhaz_fit"` objects to compare against
  `object`. Non-`"genhaz_fit"` arguments are ignored.

## Value

An analysis-of-deviance table of class `c("anova", "data.frame")` with
one row per model and columns `Df` (model degrees of freedom), `LogLik`,
`Chisq`, `Chi Df`, and `Pr(>Chisq)` (the comparison columns are `NA` for
the first model).

## Examples

``` r
if (FALSE) { # \dontrun{
fit_ph <- fit_genhaz(Surv(dat$time, dat$event), ~X,
                     data = dat, model_type = "PH", n_knots = 6)
fit_gh <- fit_genhaz(Surv(dat$time, dat$event), ~X,
                     data = dat, model_type = "GH", n_knots = 6)
anova(fit_ph, fit_gh)
} # }
```
