# Likelihood ratio test between nested GH models

Computes a chi-squared likelihood ratio test statistic comparing a
nested (restricted) model against a more general model.

## Usage

``` r
LR(fit_nested, fit_general)
```

## Arguments

- fit_nested:

  Fitted model object (the restricted model).

- fit_general:

  Fitted model object (the general model).

## Value

Named numeric vector with `LR-statistic` and `p_value`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit_ph <- fit_genhaz(Surv(dat$time, dat$event), ~X,
                     data = dat, model_type = "PH", n_knots = 6)
fit_gh <- fit_genhaz(Surv(dat$time, dat$event), ~X,
                     data = dat, model_type = "GH", n_knots = 6)
LR(fit_ph, fit_gh)
} # }
```
