# Plot estimated hazard function

Quick base-graphics plot of the estimated hazard function from a fitted
GH model at a given covariate pattern.

## Usage

``` r
plot_hazard(fit, covariate, time, xlim = NULL, ylim = NULL, ...)
```

## Arguments

- fit:

  A fitted model object from
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

- covariate:

  Numeric vector of covariate values.

- time:

  Numeric vector of time points at which to evaluate the hazard.

- xlim:

  Optional x-axis limits passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- ylim:

  Optional y-axis limits passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- ...:

  Additional graphical parameters passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns the estimated hazard values.

## Examples

``` r
if (FALSE) { # \dontrun{
t_grid <- seq(0.01, 8, length.out = 300)
plot_hazard(fit, covariate = 0, time = t_grid)
} # }
```
