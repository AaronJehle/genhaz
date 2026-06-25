# Plot a genhaz prediction

Plots the output of
[`predict.genhaz_fit()`](https://aaronjehle.github.io/genhaz/reference/predict.genhaz_fit.md)
— any of the eight prediction types. Each distinct covariate pattern
becomes a coloured line; confidence bands are drawn as dashed lines of
the same colour. `ylim` is always derived from the full range of CI
bands plus the point estimate, so nothing is clipped. Axis labels and
title are set automatically from the prediction type and can be
overridden via `ylab` / `main`.

## Usage

``` r
# S3 method for class 'genhaz_pred'
plot(
  x,
  col = NULL,
  lty = 1,
  xlab = "Time",
  ylab = NULL,
  main = NULL,
  legend = TRUE,
  ci = TRUE,
  ...
)
```

## Arguments

- x:

  A `"genhaz_pred"` object returned by
  [`predict.genhaz_fit()`](https://aaronjehle.github.io/genhaz/reference/predict.genhaz_fit.md).

- col:

  Colour vector (recycled over groups). Defaults to `2, 3, 4, ...` (R's
  standard colour sequence, skipping black).

- lty:

  Line type for point-estimate curves. Default `1`.

- xlab:

  X-axis label. Default `"Time"`.

- ylab:

  Y-axis label. Derived from the prediction type when `NULL`.

- main:

  Plot title. Derived from the prediction type when `NULL`.

- legend:

  Logical; draw a legend when multiple groups are present? Default
  `TRUE`.

- ci:

  Logical; overlay dashed confidence band lines? Default `TRUE`.

- ...:

  Additional graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `x`.

## Examples

``` r
if (FALSE) { # \dontrun{
nd <- data.frame(X = c(0, 1))
rownames(nd) <- c("X = 0", "X = 1")
t_grid <- seq(0.01, 8, length.out = 200)

pred <- predict(fit, newdata = nd, times = t_grid, type = "hazard_ratio")
plot(pred, col = "purple")
abline(h = 1, lty = 2)
} # }
```
