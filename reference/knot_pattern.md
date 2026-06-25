# Compute knot positions from event times

Places `n_knots` equally spaced (in probability) quantiles of the
log-event-times. Boundary knots are placed at the 5 % and, optionally,
the 95 % or 100 % quantile of the event-time distribution.

## Usage

``` r
knot_pattern(time, event, n_knots, limit_to_95 = TRUE)
```

## Arguments

- time:

  Numeric vector of observed times.

- event:

  Integer vector of event indicators (1 = event, 0 = censored), same
  length as `time`.

- n_knots:

  Total number of knots, including the two boundary knots.

- limit_to_95:

  Logical. If `TRUE` (default) the upper boundary knot is placed at the
  95 % quantile; if `FALSE` it is placed at the maximum observed event
  time.

## Value

Named numeric vector of length `n_knots` with knot positions on the
*log-time* scale, as required by
[`genhaz_work()`](https://aaronjehle.github.io/genhaz/reference/genhaz_work.md)
and
[`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

## Examples

``` r
set.seed(1)
t  <- rexp(200)
ev <- rbinom(200, 1, 0.8)
knot_pattern(t, ev, n_knots = 6)
#>          5%         23%         41%         59%         77%         95% 
#> -2.39966844 -1.09595798 -0.43242925 -0.01986458  0.31748801  0.99421752 
```
