# Simulate data for a named scenario

Convenience wrapper around
[`sim_mix_weib_gh()`](https://aaronjehle.github.io/genhaz/reference/sim_mix_weib_gh.md)
using the three pre-defined simulation scenarios (see
[`mixWeibSc()`](https://aaronjehle.github.io/genhaz/reference/mixWeibSc.md)).

## Usage

``` r
sim_scenario(scenario, beta1, beta2, n = 1000, tmax = 10)
```

## Arguments

- scenario:

  Integer (1, 2, or 3).

- beta1:

  AFT coefficient.

- beta2:

  PH coefficient.

- n:

  Sample size.

- tmax:

  Administrative censoring time.

## Value

Data frame with columns `time`, `X`, `event`, `T_true`.

## Examples

``` r
set.seed(1)
dat <- sim_scenario(1, beta1 = 0.5, beta2 = 0.5, n = 500)
table(dat$event)
#> 
#>   0   1 
#>  94 406 
```
