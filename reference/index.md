# Package index

## Model fitting

- [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md)
  : Fit a penalized general hazard model
- [`genhaz_work()`](https://aaronjehle.github.io/genhaz/reference/genhaz_work.md)
  : General hazard model workhorse

## S3 methods

- [`print(`*`<genhaz_fit>`*`)`](https://aaronjehle.github.io/genhaz/reference/print.genhaz_fit.md)
  : Print a fitted GH model
- [`summary(`*`<genhaz_fit>`*`)`](https://aaronjehle.github.io/genhaz/reference/summary.genhaz_fit.md)
  : Summarise a fitted GH model
- [`predict(`*`<genhaz_fit>`*`)`](https://aaronjehle.github.io/genhaz/reference/predict.genhaz_fit.md)
  : Predict from a fitted GH model
- [`plot(`*`<genhaz_fit>`*`)`](https://aaronjehle.github.io/genhaz/reference/plot.genhaz_fit.md)
  : Plot estimated curves from a fitted GH model

## Inference

- [`CI()`](https://aaronjehle.github.io/genhaz/reference/CI.md) :
  Pointwise confidence bands for hazard and survival functions
- [`waldCI()`](https://aaronjehle.github.io/genhaz/reference/waldCI.md)
  : Wald confidence interval for a single parameter
- [`waldCI_minus()`](https://aaronjehle.github.io/genhaz/reference/waldCI_minus.md)
  : Wald confidence interval for the difference of two parameters
- [`LR()`](https://aaronjehle.github.io/genhaz/reference/LR.md) :
  Likelihood ratio test between nested GH models

## Post-processing

- [`post()`](https://aaronjehle.github.io/genhaz/reference/post.md) :
  Post-estimation from a fitted general hazard model

## Simulation

- [`sim_scenario()`](https://aaronjehle.github.io/genhaz/reference/sim_scenario.md)
  : Simulate data for a named scenario
- [`mixWeibSc()`](https://aaronjehle.github.io/genhaz/reference/mixWeibSc.md)
  : Mixture Weibull model for a named simulation scenario
- [`mixWeib()`](https://aaronjehle.github.io/genhaz/reference/mixWeib.md)
  : Mixture Weibull GH model quantities
- [`sim_mix_weib_gh()`](https://aaronjehle.github.io/genhaz/reference/sim_mix_weib_gh.md)
  : Simulate data from a mixture Weibull GH model

## Utilities

- [`knot_pattern()`](https://aaronjehle.github.io/genhaz/reference/knot_pattern.md)
  : Compute knot positions from event times
- [`plot_hazard()`](https://aaronjehle.github.io/genhaz/reference/plot_hazard.md)
  : Plot estimated hazard function
- [`translate_time()`](https://aaronjehle.github.io/genhaz/reference/translate_time.md)
  : Time-varying acceleration factor between two covariate patterns
- [`translate_time2()`](https://aaronjehle.github.io/genhaz/reference/translate_time2.md)
  : Time-varying tau mapping between two covariate patterns

## Data

- [`fit_melanoma`](https://aaronjehle.github.io/genhaz/reference/fit_melanoma.md)
  : Pre-computed melanoma GH model fit
- [`mixture_weibull`](https://aaronjehle.github.io/genhaz/reference/mixture_weibull.md)
  : Mixture Weibull survival, cumulative hazard, and hazard functions

## Internal (C++ / low-level)

- [`print(`*`<summary.genhaz_fit>`*`)`](https://aaronjehle.github.io/genhaz/reference/print.summary.genhaz_fit.md)
  : Print a summary of a fitted GH model

- [`crs_cpp()`](https://aaronjehle.github.io/genhaz/reference/crs_cpp.md)
  : C++ reimplementation of survPen::crs, including derivatives.

- [`smf_cpp()`](https://aaronjehle.github.io/genhaz/reference/smf_cpp.md)
  :

  C++ reimplementation of `smf`, including derivatives.

- [`quantile_type7()`](https://aaronjehle.github.io/genhaz/reference/quantile_type7.md)
  :

  Armadillo-based reimplementation of
  [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html) for
  `type=7`.
