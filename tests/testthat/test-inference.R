local({
  surv <- Surv(.sim_dat$time, .sim_dat$event)
  .fit <- fit_genhaz(surv, ~X, data = .sim_dat, n_knots = .n_knots, lambda = 5)

  test_that("waldCI returns lower < upper", {
    ci <- waldCI(.fit, "beta2_X")
    expect_lt(ci["lower"], ci["upper"])
  })

  test_that("waldCI contains the point estimate", {
    ci <- waldCI(.fit, "beta2_X")
    expect_gte(.fit$par["beta2_X"], ci["lower"])
    expect_lte(.fit$par["beta2_X"], ci["upper"])
  })

  test_that("waldCI errors on unknown parameter", {
    expect_error(waldCI(.fit, "nonexistent_param"))
  })

  test_that("waldCI_minus returns lower < upper", {
    ci <- waldCI_minus(.fit, "beta1_X", "beta2_X")
    expect_lt(ci["lower"], ci["upper"])
  })

  test_that("LR test returns finite statistic and valid p-value", {
    # LR = 2*(negll_nested - negll_general); should be >= 0 at global optima,
    # but small-sample numerical optimization may yield small negative values.
    fit_gh0 <- fit_genhaz(Surv(.sim_dat$time, .sim_dat$event), ~X,
                          data = .sim_dat, n_knots = .n_knots, lambda = 0)
    fit_ph0 <- fit_genhaz(Surv(.sim_dat$time, .sim_dat$event), ~X,
                          data = .sim_dat, n_knots = .n_knots,
                          lambda = 0, model_type = "PH")
    lr <- LR(fit_ph0, fit_gh0)
    expect_true(is.finite(lr["LR-statistic"]))
    expect_gte(lr["p_value"], 0)
    expect_lte(lr["p_value"], 1)
  })
})
