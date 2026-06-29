local({
  surv <- Surv(.sim_dat$time, .sim_dat$event)
  .fit <- fit_genhaz(surv, ~X, data = .sim_dat, n_knots = .n_knots, lambda = 5)

  test_that("confint returns lower < upper", {
    ci <- confint(.fit, "beta2_X")
    expect_lt(ci["beta2_X", 1], ci["beta2_X", 2])
  })

  test_that("confint contains the point estimate", {
    ci <- confint(.fit, "beta2_X")
    expect_gte(.fit$par["beta2_X"], ci["beta2_X", 1])
    expect_lte(.fit$par["beta2_X"], ci["beta2_X", 2])
  })

  test_that("confint with no parm covers all parameters", {
    ci <- confint(.fit)
    expect_equal(rownames(ci), .fit$parnames)
  })

  test_that("confint errors on unknown parameter", {
    expect_error(confint(.fit, "nonexistent_param"))
  })

  test_that("confint(diff = TRUE) returns lower < upper", {
    ci <- confint(.fit, c("beta1_X", "beta2_X"), diff = TRUE)
    expect_equal(nrow(ci), 1L)
    expect_lt(ci[1, 1], ci[1, 2])
  })

  test_that("anova returns a valid deviance table", {
    # Chisq = 2*(negll_nested - negll_general); should be >= 0 at global optima,
    # but small-sample numerical optimization may yield small negative values.
    fit_gh0 <- fit_genhaz(Surv(.sim_dat$time, .sim_dat$event), ~X,
                          data = .sim_dat, n_knots = .n_knots, lambda = 0)
    fit_ph0 <- fit_genhaz(Surv(.sim_dat$time, .sim_dat$event), ~X,
                          data = .sim_dat, n_knots = .n_knots,
                          lambda = 0, model_type = "PH")
    tab <- anova(fit_ph0, fit_gh0)
    expect_s3_class(tab, "anova")
    expect_equal(nrow(tab), 2L)
    expect_true(is.finite(tab[2L, "Chisq"]))
    expect_gte(tab[2L, "Pr(>Chisq)"], 0)
    expect_lte(tab[2L, "Pr(>Chisq)"], 1)
  })

  test_that("anova errors with a single model", {
    expect_error(anova(.fit))
  })
})
