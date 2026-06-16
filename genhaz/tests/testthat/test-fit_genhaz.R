test_that("fit_genhaz returns a list with expected components", {
  surv <- Surv(.sim_dat$time, .sim_dat$event)
  fit  <- fit_genhaz(surv, ~X, data = .sim_dat, n_knots = .n_knots,
                     lambda = 5, model_type = "GH")
  expect_s3_class(fit, "genhaz_fit")
  expect_true(all(c("par", "se", "z", "p_values", "var", "AIC",
                    "edf", "lambda", "knots", "Z") %in% names(fit)))
})

test_that("fit_genhaz parameter names match formula covariates", {
  surv <- Surv(.sim_dat$time, .sim_dat$event)
  fit  <- fit_genhaz(surv, ~X, data = .sim_dat, n_knots = .n_knots, lambda = 5)
  expect_true("beta1_X" %in% names(fit$par))
  expect_true("beta2_X" %in% names(fit$par))
})

test_that("fit_genhaz se vector has no NAs", {
  surv <- Surv(.sim_dat$time, .sim_dat$event)
  fit  <- fit_genhaz(surv, ~X, data = .sim_dat, n_knots = .n_knots, lambda = 5)
  expect_true(all(!is.na(fit$se)))
})

test_that("fit_genhaz with profile = TRUE returns a finite lambda", {
  surv <- Surv(.sim_dat$time, .sim_dat$event)
  fit  <- fit_genhaz(surv, ~X, data = .sim_dat, n_knots = .n_knots,
                     profile = TRUE, tol_LCV = 0.05)
  expect_true(is.finite(fit$lambda))
  expect_gt(fit$lambda, 0)
})

test_that("fit_genhaz PH model has beta1 = 0", {
  surv <- Surv(.sim_dat$time, .sim_dat$event)
  fit  <- fit_genhaz(surv, ~X, data = .sim_dat, n_knots = .n_knots,
                     lambda = 5, model_type = "PH")
  expect_equal(unname(fit$par["beta1_X"]), 0)
})

test_that("fit_genhaz AFT model has beta2 = 0", {
  surv <- Surv(.sim_dat$time, .sim_dat$event)
  fit  <- fit_genhaz(surv, ~X, data = .sim_dat, n_knots = .n_knots,
                     lambda = 5, model_type = "AFT")
  expect_equal(unname(fit$par["beta2_X"]), 0)
})

test_that("post() returns h and H of correct length", {
  surv <- Surv(.sim_dat$time, .sim_dat$event)
  fit  <- fit_genhaz(surv, ~X, data = .sim_dat, n_knots = .n_knots, lambda = 5)
  t_grid <- seq(0.1, 5, length.out = 50)
  hv <- post(fit, X = 0, time = t_grid, res = "h")
  Hv <- post(fit, X = 0, time = t_grid, res = "H")
  expect_length(hv, 50)
  expect_length(Hv, 50)
  expect_true(all(hv > 0))
  expect_true(all(Hv > 0))
})

test_that("CI() returns a data frame with 10 columns and lower<upper", {
  surv <- Surv(.sim_dat$time, .sim_dat$event)
  fit  <- fit_genhaz(surv, ~X, data = .sim_dat, n_knots = .n_knots, lambda = 5)
  t_grid <- seq(0.1, 5, length.out = 20)
  ci  <- CI(fit, t_grid, covariate = 0)
  expect_s3_class(ci, "data.frame")
  expect_equal(ncol(ci), 10)
  expect_true(all(ci$lower_h <= ci$h + 1e-10))
  expect_true(all(ci$upper_h >= ci$h - 1e-10))
})
