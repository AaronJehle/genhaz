# Tests for predict.genhaz_fit(), focused on the time-varying acceleration
# EFFECT type = "acc_factor".
#
# Direction (S(t|x1) = S0(tau(t)), i.e. H0(tau) = H1(t)): observe row 1 at t,
# warp row 2. phi(t) = tau'(t) = h_row1(t) / h_row2(tau), and acc_factor returns
# the log-scale effect f(t) = log phi(t) (= beta1 in the constant-AFT limit),
# with a symmetric delta-method CI f +/- z*sqrt(var). time_ratio returns tau/t.
#
# Fixture: row 1 = X = 0, row 2 = X = 1 (the warped row has the higher hazard
# here, so tau stays within the model support for the simulated fit).

.fit_gen <- local({
  surv <- Surv(.sim_dat$time, .sim_dat$event)
  fit_genhaz(surv, ~X, data = .sim_dat, n_knots = .n_knots,
             lambda = 5, model_type = "GH")
})

.nd_gen <- {
  nd <- data.frame(X = c(0, 1))
  rownames(nd) <- c("row1", "row2")
  nd
}

test_that("acc_factor estimate equals f = log(h1(t)/h0(tau)) with H0(tau)=H1(t)", {
  fit   <- .fit_gen
  times <- c(0.5, 1, 2, 3)

  pr <- predict(fit, newdata = .nd_gen, times = times, type = "acc_factor")

  # manual: observe row 1 (X = 0) at t, warp row 2 (X = 1)
  f_manual <- vapply(times, function(t) {
    H1t <- post(fit, 0, t, "H")
    tau <- uniroot(function(u) post(fit, 1, u, "H") - H1t,
                   lower = 1e-8, upper = max(times) * 20, tol = 1e-10)$root
    log(post(fit, 0, t, "h") / post(fit, 1, tau, "h"))
  }, numeric(1))

  expect_equal(pr$estimate, f_manual, tolerance = 1e-5)
  # symmetric CI on the (log) effect scale
  expect_equal(pr$estimate, (pr$lower + pr$upper) / 2, tolerance = 1e-8)
})

test_that("exp(acc_factor) equals the numerical derivative of tau", {
  # exp(f(t)) = phi(t) = tau'(t). The identity holds up to the Gauss-Legendre
  # quadrature accuracy of H, so we check on an interior grid (~1e-6 agreement).
  fit   <- .fit_gen
  times <- c(0.5, 1, 2)
  eps   <- 1e-4

  tau_of <- function(t) vapply(t, function(tt)
    uniroot(function(u) post(fit, 1, u, "H") - post(fit, 0, tt, "H"),
            lower = 1e-8, upper = 60, tol = 1e-12)$root, numeric(1))

  pr      <- predict(fit, newdata = .nd_gen, times = times, type = "acc_factor")
  dtau_dt <- (tau_of(times + eps) - tau_of(times - eps)) / (2 * eps)
  expect_equal(exp(pr$estimate), dtau_dt, tolerance = 1e-3)
})

test_that("acc_factor matches log(translate_time())", {
  skip_if_not_installed("rstpm2")
  fit   <- .fit_gen
  times <- c(0.5, 1, 2, 3)

  # translate_time(cov0 = row 1 observed at t, cov1 = row 2 warped)
  pr  <- predict(fit, newdata = .nd_gen, times = times, type = "acc_factor")
  phi <- translate_time(fit, covariate0 = 0, covariate1 = 1, time = times)
  expect_equal(pr$estimate, log(as.numeric(phi)), tolerance = 1e-4)
})

test_that("acc_factor delta-method variance matches theta-perturbation FD", {
  fit   <- .fit_gen
  times <- c(0.5, 1, 2)
  z     <- qnorm(0.975)
  eps   <- 1e-5

  # f(t) = log phi(t) as a function of the full parameter vector theta
  ffun <- function(theta_vec, t) {
    f2     <- fit
    f2$par <- theta_vec
    H1t    <- post(f2, 0, t, "H")
    tau    <- uniroot(function(u) post(f2, 1, u, "H") - H1t,
                      lower = 1e-8, upper = 60, tol = 1e-12)$root
    log(post(f2, 0, t, "h") / post(f2, 1, tau, "h"))
  }

  pr <- predict(fit, newdata = .nd_gen, times = times, type = "acc_factor")
  # variance implied by the returned (symmetric) CI on the effect scale
  var_pred <- ((pr$estimate - pr$lower) / z)^2

  theta <- fit$par
  p     <- length(theta)
  var_fd <- vapply(times, function(t) {
    g <- numeric(p)
    for (j in seq_len(p)) {
      tp <- tm <- theta
      tp[j] <- theta[j] + eps
      tm[j] <- theta[j] - eps
      g[j] <- (ffun(tp, t) - ffun(tm, t)) / (2 * eps)
    }
    as.numeric(crossprod(g, fit$var %*% g))
  }, numeric(1))

  expect_equal(var_pred, var_fd, tolerance = 1e-3)
})

test_that("time_ratio warps row 2 (TR < 1 when the warped group has higher hazard)", {
  fit   <- .fit_gen
  times <- c(1, 2, 3)
  tr <- predict(fit, newdata = .nd_gen, times = times, type = "time_ratio",
                ci = FALSE)
  # row 2 (X = 1) has the higher hazard, so it reaches row 1's cumulative
  # hazard earlier: tau < t, i.e. TR < 1.
  expect_true(all(tr$estimate < 1))
  expect_true(all(is.finite(tr$estimate)))
})

test_that("acc_factor requires exactly two rows in newdata", {
  fit <- .fit_gen
  nd1 <- data.frame(X = 0); rownames(nd1) <- "only"
  expect_error(
    predict(fit, newdata = nd1, times = c(1, 2), type = "acc_factor"),
    "exactly 2 rows"
  )
})

test_that("time_ratio/acc_factor cap tau at tau_max with NA + warning", {
  fit   <- .fit_gen
  times <- c(0.5, 1, 2)

  for (ty in c("time_ratio", "acc_factor")) {
    # tiny tau_max -> no in-support solution -> all NA + warning
    expect_warning(
      pr <- predict(fit, newdata = .nd_gen, times = times, type = ty,
                    tau_max = 1e-6),
      "tau_max",
      info = ty
    )
    expect_true(all(is.na(pr$estimate)), info = ty)

    # default (in-support) tau_max -> finite estimates for these interior times
    pr2 <- predict(fit, newdata = .nd_gen, times = times, type = ty)
    expect_true(all(is.finite(pr2$estimate)), info = ty)
  }
})

test_that("ci = FALSE returns NA bounds and unchanged point estimates", {
  fit   <- .fit_gen
  times <- c(0.5, 1, 2, 3)

  for (ty in c("hazard", "survival", "cumhaz", "rmst")) {
    full <- predict(fit, newdata = .nd_gen, times = times, type = ty)
    noci <- predict(fit, newdata = .nd_gen, times = times, type = ty, ci = FALSE)
    expect_true(all(is.na(noci$lower)) && all(is.na(noci$upper)), info = ty)
    expect_equal(noci$estimate, full$estimate, tolerance = 1e-8, info = ty)
  }

  for (ty in c("surv_diff", "rmst_diff", "hazard_ratio", "time_ratio",
               "acc_factor")) {
    full <- predict(fit, newdata = .nd_gen, times = times, type = ty)
    noci <- predict(fit, newdata = .nd_gen, times = times, type = ty, ci = FALSE)
    expect_true(all(is.na(noci$lower)) && all(is.na(noci$upper)), info = ty)
    expect_equal(noci$estimate, full$estimate, tolerance = 1e-8, info = ty)
  }
})
