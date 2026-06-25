test_that("genhaz_work returns a scalar negll at theta0", {
  nll <- genhaz_work(.theta0, .sim_dat$time, .X_mat, .knots, .Z,
                     event = .sim_dat$event, lambda = 1, res = "negll",
                     model_type = "GH")
  expect_length(nll, 1)
  expect_true(is.finite(nll))
})

test_that("gradient_negll has correct length", {
  g <- genhaz_work(.theta0, .sim_dat$time, .X_mat, .knots, .Z,
                   event = .sim_dat$event, lambda = 1,
                   res = "gradient_negll", model_type = "GH")
  expect_length(g, length(.theta0))
  expect_true(all(is.finite(g)))
})

test_that("gradient_negll passes finite-difference check", {
  n_chk  <- 20
  t_chk  <- .sim_dat$time[seq_len(n_chk)]
  X_chk  <- .X_mat[seq_len(n_chk), , drop = FALSE]
  ev_chk <- .sim_dat$event[seq_len(n_chk)]
  eps    <- 1e-6

  g_ana <- genhaz_work(.theta0, t_chk, X_chk, .knots, .Z,
                       event = ev_chk, lambda = 1,
                       res = "gradient_negll", model_type = "GH")
  g_fd  <- numeric(length(.theta0))
  for (j in seq_along(.theta0)) {
    tp <- tm <- .theta0
    tp[j] <- .theta0[j] + eps
    tm[j] <- .theta0[j] - eps
    nll_p <- genhaz_work(tp, t_chk, X_chk, .knots, .Z, event = ev_chk,
                         lambda = 1, res = "negll", model_type = "GH")
    nll_m <- genhaz_work(tm, t_chk, X_chk, .knots, .Z, event = ev_chk,
                         lambda = 1, res = "negll", model_type = "GH")
    g_fd[j] <- (nll_p - nll_m) / (2 * eps)
  }
  expect_lt(max(abs(g_ana - g_fd)), 1e-4)
})

test_that("gradient_log_h passes finite-difference check", {
  n_chk <- 5
  t_chk <- .sim_dat$time[seq_len(n_chk)]
  X_chk <- .X_mat[seq_len(n_chk), , drop = FALSE]
  eps   <- 1e-6

  dlogh_ana <- genhaz_work(.theta0, t_chk, X_chk, .knots, .Z,
                           res = "gradient_log_h", model_type = "GH")
  dlogh_fd  <- matrix(0, n_chk, length(.theta0))
  for (j in seq_along(.theta0)) {
    tp <- tm <- .theta0
    tp[j] <- .theta0[j] + eps
    tm[j] <- .theta0[j] - eps
    logh_p <- genhaz_work(tp, t_chk, X_chk, .knots, .Z,
                          res = "log_h", model_type = "GH")
    logh_m <- genhaz_work(tm, t_chk, X_chk, .knots, .Z,
                          res = "log_h", model_type = "GH")
    dlogh_fd[, j] <- (logh_p - logh_m) / (2 * eps)
  }
  expect_lt(max(abs(dlogh_ana - dlogh_fd)), 1e-4)
})


test_that("gradient_H passes finite-difference check", {
  n_chk <- 5
  t_chk <- .sim_dat$time[seq_len(n_chk)]
  X_chk <- .X_mat[seq_len(n_chk), , drop = FALSE]
  eps   <- 1e-6

  dH_ana <- genhaz_work(.theta0, t_chk, X_chk, .knots, .Z,
                        res = "gradient_H", model_type = "GH")
  dH_fd  <- matrix(0, n_chk, length(.theta0))
  for (j in seq_along(.theta0)) {
    tp <- tm <- .theta0
    tp[j] <- .theta0[j] + eps
    tm[j] <- .theta0[j] - eps
    H_p <- genhaz_work(tp, t_chk, X_chk, .knots, .Z, res = "H", model_type = "GH")
    H_m <- genhaz_work(tm, t_chk, X_chk, .knots, .Z, res = "H", model_type = "GH")
    dH_fd[, j] <- (H_p - H_m) / (2 * eps)
  }
  expect_lt(max(abs(dH_ana - dH_fd)), 1e-4)
})

test_that("hessian_negll passes finite-difference check", {
  n_chk  <- 5
  t_chk  <- .sim_dat$time[seq_len(n_chk)]
  X_chk  <- .X_mat[seq_len(n_chk), , drop = FALSE]
  ev_chk <- .sim_dat$event[seq_len(n_chk)]
  eps    <- 1e-5

  H_ana <- genhaz_work(.theta0, t_chk, X_chk, .knots, .Z,
                       event = ev_chk, lambda = 1,
                       res = "hessian_negll", model_type = "GH")
  H_fd  <- matrix(0, length(.theta0), length(.theta0))
  for (j in seq_along(.theta0)) {
    tp <- tm <- .theta0
    tp[j] <- .theta0[j] + eps
    tm[j] <- .theta0[j] - eps
    gp <- genhaz_work(tp, t_chk, X_chk, .knots, .Z, event = ev_chk,
                      lambda = 1, res = "gradient_negll", model_type = "GH")
    gm <- genhaz_work(tm, t_chk, X_chk, .knots, .Z, event = ev_chk,
                      lambda = 1, res = "gradient_negll", model_type = "GH")
    H_fd[j, ] <- (gp - gm) / (2 * eps)
  }
  expect_lt(max(abs(H_ana - H_fd)), 1e-3)
})

test_that("dLCV_full passes finite difference check", {
  t_chk  <- .sim_dat$time
  X_chk  <- .X_mat[, , drop = FALSE]
  ev_chk <- .sim_dat$event
  eps    <- 1e-5

  dLCV_ana <- genhaz_work(.theta0, t_chk, X_chk, .knots, .Z,
                          event = ev_chk, lambda = exp(10),
                          res = "dLCV_full", model_type = "GH")
  
  LCV_p <- genhaz_work(.theta0, t_chk, X_chk, .knots, .Z,
                       event = ev_chk, lambda = exp(10 + eps),
                       res = "LCV", model_type = "GH")
  LCV_m <- genhaz_work(.theta0, t_chk, X_chk, .knots, .Z,
                       event = ev_chk, lambda = exp(10 - eps),
                       res = "LCV", model_type = "GH")
  
  dLCV_fd <- (LCV_p - LCV_m) / (2 * eps)
  # print(c(dLCV_ana = dLCV_ana, dLCV_fd = dLCV_fd))
  expect_lt(max(abs(dLCV_ana - dLCV_fd)), 1e-4)
})

test_that("log_h and h are self-consistent (exp(log_h) == h)", {
  lh <- genhaz_work(.theta0, .sim_dat$time[1:10], .X_mat[1:10, , drop = FALSE],
                    .knots, .Z, res = "log_h", model_type = "GH")
  hv <- genhaz_work(.theta0, .sim_dat$time[1:10], .X_mat[1:10, , drop = FALSE],
                    .knots, .Z, res = "h", model_type = "GH")
  expect_equal(exp(lh), hv, tolerance = 1e-10)
})

test_that("H is non-negative and monotone in time", {
  t_inc <- sort(runif(20, 0.1, 8))
  Hv    <- genhaz_work(.theta0, t_inc, .X_mat[seq_len(20), , drop = FALSE],
                       .knots, .Z, res = "H", model_type = "GH")
  expect_true(all(Hv >= 0))
  expect_true(all(diff(Hv) >= -1e-10))   # monotone (allow tiny numerical noise)
})

test_that("PH model forces beta1 to zero after fitting", {
  fit <- genhaz_work(.theta0, .sim_dat$time, .X_mat, .knots, .Z,
                     event = .sim_dat$event, lambda = 5,
                     res = "fit", model_type = "PH")
  idx_b1 <- (.n_knots) + 1   # position of beta1_1
  expect_equal(unname(fit$par[idx_b1]), 0)
})

test_that("AFT model forces beta2 to zero after fitting", {
  fit <- genhaz_work(.theta0, .sim_dat$time, .X_mat, .knots, .Z,
                     event = .sim_dat$event, lambda = 5,
                     res = "fit", model_type = "AFT")
  idx_b2 <- (.n_knots) + 2   # position of beta2_1
  expect_equal(unname(fit$par[idx_b2]), 0)
})
