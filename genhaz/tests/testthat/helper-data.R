# Shared fixtures for all tests.
# Called automatically by testthat before test files.

library(survival)

set.seed(20240101)

# Small dataset simulated from scenario 1)
.sim_dat <- sim_szenario(szenario = 1, beta1 = 0.5, beta2 = 0.5, n = 500)

# Shared knots and projection matrix
.n_knots  <- 8
.knots    <- knot_pattern(.sim_dat$time, .sim_dat$event, n_knots = .n_knots)
.Z        <- attr(smf_cpp(log(.sim_dat$time + 1e-10),
                          knots = .knots, intercept = FALSE, derivs = 0L), "Z")
.X_mat    <- matrix(.sim_dat$X, ncol = 1)
.theta0   <- rep(0, 1 + (.n_knots - 1) + 2)   # intercept + splines + 2 betas
