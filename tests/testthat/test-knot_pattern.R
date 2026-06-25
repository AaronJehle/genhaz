test_that("knot_pattern returns a vector of the right length", {
  t  <- c(0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  ev <- c(1,   1, 0, 1, 1, 0, 1, 1, 0, 1, 1)
  k  <- knot_pattern(t, ev, n_knots = 5)
  expect_length(k, 5)
})

test_that("knot_pattern values are on the log-time scale", {
  t  <- exp(seq(0, 2, length.out = 50))
  ev <- rep(1, 50)
  k  <- knot_pattern(t, ev, n_knots = 4)
  expect_true(all(k >= 0))   # log(t) > 0 since t > 1
})

test_that("limit_to_95 = FALSE uses the full range of event times", {
  set.seed(7)
  t  <- rexp(100)
  ev <- rep(1, 100)
  k95  <- knot_pattern(t, ev, n_knots = 4, limit_to_95 = TRUE)
  k100 <- knot_pattern(t, ev, n_knots = 4, limit_to_95 = FALSE)
  # Upper boundary knot should be at or above the 95 % one
  expect_true(max(k100) >= max(k95))
})

test_that("knot_pattern ignores censored observations", {
  t  <- c(1, 2, 3, 100)   # 100 is a huge censored time
  ev <- c(1, 1, 1, 0)
  k  <- knot_pattern(t, ev, n_knots = 3)
  # all knots must come from log(1), log(2), log(3)
  expect_true(max(k) <= log(3) + 1e-9)
})
