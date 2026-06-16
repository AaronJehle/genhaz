test_that("S_mix is between 0 and 1", {
  t <- seq(0.01, 5, length.out = 50)
  s <- mixWeib("S", t, p = 0.8, lambda1 = 0.1, lambda2 = 0.1,
               gamma1 = 3, gamma2 = 1.6, X = 0, beta1 = 0.5, beta2 = 0.5)
  expect_true(all(s >= 0 & s <= 1))
})

test_that("S_mix is non-increasing over time", {
  t <- seq(0.01, 5, length.out = 100)
  s <- mixWeib("S", t, p = 0.8, lambda1 = 0.1, lambda2 = 0.1,
               gamma1 = 3, gamma2 = 1.6, X = 0, beta1 = 0, beta2 = 0)
  expect_true(all(diff(s) <= 1e-10))
})

test_that("H = -log(S)", {
  t <- seq(0.1, 4, length.out = 30)
  S <- mixWeib("S", t, p = 0.5, lambda1 = 1, lambda2 = 1,
               gamma1 = 1.5, gamma2 = 0.5, X = 1, beta1 = 0.3, beta2 = 0.2)
  H <- mixWeib("H", t, p = 0.5, lambda1 = 1, lambda2 = 1,
               gamma1 = 1.5, gamma2 = 0.5, X = 1, beta1 = 0.3, beta2 = 0.2)
  expect_equal(-log(S), H, tolerance = 1e-10)
})

test_that("h is non-negative", {
  t <- seq(0.01, 5, length.out = 50)
  hv <- mixWeib("h", t, p = 0.8, lambda1 = 0.1, lambda2 = 0.1,
                gamma1 = 3, gamma2 = 1.6, X = 1, beta1 = 0.5, beta2 = 0.5)
  expect_true(all(hv >= 0))
})

test_that("mixWeibSz dispatches correctly for all 3 scenarios", {
  t <- 1:3
  for (sz in 1:3) {
    expect_length(mixWeibSz(sz, "h", t), 3)
    expect_length(mixWeibSz(sz, "S", t), 3)
    expect_length(mixWeibSz(sz, "H", t), 3)
  }
})

test_that("mixWeibSz stops on invalid scenario", {
  expect_error(mixWeibSz(99, "h", 1:3))
})

test_that("sim_mix_weib_gh returns data frame with correct columns", {
  set.seed(7)
  dat <- sim_mix_weib_gh(n = 50, p = 0.8, lambda1 = 0.1, lambda2 = 0.1,
                          gamma1 = 3, gamma2 = 1.6, beta1 = 0.5, beta2 = 0.5)
  expect_s3_class(dat, "data.frame")
  expect_true(all(c("time", "X", "event", "T_true") %in% names(dat)))
  expect_equal(nrow(dat), 50)
})

test_that("sim_mix_weib_gh event is 0 or 1", {
  set.seed(8)
  dat <- sim_mix_weib_gh(n = 100, p = 0.5, lambda1 = 1, lambda2 = 1,
                          gamma1 = 1.5, gamma2 = 0.5, beta1 = 0, beta2 = 0)
  expect_true(all(dat$event %in% c(0L, 1L)))
})

test_that("sim_szenario wrappers work for all 3 scenarios", {
  for (sz in 1:3) {
    set.seed(sz)
    dat <- sim_szenario(sz, beta1 = 0.3, beta2 = 0.3, n = 50)
    expect_s3_class(dat, "data.frame")
    expect_equal(nrow(dat), 50)
  }
})
