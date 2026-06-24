test_that("gauss_legendre nodes lie in (-1, 1) and weights sum to 2", {
  gl <- genhaz:::gauss_legendre(25L)

  expect_true(all(gl$x > -1 & gl$x < 1))
  expect_equal(sum(gl$w), 2, tolerance = 1e-12)
})

test_that("gauss_legendre integrates a polynomial exactly", {
  # GL with n points is exact for polynomials of degree <= 2n-1.
  # ∫₋₁¹ x^4 dx = 2/5; n=25 is more than sufficient.
  gl <- genhaz:::gauss_legendre(25L)
  expect_equal(sum(gl$w * gl$x^4), 2 / 5, tolerance = 1e-12)
})

test_that("gauss_legendre result is cached after first call", {
  genhaz:::gauss_legendre(7L)
  expect_true(exists("7", envir = genhaz:::.gl_cache, inherits = FALSE))

  gl1 <- genhaz:::gauss_legendre(7L)
  gl2 <- genhaz:::gauss_legendre(7L)
  expect_identical(gl1, gl2)
})
