#' Mixture Weibull survival, cumulative hazard, and hazard functions
#'
#' @description
#' A two-component mixture Weibull model with a general hazard (GH)
#' parametrisation. These functions are primarily intended for generating
#' synthetic survival data for simulation studies.
#'
#' The survival function is
#' \deqn{S(t \mid X) = \bigl[p \exp(-\lambda_1 t^{\gamma_1} e^{X \beta_1 \gamma_1}) +
#'   (1-p)\exp(-\lambda_2 t^{\gamma_2} e^{X \beta_1 \gamma_2})\bigr]^{\exp(\beta_2 X)}}
#'
#' @name mixture_weibull
NULL

# Internal: mixture Weibull survival function
S_mix <- function(t, p, lambda1, lambda2, gamma1, gamma2, X, beta1, beta2) {
  (p * exp(-lambda1 * t^gamma1 * exp(X * beta1 * gamma1)) +
     (1 - p) * exp(-lambda2 * t^gamma2 * exp(X * beta1 * gamma2)))^exp(beta2 * X)
}

# Internal: mixture Weibull cumulative hazard
H_mix <- function(t, p, lambda1, lambda2, gamma1, gamma2, X, beta1, beta2) {
  -log(S_mix(t, p, lambda1, lambda2, gamma1, gamma2, X, beta1, beta2))
}

# Internal: mixture Weibull AFT hazard (beta2 = 0 arm)
h_mix_aft <- function(t, p, lambda1, lambda2, gamma1, gamma2, X, beta) {
  S_t <- S_mix(t, p, lambda1, lambda2, gamma1, gamma2, X, beta, 0)
  (p * lambda1 * gamma1 * t^(gamma1 - 1) * exp(X * beta * gamma1) *
     exp(-lambda1 * t^gamma1 * exp(X * beta * gamma1)) +
     (1 - p) * lambda2 * gamma2 * t^(gamma2 - 1) * exp(X * beta * gamma2) *
     exp(-lambda2 * t^gamma2 * exp(X * beta * gamma2))) / S_t
}

# Internal: full mixture Weibull hazard
h_mix <- function(t, p, lambda1, lambda2, gamma1, gamma2, X, beta1, beta2) {
  h_mix_aft(t, p, lambda1, lambda2, gamma1, gamma2, X, beta1) * exp(beta2 * X)
}

#' Mixture Weibull GH model quantities
#'
#' Evaluates the survival function (`"S"`), cumulative hazard (`"H"`), or
#' hazard (`"h"`) of a two-component mixture Weibull GH model.
#'
#' @param res Character: `"S"`, `"H"`, or `"h"`.
#' @param t Numeric vector of time points.
#' @param p Mixture proportion (weight of the first component), in (0, 1).
#' @param lambda1,lambda2 Scale parameters of the two Weibull components.
#' @param gamma1,gamma2 Shape parameters of the two Weibull components.
#' @param X Numeric covariate value (scalar or vector).
#' @param beta1 AFT-type covariate coefficient.
#' @param beta2 PH-type covariate coefficient.
#'
#' @return Numeric vector of the requested quantity.
#' @export
#'
#' @examples
#' t <- seq(0.01, 5, length.out = 100)
#' mixWeib("h", t, p = 0.8, lambda1 = 0.1, lambda2 = 0.1,
#'         gamma1 = 3, gamma2 = 1.6, X = 0, beta1 = 0.5, beta2 = 0.5)
mixWeib <- function(res = "S", t, p, lambda1, lambda2, gamma1, gamma2,
                    X, beta1, beta2) {
  switch(res,
    S = S_mix(t, p, lambda1, lambda2, gamma1, gamma2, X, beta1, beta2),
    H = H_mix(t, p, lambda1, lambda2, gamma1, gamma2, X, beta1, beta2),
    h = h_mix(t, p, lambda1, lambda2, gamma1, gamma2, X, beta1, beta2),
    stop("'res' must be one of \"S\", \"H\", \"h\".")
  )
}

#' Mixture Weibull model for a named simulation scenario
#'
#' Convenience wrapper around [mixWeib()] for the three simulation scenarios
#' used in the accompanying study:
#' \describe{
#'   \item{Scenario 1}{Bathtub-like hazard (p=0.8, lambda1=lambda2=0.1,
#'     gamma1=3, gamma2=1.6).}
#'   \item{Scenario 2}{Unimodal hazard (p=0.5, lambda1=lambda2=1,
#'     gamma1=1.5, gamma2=0.5).}
#'   \item{Scenario 3}{Bimodal hazard (p=0.26, lambda1=0.02, lambda2=0.5,
#'     gamma1=3, gamma2=0.7).}
#' }
#'
#' @param scenario Integer (1, 2, or 3) selecting the scenario.
#' @param res Character: `"S"`, `"H"`, or `"h"`.
#' @param t Numeric vector of time points.
#' @param X Covariate value. Default 0.
#' @param beta1 AFT coefficient. Default 0.
#' @param beta2 PH coefficient. Default 0.
#'
#' @return Numeric vector.
#' @export
#'
#' @examples
#' t <- seq(0.01, 5, length.out = 100)
#' h1 <- mixWeibSc(1, "h", t, X = 1, beta1 = 0.5, beta2 = 0.5)
mixWeibSc <- function(scenario, res = "S", t, X = 0, beta1 = 0, beta2 = 0) {
  params <- switch(as.character(scenario),
    "1" = list(p = 0.8, lambda1 = 0.1, lambda2 = 0.1,  gamma1 = 3,    gamma2 = 1.6),
    "2" = list(p = 0.5, lambda1 = 1,   lambda2 = 1,    gamma1 = 1.5,  gamma2 = 0.5),
    "3" = list(p = 0.26,lambda1 = 0.02,lambda2 = 0.5,  gamma1 = 3,    gamma2 = 0.7),
    stop("'scenario' must be 1, 2, or 3.")
  )
  mixWeib(res, t, params$p, params$lambda1, params$lambda2,
          params$gamma1, params$gamma2, X, beta1, beta2)
}

#' Simulate data from a mixture Weibull GH model
#'
#' Simulates right-censored survival data from a two-component mixture Weibull
#' GH model with a binary covariate X ~ Bernoulli(0.5). Censoring times are
#' drawn from Uniform(0, `tmax`).
#'
#' @param n Sample size.
#' @param p Mixture proportion.
#' @param lambda1,lambda2 Weibull scale parameters.
#' @param gamma1,gamma2 Weibull shape parameters.
#' @param beta1 AFT coefficient for X.
#' @param beta2 PH coefficient for X.
#' @param tmax Administrative censoring time (uniform censoring upper bound).
#'
#' @return Data frame with columns `time`, `X`, `event`, and `T_true`.
#' @export
#'
#' @examples
#' set.seed(42)
#' dat <- sim_mix_weib_gh(n = 200, p = 0.8, lambda1 = 0.1, lambda2 = 0.1,
#'                        gamma1 = 3, gamma2 = 1.6, beta1 = 0.5, beta2 = 0.5)
#' head(dat)
sim_mix_weib_gh <- function(n, p, lambda1, lambda2, gamma1, gamma2,
                             beta1, beta2, tmax = 10) {
  X      <- rbinom(n, 1, 0.5)
  U      <- runif(n)
  T_true <- numeric(n)
  for (i in seq_len(n)) {
    target    <- function(t) S_mix(t, p, lambda1, lambda2, gamma1, gamma2,
                                   X[i], beta1, beta2) - U[i]
    T_true[i] <- uniroot(target, interval = c(0, 100),
                          extendInt = "yes")$root
  }
  t_cens <- runif(n, min = 0, max = tmax)
  t_obs  <- pmin(T_true, t_cens)
  event  <- as.integer(T_true <= t_obs)
  data.frame(time = t_obs, X = X, event = event, T_true = T_true)
}

#' Simulate data for a named scenario
#'
#' Convenience wrapper around [sim_mix_weib_gh()] using the three pre-defined
#' simulation scenarios (see [mixWeibSc()]).
#'
#' @param scenario Integer (1, 2, or 3).
#' @param beta1 AFT coefficient.
#' @param beta2 PH coefficient.
#' @param n Sample size.
#' @param tmax Administrative censoring time.
#'
#' @return Data frame with columns `time`, `X`, `event`, `T_true`.
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- sim_scenario(1, beta1 = 0.5, beta2 = 0.5, n = 500)
#' table(dat$event)
sim_scenario <- function(scenario, beta1, beta2, n = 1000, tmax = 10) {
  params <- switch(as.character(scenario),
    "1" = list(p = 0.8, lambda1 = 0.1, lambda2 = 0.1,  gamma1 = 3,   gamma2 = 1.6),
    "2" = list(p = 0.5, lambda1 = 1,   lambda2 = 1,    gamma1 = 1.5, gamma2 = 0.5),
    "3" = list(p = 0.26,lambda1 = 0.02,lambda2 = 0.5,  gamma1 = 3,   gamma2 = 0.7),
    stop("'scenario' must be 1, 2, or 3.")
  )
  sim_mix_weib_gh(n = n, p = params$p, lambda1 = params$lambda1,
                  lambda2 = params$lambda2, gamma1 = params$gamma1,
                  gamma2 = params$gamma2, beta1 = beta1, beta2 = beta2,
                  tmax = tmax)
}
