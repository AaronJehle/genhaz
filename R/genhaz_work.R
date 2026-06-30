# Package-level cache for Gauss-Legendre nodes/weights, keyed by n.
# Avoids recomputing inside the many genhaz_work() calls during LCV profiling.
.gl_cache <- new.env(parent = emptyenv())

# Gauss-Legendre quadrature via the Golub-Welsch algorithm.
# Returns list(x = nodes, w = weights) on [-1, 1], same format as
# pracma::gaussLegendre(n, -1, 1). Results are cached by n for the session.
#' @noRd
gauss_legendre <- function(n) {
  key <- as.character(n)
  if (!exists(key, envir = .gl_cache, inherits = FALSE)) {
    i   <- seq_len(n - 1L)
    b   <- i / sqrt(4 * i^2 - 1)
    J   <- matrix(0, n, n)
    J[cbind(i,      i + 1L)] <- b
    J[cbind(i + 1L, i     )] <- b
    eig <- eigen(J, symmetric = TRUE)
    assign(key,
           list(x = eig$values, w = 2 * eig$vectors[1L, ]^2),
           envir = .gl_cache)
  }
  get(key, envir = .gl_cache, inherits = FALSE)
}

#' General hazard model workhorse
#'
#' Calculates various quantities from a general hazard (GH) model depending on
#' the `res` argument: negative log-likelihood and its gradient/Hessian, the
#' cumulative hazard H and its derivatives, effective degrees of freedom, the
#' LCV criterion, or a fully fitted model object. Should not be called directly
#' by users; use [fit_genhaz()] and [post()] instead.
#'
#' @param theta Numeric vector of parameters: intercept, spline coefficients,
#'   then beta1 (AFT-style effects) and beta2 (PH-style effects).
#' @param time Numeric vector of observed (or evaluation) times.
#' @param X Numeric matrix of covariates (rows = observations).
#' @param knots Numeric vector of knot positions on the log-time scale,
#'   including boundary knots.
#' @param Z Projection matrix for the spline basis (from [smf_cpp()]).
#' @param event Integer event indicator vector (1 = event, 0 = censored).
#'   Required for fitting.
#' @param lambda Non-negative smoothing parameter for the spline penalty.
#' @param res Character string selecting the quantity to return. One of
#'   `"log_h"`, `"gradient_log_h"`, `"hessian_log_h"`, `"dlogh_dt"`
#'   (derivative of the log-hazard w.r.t. time), `"dh_dt"` (derivative of the
#'   hazard w.r.t. time), `"h"`,
#'   `"gradient_h"`, `"hessian_h"`, `"H"`, `"gradient_H"`, `"hessian_H"`,
#'   `"negll"`, `"negll_upen"`, `"scores"`, `"gradient_negll"`,
#'   `"hessian_negll"`, `"hessian_negll_upen"`, `"basis"`, `"edf"`, `"LCV"`,
#'   `"dLCV"`, `"dLCV_full"`, `"fit"`, or `"fit_LCV"`. The last two fit a
#'   model: `"fit"` at a fixed `lambda`, `"fit_LCV"` after profiling `lambda`
#'   via LCV. `"dLCV"` returns the first-order analytical derivative of LCV
#'   w.r.t. `log(lambda)` (evaluated at the MLE for the given `lambda`).
#'   `"dLCV_full"` returns the same derivative including the third-derivative
#'   correction term from the unpenalised Hessian's dependence on `lambda`
#'   (supports `"rc"` and `"lt_rc"`; falls back to `"dLCV"` for `"ic"`).
#' @param model_type Character vector specifying the model type for each
#'   covariate: `"GH"` (general hazard), `"PH"` (proportional hazards,
#'   beta1 = 0), `"AFT"` (accelerated failure time, beta2 = 0), or
#'   `"AH"` (additive hazards, beta1 = -beta2). If length 1, recycled to all
#'   covariates.
#' @param time2 Numeric vector for the second time variable. Used as the
#'   truncation time for `cens_type = "lt_rc"`, or as the right interval
#'   endpoint for `cens_type = "ic"`.
#' @param cens_type Censoring type: `"rc"` (right-censoring, default),
#'   `"lt_rc"` (left-truncation + right-censoring), or `"ic"` (interval
#'   censoring).
#' @param gaussified Logical; use Gauss-Legendre quadrature for integrating
#'   the cumulative hazard? Default `TRUE`. Setting to `FALSE` requires
#'   **rstpm2** to be installed.
#' @param nquad Number of quadrature points when `gaussified = TRUE`.
#'   Default 25.
#' @param control List of optimisation control parameters passed to
#'   [nlminb()]: `trace`, `x.tol`, `rel.tol`, `max_iter`.
#' @param lcv_method Character string. When `res = "fit_LCV"`, selects the
#'   strategy for optimising `log(lambda)`. One of:
#'   \describe{
#'     \item{`"full"` (default)}{Root-find on the full LCV gradient including
#'       the third-derivative correction (`dLCV_full`). Most accurate.}
#'     \item{`"approx"`}{Root-find on the first-order LCV gradient (`dLCV`),
#'       which ignores the dependence of the Hessian on `lambda`. Faster than
#'       `"full"`.}
#'     \item{`"optimize"`}{Directly minimise LCV without gradient information
#'       via `optimize()`. Equivalent to the original behaviour before analytic
#'       gradients were added.}
#'   }
#'   Ignored for all other `res` values and falls back to `lcv_gradient` for
#'   `cens_type = "ic"`.
#' @param interval Length-2 numeric vector giving the search interval for
#'   `log(lambda)` during LCV optimisation. Default `c(0, 13)`.
#' @param tol_LCV Convergence tolerance for LCV optimisation. Default `1e-3`.
#' @param ... Additional arguments (currently unused).
#'
#' @return Depends on `res`. When `res = "fit"` or `"fit_LCV"`, returns a list
#'   with components `par`, `se`, `z`, `p_values`, `AIC`, `edf`, `lambda`,
#'   `var`, and further fitting details.
#' @export
genhaz_work <- function(theta, time, X, knots, Z, event = NULL, lambda = 1,
                        res = c("log_h", "gradient_log_h", "hessian_log_h",
                                "dlogh_dt", "dh_dt",
                                "h", "gradient_h", "hessian_h",
                                "H", "gradient_H", "hessian_H",
                                "negll", "negll_upen", "scores",
                                "gradient_negll", "hessian_negll",
                                "hessian_negll_upen",
                                "basis", "edf", "LCV", "dLCV", "dLCV_full",
                                "fit", "fit_LCV"),
                        model_type = "GH", time2 = NULL, cens_type = "rc",
                        gaussified = TRUE, nquad = 25,
                        control = list(trace = 0, x.tol = 1e-12,
                                       rel.tol = 1e-12, max_iter = 100),
                        lcv_method = c("full", "approx", "optimize"),
                        interval = c(0, 13), tol_LCV = 1e-3,
                        ...) {

  res        <- match.arg(res)
  lcv_method <- match.arg(lcv_method)
  # Parameter vector layout: theta = (intercept, beta0[2..nb0], beta1[1..nb1], beta2[1..nb2])
  # nb0 = number of spline basis coefs including the intercept; nb1 = nb2 = number of covariates.
  # Total length of theta: nb0 + 2*nb1.
  nb0 <- length(knots)
  nb1 <- (length(theta) - nb0) / 2
  nb2 <- nb1

  if (!(cens_type %in% c("rc", "ic", "lt_rc"))) {
    stop("Invalid censoring type, supported are rc, ic, lt_rc")
  }

  if (is.null(dim(X))) {
    X <- matrix(X, ncol = 1)
  }

  if (nrow(X) != length(time)) {
    if (nrow(X) == 1) {
      X <- X[rep(1, length(time)), , drop = FALSE]
    } else {
      stop("Number of rows in X and length of time must be the same.")
    }
  }

  allowed_models <- c("GH", "PH", "AFT", "AH")
  if (length(model_type) == 1) {
    model_type <- rep(model_type, nb1)
  } else if (length(model_type) != nb1) {
    stop("Length of 'model_type' must be 1 or equal to the number of covariates.")
  }
  if (!all(model_type %in% allowed_models)) {
    stop("All elements of 'model_type' must be one of: 'GH', 'PH', 'AFT', 'AH'")
  }

  init   <- theta

  # Sub-model constraints are enforced by zeroing out irrelevant covariate columns.
  # PH  (beta1 = 0): X1_filtered zeros PH columns; time-shift eta1 = 0 for those.
  # AFT (beta2 = 0): X2_filtered zeros AFT columns.
  # AH  (beta1 = -beta2): X3_filtered carries the beta1 contribution to (beta1+beta2)^T x;
  #   the equality constraint is imposed as a box bound in "fit" and post-processed there.
  # eta1 = beta1^T x  (log-time shift); eta2 = (beta1+beta2)^T x  (log-hazard shift).
  is_gh  <- ifelse(model_type == "GH",  1, 0)
  is_ph  <- ifelse(model_type == "PH",  1, 0)
  is_aft <- ifelse(model_type == "AFT", 1, 0)
  is_ah  <- ifelse(model_type == "AH",  1, 0)

  is_not_gh  <- ifelse(model_type == "GH",  0, 1)
  is_not_ph  <- ifelse(model_type == "PH",  0, 1)
  is_not_aft <- ifelse(model_type == "AFT", 0, 1)
  is_not_ah  <- ifelse(model_type == "AH",  0, 1)

  X1_filtered <- X * matrix(is_not_ph,        nrow = nrow(X), ncol = nb1, byrow = TRUE)
  X2_filtered <- X * matrix(is_gh + is_ph,    nrow = nrow(X), ncol = nb2, byrow = TRUE)
  X3_filtered <- X1_filtered * matrix(is_not_ah, nrow = nrow(X), ncol = nb1, byrow = TRUE)

  filter1 <- is_gh + is_aft
  filter2 <- is_gh + is_ph

  eta1 <- function(theta) {
    beta1 <- theta[nb0 + 1:nb1]
    drop(X %*% (beta1 * is_not_ph))
  }
  eta2 <- function(theta) {
    beta1 <- theta[nb0 + 1:nb1]
    beta2 <- theta[nb0 + nb1 + 1:nb2]
    drop(X %*% (beta1 * filter1 + beta2 * filter2))
  }

  # Gauss-Legendre nodes on [-1,1] mapped to [0,t]: u = t*(x+1)/2, weight factor t/2.
  # Approximation: integral_0^t h(u) du ~= sum_i (w_i * t/2) * h(t * (x_i+1)/2).
  # Pre-replicated covariate matrices avoid repeated row-indexing inside the quadrature loops
  # for gradient_H and hessian_H_sum (used only when n_arg == n_obs).
  cc <- gauss_legendre(nquad)
  quad_scale  <- (cc$x + 1) / 2
  quad_w_half <- cc$w / 2
  n_obs <- length(time)

  .rep_quad <- rep(seq_len(n_obs), nquad)
  X1f_quad  <- X1_filtered[.rep_quad, , drop = FALSE]
  X2f_quad  <- X2_filtered[.rep_quad, , drop = FALSE]
  X3f_quad  <- X3_filtered[.rep_quad, , drop = FALSE]

  # Integrated second-derivative penalty on the spline coefs beta0:
  #   S_{ij} = integral B''_i(t) B''_j(t) dt  (i,j in 1..nb0, zero for beta1/beta2).
  # Normalised so tr(S) = nb0 - 1, matching the survPen convention.
  smf_base <- smf_cpp(log(time + 1e-10), knots = knots, Z = Z, intercept = FALSE)
  S_penal_wo_intercept <- attr(smf_base, "pen")
  S_penal_wo_intercept <- S_penal_wo_intercept * (nb0 - 1) /
    sum(diag(S_penal_wo_intercept))

  S_penal <- matrix(0,
                    nrow = nrow(S_penal_wo_intercept) + 1,
                    ncol = ncol(S_penal_wo_intercept) + 1)
  S_penal[-1, -1] <- S_penal_wo_intercept

  # Restricted-cubic-spline basis (no intercept), projected by Z.
  # For derivs > 0, smf_cpp() internally computes crs_cpp(derivs = 0) only to read
  # its knots attribute, then discards that basis and returns
  # crs_cpp(derivs, intercept = TRUE) %*% Z. We skip that wasted 0th-derivative
  # pass by calling crs_cpp() directly (it is exported), which is the single
  # biggest cost in the GH/AFT fit -- the derivative bases are evaluated at the
  # n * nquad cumulative-hazard quadrature points. The derivs = 0 path keeps
  # using smf_cpp() (already optimal). Results are identical up to BLAS rounding.
  B <- function(time_arg, derivs = 0L) {
    if (derivs == 0L) {
      smf_cpp(time_arg, knots = knots, Z = Z, intercept = FALSE)
    } else {
      unclass(crs_cpp(time_arg, knots = knots, derivs = derivs,
                      intercept = TRUE)) %*% Z
    }
  }

  # --- Basis caching ---------------------------------------------------------
  # The spline basis B(log(t) + eta1) depends on theta ONLY through eta1, i.e.
  # only through beta1. During a penalised fit, nlminb evaluates the objective,
  # gradient and Hessian at the same theta (hence the same beta1), and for
  # PH / baseline models beta1 is box-constrained so eta1 never changes at all.
  # We therefore memoise the basis matrices keyed on beta1, recomputing only
  # when beta1 changes. The fixed (theta-independent) log-time grids and the
  # quadrature weight vectors are built once here.
  #
  # has_time_effect is FALSE when no covariate shifts the log-time scale (pure
  # PH / baseline). In that case eta1 == 0 and every first/second-derivative
  # basis enters the gradient and Hessian only multiplied by X1_filtered == 0,
  # so those basis evaluations can be skipped entirely.
  has_time_effect <- any(is_not_ph == 1)

  lt_obs_t1  <- log(time + 1e-10)
  lt_quad_t1 <- log(c(outer(time, quad_scale)) + 1e-10)
  w_quad_t1  <- c(outer(time, quad_w_half))
  if (!is.null(time2)) {
    lt_quad_t2 <- log(c(outer(time2, quad_scale)) + 1e-10)
    w_quad_t2  <- c(outer(time2, quad_w_half))
  }

  .bcache <- new.env(parent = emptyenv())
  .bcache$beta1 <- NULL
  # Cached basis on a fixed grid. kind = "obs" (basis at time_arg itself) or
  # "quad" (time_arg expanded over the GL nodes, observation index fastest).
  # Results are keyed on (kind, grid, derivs) and dropped whenever beta1 moves.
  # Grids that are neither `time` nor `time2` (not produced by the fit/predict
  # paths) fall back to a direct, uncached evaluation.
  cbasis <- function(theta, time_arg, derivs, kind) {
    beta1 <- theta[nb0 + seq_len(nb1)]
    if (!identical(beta1, .bcache$beta1)) {
      .bcache$beta1 <- beta1
      .bcache$e1    <- eta1(theta)
      .bcache$store <- new.env(parent = emptyenv())
    }
    role <- if (identical(time_arg, time)) {
      "t1"
    } else if (!is.null(time2) && identical(time_arg, time2)) {
      "t2"
    } else {
      e1  <- eta1(theta)
      arg <- if (kind == "obs") log(time_arg + 1e-10) + e1
             else log(c(outer(time_arg, quad_scale)) + 1e-10) + rep(e1, nquad)
      return(B(arg, derivs = derivs))
    }
    key <- paste0(kind, role, derivs)
    val <- .bcache$store[[key]]
    if (!is.null(val)) return(val)
    e1  <- .bcache$e1
    arg <- if (kind == "obs") {
      (if (role == "t1") lt_obs_t1 else log(time2 + 1e-10)) + e1
    } else {
      (if (role == "t1") lt_quad_t1 else lt_quad_t2) + rep(e1, nquad)
    }
    val <- B(arg, derivs = derivs)
    .bcache$store[[key]] <- val
    val
  }

  # Log-time GH model: log h(t|x) = B(log(t) + eta1)^T beta0 + intercept + eta2
  # B(.) is the restricted cubic spline basis evaluated on the shifted log-time scale;
  # theta[1] is the intercept, theta[2:nb0] are the spline coefs.
  log_h <- function(theta, time_arg) {
    beta0 <- theta[2:nb0]
    drop(cbasis(theta, time_arg, 0L, "obs") %*% beta0 + theta[1] + eta2(theta))
  }
  h_fn <- function(theta, time_arg) exp(log_h(theta, time_arg))

  # Time-derivative of the log-hazard (and hazard) at fixed theta and X.
  # log h(t|x) = B(log(t) + eta1)^T beta0 + intercept + eta2, so
  #   d log h / dt = (B'(log(t) + eta1)^T beta0) * d/dt[log(t + eps)]
  #               = (B'(log(t) + eta1)^T beta0) / (t + eps).
  # d h / dt = h(t|x) * d log h / dt. Needed for the delta-method CI of the
  # time-varying acceleration factor (see R/generics.R, type = "acc_factor").
  dlogh_dt <- function(theta, time_arg) {
    beta0   <- theta[2:nb0]
    B_prime <- cbasis(theta, time_arg, 1L, "obs")
    drop(B_prime %*% beta0) / (time_arg + 1e-10)
  }
  dh_dt <- function(theta, time_arg) {
    h_fn(theta, time_arg) * dlogh_dt(theta, time_arg)
  }

  Quadrature <- function(f, theta, time_arg) {
    Integral <- 0
    for (i in seq_len(nquad)) {
      wi <- cc$w[i] * time_arg / 2
      ti <- (cc$x[i] + 1) * time_arg / 2
      Integral <- Integral + wi * f(theta, ti)
    }
    Integral
  }

  # Gradient of log h w.r.t. theta: returns (n x p) matrix with column blocks
  #   intercept: 1
  #   d/d_beta0: B(log(t) + eta1)
  #   d/d_beta1: B'(log(t) + eta1)^T beta0 * X1_filtered + X3_filtered
  #   d/d_beta2: X2_filtered
  # Derivs.txt "Derivatives of the log-hazard -- log-time formulation".
  gradient_log_h <- function(theta, time_arg) {
    beta0 <- theta[2:nb0]
    db0   <- cbasis(theta, time_arg, 0L, "obs")
    db1   <- if (has_time_effect) {
      drop(cbasis(theta, time_arg, 1L, "obs") %*% beta0) * X1_filtered + X3_filtered
    } else {
      X3_filtered                       # == 0 when no covariate shifts log-time
    }
    db2   <- X2_filtered
    cbind(1, db0, db1, db2)
  }
  gradient_h <- function(theta, time_arg) {
    h_fn(theta, time_arg) * gradient_log_h(theta, time_arg)
  }

  # Hessian of log h: (n x p x p) array; per Derivs.txt only two blocks are non-zero:
  #   (beta0, beta1): B'(log(t) + eta1) * X1_filtered
  #   (beta1, beta1): B''(log(t) + eta1)^T beta0 * X1f * X1f^T
  # All other blocks (beta0/beta0, beta0/beta2, beta1/beta2, beta2/beta2) are zero.
  hessian_log_h <- function(theta, time_arg) {
    beta0   <- theta[2:nb0]
    Hessian  <- array(0, c(length(time_arg), length(theta), length(theta)))
    if (!has_time_effect) return(Hessian)   # all non-zero blocks vanish
    B_prime  <- cbasis(theta, time_arg, 1L, "obs")
    B_prime2 <- cbasis(theta, time_arg, 2L, "obs")
    index0   <- 2:nb0; index1 <- nb0 + 1:nb1
    for (k in seq_len(nb1)) {
      Hessian[, index0, index1[k]] <-
        Hessian[, index1[k], index0] <- B_prime * X1_filtered[, k]
    }
    for (k in seq_len(nb1)) {
      Hessian[, index1, index1[k]] <-
        drop(B_prime2 %*% beta0) * X1_filtered * X1_filtered[, k]
    }
    Hessian
  }
  hessian_h <- function(theta, time_arg) {
    hv     <- h_fn(theta, time_arg)
    dlogh  <- gradient_log_h(theta, time_arg)
    d2logh <- hessian_log_h(theta, time_arg)
    m      <- ncol(dlogh)
    dlogh2 <- dlogh[, rep(seq_len(m), each = m)] *
              dlogh[, rep(seq_len(m), times = m)]
    dlogh2 <- array(dlogh2, dim(d2logh))
    hv * (dlogh2 + d2logh)
  }

  # H(t|x) = integral_0^t h(u|x) du, approximated by GL quadrature.
  # outer(time, quad_scale) expands all n_obs x nquad evaluation points at once;
  # rowSums collapses back to n_obs values.
  H_fn <- function(theta, time_arg) {
    n_arg  <- length(time_arg)
    beta0  <- theta[2:nb0]
    eta2_v <- rep(eta2(theta), nquad)
    w_q    <- if (identical(time_arg, time)) w_quad_t1
              else if (!is.null(time2) && identical(time_arg, time2)) w_quad_t2
              else c(outer(time_arg, quad_w_half))
    B_all  <- cbasis(theta, time_arg, 0L, "quad")
    h_all  <- exp(drop(B_all %*% beta0) + theta[1] + eta2_v)
    rowSums(matrix(w_q * h_all, n_arg, nquad))
  }

  # dH/d_theta = integral_0^t d(h(u|x))/d_theta du  (Leibniz / differentiate under integral).
  # Loops over quadrature points to keep the (n x p) output in memory.
  gradient_H <- function(theta, time_arg) {
    n_arg  <- length(time_arg)
    beta0  <- theta[2:nb0]
    eta2_v <- rep(eta2(theta), nquad)
    if (identical(time_arg, time)) {
      w_q <- w_quad_t1; X1f <- X1f_quad; X2f <- X2f_quad; X3f <- X3f_quad
    } else if (!is.null(time2) && identical(time_arg, time2)) {
      w_q <- w_quad_t2; X1f <- X1f_quad; X2f <- X2f_quad; X3f <- X3f_quad
    } else {
      w_q <- c(outer(time_arg, quad_w_half))
      ri  <- rep(seq_len(n_arg), nquad)
      X1f <- X1_filtered[ri, , drop = FALSE]
      X2f <- X2_filtered[ri, , drop = FALSE]
      X3f <- X3_filtered[ri, , drop = FALSE]
    }
    B_all <- cbasis(theta, time_arg, 0L, "quad")
    wh    <- w_q * exp(drop(B_all %*% beta0) + theta[1] + eta2_v)
    db1   <- if (has_time_effect) {
      drop(cbasis(theta, time_arg, 1L, "quad") %*% beta0) * X1f + X3f
    } else {
      X3f
    }
    wint  <- wh * cbind(1, B_all, db1, X2f)
    # Sum the nquad blocks (observation index fastest) back to one row per obs.
    result <- wint[seq_len(n_arg), , drop = FALSE]
    for (q in seq_len(nquad)[-1L])
      result <- result + wint[((q - 1L) * n_arg + 1L):(q * n_arg), , drop = FALSE]
    result
  }

  # d2H/d_theta2 = integral_0^t d2(h(u|x))/d_theta2 du  (Leibniz rule).
  # Returns (n x p x p) array; use hessian_H_sum for the fitting path (avoids the large array).
  hessian_H <- function(theta, time_arg) {
    n_arg   <- length(time_arg)
    N_q     <- n_arg * nquad
    w_all_q <- if (identical(time_arg, time)) w_quad_t1
               else if (!is.null(time2) && identical(time_arg, time2)) w_quad_t2
               else c(outer(time_arg, quad_w_half))
    beta0   <- theta[2:nb0]
    eta2_v  <- rep(eta2(theta), nquad)
    B_all   <- cbasis(theta, time_arg, 0L, "quad")
    Bp_all  <- cbasis(theta, time_arg, 1L, "quad")
    Bp2_all <- cbasis(theta, time_arg, 2L, "quad")
    h_all   <- exp(drop(B_all %*% beta0) + theta[1] + eta2_v)
    if (n_arg == n_obs) {
      X1f_q <- X1f_quad; X2f_q <- X2f_quad; X3f_q <- X3f_quad
    } else {
      ri    <- rep(seq_len(n_arg), nquad)
      X1f_q <- X1_filtered[ri, , drop = FALSE]
      X2f_q <- X2_filtered[ri, , drop = FALSE]
      X3f_q <- X3_filtered[ri, , drop = FALSE]
    }
    p      <- length(theta)
    index0 <- 2:nb0; index1 <- nb0 + 1:nb1
    glh    <- cbind(1, B_all, drop(Bp_all %*% beta0) * X1f_q + X3f_q, X2f_q)
    d2lh   <- array(0, c(N_q, p, p))
    for (k in seq_len(nb1)) {
      d2lh[, index0, index1[k]] <-
        d2lh[, index1[k], index0] <- Bp_all * X1f_q[, k]
    }
    Bp2_b0 <- drop(Bp2_all %*% beta0)
    for (k in seq_len(nb1)) {
      d2lh[, index1, index1[k]] <- Bp2_b0 * X1f_q * X1f_q[, k]
    }
    m_val  <- p
    glh2   <- glh[, rep(seq_len(m_val), each = m_val)] *
               glh[, rep(seq_len(m_val), times = m_val)]
    comp   <- glh2 + matrix(d2lh, N_q, p * p)
    wcomp  <- (w_all_q * h_all) * comp
    result_2d <- Reduce("+", lapply(seq_len(nquad), function(q)
      wcomp[((q - 1L) * n_arg + 1L):(q * n_arg), , drop = FALSE]))
    array(result_2d, c(n_arg, p, p))
  }

  # Memory-efficient variant of hessian_H: accumulates directly into a (p x p) result
  # by summing over all observations and quadrature points, avoiding the N_q x p x p array.
  # Used in hessian_negll for "rc" and "lt_rc" -- the common path during model fitting.
  hessian_H_sum <- function(theta, time_arg) {
    n_arg  <- length(time_arg)
    beta0  <- theta[2:nb0]
    eta2_v <- rep(eta2(theta), nquad)
    p      <- length(theta)
    index0 <- 2:nb0; index1 <- nb0 + seq_len(nb1)
    if (identical(time_arg, time)) {
      w_q <- w_quad_t1; X1f <- X1f_quad; X2f <- X2f_quad; X3f <- X3f_quad
    } else if (!is.null(time2) && identical(time_arg, time2)) {
      w_q <- w_quad_t2; X1f <- X1f_quad; X2f <- X2f_quad; X3f <- X3f_quad
    } else {
      w_q <- c(outer(time_arg, quad_w_half))
      ri  <- rep(seq_len(n_arg), nquad)
      X1f <- X1_filtered[ri, , drop = FALSE]
      X2f <- X2_filtered[ri, , drop = FALSE]
      X3f <- X3_filtered[ri, , drop = FALSE]
    }
    B_all <- cbasis(theta, time_arg, 0L, "quad")
    wh    <- w_q * exp(drop(B_all %*% beta0) + theta[1] + eta2_v)
    if (!has_time_effect) {
      # X1_filtered == 0: the d2(log h) blocks vanish and the beta1 columns of
      # the score are zero, so only the score outer-product remains.
      return(crossprod(sqrt(wh) * cbind(1, B_all, X3f, X2f)))
    }
    Bp_all  <- cbasis(theta, time_arg, 1L, "quad")
    Bp2_all <- cbasis(theta, time_arg, 2L, "quad")
    glh     <- cbind(1, B_all, drop(Bp_all %*% beta0) * X1f + X3f, X2f)
    # Score outer-product summed over all observations and quadrature points.
    result  <- crossprod(sqrt(wh) * glh)
    # d2(log h) contributions, summed the same way (crossprod over N_q rows).
    blk01   <- crossprod(Bp_all, wh * X1f)                       # (nb0-1) x nb1
    result[index0, index1] <- result[index0, index1] + blk01
    result[index1, index0] <- result[index1, index0] + t(blk01)
    result[index1, index1] <- result[index1, index1] +
      crossprod(X1f, (wh * drop(Bp2_all %*% beta0)) * X1f)
    result
  }

  # Penalised negative log-likelihood: -l(theta) + lambda/2 * theta^T S theta.
  # Individual log-lik: delta * log h(t) - H(t)  (right-censoring).
  # Left-truncation adds H(t_trunc); interval censoring uses log[S(t1) - S(t2)].
  negll <- function(theta, time_arg, penalized = TRUE, ...) {
    logh <- log_h(theta, time_arg)
    H1   <- H_fn(theta, time_arg)
    penal <- 0
    if (penalized) {
      penal <- lambda / 2 * t(theta[1:nb0]) %*% S_penal %*% theta[1:nb0]
    }
    if (cens_type == "rc") {
      return(-sum(logh * event - H1) + penal)
    }
    H2 <- H_fn(theta, time2)
    if (cens_type == "lt_rc") {
      return(-sum(logh * event - H1 + H2) + penal)
    }
    if (cens_type == "ic") {
      return(-sum(log(exp(-H1) - exp(-H2))) + penal)
    }
  }

  gradient_negll <- function(theta, time_arg, penalized = TRUE, ...) {
    dlogh <- gradient_log_h(theta, time_arg)
    dH    <- gradient_H(theta, time_arg)
    penal <- 0
    if (penalized) {
      penal <- c(lambda * S_penal %*% theta[1:nb0], rep(0, 2 * nb1))
    }
    if (cens_type == "rc") {
      return(-colSums(dlogh * event - dH) + penal)
    }
    dH2 <- gradient_H(theta, time2)
    if (cens_type == "lt_rc") {
      return(-colSums(dlogh * event - dH + dH2) + penal)
    }
    H1 <- H_fn(theta, time_arg)
    H2 <- H_fn(theta, time2)
    if (cens_type == "ic") {
      s <- (exp(-H2) * dH2 - exp(-H1) * dH) / (exp(-H1) - exp(-H2))
      return(-colSums(s) + penal)
    }
  }

  scores <- function(theta, time_arg) {
    dlogh <- gradient_log_h(theta, time_arg)
    dH    <- gradient_H(theta, time_arg)
    if (cens_type == "rc") {
      return(dlogh * event - dH)
    }
    dH2 <- gradient_H(theta, time2)
    if (cens_type == "lt_rc") {
      return(dlogh * event - dH + dH2)
    }
    H1 <- H_fn(theta, time_arg)
    H2 <- H_fn(theta, time2)
    if (cens_type == "ic") {
      s <- (exp(-H2) * dH2 - exp(-H1) * dH) / (exp(-H1) - exp(-H2))
      return(s)
    }
  }

  hessian_negll <- function(theta, time_arg, penalized = TRUE, ...) {
    p      <- length(theta)
    index0 <- 2:nb0; index1 <- nb0 + seq_len(nb1)
    penal  <- matrix(0, p, p)
    if (penalized) {
      penal <- rbind(
        cbind(lambda * S_penal, matrix(0, nb0, 2 * nb1)),
        matrix(0, 2 * nb1, nb0 + 2 * nb1)
      )
    }
    if (cens_type %in% c("rc", "lt_rc")) {
      d2H_sum <- hessian_H_sum(theta, time_arg)
      d2lh_ev <- matrix(0, p, p)
      if (has_time_effect) {
        beta0   <- theta[2:nb0]
        Bp_obs  <- cbasis(theta, time_arg, 1L, "obs")
        Bp2_obs <- cbasis(theta, time_arg, 2L, "obs")
        blk01   <- crossprod(Bp_obs, event * X1_filtered)         # (nb0-1) x nb1
        d2lh_ev[index0, index1] <- blk01
        d2lh_ev[index1, index0] <- t(blk01)
        d2lh_ev[index1, index1] <-
          crossprod(X1_filtered, (event * drop(Bp2_obs %*% beta0)) * X1_filtered)
      }
      if (cens_type == "rc") return(-d2lh_ev + d2H_sum + penal)
      d2H2_sum <- hessian_H_sum(theta, time2)
      return(-d2lh_ev + d2H_sum - d2H2_sum + penal)
    }
    # Interval censoring: array-based path. The penalty is added here just as for
    # rc/lt_rc (negll and gradient_negll already penalise the ic case), so that
    # the penalised Hessian -- and hence fit$var, edf, LCV and dLCV -- are
    # consistent across all censoring types.
    d2H  <- hessian_H(theta, time_arg)
    d2H2 <- hessian_H(theta, time2)
    H1   <- H_fn(theta, time_arg)
    H2   <- H_fn(theta, time2)
    dH   <- gradient_H(theta, time_arg)
    dH2  <- gradient_H(theta, time2)
    S1   <- exp(-H1); S2 <- exp(-H2)
    v    <- (S2 * dH2 - S1 * dH)
    a    <- array(rep(v, times = p), dim = c(length(time_arg), p, p))
    b    <- aperm(a, c(1, 3, 2))
    p1   <- a * b / (S1 - S2)^2
    v    <- dH
    a    <- array(rep(v, times = p), dim = c(length(time_arg), p, p))
    b    <- aperm(a, c(1, 3, 2))
    dH_2 <- a * b
    v    <- dH2
    a    <- array(rep(v, times = p), dim = c(length(time_arg), p, p))
    b    <- aperm(a, c(1, 3, 2))
    dH2_2 <- a * b
    # With g = S1 - S2, the per-observation negll Hessian is -d^2 log(g)/dtheta^2
    #   = -g''/g + g' g'^T / g^2 = p2 + p1,
    # where p2 = -g''/g and p1 = g' g'^T / g^2. (The previous code used
    # -(p2 - p1), flipping the sign of the second-derivative term p2, which did
    # not match the finite-difference Hessian of gradient_negll.)
    p2   <- (S1 * d2H - S1 * dH_2 - S2 * d2H2 + S2 * dH2_2) / (S1 - S2)
    apply(p1 + p2, 2:3, sum) + penal
  }

  # --- res = "H" dispatch (before edf/fit) ---
  if (res == "H") {
    if (!gaussified) {
      inner <- function(t) h_fn(theta, t)
      return(mapply(function(a, b) integrate(inner, a, b)$value,
                    rep(0, length(time)), time))
    }
    return(H_fn(theta, time))
  }

  # Effective degrees of freedom: edf = tr(H_pen^{-1} * H_unpen) - #{constrained params}.
  # modvec adds 1 on diagonal for constrained params (beta1=0 for PH; beta2=0 for AFT/AH)
  # so H_pen + diag(modvec) stays invertible; those positions contribute 0 to the trace.
  # Only the unpenalised Hessian is evaluated; the penalised one differs from it
  # by exactly lambda * S (the penalty touches the nb0 block only and is now applied
  # for every censoring type), so H_pen = H_upen + lambda * S_full -- no second
  # Hessian evaluation needed.
  edf_fn <- function(theta_lambda, roh) {
    modvec <- c(rep(0, nb0), is_ph, is_aft + is_ah)
    tryCatch({
      H_upen <- genhaz_work(theta = theta_lambda, time = time, X = X,
                            knots = knots, Z = Z, event = event,
                            lambda = exp(roh), cens_type = cens_type,
                            model_type = model_type, time2 = time2,
                            res = "hessian_negll_upen", control = control)
      S_full <- matrix(0, length(theta_lambda), length(theta_lambda))
      S_full[seq_len(nb0), seq_len(nb0)] <- S_penal
      H_pen  <- H_upen + exp(roh) * S_full
      sum(diag(solve(H_pen + diag(modvec)) %*% H_upen)) - sum(is_not_gh)
    }, error = function(e) {
      message("Error computing edf for roh = ", roh, ": ", conditionMessage(e))
      Inf
    })
  }

  # First-order analytical gradient dLCV/d_rho (rho = log lambda) via implicit differentiation.
  # Uses d(theta_hat)/d_rho = H_pen^{-1} lambda S theta (implicit function theorem).
  # Ignores the dependence of H_unpen on theta_hat (i.e., treats third-derivative terms as zero).
  # Select via lcv_method = "approx"; faster than dLCV_full.
  lcv_gradient <- function(fit_lambda, roh) {
    lam    <- exp(roh)
    theta  <- fit_lambda$par
    var_m  <- fit_lambda$var
    S_m    <- fit_lambda$S
    nb0_f  <- fit_lambda$nb0
    nb1_f  <- fit_lambda$nb1
    n_f    <- fit_lambda$n
    s_theta <- c(drop(S_m %*% theta[1:nb0_f]), rep(0, 2 * nb1_f))
    lam_s   <- lam * s_theta
    d_unpen <- drop(lam_s %*% (var_m %*% lam_s))
    S_full  <- rbind(cbind(S_m, matrix(0, nb0_f, 2 * nb1_f)),
                     matrix(0, 2 * nb1_f, nb0_f + 2 * nb1_f))
    H_npen  <- -fit_lambda$hessian
    H_nupen <- H_npen - lam * S_full
    d_edf   <- -lam * sum(diag(var_m %*% (S_full %*% (var_m %*% H_nupen))))
    d_unpen + log(n_f) / 2 * d_edf
  }
  
  # LCV criterion to minimise over rho = log(lambda):
  #   LCV(rho) = -l_unpen(theta_hat(rho)) + edf(rho) * log(n) / 2
  # theta_hat is re-optimised at each rho; init is updated in-place to warm-start the next call.
  LCV_fn <- function(roh) {
    lam        <- exp(roh)
    fit_lambda <- tryCatch(
      genhaz_work(init, time, X, knots, Z,
                  event = event, lambda = lam,
                  res = "fit", model_type = model_type,
                  time2 = time2, cens_type = cens_type,
                  gaussified = gaussified, nquad = nquad,
                  control = control),
      error = function(e) NULL
    )
    if (is.null(fit_lambda)) return(Inf)
    theta_lambda <- fit_lambda$par
    init <<- theta_lambda
    edf_val  <- fit_lambda$edf
    penalty  <- (lam / 2) * t(theta_lambda[1:nb0]) %*% S_penal %*% theta_lambda[1:nb0]
    unpen_negll <- fit_lambda$negll - as.numeric(penalty)
    unpen_negll + edf_val * log(length(time)) / 2
  }
  
  # Evaluates lcv_gradient (first-order approximation) at current theta for the given rho.
  dLCV_eval <- function(roh) {
    fl <- genhaz_work(theta, time, X, knots, Z,
                      event = event, lambda = exp(roh), res = "fit",
                      model_type = model_type, time2 = time2,
                      cens_type = cens_type, gaussified = gaussified,
                      nquad = nquad, control = control)
    lcv_gradient(fl, roh)
  }

  # Contracted third derivative of ln h, summed with weights.
  # Returns p×p matrix: Σₚ weights[p] · Σₗ v[l] · d³(ln h_p)/(dθⱼ dθₖ dθₗ).
  # Only two blocks are non-zero (log-time formulation, Derivs.txt):
  #   (β₀, β₁, β₁): B''(ln t + η₁)ₐ · X1f[,b] · (X1f %*% v1)
  #   (β₁, β₁, β₁): B'''(ln t + η₁)ᵀβ₀ · X1f[,a] · X1f[,b] · (X1f %*% v1)
  third_deriv_log_h_contract <- function(theta, time_arg, v,
                                         weights = rep(1, length(time_arg))) {
    beta0   <- theta[2:nb0]
    e1      <- eta1(theta)
    log_te  <- log(time_arg + 1e-10) + e1
    Bp2     <- B(log_te, derivs = 2L)
    Bp3     <- B(log_te, derivs = 3L)
    Bp3_b0  <- drop(Bp3 %*% beta0)
    v0      <- v[2:nb0]
    v1      <- v[nb0 + seq_len(nb1)]
    cv1     <- drop(X1_filtered %*% v1)
    cv0_b2  <- drop(Bp2 %*% v0)   # B''(log_te) %*% v0: contribution via d³lnh/(β₁,β₁,β₀)
    p_dim   <- length(theta)
    index0  <- 2:nb0
    index1  <- nb0 + seq_len(nb1)
    T_mat   <- matrix(0, p_dim, p_dim)
    w_cv1       <- weights * cv1
    contrib_0_1 <- crossprod(Bp2, w_cv1 * X1_filtered)
    T_mat[index0, index1] <- contrib_0_1
    T_mat[index1, index0] <- t(contrib_0_1)
    T_mat[index1, index1] <- crossprod(X1_filtered,
                                       (weights * (Bp3_b0 * cv1 + cv0_b2)) * X1_filtered)
    T_mat
  }

  # Contracted third derivative of H = ∫₀ᵗ h du, via Gauss-Legendre quadrature.
  # Returns p×p matrix: Σₚ ∫₀^{tₚ} Σₗ v[l] · d³h(u)/(dθⱼ dθₖ dθₗ) du.
  # Uses chain rule: d³h/dθ³·v = h·(d3lh_v + gv·d2lh + sym_outer(d2lh_v,glh) + gv·outer(glh,glh))
  # No n×p×p arrays: all accumulation via p×p crossprod operations.
  compute_third_deriv_H_contract <- function(theta, time_arg, v) {
    n_arg  <- length(time_arg)
    beta0  <- theta[2:nb0]
    e1     <- eta1(theta)
    e2     <- eta2(theta)
    p_dim  <- length(theta)
    index0 <- 2:nb0
    index1 <- nb0 + seq_len(nb1)
    v0     <- v[index0]
    v1     <- v[index1]
    X1f    <- if (n_arg == n_obs) X1_filtered else X1_filtered[seq_len(n_arg), , drop = FALSE]
    X2f    <- if (n_arg == n_obs) X2_filtered else X2_filtered[seq_len(n_arg), , drop = FALSE]
    X3f    <- if (n_arg == n_obs) X3_filtered else X3_filtered[seq_len(n_arg), , drop = FALSE]
    cv1    <- drop(X1f %*% v1)
    T_H    <- matrix(0, p_dim, p_dim)
    for (q in seq_len(nquad)) {
      ti      <- time_arg * quad_scale[q]
      wi      <- time_arg * quad_w_half[q]
      log_tei <- log(ti + 1e-10) + e1
      Bi      <- B(log_tei)
      Bpi     <- B(log_tei, derivs = 1L)
      Bp2i    <- B(log_tei, derivs = 2L)
      Bp3i    <- B(log_tei, derivs = 3L)
      Bp_b0i  <- drop(Bpi  %*% beta0)
      Bp2_b0i <- drop(Bp2i %*% beta0)
      Bp3_b0i <- drop(Bp3i %*% beta0)
      hi      <- exp(drop(Bi %*% beta0) + theta[1] + e2)
      whi     <- wi * hi
      glhi    <- cbind(1, Bi, Bp_b0i * X1f + X3f, X2f)
      gv      <- drop(glhi %*% v)
      cv0      <- drop(Bpi  %*% v0)
      cv0_b2   <- drop(Bp2i %*% v0)  # B''(log_te) %*% v0: contribution via d³lnh/(β₁,β₁,β₀)
      # Hessian of ln h times v (n × p_dim, non-zero only at index0 and index1)
      d2lh_v  <- matrix(0, n_arg, p_dim)
      d2lh_v[, index0] <- Bpi * cv1
      d2lh_v[, index1] <- (Bp2_b0i * cv1 + cv0) * X1f
      # Sparse: (d3lh_v + gv·d2lh) contributions to (index0,index1) and (index1,index1)
      coef_0_1    <- gv * Bpi + cv1 * Bp2i
      contrib_0_1 <- crossprod(coef_0_1, whi * X1f)
      T_H[index0, index1] <- T_H[index0, index1] + contrib_0_1
      T_H[index1, index0] <- T_H[index1, index0] + t(contrib_0_1)
      T_H[index1, index1] <- T_H[index1, index1] +
        crossprod(X1f, (whi * (gv * Bp2_b0i + cv1 * Bp3_b0i + cv0_b2)) * X1f)
      # Dense: sym_outer(d2lh_v, glh) + gv·outer(glh, glh)
      T_H <- T_H +
        crossprod(d2lh_v, whi * glhi) +
        crossprod(glhi,   whi * d2lh_v) +
        crossprod(glhi,  (whi * gv) * glhi)
    }
    T_H
  }

  # Full LCV gradient including third-derivative correction.
  # Extends lcv_gradient() by adding Trace(var_m·d3l_v·(var_m·H_nupen − I)).
  # Supports cens_type "rc" and "lt_rc"; falls back to lcv_gradient for "ic".
  dLCV_full <- function(fit_lambda, roh) {
    if (cens_type == "ic") return(lcv_gradient(fit_lambda, roh))
    lam    <- exp(roh)
    theta_ <- fit_lambda$par
    var_m  <- fit_lambda$var
    S_m    <- fit_lambda$S
    nb0_f  <- fit_lambda$nb0
    nb1_f  <- fit_lambda$nb1
    n_f    <- fit_lambda$n
    p      <- length(theta_)
    S_full <- rbind(cbind(S_m, matrix(0, nb0_f, 2 * nb1_f)),
                    matrix(0, 2 * nb1_f, nb0_f + 2 * nb1_f))
    s_th   <- c(drop(S_m %*% theta_[1:nb0_f]), rep(0, 2 * nb1_f))
    lam_s  <- lam * s_th
    # dθ̂/dρ = H_pen⁻¹ λSθ = −var_m · lam_s
    dth    <- -drop(var_m %*% lam_s)
    d_unpen <- drop(lam_s %*% (var_m %*% lam_s))
    H_npen  <- -fit_lambda$hessian
    H_nupen <- H_npen - lam * S_full
    d_edf_current <- -lam * sum(diag(var_m %*% (S_full %*% (var_m %*% H_nupen))))
    # d3l_v = contracted third deriv of l with dθ̂/dρ
    T_obs <- third_deriv_log_h_contract(theta_, time, dth, weights = event)
    T_H1  <- compute_third_deriv_H_contract(theta_, time, dth)
    d3l_v <- T_obs - T_H1
    if (cens_type == "lt_rc") {
      d3l_v <- d3l_v + compute_third_deriv_H_contract(theta_, time2, dth)
    }
    # Correction: Trace(var_m · d3l_v · (var_m · H_nupen − I))
    var_d3l  <- var_m %*% d3l_v
    var_H_nu <- var_m %*% H_nupen
    d_edf_correction <- sum(var_d3l * t(var_H_nu - diag(p)))
    d_edf_full <- d_edf_current + d_edf_correction
    d_unpen + log(n_f) / 2 * d_edf_full
  }

  dLCV_full_eval <- function(roh) {
    fl <- genhaz_work(theta, time, X, knots, Z,
                      event = event, lambda = exp(roh), res = "fit",
                      model_type = model_type, time2 = time2,
                      cens_type = cens_type, gaussified = gaussified,
                      nquad = nquad, control = control)
    dLCV_full(fl, roh)
  }

  # --- res = "fit": optimise theta at fixed lambda, compute SE / edf / AIC / LCV ---
  # Sub-model constraints: lower = upper = 0 for beta1 (PH) or beta2 (AFT/AH) as box bounds.
  # AH post-processing sets beta2 <- -beta1 after optimisation.
  # Returns: par, se, z, p_values, var, edf, AIC, LCV, lambda, knots, Z, S, n, nb0, nb1.
  if (res == "fit") {
    lower <- rep(-Inf, length(theta))
    upper <- rep(Inf, length(theta))
    for (i in seq_len(nb1)) {
      if (model_type[i] == "PH")  lower[nb0 + i] <- upper[nb0 + i] <- 0
      if (model_type[i] %in% c("AH", "AFT"))
        lower[nb0 + nb1 + i] <- upper[nb0 + nb1 + i] <- 0
    }
    if (!exists("init", inherits = FALSE) || is.null(init)) init <- theta
    fit <- nlminb(init, negll, gradient_negll, hessian_negll,
                  lower = lower, upper = upper, time = time, penalized = TRUE)
    fit$par <- fit$par * c(rep(1, nb0), is_not_ph, is_ph + is_gh) +
      c(rep(0, nb0 + nb1), -is_ah * fit$par[(nb0 + 1):(nb0 + nb1)])
    fit$parnames <- c("intercept", paste0("s", 1:(nb0 - 1)),
                      paste0("beta1_", seq_len(nb1)),
                      paste0("beta2_", seq_len(nb2)))
    names(fit$par) <- fit$parnames
    fit$hessian <- -hessian_negll(fit$par, time)
    rownames(fit$hessian) <- colnames(fit$hessian) <- fit$parnames
    modvec <- c(rep(0, nb0), is_ph, is_aft + is_ah)
    h_pen_inv <- tryCatch(
      solve(-fit$hessian + diag(modvec)) - diag(modvec),
      error = function(e) {
        warning("Non-invertible Hessian: ", conditionMessage(e),
                " -- standard errors set to NA")
        matrix(NA_real_, nrow = length(fit$par), ncol = length(fit$par))
      }
    )
    fit$var      <- h_pen_inv
    rownames(fit$var) <- colnames(fit$var) <- fit$parnames
    fit$se       <- sqrt(diag(fit$var))
    fit$z        <- fit$par / fit$se
    fit$p_values <- 2 * (1 - pnorm(abs(fit$z)))
    fit$df       <- length(theta) - sum(is_not_gh)
    fit$AIC      <- 2 * fit$objective + fit$df
    fit$model_type <- model_type
    fit$knots    <- knots
    fit$Z        <- Z
    fit$S        <- S_penal
    fit$negll    <- fit$objective
    # edf = tr(H_pen^{-1} H_upen) - #{constrained}. Reuse the Hessian and its
    # (penalised) inverse already computed above instead of re-evaluating the
    # Hessian twice via edf_fn(): H_pen = -fit$hessian, H_upen = H_pen - lambda*S,
    # and solve(H_pen + diag(modvec)) == h_pen_inv + diag(modvec).
    if (anyNA(h_pen_inv)) {
      fit$edf <- Inf
    } else {
      # H_pen = -fit$hessian; H_upen = H_pen - lambda*S (the penalty is now applied
      # for every censoring type). solve(H_pen + diag(modvec)) == h_pen_inv + diag(modvec).
      S_full <- matrix(0, length(theta), length(theta))
      S_full[seq_len(nb0), seq_len(nb0)] <- S_penal
      H_upen  <- -fit$hessian - lambda * S_full
      fit$edf <- sum(diag((h_pen_inv + diag(modvec)) %*% H_upen)) - sum(is_not_gh)
    }
    fit$lambda   <- lambda
    fit$penalty  <- (lambda / 2) * t(fit$par[1:nb0]) %*% S_penal %*% fit$par[1:nb0]
    fit$unpen_negll <- fit$negll - as.numeric(fit$penalty)
    fit$LCV      <- fit$unpen_negll + fit$edf * log(length(time)) / 2
    fit$n        <- length(time)
    fit$nb0      <- nb0
    fit$nb1      <- nb1
    fit$X_center <- attr(X, "scaled:center")
    names(fit$se) <- fit$parnames
    return(fit)
  }

  # --- res = "fit_LCV" ---
  if (res == "fit_LCV") {
    init <<- theta
    if (lcv_method == "optimize") {
      roh_opt <- optimize(LCV_fn, interval = interval, tol = tol_LCV)$minimum
    } else {
      dlcv_fn <- if (lcv_method == "full") dLCV_full else lcv_gradient
      dLCV_fn <- function(roh) {
        fl <- genhaz_work(init, time, X, knots, Z,
                          event = event, lambda = exp(roh), res = "fit",
                          model_type = model_type, time2 = time2,
                          cens_type = cens_type, gaussified = gaussified,
                          nquad = nquad, control = control)
        init <<- fl$par
        dlcv_fn(fl, roh)
      }
      roh_opt <- tryCatch(
        uniroot(dLCV_fn, interval = interval, tol = tol_LCV,
                extendInt = "yes")$root,
        error = function(e) optimize(LCV_fn, interval = interval,
                                     tol = tol_LCV)$minimum
      )
    }
    fit <- genhaz_work(init, time, X, knots, Z,
                       event = event, lambda = exp(roh_opt), res = "fit",
                       model_type = model_type, time2 = time2,
                       cens_type = cens_type, gaussified = gaussified,
                       nquad = nquad, control = control)
    fit$lambda  <- exp(roh_opt)
    fit$roh_opt <- roh_opt
    return(fit)
  }

  # --- remaining res values ---
  out <- switch(res,
    basis          = B(log(time + 1e-10) + eta1(theta)),
    log_h          = log_h(theta, time),
    gradient_log_h = gradient_log_h(theta, time),
    hessian_log_h  = hessian_log_h(theta, time),
    dlogh_dt       = dlogh_dt(theta, time),
    dh_dt          = dh_dt(theta, time),
    h              = h_fn(theta, time),
    gradient_h     = gradient_h(theta, time),
    hessian_h      = hessian_h(theta, time),
    negll          = negll(theta, time, penalized = TRUE),
    negll_upen     = negll(theta, time, penalized = FALSE),
    gradient_negll = gradient_negll(theta, time, penalized = TRUE),
    scores         = scores(theta, time),
    hessian_negll  = hessian_negll(theta, time, penalized = TRUE),
    hessian_negll_upen = hessian_negll(theta, time, penalized = FALSE),
    gradient_H     = gradient_H(theta, time),
    hessian_H      = hessian_H(theta, time),
    LCV            = LCV_fn(log(lambda)),
    dLCV           = dLCV_eval(log(lambda)),
    dLCV_full      = dLCV_full_eval(log(lambda)),
    edf            = edf_fn(theta, log(lambda)),
    NULL
  )
  if (!is.null(out)) return(out)
  stop("res not matched: ", res)
}
