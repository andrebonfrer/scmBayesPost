# =============================================================================
#  scmBayesPost — Gibbs samplers (list-based, no X_block materialisation)
#
#  All samplers operate on X_list / y_list / w_list from gdata directly.
#  X_block (157M+ rows) is never constructed. Per-unit [K x K] systems
#  are solved independently, collapsing to scalar arithmetic for K=1.
#
#  Contents
#  --------
#  resolve_sampler_control()
#  gibbs_postscm()                      dispatcher
#  gibbs_sampling_simple()              no first stage, no moderators
#  gibbs_sampling_moderators()          second-stage moderators only
#  gibbs_sampling_selection()           Bayesian probit first stage
#  gibbs_sampling_selection_moderators() first stage + moderators
#  .precompute_unit_stats()             precompute XtWX, XtWy per unit
#  .sample_beta_units()                 sample beta unit by unit
#  .compute_wss()                       weighted sum of squares
#  .sample_z_star()                     fast truncated normal sampler
# =============================================================================


# -----------------------------------------------------------------------------
#' Resolve sampler control parameters
#' @keywords internal
resolve_sampler_control <- function(control = NULL, has_Z = FALSE) {
  defaults <- list(
    Sigma_gamma_prior   = 10,
    mu_gamma_prior      = NULL,
    a_sigma_alpha_prior = 2,
    b_sigma_alpha_prior = 2,
    a_sigma_tau_prior   = 2,
    b_sigma_tau_prior   = 2,
    mu_delta_prior      = NULL,
    Sigma_delta_prior   = NULL,
    sigma2_rho_prior    = 10
  )
  if (is.null(control)) return(defaults)
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0)
    stop(sprintf("Unknown control parameter(s): %s",
                 paste(unknown, collapse = ", ")), call. = FALSE)
  utils::modifyList(defaults, control)
}


# -----------------------------------------------------------------------------
#' Precompute per-unit sufficient statistics
#'
#' Computes X_j'W_j X_j and X_j'W_j y_j for each unit j. Called once before
#' the Gibbs loop (XtWX never changes; XtWy is recomputed each iteration
#' via .compute_XtWy_units when y_tilde changes).
#'
#' @param X_list List of J0 matrices \[n_j x K\].
#' @param w_list List of J0 weight vectors \[n_j\].
#' @param y_list List of J0 outcome vectors \[n_j\]. NULL for XtWX only.
#' @param K      Integer. Covariates per unit.
#' @param J0     Integer. Number of treated units.
#'
#' @return List with XtWX (list of J0 \[K x K\] matrices) and optionally
#'   XtWy (list of J0 \[K\] vectors).
#' @keywords internal
.precompute_unit_stats <- function(X_list, w_list, y_list = NULL, K, J0) {
  XtWX <- vector("list", J0)
  XtWy <- if (!is.null(y_list)) vector("list", J0) else NULL

  for (j in seq_len(J0)) {
    Xj <- X_list[[j]]
    wj <- w_list[[j]]
    if (K == 1L) {
      xj        <- as.numeric(Xj)
      XtWX[[j]] <- sum(wj * xj^2)           # scalar
      if (!is.null(y_list))
        XtWy[[j]] <- sum(wj * xj * y_list[[j]])  # scalar
    } else {
      WXj       <- wj * Xj
      XtWX[[j]] <- crossprod(Xj, WXj)       # [K x K]
      if (!is.null(y_list))
        XtWy[[j]] <- as.numeric(crossprod(Xj, wj * y_list[[j]]))  # [K]
    }
  }
  list(XtWX = XtWX, XtWy = XtWy)
}


# -----------------------------------------------------------------------------
#' Compute XtWy per unit for current y_tilde
#' @keywords internal
.compute_XtWy_units <- function(X_list, w_list, y_tilde_list, K, J0) {
  XtWy <- vector("list", J0)
  for (j in seq_len(J0)) {
    Xj <- X_list[[j]]
    wj <- w_list[[j]]
    yj <- y_tilde_list[[j]]
    XtWy[[j]] <- if (K == 1L) sum(wj * as.numeric(Xj) * yj)
    else as.numeric(crossprod(Xj, wj * yj))
  }
  XtWy
}


# -----------------------------------------------------------------------------
#' Sample beta unit by unit (block-diagonal posterior)
#'
#' For K=1: each unit's posterior is scalar N(m_j, v_j) — no matrix ops.
#' For K>1: Cholesky of \[K x K\] per unit.
#'
#' @param XtWX      List of J0 \[K x K\] matrices (precomputed).
#' @param XtWy      List of J0 \[K\] vectors (current iteration).
#' @param tau       Numeric vector length K: current dispersion SDs.
#' @param sigma2    Scalar: current observation noise variance.
#' @param K         Integer.
#' @param J0        Integer.
#' @param prior_mean Optional list of J0 \[K\] vectors. Default NULL (zero).
#'
#' @return Numeric vector length K*J0.
#' @keywords internal
.sample_beta_units <- function(XtWX, XtWy, tau, sigma2, K, J0,
                               prior_mean = NULL) {
  tau_inv2 <- 1 / tau^2    # precision
  beta_out <- numeric(K * J0)

  for (j in seq_len(J0)) {
    pm_j <- if (is.null(prior_mean)) rep(0, K) else prior_mean[[j]]

    if (K == 1L) {
      # Scalar fast path — no matrix operations at all
      v_j <- 1 / (XtWX[[j]] / sigma2 + tau_inv2)
      m_j <- v_j * (XtWy[[j]] / sigma2 + tau_inv2 * pm_j)
      beta_out[j] <- stats::rnorm(1L, mean = m_j, sd = sqrt(v_j))
    } else {
      A_j  <- XtWX[[j]] / sigma2 + diag(tau_inv2, K)
      b_j  <- XtWy[[j]] / sigma2 + tau_inv2 * pm_j
      cA_j <- chol(A_j)
      Ai_j <- chol2inv(cA_j)
      mu_j <- as.numeric(Ai_j %*% b_j)
      beta_out[((j-1L)*K+1L):(j*K)] <-
        mu_j + as.numeric(backsolve(cA_j, stats::rnorm(K)))
    }
  }
  beta_out
}


# -----------------------------------------------------------------------------
#' Weighted sum of squared residuals (for sigma2 update)
#' @keywords internal
.compute_wss <- function(X_list, w_list, y_list, beta_bd, K, J0,
                         nu_list = NULL, rho = 0) {
  wss <- 0
  for (j in seq_len(J0)) {
    bj  <- beta_bd[((j-1L)*K+1L):(j*K)]
    Xj  <- X_list[[j]]
    rj  <- y_list[[j]] - if (K == 1L) as.numeric(Xj) * bj
    else as.numeric(Xj %*% bj)
    if (!is.null(nu_list)) rj <- rj - rho * nu_list[[j]]
    wss <- wss + sum(w_list[[j]] * rj^2)
  }
  wss
}


# -----------------------------------------------------------------------------
#' Fast truncated normal sampler (inverse-CDF, no ifelse dispatch)
#' @keywords internal
.sample_z_star <- function(eta, d_vec, treat_idx, ctrl_idx) {
  z <- numeric(length(eta))
  if (length(treat_idx) > 0L) {
    mu_t <- eta[treat_idx]
    p_lo <- stats::pnorm(-mu_t)
    u    <- stats::runif(length(treat_idx), min = p_lo, max = 1)
    z[treat_idx] <- mu_t + stats::qnorm(u)
  }
  if (length(ctrl_idx) > 0L) {
    mu_c <- eta[ctrl_idx]
    p_hi <- stats::pnorm(-mu_c)
    u    <- stats::runif(length(ctrl_idx), min = 0, max = p_hi)
    z[ctrl_idx] <- mu_c + stats::qnorm(u)
  }
  z
}


# =============================================================================
#' Gibbs sampler dispatcher
#'
#' @param gdata   Output from prepare_data_general().
#' @param n_iter  Integer. Total Gibbs iterations.
#' @param burn_in Integer. Burn-in iterations discarded.
#' @param control Optional list of prior/control parameters.
#'
#' @return List of posterior draw matrices.
#' @export
gibbs_postscm <- function(gdata,
                          n_iter  = 1000,
                          burn_in = 500,
                          control = NULL) {

  if (is.null(gdata$cov))
    stop("gdata$cov missing.", call. = FALSE)
  if (is.null(gdata$cov$intX))
    stop("gdata$cov$intX not set.", call. = FALSE)

  has_Z  <- !is.null(gdata$Z_block)
  has_fs <- !is.null(gdata$first_stage) &&
    isTRUE(gdata$cov$first_stage == "selection_probit_bayes")
  M      <- if (!is.null(gdata$cov$M)) gdata$cov$M else 1L

  ctrl <- resolve_sampler_control(control = control, has_Z = has_Z)

  if (M == 1L && has_fs && has_Z)
    return(gibbs_sampling_selection_moderators(gdata, n_iter, burn_in, ctrl))
  if (M == 1L && has_fs && !has_Z)
    return(gibbs_sampling_selection(gdata, n_iter, burn_in, ctrl))
  if (M == 1L && !has_Z)
    return(gibbs_sampling_simple(gdata, n_iter, burn_in, ctrl))
  if (M == 1L && has_Z)
    return(gibbs_sampling_moderators(gdata, n_iter, burn_in, ctrl))
  if (M > 1L)
    stop("Multi-outcome model not yet implemented.", call. = FALSE)

  stop("Unhandled model configuration.", call. = FALSE)
}


# =============================================================================
#' Gibbs sampler — no first stage, no moderators
#' @keywords internal
gibbs_sampling_simple <- function(gdata, n_iter = 1000, burn_in = 500,
                                  control = NULL) {
  ctrl <- resolve_sampler_control(control, has_Z = FALSE)

  X_list <- gdata$X_list
  y_list <- gdata$y_list
  w_list <- gdata$w_list
  K      <- length(gdata$cov$Xcols)
  J0     <- gdata$cov$J0

  a_sa <- ctrl$a_sigma_alpha_prior; b_sa <- ctrl$b_sigma_alpha_prior
  a_st <- ctrl$a_sigma_tau_prior;   b_st <- ctrl$b_sigma_tau_prior

  # Precompute XtWX once (fixed); XtWy also fixed (y never changes here)
  pre    <- .precompute_unit_stats(X_list, w_list, y_list, K, J0)
  XtWX   <- pre$XtWX
  XtWy0  <- pre$XtWy   # fixed — no rho/nu in this sampler
  N_total <- sum(sapply(y_list, length))

  beta_bd <- numeric(K * J0)
  sigma2  <- 1
  tau     <- rep(1, K)

  n_save         <- n_iter - burn_in
  beta_samples   <- matrix(0, n_save, K * J0)
  sigma2_samples <- numeric(n_save)
  tau_samples    <- matrix(0, n_save, K)

  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)

  for (iter in seq_len(n_iter)) {

    # beta
    beta_bd <- .sample_beta_units(XtWX, XtWy0, tau, sigma2, K, J0)

    # sigma2
    wss    <- .compute_wss(X_list, w_list, y_list, beta_bd, K, J0)
    sigma2 <- 1 / stats::rgamma(1, shape = a_sa + N_total / 2,
                                rate  = b_sa + wss / 2)

    # tau
    bmat <- matrix(beta_bd, ncol = K, byrow = TRUE)
    for (k in seq_len(K))
      tau[k] <- sqrt(1 / stats::rgamma(1,
                                       shape = a_st + J0 / 2,
                                       rate  = b_st + sum(bmat[, k]^2) / 2))

    if (iter > burn_in) {
      s <- iter - burn_in
      beta_samples[s, ] <- beta_bd
      sigma2_samples[s] <- sigma2
      tau_samples[s, ]  <- tau
    }
    utils::setTxtProgressBar(pb, iter)
  }
  close(pb)
  colnames(tau_samples) <- gdata$cov$Xcols

  list(beta_samples = beta_samples, sigma2_samples = sigma2_samples,
       tau_samples  = tau_samples)
}


# =============================================================================
#' Gibbs sampler — second-stage moderators, no first stage
#' @keywords internal
gibbs_sampling_moderators <- function(gdata, n_iter = 1000, burn_in = 500,
                                      control = NULL) {
  ctrl <- resolve_sampler_control(control, has_Z = TRUE)

  X_list <- gdata$X_list
  y_list <- gdata$y_list
  w_list <- gdata$w_list
  Z      <- gdata$Z_block
  K      <- length(gdata$cov$Xcols)
  J0     <- gdata$cov$J0
  G      <- ncol(Z)
  k_tr   <- gdata$cov$intX

  a_sa <- ctrl$a_sigma_alpha_prior; b_sa <- ctrl$b_sigma_alpha_prior
  a_st <- ctrl$a_sigma_tau_prior;   b_st <- ctrl$b_sigma_tau_prior

  mu_gamma_prior <- if (is.null(ctrl$mu_gamma_prior)) rep(0, G)
  else as.numeric(ctrl$mu_gamma_prior)
  Sigma_gamma_prior_inv <- if (length(ctrl$Sigma_gamma_prior) == 1)
    diag(1 / ctrl$Sigma_gamma_prior, G)
  else diag(1 / as.numeric(ctrl$Sigma_gamma_prior), G)

  pre   <- .precompute_unit_stats(X_list, w_list, NULL, K, J0)
  XtWX  <- pre$XtWX
  N_total <- sum(sapply(y_list, length))

  beta_bd <- numeric(K * J0)
  gamma   <- rep(0, G)
  sigma2  <- 1
  tau     <- rep(1, K)

  n_save         <- n_iter - burn_in
  beta_samples   <- matrix(0, n_save, K * J0)
  gamma_samples  <- matrix(0, n_save, G)
  sigma2_samples <- numeric(n_save)
  tau_samples    <- matrix(0, n_save, K)

  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)

  for (iter in seq_len(n_iter)) {

    # Moderator prior mean: per-unit, only treatment covariate k_tr gets z_j'gamma
    z_gamma <- as.numeric(Z %*% gamma)   # [J0]
    pm_list <- lapply(seq_len(J0), function(j) {
      pm <- rep(0, K); pm[k_tr] <- z_gamma[j]; pm
    })

    # beta with moderator prior mean
    XtWy <- .compute_XtWy_units(X_list, w_list, y_list, K, J0)
    beta_bd <- .sample_beta_units(XtWX, XtWy, tau, sigma2, K, J0,
                                  prior_mean = pm_list)

    # gamma
    bmat    <- matrix(beta_bd, ncol = K, byrow = TRUE)
    beta_tr <- bmat[, k_tr]
    V_gamma <- solve(crossprod(Z) / tau[k_tr]^2 + Sigma_gamma_prior_inv)
    rhs     <- as.numeric(crossprod(Z, beta_tr)) / tau[k_tr]^2 +
      as.numeric(Sigma_gamma_prior_inv %*% mu_gamma_prior)
    gamma   <- rMVNormCovariance(1, mu = as.numeric(V_gamma %*% rhs),
                                 Sigma = V_gamma)

    # sigma2
    wss    <- .compute_wss(X_list, w_list, y_list, beta_bd, K, J0)
    sigma2 <- 1 / stats::rgamma(1, shape = a_sa + N_total / 2,
                                rate  = b_sa + wss / 2)

    # tau — moderator-adjusted residuals for treatment covariate
    z_gamma <- as.numeric(Z %*% gamma)
    for (k in seq_len(K)) {
      centered <- bmat[, k] - if (k == k_tr) z_gamma else 0
      tau[k]   <- sqrt(1 / stats::rgamma(1,
                                         shape = a_st + J0 / 2,
                                         rate  = b_st + sum(centered^2) / 2))
    }

    if (iter > burn_in) {
      s <- iter - burn_in
      beta_samples[s, ]  <- beta_bd
      gamma_samples[s, ] <- as.numeric(gamma)
      sigma2_samples[s]  <- sigma2
      tau_samples[s, ]   <- tau
    }
    utils::setTxtProgressBar(pb, iter)
  }
  close(pb)
  colnames(gamma_samples) <- colnames(Z)
  colnames(tau_samples)   <- gdata$cov$Xcols

  list(beta_samples = beta_samples, gamma_samples = gamma_samples,
       sigma2_samples = sigma2_samples, tau_samples = tau_samples)
}


# =============================================================================
#' Gibbs sampler — Bayesian probit first stage, no moderators
#' @keywords internal
gibbs_sampling_selection <- function(gdata, n_iter = 1000, burn_in = 500,
                                     control = NULL) {
  ctrl <- resolve_sampler_control(control, has_Z = FALSE)

  X_list       <- gdata$X_list
  y_list       <- gdata$y_list
  w_list       <- gdata$w_list
  row_idx_list <- gdata$row_idx_list
  K            <- length(gdata$cov$Xcols)
  J0           <- gdata$cov$J0

  fs    <- gdata$first_stage
  X_fs  <- fs$X_fs
  d_vec <- as.numeric(fs$d)
  n_obs <- length(d_vec)
  p_fs  <- ncol(X_fs)

  a_sa <- ctrl$a_sigma_alpha_prior; b_sa <- ctrl$b_sigma_alpha_prior
  a_st <- ctrl$a_sigma_tau_prior;   b_st <- ctrl$b_sigma_tau_prior
  s2rp <- if (!is.null(ctrl$sigma2_rho_prior)) ctrl$sigma2_rho_prior else 10

  mu_delta_prior <- if (is.null(ctrl$mu_delta_prior)) rep(0, p_fs)
  else as.numeric(ctrl$mu_delta_prior)
  Sigma_delta_prior_inv <- if (is.null(ctrl$Sigma_delta_prior))
    diag(1 / 10, p_fs)
  else {
    sdp <- as.numeric(ctrl$Sigma_delta_prior)
    if (length(sdp) == 1) diag(1 / sdp, p_fs) else diag(1 / sdp)
  }

  # Precompute fixed quantities
  XfsXfs       <- crossprod(X_fs)
  A_delta_base <- XfsXfs + Sigma_delta_prior_inv
  pre          <- .precompute_unit_stats(X_list, w_list, NULL, K, J0)
  XtWX         <- pre$XtWX
  N_total      <- sum(sapply(y_list, length))
  treat_idx    <- which(d_vec == 1L)
  ctrl_idx     <- which(d_vec == 0L)

  # Precompute per-unit weight sums for rho update
  w_unit_sums <- sapply(w_list, sum)

  beta_bd <- numeric(K * J0)
  rho     <- 0
  delta   <- if (!is.null(fs$delta0)) as.numeric(fs$delta0) else rep(0, p_fs)
  sigma2  <- 1
  tau     <- rep(1, K)
  z_star  <- ifelse(d_vec == 1, 0.5, -0.5)

  n_save         <- n_iter - burn_in
  beta_samples   <- matrix(0, n_save, K * J0)
  rho_samples    <- numeric(n_save)
  delta_samples  <- matrix(0, n_save, p_fs)
  sigma2_samples <- numeric(n_save)
  tau_samples    <- matrix(0, n_save, K)

  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)

  for (iter in seq_len(n_iter)) {

    # BLOCK 1: z*
    eta    <- as.numeric(X_fs %*% delta)
    z_star <- .sample_z_star(eta, d_vec, treat_idx, ctrl_idx)

    # BLOCK 2: delta
    b_delta     <- as.numeric(crossprod(X_fs, z_star)) +
      as.numeric(Sigma_delta_prior_inv %*% mu_delta_prior)
    cA_d        <- chol(A_delta_base)
    A_delta_inv <- chol2inv(cA_d)
    delta <- rMVNormCovariance(1,
                               mu    = as.numeric(A_delta_inv %*% b_delta),
                               Sigma = A_delta_inv)

    # BLOCK 3: nu_hat — split by unit using row_idx_list (no stacking)
    nu_hat   <- z_star - as.numeric(X_fs %*% delta)
    nu_list  <- lapply(row_idx_list, function(idx) nu_hat[idx])

    # BLOCK 4a: beta — y_tilde_j = y_j - rho * nu_j
    y_tilde_list <- lapply(seq_len(J0), function(j)
      y_list[[j]] - rho * nu_list[[j]])
    XtWy <- .compute_XtWy_units(X_list, w_list, y_tilde_list, K, J0)
    beta_bd <- .sample_beta_units(XtWX, XtWy, tau, sigma2, K, J0)

    # BLOCK 4b: rho (scalar)
    # nuWnu = sum_j sum_i w_ij * nu_ij^2
    # nuWy2 = sum_j sum_i w_ij * nu_ij * (y_ij - x_ij' beta_j)
    nuWnu <- 0; nuWy2 <- 0
    for (j in seq_len(J0)) {
      bj   <- beta_bd[((j-1L)*K+1L):(j*K)]
      Xj   <- X_list[[j]]
      r2j  <- y_list[[j]] - if (K==1L) as.numeric(Xj)*bj
      else as.numeric(Xj %*% bj)
      wj   <- w_list[[j]]
      nuj  <- nu_list[[j]]
      nuWnu <- nuWnu + sum(wj * nuj^2)
      nuWy2 <- nuWy2 + sum(wj * nuj * r2j)
    }
    v_rho <- 1 / (nuWnu / sigma2 + 1 / s2rp)
    rho   <- stats::rnorm(1, mean = v_rho * nuWy2 / sigma2, sd = sqrt(v_rho))

    # BLOCK 5: sigma2
    wss    <- .compute_wss(X_list, w_list, y_list, beta_bd, K, J0,
                           nu_list = nu_list, rho = rho)
    sigma2 <- 1 / stats::rgamma(1, shape = a_sa + N_total / 2,
                                rate  = b_sa + wss / 2)

    # BLOCK 6: tau
    bmat <- matrix(beta_bd, ncol = K, byrow = TRUE)
    for (k in seq_len(K))
      tau[k] <- sqrt(1 / stats::rgamma(1,
                                       shape = a_st + J0 / 2,
                                       rate  = b_st + sum(bmat[, k]^2) / 2))

    if (iter > burn_in) {
      s <- iter - burn_in
      beta_samples[s, ]  <- beta_bd
      rho_samples[s]     <- rho
      delta_samples[s, ] <- as.numeric(delta)
      sigma2_samples[s]  <- sigma2
      tau_samples[s, ]   <- tau
    }
    utils::setTxtProgressBar(pb, iter)
  }
  close(pb)
  colnames(delta_samples) <- colnames(X_fs)
  colnames(tau_samples)   <- gdata$cov$Xcols

  list(beta_samples = beta_samples, rho_samples = rho_samples,
       delta_samples = delta_samples, sigma2_samples = sigma2_samples,
       tau_samples = tau_samples)
}


# =============================================================================
#' Gibbs sampler — Bayesian probit first stage + second-stage moderators
#' @keywords internal
gibbs_sampling_selection_moderators <- function(gdata, n_iter = 1000,
                                                burn_in = 500,
                                                control = NULL) {
  ctrl <- resolve_sampler_control(control, has_Z = TRUE)

  X_list       <- gdata$X_list
  y_list       <- gdata$y_list
  w_list       <- gdata$w_list
  row_idx_list <- gdata$row_idx_list
  Z            <- gdata$Z_block
  K            <- length(gdata$cov$Xcols)
  J0           <- gdata$cov$J0
  G            <- ncol(Z)
  k_tr         <- gdata$cov$intX

  fs    <- gdata$first_stage
  X_fs  <- fs$X_fs
  d_vec <- as.numeric(fs$d)
  n_obs <- length(d_vec)
  p_fs  <- ncol(X_fs)

  a_sa <- ctrl$a_sigma_alpha_prior; b_sa <- ctrl$b_sigma_alpha_prior
  a_st <- ctrl$a_sigma_tau_prior;   b_st <- ctrl$b_sigma_tau_prior
  s2rp <- if (!is.null(ctrl$sigma2_rho_prior)) ctrl$sigma2_rho_prior else 10

  mu_gamma_prior <- if (is.null(ctrl$mu_gamma_prior)) rep(0, G)
  else as.numeric(ctrl$mu_gamma_prior)
  Sigma_gamma_prior_inv <- if (length(ctrl$Sigma_gamma_prior) == 1)
    diag(1 / ctrl$Sigma_gamma_prior, G)
  else diag(1 / as.numeric(ctrl$Sigma_gamma_prior), G)

  mu_delta_prior <- if (is.null(ctrl$mu_delta_prior)) rep(0, p_fs)
  else as.numeric(ctrl$mu_delta_prior)
  Sigma_delta_prior_inv <- if (is.null(ctrl$Sigma_delta_prior))
    diag(1 / 10, p_fs)
  else {
    sdp <- as.numeric(ctrl$Sigma_delta_prior)
    if (length(sdp) == 1) diag(1 / sdp, p_fs) else diag(1 / sdp)
  }

  XfsXfs       <- crossprod(X_fs)
  A_delta_base <- XfsXfs + Sigma_delta_prior_inv
  pre          <- .precompute_unit_stats(X_list, w_list, NULL, K, J0)
  XtWX         <- pre$XtWX
  N_total      <- sum(sapply(y_list, length))
  treat_idx    <- which(d_vec == 1L)
  ctrl_idx     <- which(d_vec == 0L)

  beta_bd <- numeric(K * J0)
  gamma   <- rep(0, G)
  rho     <- 0
  delta   <- if (!is.null(fs$delta0)) as.numeric(fs$delta0) else rep(0, p_fs)
  sigma2  <- 1
  tau     <- rep(1, K)
  z_star  <- ifelse(d_vec == 1, 0.5, -0.5)

  n_save         <- n_iter - burn_in
  beta_samples   <- matrix(0, n_save, K * J0)
  gamma_samples  <- matrix(0, n_save, G)
  rho_samples    <- numeric(n_save)
  delta_samples  <- matrix(0, n_save, p_fs)
  sigma2_samples <- numeric(n_save)
  tau_samples    <- matrix(0, n_save, K)

  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)

  for (iter in seq_len(n_iter)) {

    # BLOCK 1: z*
    eta    <- as.numeric(X_fs %*% delta)
    z_star <- .sample_z_star(eta, d_vec, treat_idx, ctrl_idx)

    # BLOCK 2: delta
    b_delta     <- as.numeric(crossprod(X_fs, z_star)) +
      as.numeric(Sigma_delta_prior_inv %*% mu_delta_prior)
    cA_d        <- chol(A_delta_base)
    A_delta_inv <- chol2inv(cA_d)
    delta <- rMVNormCovariance(1,
                               mu    = as.numeric(A_delta_inv %*% b_delta),
                               Sigma = A_delta_inv)

    # BLOCK 3: nu_hat split by unit
    nu_hat  <- z_star - as.numeric(X_fs %*% delta)
    nu_list <- lapply(row_idx_list, function(idx) nu_hat[idx])

    # BLOCK 4a: beta with moderator prior mean + rho partial-out
    z_gamma  <- as.numeric(Z %*% gamma)
    pm_list  <- lapply(seq_len(J0), function(j) {
      pm <- rep(0, K); pm[k_tr] <- z_gamma[j]; pm
    })
    y_tilde_list <- lapply(seq_len(J0), function(j)
      y_list[[j]] - rho * nu_list[[j]])
    XtWy    <- .compute_XtWy_units(X_list, w_list, y_tilde_list, K, J0)
    beta_bd <- .sample_beta_units(XtWX, XtWy, tau, sigma2, K, J0,
                                  prior_mean = pm_list)

    # BLOCK 4b: gamma
    bmat    <- matrix(beta_bd, ncol = K, byrow = TRUE)
    beta_tr <- bmat[, k_tr]
    V_gamma <- solve(crossprod(Z) / tau[k_tr]^2 + Sigma_gamma_prior_inv)
    rhs     <- as.numeric(crossprod(Z, beta_tr)) / tau[k_tr]^2 +
      as.numeric(Sigma_gamma_prior_inv %*% mu_gamma_prior)
    gamma   <- rMVNormCovariance(1, mu = as.numeric(V_gamma %*% rhs),
                                 Sigma = V_gamma)

    # BLOCK 4c: rho (scalar)
    nuWnu <- 0; nuWy2 <- 0
    z_gamma2 <- as.numeric(Z %*% gamma)
    for (j in seq_len(J0)) {
      bj   <- beta_bd[((j-1L)*K+1L):(j*K)]
      Xj   <- X_list[[j]]
      r2j  <- y_list[[j]] - if (K==1L) as.numeric(Xj)*bj
      else as.numeric(Xj %*% bj)
      wj   <- w_list[[j]]
      nuj  <- nu_list[[j]]
      nuWnu <- nuWnu + sum(wj * nuj^2)
      nuWy2 <- nuWy2 + sum(wj * nuj * r2j)
    }
    v_rho <- 1 / (nuWnu / sigma2 + 1 / s2rp)
    rho   <- stats::rnorm(1, mean = v_rho * nuWy2 / sigma2, sd = sqrt(v_rho))

    # BLOCK 5: sigma2
    wss    <- .compute_wss(X_list, w_list, y_list, beta_bd, K, J0,
                           nu_list = nu_list, rho = rho)
    sigma2 <- 1 / stats::rgamma(1, shape = a_sa + N_total / 2,
                                rate  = b_sa + wss / 2)

    # BLOCK 6: tau — moderator-adjusted for treatment covariate
    z_gamma3 <- as.numeric(Z %*% gamma)
    for (k in seq_len(K)) {
      centered <- bmat[, k] - if (k == k_tr) z_gamma3 else 0
      tau[k]   <- sqrt(1 / stats::rgamma(1,
                                         shape = a_st + J0 / 2,
                                         rate  = b_st + sum(centered^2) / 2))
    }

    if (iter > burn_in) {
      s <- iter - burn_in
      beta_samples[s, ]  <- beta_bd
      gamma_samples[s, ] <- as.numeric(gamma)
      rho_samples[s]     <- rho
      delta_samples[s, ] <- as.numeric(delta)
      sigma2_samples[s]  <- sigma2
      tau_samples[s, ]   <- tau
    }
    utils::setTxtProgressBar(pb, iter)
  }
  close(pb)
  colnames(gamma_samples) <- colnames(Z)
  colnames(delta_samples) <- colnames(X_fs)
  colnames(tau_samples)   <- gdata$cov$Xcols

  list(beta_samples = beta_samples, gamma_samples = gamma_samples,
       rho_samples = rho_samples, delta_samples = delta_samples,
       sigma2_samples = sigma2_samples, tau_samples = tau_samples)
}


# =============================================================================
#  Standalone Gibbs utilities (unchanged from original)
# =============================================================================

#' Gibbs draw for binomial probit with sigma^2
#'
#' Performs one Gibbs sampling iteration for a binomial probit model with
#' latent variable augmentation (Albert and Chib 1993), including an update
#' for the variance parameter \eqn{\sigma^2}.
#'
#' @param X         Matrix of predictors \eqn{\lbrack n \times p \rbrack}.
#' @param y         Binary response vector (0/1) of length \eqn{n}.
#' @param beta      Current regression coefficients, length \eqn{p}.
#' @param Sigma_prior Prior covariance matrix \eqn{\lbrack p \times p \rbrack}
#'   for \eqn{\beta}.
#' @param mu_prior  Prior mean vector, length \eqn{p}.
#' @param sigma2    Current value of \eqn{\sigma^2}. Default 1.
#' @param alpha_0   Shape parameter of the Inv-Gamma prior on \eqn{\sigma^2}.
#'   Default 2.
#' @param beta_0    Rate parameter of the Inv-Gamma prior on \eqn{\sigma^2}.
#'   Default 2.
#'
#' @return A list with elements \code{beta}, \code{sigma2}, and \code{z}.
#' @importFrom truncnorm rtruncnorm
#' @export
gibbs_binomial_probit <- function(X, y, beta, Sigma_prior,
                                  mu_prior, sigma2 = 1,
                                  alpha_0 = 2, beta_0 = 2) {
  n <- length(y); p <- ncol(X)
  z <- rep(0, n)
  for (i in seq_len(n)) {
    mu_i <- X[i, ] %*% beta
    z[i] <- if (y[i] == 1)
      truncnorm::rtruncnorm(1, a = 0,   mean = mu_i, sd = sqrt(sigma2))
    else
      truncnorm::rtruncnorm(1, b = 0,   mean = mu_i, sd = sqrt(sigma2))
  }
  Sp_inv <- solve(Sigma_prior)
  Sp     <- solve(t(X) %*% X / sigma2 + Sp_inv)
  mp     <- Sp %*% (t(X) %*% z / sigma2 + Sp_inv %*% mu_prior)
  bn     <- MASS::mvrnorm(1, mp, Sp)
  as_    <- alpha_0 + n / 2
  bs_    <- beta_0 + 0.5 * sum((z - X %*% bn)^2)
  list(beta = bn, sigma2 = 1 / stats::rgamma(1, shape = as_, rate = bs_), z = z)
}

#' Single Gibbs draw for OLS model
#'
#' Performs one iteration of a Gibbs sampler for a simple OLS regression
#' model under a flat prior on \eqn{\beta} and an improper Jeffreys prior
#' on \eqn{\sigma^2}.
#'
#' @param y      Response vector \eqn{\lbrack n \times 1 \rbrack}.
#' @param X      Design matrix \eqn{\lbrack n \times p \rbrack}, including intercept.
#' @param beta   Current coefficient vector, length \eqn{p}.
#' @param sigma2 Current error variance.
#'
#' @return A list with elements \code{betaZ}, \code{sigmaZ2},
#'   and \code{residuals}.
#' @importFrom MASS mvrnorm
#' @export
gibbs_sampler_one_draw <- function(y, X, beta, sigma2) {
  n <- nrow(X)
  XtX_inv <- solve(t(X) %*% X)
  bs  <- MASS::mvrnorm(1, XtX_inv %*% t(X) %*% y, sigma2 * XtX_inv)
  res <- y - X %*% bs
  list(betaZ = bs,
       sigmaZ2 = 1 / stats::rgamma(1, shape = n/2, rate = sum(res^2)/2),
       residuals = res)
}

#' Gibbs draw for multivariate regression with covariance option
#'
#' Performs a single Gibbs sampling draw for a multivariate regression model.
#' Each outcome has its own predictor matrix. Supports diagonal (independent)
#' or dense (correlated) residual covariance structures.
#'
#' @param y_list         List of outcome vectors, one per dependent variable.
#' @param X_list         List of predictor matrices, one per outcome.
#' @param sigma_param    For \code{"diagonal"}: list of initial
#'   \eqn{\sigma^2_j} values. For \code{"dense"}: initial covariance
#'   matrix \eqn{\Sigma}.
#' @param covariance_type \code{"diagonal"} (default) or \code{"dense"}.
#'
#' @return A list with elements \code{beta_samples}, \code{sigma2_samples},
#'   \code{residuals}, \code{predicted}, and \code{covariance}.
#' @export
gibbs_sampler_one_draw_cov <- function(y_list, X_list, sigma_param,
                                       covariance_type = "diagonal") {
  p <- length(y_list); n <- nrow(X_list[[1]])
  beta_samples <- vector("list", p)
  resid_mat <- pred_mat <- matrix(0, n, p)
  sigma2_list <- if (covariance_type == "diagonal") sigma_param else NULL
  Sigma_inv   <- if (covariance_type == "dense") solve(sigma_param) else NULL
  for (j in seq_len(p)) {
    Xj <- X_list[[j]]; yj <- y_list[[j]]
    Ai <- solve(t(Xj) %*% Xj)
    bv <- if (covariance_type == "diagonal") sigma2_list[[j]] * Ai
    else Sigma_inv[j,j] * Ai
    beta_samples[[j]] <- MASS::mvrnorm(1, Ai %*% t(Xj) %*% yj, bv)
    rj <- yj - Xj %*% beta_samples[[j]]
    resid_mat[,j] <- rj; pred_mat[,j] <- Xj %*% beta_samples[[j]]
    if (covariance_type == "diagonal")
      sigma2_list[j] <- 1 / stats::rgamma(1, shape=n/2, rate=sum(rj^2)/2)
  }
  ss <- if (covariance_type == "dense") {
    MCMCpack::riwish(n + p + 1, t(resid_mat) %*% resid_mat)
  } else sigma2_list
  list(beta_samples = beta_samples, sigma2_samples = ss,
       residuals = as.data.frame(resid_mat),
       predicted = as.data.frame(pred_mat),
       covariance = if (covariance_type == "dense") ss else NULL)
}

utils::globalVariables(c("estimate","term","lower","upper","sig",".grp"))
