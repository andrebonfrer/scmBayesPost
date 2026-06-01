# =============================================================================
#  scmBayesPost — Gibbs samplers (all variants) + utilities
#
#  Contents
#  --------
#  resolve_sampler_control()        internal helper: merge control with defaults
#  gibbs_postscm()                  dispatcher
#  gibbs_sampling_simple()          single outcome, no moderators, no first stage
#  gibbs_sampling_moderators()      single outcome, with second-stage moderators
#  gibbs_sampling_selection()                 single outcome, Bayesian probit first stage
#  gibbs_sampling_selection_moderators()  first stage + second-stage moderators
#  .build_nu_hat_stacked()          internal helper for selection sampler
#  gibbs_binomial_probit()          standalone Albert-Chib probit Gibbs draw
#  gibbs_sampler_one_draw()         standalone OLS Gibbs draw
#  gibbs_sampler_one_draw_cov()     standalone multivariate regression Gibbs draw
# =============================================================================


# -----------------------------------------------------------------------------
#' Resolve sampler control parameters
#'
#' Internal helper that merges user-supplied sampler controls with defaults.
#' Supports parameters for all sampler variants (simple, moderators, selection).
#'
#' @param control Optional named list of control parameters.
#' @param has_Z   Logical; whether a second-stage moderator model is present.
#'
#' @return A named list of resolved control parameters.
#' @keywords internal
resolve_sampler_control <- function(control = NULL, has_Z = FALSE) {

  defaults <- list(
    # --- second-stage gamma prior (used only when Z_block present)
    Sigma_gamma_prior   = 10,
    mu_gamma_prior      = NULL,
    # --- observation-noise variance prior (Inv-Gamma)
    a_sigma_alpha_prior = 2,
    b_sigma_alpha_prior = 2,
    # --- coefficient-dispersion prior (Inv-Gamma on tau^2)
    a_sigma_tau_prior   = 2,
    b_sigma_tau_prior   = 2,
    # --- first-stage probit prior (used only when first_stage = "selection_probit_bayes")
    mu_delta_prior      = NULL,   # default: zero vector (set inside sampler)
    Sigma_delta_prior   = NULL    # default: scalar 10  (set inside sampler)
  )

  if (is.null(control)) return(defaults)

  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0) {
    stop(
      sprintf("Unknown control parameter(s): %s",
              paste(unknown, collapse = ", ")),
      call. = FALSE
    )
  }

  utils::modifyList(defaults, control)
}


# -----------------------------------------------------------------------------
#' Gibbs sampler dispatcher for post-SCM models
#'
#' Dispatches to the appropriate internal Gibbs sampler depending on whether
#' the prepared data object includes a Bayesian probit first stage, moderators,
#' and whether the model is single- or multi-outcome.
#'
#' @param gdata   Output from \code{prepare_data_general()}.
#' @param n_iter  Integer. Total number of Gibbs iterations.
#' @param burn_in Integer. Number of initial iterations discarded.
#' @param control Optional list of sampler control parameters and priors.
#'   Supported entries:
#'   \describe{
#'     \item{\code{Sigma_gamma_prior}}{Prior variance for second-stage
#'       \eqn{\gamma} coefficients. Scalar or vector of length
#'       \code{ncol(Z_block)}.}
#'     \item{\code{mu_gamma_prior}}{Prior mean vector for \eqn{\gamma}.
#'       Defaults to zero.}
#'     \item{\code{a_sigma_alpha_prior}, \code{b_sigma_alpha_prior}}{
#'       Shape/rate for the observation-noise variance prior.}
#'     \item{\code{a_sigma_tau_prior}, \code{b_sigma_tau_prior}}{
#'       Shape/rate for the coefficient-dispersion prior.}
#'     \item{\code{mu_delta_prior}}{Prior mean for first-stage probit
#'       coefficients \eqn{\delta}. Numeric vector of length
#'       \code{ncol(X_fs)}. Defaults to zero.}
#'     \item{\code{Sigma_delta_prior}}{Prior variance for \eqn{\delta}.
#'       Scalar or vector of length \code{ncol(X_fs)}. Default 10.}
#'   }
#'
#' @return A list of posterior draws.
#' @export
gibbs_postscm <- function(gdata,
                          n_iter  = 1000,
                          burn_in = 500,
                          control = NULL) {

  if (is.null(gdata$cov))
    stop("gdata$cov missing", call. = FALSE)
  if (is.null(gdata$cov$intX))
    stop("Treatment coefficient index (intX) not set in gdata.", call. = FALSE)

  has_Z  <- !is.null(gdata$Z_block)
  has_fs <- !is.null(gdata$first_stage) &&
    isTRUE(gdata$cov$first_stage == "selection_probit_bayes")
  M      <- if (!is.null(gdata$cov$M)) gdata$cov$M else 1

  ctrl <- resolve_sampler_control(control = control, has_Z = has_Z)

  # ---------- single outcome, Bayesian selection first stage ----------

  if (M == 1 && has_fs && has_Z) {
    return(
      gibbs_sampling_selection_moderators(
        gdata   = gdata,
        n_iter  = n_iter,
        burn_in = burn_in,
        control = ctrl
      )
    )
  }

  if (M == 1 && has_fs && !has_Z) {
    return(
      gibbs_sampling_selection(
        gdata   = gdata,
        n_iter  = n_iter,
        burn_in = burn_in,
        control = ctrl
      )
    )
  }

  # ---------- single outcome, no first stage ----------
  if (M == 1 && !has_Z) {
    return(
      gibbs_sampling_simple(
        gdata   = gdata,
        n_iter  = n_iter,
        burn_in = burn_in,
        control = ctrl
      )
    )
  }

  if (M == 1 && has_Z) {
    return(
      gibbs_sampling_moderators(
        gdata   = gdata,
        n_iter  = n_iter,
        burn_in = burn_in,
        control = ctrl
      )
    )
  }

  # ---------- multi outcome ----------
  if (M > 1 && !has_Z)
    stop("Multi-outcome model without moderators not implemented yet.",
         call. = FALSE)
  if (M > 1 && has_Z)
    stop("Multi-outcome moderator model not implemented yet.", call. = FALSE)

  stop("Unhandled model configuration in gibbs_postscm().", call. = FALSE)
}


# -----------------------------------------------------------------------------
#' Gibbs sampler for post-SCM model without moderators
#'
#' Internal sampler used when no second-stage moderators (f.Z),
#' no instruments, and no selection equations are present.
#'
#' @keywords internal
gibbs_sampling_simple <- function(gdata,
                                  n_iter  = 1000,
                                  burn_in = 500,
                                  control = NULL) {

  ctrl <- resolve_sampler_control(control, has_Z = FALSE)

  y_block <- gdata$Y_block
  X_block <- gdata$X_block
  W_block <- gdata$W

  if (is.null(gdata$cov$intX))
    stop("gdata$cov$intX is NULL.", call. = FALSE)

  K  <- length(gdata$cov$Xcols)
  J0 <- gdata$cov$J0

  # priors
  a_sigma_alpha_prior <- ctrl$a_sigma_alpha_prior
  b_sigma_alpha_prior <- ctrl$b_sigma_alpha_prior
  a_sigma_tau_prior   <- ctrl$a_sigma_tau_prior
  b_sigma_tau_prior   <- ctrl$b_sigma_tau_prior

  # initial values
  beta   <- matrix(0, K * J0, 1)
  sigma2 <- 1
  tau    <- rep(1, K)

  # precompute
  XtW  <- Matrix::t(X_block) %*% W_block
  XtWX <- XtW %*% X_block
  XtWy <- XtW %*% y_block

  # storage
  n_save         <- n_iter - burn_in
  beta_samples   <- matrix(0, n_save, K * J0)
  sigma2_samples <- numeric(n_save)
  tau_samples    <- matrix(0, n_save, K)

  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)

  for (iter in seq_len(n_iter)) {

    # ---- beta update
    Sigma_beta_prior_inv <- Matrix::kronecker(diag(J0), diag(1 / tau))
    A     <- XtWX / sigma2 + Sigma_beta_prior_inv
    b     <- XtWy / sigma2
    cA    <- Matrix::chol(A)
    A_inv <- Matrix::chol2inv(cA)
    mu_beta <- A_inv %*% b
    beta  <- rMVNormCovariance(1, mu = as.numeric(mu_beta), Sigma = A_inv)

    # ---- sigma2 update
    residuals  <- y_block - X_block %*% beta
    alpha      <- a_sigma_alpha_prior + length(y_block) / 2
    beta_param <- b_sigma_alpha_prior +
      as.numeric(Matrix::t(residuals) %*% W_block %*% residuals) / 2
    sigma2 <- 1 / stats::rgamma(1, shape = alpha, rate = beta_param)

    # ---- tau update
    beta_matrix <- matrix(beta, ncol = K, byrow = TRUE)
    for (k in seq_len(K)) {
      tau[k] <- sqrt(1 / stats::rgamma(
        1,
        shape = a_sigma_tau_prior + J0 / 2,
        rate  = b_sigma_tau_prior + sum(beta_matrix[, k]^2) / 2
      ))
    }

    # ---- store
    if (iter > burn_in) {
      s <- iter - burn_in
      beta_samples[s, ]   <- as.numeric(beta)
      sigma2_samples[s]   <- sigma2
      tau_samples[s, ]    <- tau
    }

    utils::setTxtProgressBar(pb, iter)
  }

  close(pb)
  colnames(tau_samples) <- gdata$cov$Xcols

  list(
    beta_samples   = beta_samples,
    sigma2_samples = sigma2_samples,
    tau_samples    = tau_samples
  )
}


# -----------------------------------------------------------------------------
#' Gibbs sampler for single-outcome post-SCM model with moderators
#'
#' Internal sampler for the generalized scmBayesPost object when a second-stage
#' moderator equation is present, but no IV and no selection equation are used.
#'
#' @keywords internal
gibbs_sampling_moderators <- function(gdata,
                                      n_iter  = 1000,
                                      burn_in = 500,
                                      control = NULL) {

  ctrl <- resolve_sampler_control(control, has_Z = TRUE)

  y_block <- gdata$Y_block
  X_block <- gdata$X_block
  W_block <- gdata$W
  Z       <- gdata$Z_block

  if (is.null(Z))              stop("gdata$Z_block is NULL.", call. = FALSE)
  if (is.null(gdata$cov$intX)) stop("gdata$cov$intX is NULL.", call. = FALSE)

  K    <- length(gdata$cov$Xcols)
  J0   <- gdata$cov$J0
  G    <- ncol(Z)
  k_tr <- gdata$cov$intX

  # priors
  Sigma_gamma_prior <- ctrl$Sigma_gamma_prior

  mu_gamma_prior <- if (is.null(ctrl$mu_gamma_prior)) {
    rep(0, G)
  } else {
    mg <- as.numeric(ctrl$mu_gamma_prior)
    if (length(mg) != G)
      stop(sprintf("control$mu_gamma_prior must have length %d, got %d.",
                   G, length(mg)), call. = FALSE)
    mg
  }

  a_sigma_alpha_prior <- ctrl$a_sigma_alpha_prior
  b_sigma_alpha_prior <- ctrl$b_sigma_alpha_prior
  a_sigma_tau_prior   <- ctrl$a_sigma_tau_prior
  b_sigma_tau_prior   <- ctrl$b_sigma_tau_prior

  # gamma prior precision
  Sigma_gamma_prior_inv <- if (length(Sigma_gamma_prior) == 1) {
    diag(1 / Sigma_gamma_prior, G)
  } else {
    sg <- as.numeric(Sigma_gamma_prior)
    if (length(sg) != G)
      stop(sprintf("control$Sigma_gamma_prior must have length 1 or %d, got %d.",
                   G, length(sg)), call. = FALSE)
    diag(1 / sg, G)
  }

  # initial values
  beta   <- matrix(0, K * J0, 1)
  gamma  <- rep(0, G)
  sigma2 <- 1
  tau    <- rep(1, K)

  # precompute
  XtW  <- Matrix::t(X_block) %*% W_block
  XtWX <- XtW %*% X_block
  XtWy <- XtW %*% y_block

  # storage
  n_save         <- n_iter - burn_in
  beta_samples   <- matrix(0, n_save, K * J0)
  gamma_samples  <- matrix(0, n_save, G)
  sigma2_samples <- numeric(n_save)
  tau_samples    <- matrix(0, n_save, K)

  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)

  for (iter in seq_len(n_iter)) {

    # ---- beta update
    dZ     <- rep(0, K); dZ[k_tr] <- 1
    Z_star <- Matrix::kronecker(Z, dZ)
    Sigma_beta_prior_inv <- Matrix::kronecker(diag(J0), diag(1 / tau))
    q      <- Z_star %*% gamma
    A      <- XtWX / sigma2 + Sigma_beta_prior_inv
    b      <- XtWy / sigma2 + Sigma_beta_prior_inv %*% q
    cA     <- Matrix::chol(A)
    A_inv  <- Matrix::chol2inv(cA)
    mu_beta <- A_inv %*% b
    beta   <- rMVNormCovariance(1, mu = as.numeric(mu_beta), Sigma = A_inv)

    # ---- gamma update
    beta_matrix <- matrix(beta, ncol = K, byrow = TRUE)
    beta_tr     <- as.numeric(beta_matrix[, k_tr, drop = TRUE])
    V_gamma     <- solve(crossprod(Z) / tau[k_tr]^2 + Sigma_gamma_prior_inv)
    rhs         <- as.numeric(crossprod(Z, beta_tr)) / tau[k_tr]^2 +
      as.numeric(Sigma_gamma_prior_inv %*% mu_gamma_prior)
    m_gamma     <- V_gamma %*% rhs
    gamma       <- rMVNormCovariance(1, mu = as.numeric(m_gamma), Sigma = V_gamma)

    # ---- sigma2 update
    residuals  <- y_block - X_block %*% beta
    alpha      <- a_sigma_alpha_prior + length(y_block) / 2
    beta_param <- b_sigma_alpha_prior +
      as.numeric(Matrix::t(residuals) %*% W_block %*% residuals) / 2
    sigma2 <- 1 / stats::rgamma(1, shape = alpha, rate = beta_param)

    # ---- tau update
    for (k in seq_len(K)) {
      tau[k] <- sqrt(1 / stats::rgamma(
        1,
        shape = a_sigma_tau_prior + J0 / 2,
        rate  = b_sigma_tau_prior +
          sum((beta_matrix[, k] -
                 Z_star[((1:J0) - 1) * K + k, , drop = FALSE] %*% gamma)^2) / 2
      ))
    }

    # ---- store
    if (iter > burn_in) {
      s <- iter - burn_in
      beta_samples[s, ]  <- as.numeric(beta)
      gamma_samples[s, ] <- as.numeric(gamma)
      sigma2_samples[s]  <- sigma2
      tau_samples[s, ]   <- tau
    }

    utils::setTxtProgressBar(pb, iter)
  }

  close(pb)
  colnames(gamma_samples) <- colnames(Z)
  colnames(tau_samples)   <- gdata$cov$Xcols

  list(
    beta_samples   = beta_samples,
    gamma_samples  = gamma_samples,
    sigma2_samples = sigma2_samples,
    tau_samples    = tau_samples
  )
}


# -----------------------------------------------------------------------------
#' Gibbs sampler for post-SCM model with Bayesian probit first stage
#'
#' Internal sampler. Uses Albert-Chib (1993) data augmentation to sample
#' the latent utilities \eqn{z^*} and first-stage probit coefficients
#' \eqn{\delta} jointly with the outcome equation parameters.
#'
#' The joint model is:
#' \deqn{z^*_i = x_{fs,i}^\top \delta + \nu_i, \quad \nu_i \sim N(0,1)}
#' \deqn{d_i = \mathbf{1}[z^*_i > 0]}
#' \deqn{y_i = x_i^\top \beta + \rho\,\hat\nu_i + \varepsilon_i,
#'       \quad \varepsilon_i \sim N(0,\sigma^2)}
#'
#' \eqn{\rho} is a single global coefficient on \eqn{\hat\nu} and is
#' sampled separately from the block-diagonal \eqn{\beta} to avoid
#' dimension mismatches with the Kronecker prior structure.
#'
#' Invoked when \code{gdata$cov$first_stage == "selection_probit_bayes"}.
#'
#' @keywords internal
gibbs_sampling_selection <- function(gdata,
                                     n_iter  = 1000,
                                     burn_in = 500,
                                     control = NULL) {

  ctrl <- resolve_sampler_control(control, has_Z = FALSE)

  # ---- unpack outcome equation objects
  y_block <- gdata$Y_block
  X_block <- gdata$X_block
  W_block <- gdata$W
  K  <- length(gdata$cov$Xcols)   # covariates per unit (outcome equation)
  J0 <- gdata$cov$J0              # number of treated units

  # ---- unpack first-stage objects
  fs    <- gdata$first_stage
  X_fs  <- fs$X_fs
  d_vec <- as.numeric(fs$d)
  n_obs <- length(d_vec)
  p_fs  <- ncol(X_fs)

  if (is.null(X_fs))
    stop(paste0("gdata$first_stage$X_fs is NULL. ",
                "Re-run prepare_data_general with ",
                "first_stage = 'selection_probit_bayes'."),
         call. = FALSE)
  if (nrow(X_fs) != n_obs)
    stop("nrow(X_fs) != length(d). Check first-stage data alignment.",
         call. = FALSE)

  # ---- outcome equation priors
  a_sigma_alpha_prior <- ctrl$a_sigma_alpha_prior
  b_sigma_alpha_prior <- ctrl$b_sigma_alpha_prior
  a_sigma_tau_prior   <- ctrl$a_sigma_tau_prior
  b_sigma_tau_prior   <- ctrl$b_sigma_tau_prior

  # rho prior: N(0, sigma2_rho_prior) — separate from block-diagonal beta prior
  sigma2_rho_prior <- if (!is.null(ctrl$sigma2_rho_prior))
    as.numeric(ctrl$sigma2_rho_prior) else 10

  # ---- first-stage probit priors
  mu_delta_prior <- if (is.null(ctrl$mu_delta_prior)) {
    rep(0, p_fs)
  } else {
    md <- as.numeric(ctrl$mu_delta_prior)
    if (length(md) != p_fs)
      stop(sprintf("control$mu_delta_prior must have length %d, got %d.",
                   p_fs, length(md)), call. = FALSE)
    md
  }

  Sigma_delta_prior_inv <- if (is.null(ctrl$Sigma_delta_prior)) {
    diag(1 / 10, p_fs)
  } else {
    sdp <- as.numeric(ctrl$Sigma_delta_prior)
    if (length(sdp) == 1) {
      diag(1 / sdp, p_fs)
    } else {
      if (length(sdp) != p_fs)
        stop(sprintf(
          "control$Sigma_delta_prior must have length 1 or %d, got %d.",
          p_fs, length(sdp)), call. = FALSE)
      diag(1 / sdp)
    }
  }

  # ---- precompute first-stage cross-products (fixed across iterations)
  XfsXfs       <- crossprod(X_fs)
  A_delta_base <- XfsXfs + Sigma_delta_prior_inv

  # ---- dimension note -------------------------------------------------
  # X_block is [N_stacked x K*J0] block-diagonal (K covariates per unit).
  # nu_hat_stacked is [N_stacked x 1] — a SINGLE global column, NOT one
  # column per unit. The coefficient rho is therefore a scalar, sampled
  # separately from the K*J0 block-diagonal beta vector to avoid the
  # Kronecker dimension mismatch that arises when treating rho as part
  # of a (K+1)*J0 structure.
  # ---------------------------------------------------------------------

  # ---- initial values
  beta_bd <- matrix(0, K * J0, 1)   # unit-specific outcome coefficients
  rho     <- 0                       # scalar global selection correction
  delta   <- if (!is.null(fs$delta0)) as.numeric(fs$delta0) else rep(0, p_fs)
  sigma2  <- 1
  tau     <- rep(1, K)               # K dispersion params for beta_bd only

  # initialise z_star consistently with observed treatment
  z_star <- ifelse(d_vec == 1, 0.5, -0.5)

  # ---- storage
  n_save         <- n_iter - burn_in
  beta_samples   <- matrix(0, n_save, K * J0)
  rho_samples    <- numeric(n_save)
  delta_samples  <- matrix(0, n_save, p_fs)
  sigma2_samples <- numeric(n_save)
  tau_samples    <- matrix(0, n_save, K)

  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)

  for (iter in seq_len(n_iter)) {

    # ==================================================================
    # BLOCK 1: Sample z*_i  (Albert-Chib latent utility augmentation)
    # z*_i | delta, d_i = 1 ~ TN_{[0,  Inf)}(x_fs_i' delta, 1)
    # z*_i | delta, d_i = 0 ~ TN_{(-Inf, 0)}(x_fs_i' delta, 1)
    # ==================================================================
    eta <- as.numeric(X_fs %*% delta)

    z_star <- ifelse(
      d_vec == 1,
      truncnorm::rtruncnorm(n_obs, a =    0, b =  Inf, mean = eta, sd = 1),
      truncnorm::rtruncnorm(n_obs, a = -Inf, b =    0, mean = eta, sd = 1)
    )

    # ==================================================================
    # BLOCK 2: Sample delta | z*, X_fs
    # Posterior: N(A_delta_inv %*% b_delta, A_delta_inv)
    # A_delta = X_fs'X_fs + Sigma_delta_prior_inv   (precomputed)
    # b_delta = X_fs'z*   + Sigma_delta_prior_inv %*% mu_delta_prior
    # ==================================================================
    b_delta     <- as.numeric(crossprod(X_fs, z_star)) +
      as.numeric(Sigma_delta_prior_inv %*% mu_delta_prior)
    cA_delta    <- chol(A_delta_base)
    A_delta_inv <- chol2inv(cA_delta)
    mu_delta    <- as.numeric(A_delta_inv %*% b_delta)
    delta       <- rMVNormCovariance(1, mu = mu_delta, Sigma = A_delta_inv)

    # ==================================================================
    # BLOCK 3: Build nu_hat (stacked to match X_block row ordering)
    # nu_hat_i = z*_i - x_fs_i' delta
    # ==================================================================
    nu_hat <- z_star - as.numeric(X_fs %*% delta)

    nu_hat_stacked <- .build_nu_hat_stacked(
      nu_hat   = nu_hat,
      X_idlist = gdata$X_idlist,
      dta_id   = gdata$dtaidx[["id"]],
      dta_wID  = gdata$dtaidx[["wID"]]
    )

    nu_col <- Matrix::Matrix(nu_hat_stacked, ncol = 1)  # [N_stacked x 1]

    # ==================================================================
    # BLOCK 4a: Sample beta_bd | y, X_block, rho, nu_col, sigma2, tau
    #
    # Partial out rho: y_tilde = y - nu_col * rho
    # beta_bd | y_tilde ~ N(A_bd_inv b_bd, A_bd_inv)
    # with block-diagonal prior precision diag(J0) x diag(1/tau)
    # ==================================================================
    y_tilde <- y_block - nu_col * rho

    Sigma_beta_prior_inv <- Matrix::kronecker(diag(J0), diag(1 / tau))

    XtW_bd  <- Matrix::t(X_block) %*% W_block
    XtWX_bd <- XtW_bd %*% X_block
    XtWy_bd <- XtW_bd %*% y_tilde

    A_bd  <- XtWX_bd / sigma2 + Sigma_beta_prior_inv
    b_bd  <- XtWy_bd / sigma2

    cA_bd    <- Matrix::chol(A_bd)
    A_bd_inv <- Matrix::chol2inv(cA_bd)
    mu_bd    <- A_bd_inv %*% b_bd

    beta_bd <- rMVNormCovariance(1,
                                 mu    = as.numeric(mu_bd),
                                 Sigma = A_bd_inv)

    # ==================================================================
    # BLOCK 4b: Sample rho (scalar) | y, X_block, beta_bd, nu_col, sigma2
    #
    # Partial out beta_bd: y_tilde2 = y - X_block %*% beta_bd
    # rho | y_tilde2 ~ N(m_rho, v_rho)
    # v_rho = 1 / (nu'Wnu/sigma2 + 1/sigma2_rho_prior)
    # m_rho = v_rho * nu'W y_tilde2 / sigma2
    # ==================================================================
    y_tilde2 <- y_block - X_block %*% beta_bd

    nuWnu <- as.numeric(Matrix::t(nu_col) %*% W_block %*% nu_col)
    nuWy2 <- as.numeric(Matrix::t(nu_col) %*% W_block %*% y_tilde2)

    v_rho <- 1 / (nuWnu / sigma2 + 1 / sigma2_rho_prior)
    m_rho <- v_rho * nuWy2 / sigma2
    rho   <- stats::rnorm(1, mean = m_rho, sd = sqrt(v_rho))

    # ==================================================================
    # BLOCK 5: Sample sigma2 | y, X_block, beta_bd, rho, nu_col
    # ==================================================================
    residuals <- y_block - X_block %*% beta_bd - nu_col * rho
    alpha_s   <- a_sigma_alpha_prior + length(y_block) / 2
    beta_s    <- b_sigma_alpha_prior +
      as.numeric(Matrix::t(residuals) %*% W_block %*% residuals) / 2
    sigma2 <- 1 / stats::rgamma(1, shape = alpha_s, rate = beta_s)

    # ==================================================================
    # BLOCK 6: Sample tau_k | beta_bd  (K dispersion params, block-diag only)
    # ==================================================================
    beta_matrix <- matrix(beta_bd, ncol = K, byrow = TRUE)
    for (k in seq_len(K)) {
      tau[k] <- sqrt(1 / stats::rgamma(
        1,
        shape = a_sigma_tau_prior + J0 / 2,
        rate  = b_sigma_tau_prior + sum(beta_matrix[, k]^2) / 2
      ))
    }

    # ---- store
    if (iter > burn_in) {
      s <- iter - burn_in
      beta_samples[s, ]  <- as.numeric(beta_bd)
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

  list(
    beta_samples   = beta_samples,
    rho_samples    = rho_samples,
    delta_samples  = delta_samples,
    sigma2_samples = sigma2_samples,
    tau_samples    = tau_samples
  )
}


# -----------------------------------------------------------------------------
#' Gibbs sampler for post-SCM model with Bayesian probit first stage
#' AND second-stage moderators
#'
#' Combines Albert-Chib (1993) latent utility augmentation for the selection
#' first stage with the hierarchical moderator equation for treatment effect
#' heterogeneity. rho (selection correction) is sampled as a scalar separate
#' from the block-diagonal beta, avoiding Kronecker dimension mismatches.
#'
#' @keywords internal
gibbs_sampling_selection_moderators <- function(gdata,
                                                n_iter  = 1000,
                                                burn_in = 500,
                                                control = NULL) {

  ctrl <- resolve_sampler_control(control, has_Z = TRUE)

  # ---- unpack outcome equation objects
  y_block <- gdata$Y_block
  X_block <- gdata$X_block
  W_block <- gdata$W
  Z       <- gdata$Z_block

  if (is.null(Z))
    stop("gdata$Z_block is NULL.", call. = FALSE)

  K    <- length(gdata$cov$Xcols)
  J0   <- gdata$cov$J0
  G    <- ncol(Z)
  k_tr <- gdata$cov$intX

  # ---- unpack first-stage objects
  fs    <- gdata$first_stage
  X_fs  <- fs$X_fs
  d_vec <- as.numeric(fs$d)
  n_obs <- length(d_vec)
  p_fs  <- ncol(X_fs)

  if (is.null(X_fs))
    stop("gdata$first_stage$X_fs is NULL.", call. = FALSE)
  if (nrow(X_fs) != n_obs)
    stop("nrow(X_fs) != length(d).", call. = FALSE)

  # ---- outcome equation priors
  a_sigma_alpha_prior <- ctrl$a_sigma_alpha_prior
  b_sigma_alpha_prior <- ctrl$b_sigma_alpha_prior
  a_sigma_tau_prior   <- ctrl$a_sigma_tau_prior
  b_sigma_tau_prior   <- ctrl$b_sigma_tau_prior

  sigma2_rho_prior <- if (!is.null(ctrl$sigma2_rho_prior))
    as.numeric(ctrl$sigma2_rho_prior) else 10

  # ---- moderator (gamma) priors
  Sigma_gamma_prior <- ctrl$Sigma_gamma_prior
  mu_gamma_prior <- if (is.null(ctrl$mu_gamma_prior)) {
    rep(0, G)
  } else {
    mg <- as.numeric(ctrl$mu_gamma_prior)
    if (length(mg) != G)
      stop(sprintf("mu_gamma_prior must have length %d, got %d.", G, length(mg)),
           call. = FALSE)
    mg
  }
  Sigma_gamma_prior_inv <- if (length(Sigma_gamma_prior) == 1) {
    diag(1 / Sigma_gamma_prior, G)
  } else {
    sg <- as.numeric(Sigma_gamma_prior)
    if (length(sg) != G)
      stop(sprintf("Sigma_gamma_prior must have length 1 or %d, got %d.",
                   G, length(sg)), call. = FALSE)
    diag(1 / sg, G)
  }

  # ---- first-stage probit priors
  mu_delta_prior <- if (is.null(ctrl$mu_delta_prior)) {
    rep(0, p_fs)
  } else {
    md <- as.numeric(ctrl$mu_delta_prior)
    if (length(md) != p_fs)
      stop(sprintf("mu_delta_prior must have length %d, got %d.", p_fs, length(md)),
           call. = FALSE)
    md
  }
  Sigma_delta_prior_inv <- if (is.null(ctrl$Sigma_delta_prior)) {
    diag(1 / 10, p_fs)
  } else {
    sdp <- as.numeric(ctrl$Sigma_delta_prior)
    if (length(sdp) == 1) diag(1 / sdp, p_fs)
    else {
      if (length(sdp) != p_fs)
        stop(sprintf("Sigma_delta_prior must have length 1 or %d, got %d.",
                     p_fs, length(sdp)), call. = FALSE)
      diag(1 / sdp)
    }
  }

  # ---- precompute fixed cross-products
  XfsXfs       <- crossprod(X_fs)
  A_delta_base <- XfsXfs + Sigma_delta_prior_inv

  # ---- initial values
  beta_bd <- matrix(0, K * J0, 1)
  gamma   <- rep(0, G)
  rho     <- 0
  delta   <- if (!is.null(fs$delta0)) as.numeric(fs$delta0) else rep(0, p_fs)
  sigma2  <- 1
  tau     <- rep(1, K)

  z_star <- ifelse(d_vec == 1, 0.5, -0.5)

  # ---- storage
  n_save         <- n_iter - burn_in
  beta_samples   <- matrix(0, n_save, K * J0)
  gamma_samples  <- matrix(0, n_save, G)
  rho_samples    <- numeric(n_save)
  delta_samples  <- matrix(0, n_save, p_fs)
  sigma2_samples <- numeric(n_save)
  tau_samples    <- matrix(0, n_save, K)

  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)

  for (iter in seq_len(n_iter)) {

    # ==================================================================
    # BLOCK 1: Sample z*_i  (Albert-Chib)
    # ==================================================================
    eta <- as.numeric(X_fs %*% delta)
    z_star <- ifelse(
      d_vec == 1,
      truncnorm::rtruncnorm(n_obs, a =    0, b =  Inf, mean = eta, sd = 1),
      truncnorm::rtruncnorm(n_obs, a = -Inf, b =    0, mean = eta, sd = 1)
    )

    # ==================================================================
    # BLOCK 2: Sample delta | z*, X_fs
    # ==================================================================
    b_delta     <- as.numeric(crossprod(X_fs, z_star)) +
      as.numeric(Sigma_delta_prior_inv %*% mu_delta_prior)
    cA_delta    <- chol(A_delta_base)
    A_delta_inv <- chol2inv(cA_delta)
    mu_delta    <- as.numeric(A_delta_inv %*% b_delta)
    delta       <- rMVNormCovariance(1, mu = mu_delta, Sigma = A_delta_inv)

    # ==================================================================
    # BLOCK 3: Build nu_hat and stack to match X_block rows
    # ==================================================================
    nu_hat <- z_star - as.numeric(X_fs %*% delta)
    nu_hat_stacked <- .build_nu_hat_stacked(
      nu_hat   = nu_hat,
      X_idlist = gdata$X_idlist,
      dta_id   = gdata$dtaidx[["id"]],
      dta_wID  = gdata$dtaidx[["wID"]]
    )
    nu_col <- Matrix::Matrix(nu_hat_stacked, ncol = 1)

    # ==================================================================
    # BLOCK 4a: Sample beta_bd | y, rho, gamma, sigma2, tau
    # Moderator prior mean q enters via Sigma_beta_prior_inv %*% q
    # y_tilde = y - nu_col * rho  (partial out selection correction)
    # ==================================================================
    dZ     <- rep(0, K); dZ[k_tr] <- 1
    Z_star <- Matrix::kronecker(Z, dZ)         # [K*J0 x G] selector

    Sigma_beta_prior_inv <- Matrix::kronecker(diag(J0), diag(1 / tau))
    q       <- Z_star %*% gamma               # moderator-implied prior mean

    y_tilde <- y_block - nu_col * rho

    XtW_bd  <- Matrix::t(X_block) %*% W_block
    XtWX_bd <- XtW_bd %*% X_block
    XtWy_bd <- XtW_bd %*% y_tilde

    A_bd    <- XtWX_bd / sigma2 + Sigma_beta_prior_inv
    b_bd    <- XtWy_bd / sigma2 + Sigma_beta_prior_inv %*% q

    cA_bd    <- Matrix::chol(A_bd)
    A_bd_inv <- Matrix::chol2inv(cA_bd)
    mu_bd    <- A_bd_inv %*% b_bd
    beta_bd  <- rMVNormCovariance(1, mu = as.numeric(mu_bd), Sigma = A_bd_inv)

    # ==================================================================
    # BLOCK 4b: Sample gamma | beta_bd, tau, Z
    # Identical to gibbs_sampling_moderators gamma block
    # ==================================================================
    beta_matrix <- matrix(beta_bd, ncol = K, byrow = TRUE)
    beta_tr     <- as.numeric(beta_matrix[, k_tr, drop = TRUE])

    V_gamma <- solve(crossprod(Z) / tau[k_tr]^2 + Sigma_gamma_prior_inv)
    rhs     <- as.numeric(crossprod(Z, beta_tr)) / tau[k_tr]^2 +
      as.numeric(Sigma_gamma_prior_inv %*% mu_gamma_prior)
    m_gamma <- V_gamma %*% rhs
    gamma   <- rMVNormCovariance(1, mu = as.numeric(m_gamma), Sigma = V_gamma)

    # ==================================================================
    # BLOCK 4c: Sample rho (scalar) | y, beta_bd, nu_col, sigma2
    # ==================================================================
    y_tilde2 <- y_block - X_block %*% beta_bd
    nuWnu    <- as.numeric(Matrix::t(nu_col) %*% W_block %*% nu_col)
    nuWy2    <- as.numeric(Matrix::t(nu_col) %*% W_block %*% y_tilde2)
    v_rho    <- 1 / (nuWnu / sigma2 + 1 / sigma2_rho_prior)
    m_rho    <- v_rho * nuWy2 / sigma2
    rho      <- stats::rnorm(1, mean = m_rho, sd = sqrt(v_rho))

    # ==================================================================
    # BLOCK 5: Sample sigma2
    # ==================================================================
    residuals <- y_block - X_block %*% beta_bd - nu_col * rho
    alpha_s   <- a_sigma_alpha_prior + length(y_block) / 2
    beta_s    <- b_sigma_alpha_prior +
      as.numeric(Matrix::t(residuals) %*% W_block %*% residuals) / 2
    sigma2 <- 1 / stats::rgamma(1, shape = alpha_s, rate = beta_s)

    # ==================================================================
    # BLOCK 6: Sample tau_k | beta_bd, gamma  (moderator-adjusted residuals)
    # ==================================================================
    for (k in seq_len(K)) {
      tau[k] <- sqrt(1 / stats::rgamma(
        1,
        shape = a_sigma_tau_prior + J0 / 2,
        rate  = b_sigma_tau_prior +
          sum((beta_matrix[, k] -
                 Z_star[((1:J0) - 1) * K + k, , drop = FALSE] %*% gamma)^2) / 2
      ))
    }

    # ---- store
    if (iter > burn_in) {
      s <- iter - burn_in
      beta_samples[s, ]  <- as.numeric(beta_bd)
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

  list(
    beta_samples   = beta_samples,
    gamma_samples  = gamma_samples,
    rho_samples    = rho_samples,
    delta_samples  = delta_samples,
    sigma2_samples = sigma2_samples,
    tau_samples    = tau_samples
  )
}


# -----------------------------------------------------------------------------
#' Build stacked nu_hat vector aligned with X_block row ordering
#'
#' X_block is built by stacking J0 pseudo-panels. Each pseudo-panel contains
#' rows for (treated unit + all controls) in the order defined by pseudo_ids.
#' This helper replicates nu_hat (which is observation-level, length n_obs)
#' into the same stacked ordering so it can be appended as a column to X_block.
#'
#' @param nu_hat   Numeric vector of length n_obs.
#' @param X_idlist data.table with columns id and wID, nrow = nrow(X_block),
#'   as returned by prepare_data_general.
#' @param dta_id   Character vector of unit ids aligned with nu_hat
#'   (i.e. as.character(dta\[\[id_col\]\]) in original row order).
#' @param dta_wID  Numeric/integer vector of time ids aligned with nu_hat.
#'
#' @return Numeric vector of length nrow(X_block).
#' @keywords internal
.build_nu_hat_stacked <- function(nu_hat, X_idlist, dta_id, dta_wID) {

  n_orig <- length(dta_id)
  if (length(nu_hat) != n_orig)
    stop(sprintf(
      ".build_nu_hat_stacked: length(nu_hat) = %d != length(dta_id) = %d",
      length(nu_hat), n_orig
    ), call. = FALSE)
  if (length(dta_wID) != n_orig)
    stop(sprintf(
      ".build_nu_hat_stacked: length(dta_wID) = %d != length(dta_id) = %d",
      length(dta_wID), n_orig
    ), call. = FALSE)

  lookup <- data.table::data.table(
    id     = as.character(dta_id),
    wID    = dta_wID,
    nu_hat = nu_hat
  )
  data.table::setkeyv(lookup, c("id", "wID"))

  query <- data.table::copy(X_idlist)
  query[, id := as.character(id)]
  data.table::setkeyv(query, c("id", "wID"))

  merged <- lookup[query, on = c("id", "wID")]

  if (anyNA(merged$nu_hat))
    stop(paste0(
      ".build_nu_hat_stacked: ", sum(is.na(merged$nu_hat)),
      " unmatched rows. Check that X_idlist and dta_id share the same ",
      "(id, wID) keys."
    ), call. = FALSE)

  merged$nu_hat
}


# =============================================================================
#  Standalone Gibbs utilities
#  (kept as-is from original; no dependency on prepare_data_general)
# =============================================================================

#' Gibbs Sampling for a Binomial Probit Model with sigma^2
#'
#' Performs one Gibbs sampling iteration for a binomial probit model with
#' latent variable augmentation (Albert & Chib 1993). Includes an update
#' step for the variance parameter \eqn{\sigma^2}.
#'
#' The model is:
#' \deqn{y_i = 1 \text{ if } z_i > 0, \quad z_i = X_i\beta + \varepsilon_i,
#'       \quad \varepsilon_i \sim N(0,\sigma^2)}
#'
#' @param X        Matrix of predictors \eqn{[n \times p]}.
#' @param y        Binary response vector \eqn{(0/1)} of length \eqn{n}.
#' @param beta     Current regression coefficients, length \eqn{p}.
#' @param Sigma_prior Prior covariance matrix \eqn{[p \times p]} for \eqn{\beta}.
#' @param mu_prior Prior mean vector, length \eqn{p}.
#' @param sigma2   Current value of \eqn{\sigma^2}. Default 1.
#' @param alpha_0  Shape parameter of the Inv-Gamma prior on \eqn{\sigma^2}. Default 2.
#' @param beta_0   Rate parameter of the Inv-Gamma prior on \eqn{\sigma^2}. Default 2.
#'
#' @return A list with elements \code{beta}, \code{sigma2}, and \code{z}.
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(100 * 3), 100, 3)
#' y <- rbinom(100, 1, 0.5)
#' result <- gibbs_binomial_probit(X, y, rep(0, 3), diag(3), rep(0, 3))
#' print(result$beta)
#'
#' @importFrom truncnorm rtruncnorm
#' @export
gibbs_binomial_probit <- function(X, y, beta, Sigma_prior,
                                  mu_prior, sigma2 = 1,
                                  alpha_0 = 2, beta_0 = 2) {

  n <- length(y)
  p <- ncol(X)

  # Step 1: sample latent z from truncated normal
  z <- rep(0, n)
  for (i in seq_len(n)) {
    mu_i <- X[i, ] %*% beta
    if (y[i] == 1) {
      z[i] <- truncnorm::rtruncnorm(1, a = 0,    mean = mu_i, sd = sqrt(sigma2))
    } else {
      z[i] <- truncnorm::rtruncnorm(1, b = 0,    mean = mu_i, sd = sqrt(sigma2))
    }
  }

  # Step 2: sample beta
  Sigma_prior_inv  <- solve(Sigma_prior)
  Sigma_posterior  <- solve(t(X) %*% X / sigma2 + Sigma_prior_inv)
  mu_posterior     <- Sigma_posterior %*%
    (t(X) %*% z / sigma2 + Sigma_prior_inv %*% mu_prior)
  beta_new <- MASS::mvrnorm(1, mu_posterior, Sigma_posterior)

  # Step 3: sample sigma2
  alpha_sigma <- alpha_0 + n / 2
  beta_sigma  <- beta_0 + 0.5 * sum((z - X %*% beta_new)^2)
  sigma2_new  <- 1 / stats::rgamma(1, shape = alpha_sigma, rate = beta_sigma)

  list(beta = beta_new, sigma2 = sigma2_new, z = z)
}


#' Perform one iteration of the Gibbs sampler for an OLS model
#'
#' Samples regression coefficients \eqn{\beta} and error variance
#' \eqn{\sigma^2} from their full conditionals under a flat prior on
#' \eqn{\beta} and an improper Jeffreys prior on \eqn{\sigma^2}.
#'
#' @param y      Response vector \eqn{[n \times 1]}.
#' @param X      Design matrix \eqn{[n \times p]}, including intercept.
#' @param beta   Current coefficient vector, length \eqn{p}.
#' @param sigma2 Current error variance.
#'
#' @return A list with elements \code{betaZ}, \code{sigmaZ2}, \code{residuals}.
#'
#' @examples
#' set.seed(123)
#' X <- cbind(1, rnorm(100))
#' y <- X %*% c(2, 3) + rnorm(100)
#' out <- gibbs_sampler_one_draw(y, X, rep(0, 2), 1)
#' print(out$betaZ)
#'
#' @importFrom MASS mvrnorm
#' @export
gibbs_sampler_one_draw <- function(y, X, beta, sigma2) {

  n <- nrow(X)
  p <- ncol(X)

  XtX_inv   <- solve(t(X) %*% X)
  beta_mean <- XtX_inv %*% t(X) %*% y
  beta_var  <- sigma2 * XtX_inv
  beta_sample <- MASS::mvrnorm(1, beta_mean, beta_var)

  residuals    <- y - X %*% beta_sample
  alpha        <- n / 2
  beta_param   <- sum(residuals^2) / 2
  sigma2_sample <- 1 / stats::rgamma(1, shape = alpha, rate = beta_param)

  list(betaZ      = beta_sample,
       sigmaZ2    = sigma2_sample,
       residuals  = residuals)
}


#' Gibbs Sampler for Multivariate Regression with Covariance Option
#'
#' Performs a single Gibbs sampling draw for a multivariate regression model.
#' Each outcome has its own predictor matrix. Supports diagonal (independent)
#' or dense (correlated) residual covariance structures.
#'
#' @param y_list         List of outcome vectors, one per dependent variable.
#' @param X_list         List of predictor matrices, one per outcome.
#' @param sigma_param    For \code{"diagonal"}: list of initial \eqn{\sigma^2_j}
#'   values. For \code{"dense"}: initial covariance matrix \eqn{\Sigma}.
#' @param covariance_type \code{"diagonal"} or \code{"dense"}.
#'
#' @return A list with elements \code{beta_samples}, \code{sigma2_samples},
#'   \code{residuals}, \code{predicted}, and \code{covariance}.
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' X1 <- cbind(1, rnorm(n)); X2 <- cbind(1, rnorm(n), rnorm(n))
#' y1 <- X1 %*% c(1, 0.5) + rnorm(n)
#' y2 <- X2 %*% c(1, -0.3, 0.2) + rnorm(n)
#' out <- gibbs_sampler_one_draw_cov(list(y1, y2), list(X1, X2),
#'                                   list(1, 1), "diagonal")
#' print(out$beta_samples)
#'
#' @export
gibbs_sampler_one_draw_cov <- function(y_list, X_list, sigma_param,
                                       covariance_type = "diagonal") {

  p <- length(y_list)
  n <- nrow(X_list[[1]])

  beta_samples     <- vector("list", p)
  residuals_matrix <- matrix(0, nrow = n, ncol = p)
  predicted_matrix <- matrix(0, nrow = n, ncol = p)

  if (covariance_type == "diagonal") {
    sigma2_list <- sigma_param
  } else if (covariance_type == "dense") {
    Sigma     <- sigma_param
    Sigma_inv <- solve(Sigma)
  } else {
    stop("Invalid covariance_type. Use 'diagonal' or 'dense'.")
  }

  for (j in seq_len(p)) {
    X_j       <- X_list[[j]]
    y_j       <- y_list[[j]]
    XtX_inv_j <- solve(t(X_j) %*% X_j)
    beta_mean_j <- XtX_inv_j %*% t(X_j) %*% y_j

    if (covariance_type == "diagonal") {
      beta_var_j <- sigma2_list[[j]] * XtX_inv_j
    } else {
      beta_var_j <- Sigma_inv[j, j] * XtX_inv_j
    }

    beta_samples[[j]]    <- MASS::mvrnorm(1, beta_mean_j, beta_var_j)
    residuals_j          <- y_j - X_j %*% beta_samples[[j]]
    residuals_matrix[, j] <- residuals_j
    predicted_matrix[, j] <- X_j %*% beta_samples[[j]]

    if (covariance_type == "diagonal") {
      alpha        <- n / 2
      beta_param   <- sum(residuals_j^2) / 2
      sigma2_list[j] <- 1 / stats::rgamma(1, shape = alpha, rate = beta_param)
    }
  }

  if (covariance_type == "dense") {
    nu           <- n + p + 1
    S            <- t(residuals_matrix) %*% residuals_matrix
    sigma_sample <- MCMCpack::riwish(nu, S)
  } else {
    sigma_sample <- sigma2_list
  }

  list(
    beta_samples   = beta_samples,
    sigma2_samples = sigma_sample,
    residuals      = as.data.frame(residuals_matrix),
    predicted      = as.data.frame(predicted_matrix),
    covariance     = if (covariance_type == "dense") sigma_sample else NULL
  )
}
