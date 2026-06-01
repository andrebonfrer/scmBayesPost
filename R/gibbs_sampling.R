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
#' Sample beta block-diagonally (fast path for large J0)
#'
#' When X_block is block-diagonal with J0 blocks of size \[T_j x K\] and the
#' prior precision is block-diagonal diag(J0) x diag(1/tau), the posterior
#' precision A_bd = XtWX_bd/sigma2 + Sigma_beta_prior_inv is also block-
#' diagonal. Factorising each \[K x K\] block independently reduces cost from
#' O((K*J0)^3) to O(J0 * K^3). For K=1 each block is a scalar — the Cholesky
#' collapses to a square root and the MVN draw to a single normal draw.
#'
#' @param XtWX_blocks List of J0 matrices, each \[K x K\]: the j-th block of
#'   X_block'W X_block.
#' @param XtWy_blocks List of J0 vectors, each length K: the j-th block of
#'   X_block'W y_tilde.
#' @param tau         Numeric vector length K: current dispersion SDs.
#' @param sigma2      Scalar: current observation noise variance.
#' @param K           Integer: covariates per unit.
#' @param prior_mean  Optional list of J0 vectors (prior mean per unit block).
#'   Default NULL gives zero prior mean.
#'
#' @return Numeric vector of length K*J0: sampled beta in stacked unit order.
#' @keywords internal
.sample_beta_bd_blockdiag <- function(XtWX_blocks,
                                      XtWy_blocks,
                                      tau,
                                      sigma2,
                                      K,
                                      prior_mean = NULL) {

  J0          <- length(XtWX_blocks)
  tau_inv     <- 1 / tau^2          # precision (1/tau^2 not 1/tau)
  beta_out    <- numeric(K * J0)

  for (j in seq_len(J0)) {
    XtWX_j <- XtWX_blocks[[j]]     # [K x K]
    XtWy_j <- XtWy_blocks[[j]]     # [K]
    pm_j   <- if (is.null(prior_mean)) rep(0, K) else prior_mean[[j]]

    # Posterior precision: A_j = XtWX_j/sigma2 + diag(tau_inv)
    A_j <- XtWX_j / sigma2 + diag(tau_inv, K)

    if (K == 1L) {
      # Scalar fast path: avoid Cholesky entirely
      v_j   <- 1 / as.numeric(A_j)
      m_j   <- v_j * (as.numeric(XtWy_j) / sigma2 +
                        tau_inv * pm_j)
      beta_out[j] <- stats::rnorm(1, mean = m_j, sd = sqrt(v_j))
    } else {
      # General K > 1: Cholesky of [K x K] block
      b_j   <- XtWy_j / sigma2 + tau_inv * pm_j
      cA_j  <- chol(A_j)
      Ai_j  <- chol2inv(cA_j)
      mu_j  <- as.numeric(Ai_j %*% b_j)
      # Draw from N(mu_j, Ai_j) via Cholesky backsolve
      z     <- stats::rnorm(K)
      beta_out[((j - 1L) * K + 1L):(j * K)] <-
        mu_j + as.numeric(backsolve(cA_j, z))
    }
  }

  beta_out
}


# -----------------------------------------------------------------------------
#' Precompute per-unit XtWX and XtWy blocks from block-diagonal X_block
#'
#' Extracts the J0 diagonal blocks of X_block'W X_block and X_block'W y,
#' each of size \[K x K\] and \[K x 1\] respectively. Called once before the
#' Gibbs loop; XtWy_blocks is recomputed per iteration (since y_tilde
#' changes) but XtWX_blocks is fixed.
#'
#' @param X_block sparse \[N_stacked x K*J0\] block-diagonal Matrix.
#' @param W_block sparse diagonal \[N_stacked x N_stacked\] weight Matrix.
#' @param K       Integer. Covariates per unit.
#' @param J0      Integer. Number of treated units.
#'
#' @return List with elements \code{XtWX_blocks} and \code{row_ranges},
#'   where \code{row_ranges[[j]]} gives the row indices of block j in
#'   X_block (needed to extract XtWy_blocks per iteration).
#' @keywords internal
.precompute_XtWX_blocks <- function(X_block, W_block, K, J0) {

  # X_block column ranges: block j occupies cols (j-1)*K+1 : j*K
  XtWX_blocks <- vector("list", J0)

  # For block-diagonal X_block, block j's rows are those where only
  # columns (j-1)*K+1:j*K are non-zero. We infer row ranges from the
  # structure: each unit block has the same number of rows T_j.
  # Total rows N = sum(T_j); if balanced T_j = N/J0.
  N_stacked <- nrow(X_block)

  # Extract column ranges and corresponding row blocks
  # For a proper block-diagonal matrix, column block j is non-zero only
  # in row block j. We use Matrix::which to find non-zero rows per col block.
  col_starts <- (seq_len(J0) - 1L) * K + 1L

  row_ranges <- vector("list", J0)

  for (j in seq_len(J0)) {
    cols_j <- col_starts[j]:(col_starts[j] + K - 1L)
    # Non-zero rows for this column block
    rows_j <- unique(Matrix::which(X_block[, cols_j, drop = FALSE] != 0,
                                   arr.ind = TRUE)[, "row"])
    rows_j <- sort(rows_j)
    row_ranges[[j]] <- rows_j

    Xj <- as.matrix(X_block[rows_j, cols_j, drop = FALSE])
    Wj <- Matrix::diag(W_block)[rows_j]
    XtWX_blocks[[j]] <- crossprod(Xj, Wj * Xj)   # [K x K]
  }

  list(XtWX_blocks = XtWX_blocks, row_ranges = row_ranges)
}


# -----------------------------------------------------------------------------
#' Compute per-unit XtWy blocks (called each iteration)
#'
#' @param X_block   sparse block-diagonal design matrix.
#' @param W_block   diagonal weight matrix.
#' @param y_tilde   Numeric vector \[N_stacked\]: current partially-out outcome.
#' @param row_ranges List of J0 integer vectors from .precompute_XtWX_blocks.
#' @param K         Integer. Covariates per unit.
#' @param J0        Integer. Number of treated units.
#'
#' @return List of J0 numeric vectors each of length K.
#' @keywords internal
.compute_XtWy_blocks <- function(X_block, W_block, y_tilde,
                                 row_ranges, K, J0) {

  col_starts    <- (seq_len(J0) - 1L) * K + 1L
  W_diag        <- Matrix::diag(W_block)
  y_vec         <- as.numeric(y_tilde)
  XtWy_blocks   <- vector("list", J0)

  for (j in seq_len(J0)) {
    rows_j  <- row_ranges[[j]]
    cols_j  <- col_starts[j]:(col_starts[j] + K - 1L)
    Xj      <- as.matrix(X_block[rows_j, cols_j, drop = FALSE])
    Wj      <- W_diag[rows_j]
    XtWy_blocks[[j]] <- as.numeric(crossprod(Xj, Wj * y_vec[rows_j]))
  }

  XtWy_blocks
}


# -----------------------------------------------------------------------------
#' Fast block-diagonal matrix-vector product (replaces X_block %*% beta_bd)
#'
#' For block-diagonal X_block and coefficient vector beta_bd, computes
#' X_block %*% beta_bd by iterating over J0 unit blocks rather than
#' performing the full sparse matrix-vector product. For K=1 each block
#' reduces to a scalar multiply-and-assign, avoiding sparse matrix overhead.
#'
#' @param beta_bd   Numeric vector length K*J0.
#' @param row_ranges List of J0 integer vectors: row indices per unit block.
#' @param X_block   Sparse block-diagonal Matrix \[N_stacked x K*J0\].
#' @param K         Integer. Covariates per unit.
#' @param J0        Integer. Number of treated units.
#' @param N_stacked Integer. Total rows in X_block.
#'
#' @return Numeric vector length N_stacked.
#' @keywords internal
.Xbeta_blockdiag <- function(beta_bd, row_ranges, X_block, K, J0, N_stacked) {
  out        <- numeric(N_stacked)
  col_starts <- (seq_len(J0) - 1L) * K + 1L

  if (K == 1L) {
    # Fast path: scalar multiply per unit block
    for (j in seq_len(J0)) {
      out[row_ranges[[j]]] <- beta_bd[j]   # X_j is all-ones in treated col
      # General: multiply by X slice — for K=1 X_j %*% beta_j = beta_j * x_j
      # but x_j is a column of X_block. Fetch it:
      xj <- as.numeric(X_block[row_ranges[[j]], col_starts[j], drop = TRUE])
      out[row_ranges[[j]]] <- xj * beta_bd[j]
    }
  } else {
    for (j in seq_len(J0)) {
      rows_j <- row_ranges[[j]]
      cols_j <- col_starts[j]:(col_starts[j] + K - 1L)
      Xj     <- as.matrix(X_block[rows_j, cols_j, drop = FALSE])
      bj     <- beta_bd[cols_j]
      out[rows_j] <- as.numeric(Xj %*% bj)
    }
  }
  out
}


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
#' @section Performance:
#' For large datasets with many treated units, install \pkg{RhpcBLASctl}
#' to enable multi-threaded BLAS operations:
#' \code{install.packages("RhpcBLASctl")}. The package will use it
#' automatically if available.
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

  y_block <- gdata$Y_block
  X_block <- gdata$X_block
  W_block <- gdata$W
  K  <- length(gdata$cov$Xcols)
  J0 <- gdata$cov$J0

  # performance package
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    n_cores <- parallel::detectCores(logical = FALSE)
    if (RhpcBLASctl::blas_get_num_procs() < n_cores) {
      RhpcBLASctl::blas_set_num_threads(n_cores)
        message(sprintf(
          "RhpcBLASctl: BLAS threads set to %d for matrix operations.", n_cores
        ))
    }
  } else {
    message(paste0(
      "Install RhpcBLASctl for faster BLAS threading: ",
      "install.packages('RhpcBLASctl')"
    ))
  }

  fs    <- gdata$first_stage
  X_fs  <- fs$X_fs
  d_vec <- as.numeric(fs$d)
  n_obs <- length(d_vec)
  p_fs  <- ncol(X_fs)

  if (is.null(X_fs))
    stop("gdata$first_stage$X_fs is NULL.", call. = FALSE)
  if (nrow(X_fs) != n_obs)
    stop("nrow(X_fs) != length(d).", call. = FALSE)

  a_sigma_alpha_prior <- ctrl$a_sigma_alpha_prior
  b_sigma_alpha_prior <- ctrl$b_sigma_alpha_prior
  a_sigma_tau_prior   <- ctrl$a_sigma_tau_prior
  b_sigma_tau_prior   <- ctrl$b_sigma_tau_prior
  sigma2_rho_prior    <- if (!is.null(ctrl$sigma2_rho_prior))
    as.numeric(ctrl$sigma2_rho_prior) else 10

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

  # precompute first-stage cross-products
  XfsXfs       <- crossprod(X_fs)
  A_delta_base <- XfsXfs + Sigma_delta_prior_inv

  # precompute nu_hat alignment index (replaces per-iter data.table merge)
  nu_row_idx <- .build_nu_hat_index(
    X_idlist = gdata$X_idlist,
    dta_id   = gdata$dtaidx[["id"]],
    dta_wID  = gdata$dtaidx[["wID"]]
  )

  # precompute treated/control split for fast z* sampling
  treat_idx <- which(d_vec == 1L)
  ctrl_idx  <- which(d_vec == 0L)

  # BLOCK-DIAGONAL OPTIMISATION: precompute XtWX per-unit blocks once.
  # A_bd = XtWX_bd/sigma2 + diag(J0) x diag(1/tau) is block-diagonal.
  # Factorising J0 blocks of [K x K] costs O(J0 * K^3) vs O((K*J0)^3)
  # for the full matrix — a J0^2 = 1879^2 ~ 3.5M fold improvement.
  bd_pre     <- .precompute_XtWX_blocks(X_block, W_block, K, J0)
  XtWX_blocks <- bd_pre$XtWX_blocks
  row_ranges  <- bd_pre$row_ranges

  # precompute W_diag, N_stacked, and X column slices for fast Xbeta
  W_diag    <- Matrix::diag(W_block)
  N_stacked <- nrow(X_block)

  # For K=1: cache the single column slice per unit (avoids per-iter sparse extract)
  if (K == 1L) {
    X_cols_cache <- lapply(seq_len(J0), function(j)
      as.numeric(X_block[row_ranges[[j]], j, drop = TRUE])
    )
  } else {
    X_cols_cache <- NULL
  }

  # ---- initial values
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

    # BLOCK 1: Sample z*  (fast inverse-CDF, no ifelse dispatch)
    eta    <- as.numeric(X_fs %*% delta)
    z_star <- .sample_z_star(eta, d_vec, treat_idx, ctrl_idx)

    # BLOCK 2: Sample delta | z*, X_fs
    b_delta     <- as.numeric(crossprod(X_fs, z_star)) +
      as.numeric(Sigma_delta_prior_inv %*% mu_delta_prior)
    cA_delta    <- chol(A_delta_base)
    A_delta_inv <- chol2inv(cA_delta)
    mu_delta    <- as.numeric(A_delta_inv %*% b_delta)
    delta       <- rMVNormCovariance(1, mu = mu_delta, Sigma = A_delta_inv)

    # BLOCK 3: Build nu_hat (integer index lookup, no merge)
    nu_hat         <- z_star - as.numeric(X_fs %*% delta)
    nu_hat_stacked <- nu_hat[nu_row_idx]

    # BLOCK 4a: Sample beta_bd block-diagonally (KEY optimisation)
    # y_tilde = y - rho * nu_hat_stacked (partial out selection correction)
    y_tilde_vec <- as.numeric(y_block) - rho * nu_hat_stacked
    XtWy_blocks <- .compute_XtWy_blocks(X_block, W_block, y_tilde_vec,
                                        row_ranges, K, J0)
    beta_bd <- .sample_beta_bd_blockdiag(XtWX_blocks, XtWy_blocks,
                                         tau, sigma2, K)

    # BLOCK 4b: Sample rho (scalar)
    # y_tilde2 = y - X_block %*% beta_bd (fast block-diagonal product)
    if (K == 1L) {
      Xb <- numeric(N_stacked)
      for (j in seq_len(J0)) Xb[row_ranges[[j]]] <- X_cols_cache[[j]] * beta_bd[j]
    } else {
      Xb <- .Xbeta_blockdiag(beta_bd, row_ranges, X_block, K, J0, N_stacked)
    }
    y_tilde2_vec <- as.numeric(y_block) - Xb
    nuWnu <- sum(W_diag * nu_hat_stacked^2)
    nuWy2 <- sum(W_diag * nu_hat_stacked * y_tilde2_vec)
    v_rho <- 1 / (nuWnu / sigma2 + 1 / sigma2_rho_prior)
    m_rho <- v_rho * nuWy2 / sigma2
    rho   <- stats::rnorm(1, mean = m_rho, sd = sqrt(v_rho))

    # BLOCK 5: Sample sigma2
    res_vec <- y_tilde2_vec - rho * nu_hat_stacked
    alpha_s <- a_sigma_alpha_prior + length(y_block) / 2
    beta_s  <- b_sigma_alpha_prior + sum(W_diag * res_vec^2) / 2
    sigma2  <- 1 / stats::rgamma(1, shape = alpha_s, rate = beta_s)

    # BLOCK 6: Sample tau_k
    beta_matrix <- matrix(beta_bd, ncol = K, byrow = TRUE)
    for (k in seq_len(K)) {
      tau[k] <- sqrt(1 / stats::rgamma(
        1,
        shape = a_sigma_tau_prior + J0 / 2,
        rate  = b_sigma_tau_prior + sum(beta_matrix[, k]^2) / 2
      ))
    }

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

  # Performance improvement
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    n_cores <- parallel::detectCores(logical = FALSE)
    if (RhpcBLASctl::blas_get_num_procs() < n_cores) {
      RhpcBLASctl::blas_set_num_threads(n_cores)
        message(sprintf(
          "RhpcBLASctl: BLAS threads set to %d for matrix operations.", n_cores
        ))
    }
  } else {
    message(paste0(
      "Install RhpcBLASctl for faster BLAS threading: ",
      "install.packages('RhpcBLASctl')"
    ))
  }

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

  # precompute first-stage cross-products
  XfsXfs       <- crossprod(X_fs)
  A_delta_base <- XfsXfs + Sigma_delta_prior_inv

  # precompute Z_star (fixed: Z and K don't change across iterations)
  dZ     <- rep(0, K); dZ[k_tr] <- 1
  Z_star <- Matrix::kronecker(Z, dZ)   # [K*J0 x G]

  # precompute nu_hat alignment index
  nu_row_idx <- .build_nu_hat_index(
    X_idlist = gdata$X_idlist,
    dta_id   = gdata$dtaidx[["id"]],
    dta_wID  = gdata$dtaidx[["wID"]]
  )

  # precompute treated/control split for fast z* sampling
  treat_idx <- which(d_vec == 1L)
  ctrl_idx  <- which(d_vec == 0L)

  # BLOCK-DIAGONAL OPTIMISATION: precompute XtWX per-unit blocks once
  bd_pre      <- .precompute_XtWX_blocks(X_block, W_block, K, J0)
  XtWX_blocks <- bd_pre$XtWX_blocks
  row_ranges  <- bd_pre$row_ranges
  W_diag      <- Matrix::diag(W_block)
  N_stacked   <- nrow(X_block)

  # Cache X column slices for fast block-diagonal Xbeta (K=1 fast path)
  if (K == 1L) {
    X_cols_cache <- lapply(seq_len(J0), function(j)
      as.numeric(X_block[row_ranges[[j]], j, drop = TRUE])
    )
  } else {
    X_cols_cache <- NULL
  }

  # Precompute per-unit gamma prior mean contribution from Z_star.
  # Z_star %*% gamma gives the moderator-implied prior mean for all K*J0
  # coefficients. For the block-diagonal sampler we need it split per unit.
  # col_starts_j gives which row of Z_star is the treatment row for unit j.
  # Since k_tr is the treatment covariate index within each K-block:
  # Z_star row for unit j, covariate k_tr is (j-1)*K + k_tr
  tr_rows <- (seq_len(J0) - 1L) * K + k_tr   # length J0

  # ---- initial values
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

    # BLOCK 1: Sample z*  (fast inverse-CDF)
    eta    <- as.numeric(X_fs %*% delta)
    z_star <- .sample_z_star(eta, d_vec, treat_idx, ctrl_idx)

    # BLOCK 2: Sample delta | z*, X_fs
    b_delta     <- as.numeric(crossprod(X_fs, z_star)) +
      as.numeric(Sigma_delta_prior_inv %*% mu_delta_prior)
    cA_delta    <- chol(A_delta_base)
    A_delta_inv <- chol2inv(cA_delta)
    mu_delta    <- as.numeric(A_delta_inv %*% b_delta)
    delta       <- rMVNormCovariance(1, mu = mu_delta, Sigma = A_delta_inv)

    # BLOCK 3: Build nu_hat (integer index lookup)
    nu_hat         <- z_star - as.numeric(X_fs %*% delta)
    nu_hat_stacked <- nu_hat[nu_row_idx]

    # BLOCK 4a: Sample beta_bd block-diagonally with moderator prior mean
    # Prior mean for unit j: zero except at treatment covariate k_tr where
    # it equals z_j' gamma (the moderator-implied shrinkage target).
    gamma_prior_means <- vector("list", J0)
    z_gamma <- as.numeric(Z %*% gamma)   # [J0] moderator contributions
    for (j in seq_len(J0)) {
      pm_j <- rep(0, K)
      pm_j[k_tr] <- z_gamma[j]
      gamma_prior_means[[j]] <- pm_j
    }

    y_tilde_vec <- as.numeric(y_block) - rho * nu_hat_stacked
    XtWy_blocks <- .compute_XtWy_blocks(X_block, W_block, y_tilde_vec,
                                        row_ranges, K, J0)
    beta_bd <- .sample_beta_bd_blockdiag(XtWX_blocks, XtWy_blocks,
                                         tau, sigma2, K,
                                         prior_mean = gamma_prior_means)

    # ==================================================================
    # BLOCK 4b: Sample gamma | beta_bd, tau, Z
    beta_matrix <- matrix(beta_bd, ncol = K, byrow = TRUE)
    beta_tr     <- as.numeric(beta_matrix[, k_tr, drop = TRUE])
    V_gamma <- solve(crossprod(Z) / tau[k_tr]^2 + Sigma_gamma_prior_inv)
    rhs     <- as.numeric(crossprod(Z, beta_tr)) / tau[k_tr]^2 +
      as.numeric(Sigma_gamma_prior_inv %*% mu_gamma_prior)
    m_gamma <- V_gamma %*% rhs
    gamma   <- rMVNormCovariance(1, mu = as.numeric(m_gamma), Sigma = V_gamma)

    # BLOCK 4c: Sample rho (scalar) — fast block-diagonal product
    if (K == 1L) {
      Xb <- numeric(N_stacked)
      for (j in seq_len(J0)) Xb[row_ranges[[j]]] <- X_cols_cache[[j]] * beta_bd[j]
    } else {
      Xb <- .Xbeta_blockdiag(beta_bd, row_ranges, X_block, K, J0, N_stacked)
    }
    y_tilde2_vec <- as.numeric(y_block) - Xb
    nuWnu <- sum(W_diag * nu_hat_stacked^2)
    nuWy2 <- sum(W_diag * nu_hat_stacked * y_tilde2_vec)
    v_rho <- 1 / (nuWnu / sigma2 + 1 / sigma2_rho_prior)
    m_rho <- v_rho * nuWy2 / sigma2
    rho   <- stats::rnorm(1, mean = m_rho, sd = sqrt(v_rho))

    # BLOCK 5: Sample sigma2 — use vector ops
    res_vec <- y_tilde2_vec - rho * nu_hat_stacked
    alpha_s <- a_sigma_alpha_prior + length(y_block) / 2
    beta_s  <- b_sigma_alpha_prior + sum(W_diag * res_vec^2) / 2
    sigma2  <- 1 / stats::rgamma(1, shape = alpha_s, rate = beta_s)

    # BLOCK 6: Sample tau_k | beta_bd, gamma (moderator-adjusted residuals)
    for (k in seq_len(K)) {
      tau[k] <- sqrt(1 / stats::rgamma(
        1,
        shape = a_sigma_tau_prior + J0 / 2,
        rate  = b_sigma_tau_prior +
          sum((beta_matrix[, k] -
                 Z_star[((1:J0) - 1L) * K + k, , drop = FALSE] %*% gamma)^2) / 2
      ))
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
#' Fast truncated normal sampler via inverse-CDF
#'
#' Avoids the ifelse dispatch overhead of rtruncnorm for large vectors by
#' splitting treated/untreated indices and sampling each group separately.
#' Roughly 2-3x faster than truncnorm::rtruncnorm for n > 100k.
#'
#' @param eta       Numeric vector of linear predictors (length n_obs).
#' @param d_vec     Integer/logical vector of treatment indicators (length n_obs).
#' @param treat_idx Integer vector of indices where d_vec == 1 (precomputed).
#' @param ctrl_idx  Integer vector of indices where d_vec == 0 (precomputed).
#' @return Numeric vector z_star of length n_obs.
#' @keywords internal
.sample_z_star <- function(eta, d_vec, treat_idx, ctrl_idx) {
  n     <- length(eta)
  z     <- numeric(n)

  # Treated: TN_{[0, Inf)}(eta_i, 1) via inverse CDF
  if (length(treat_idx) > 0L) {
    mu_t  <- eta[treat_idx]
    p_lo  <- stats::pnorm(-mu_t)          # P(Z < 0) = P(z* < 0 | treated)
    u     <- stats::runif(length(treat_idx), min = p_lo, max = 1)
    z[treat_idx] <- mu_t + stats::qnorm(u)
  }

  # Untreated: TN_{(-Inf, 0)}(eta_i, 1) via inverse CDF
  if (length(ctrl_idx) > 0L) {
    mu_c  <- eta[ctrl_idx]
    p_hi  <- stats::pnorm(-mu_c)          # P(Z < 0)
    u     <- stats::runif(length(ctrl_idx), min = 0, max = p_hi)
    z[ctrl_idx] <- mu_c + stats::qnorm(u)
  }

  z
}


# -----------------------------------------------------------------------------
#' Build integer index mapping X_block rows to original dta rows
#'
#' Computes once (before the Gibbs loop) the integer vector that maps each
#' row of X_block (stacked pseudo-panel order) back to its corresponding row
#' in the original dta. Inside the loop, nu_hat_stacked = nu_hat\[nu_row_idx\]
#' replaces the expensive data.table merge in .build_nu_hat_stacked.
#'
#' @param X_idlist data.table with columns id and wID (nrow = nrow(X_block)).
#' @param dta_id   Character vector of unit ids in original dta row order.
#' @param dta_wID  Numeric/integer vector of wID in original dta row order.
#'
#' @return Integer vector of length nrow(X_block). Entry i gives the row
#'   index in the original dta corresponding to row i of X_block.
#' @keywords internal
.build_nu_hat_index <- function(X_idlist, dta_id, dta_wID) {

  # Build lookup: (id, wID) -> original row index
  lookup <- data.table::data.table(
    id      = as.character(dta_id),
    wID     = dta_wID,
    row_idx = seq_along(dta_id)
  )
  data.table::setkeyv(lookup, c("id", "wID"))

  query <- data.table::copy(X_idlist)
  query[, id := as.character(id)]
  data.table::setkeyv(query, c("id", "wID"))

  merged <- lookup[query, on = c("id", "wID")]

  if (anyNA(merged$row_idx))
    stop(paste0(
      ".build_nu_hat_index: ", sum(is.na(merged$row_idx)),
      " unmatched rows. Check that X_idlist and dta share the same (id, wID) keys."
    ), call. = FALSE)

  as.integer(merged$row_idx)
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
