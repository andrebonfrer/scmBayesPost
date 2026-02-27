#' Gibbs sampling dispatcher for post-SCM models
#'
#' Runs an appropriate Gibbs sampler depending on the complexity level encoded in gdata$cov.
#'
#' @param gdata Output from prepare_data_general().
#' @param n_iter Number of iterations.
#' @param burn_in Burn-in iterations.
#' @param run_selection_gibbs Logical; run selection equation Gibbs block (only relevant for full sampler).
#' @param do.CF Use control function when IV blocks exist.
#' @param selection.method Character; "2SRI" or "GR" (only relevant if run_selection_gibbs=TRUE).
#' @param Z_cov_dense Dense covariance for IV first stage (requires MCMCpack if TRUE).
#' @param CF.interactions Expand CF residual interactions.
#'
#' @return Posterior draws.
#' @export
gibbs_postscm <- function(gdata,
                          n_iter = 1000,
                          burn_in = 500,
                          run_selection_gibbs = FALSE,
                          do.CF = TRUE,
                          selection.method = c("2SRI","GR"),
                          Z_cov_dense = FALSE,
                          CF.interactions = TRUE) {

  selection.method <- match.arg(selection.method)

  lvl <- gdata$cov$second_stage

  if (is.null(lvl) || lvl == "none") {
    # Minimal branch ignores selection/IV blocks by design
    y <- gdata$Y_block
    X <- gdata$X_block
    W <- gdata$W

    XtW <- Matrix::t(X) %*% W
    XtWX <- XtW %*% X
    XtWy <- XtW %*% y

    p <- ncol(X)
    beta <- matrix(0, p, 1)
    sigma2 <- 1

    beta_samples <- matrix(0, n_iter - burn_in, p)
    sigma2_samples <- numeric(n_iter - burn_in)

    pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
    for (iter in 1:n_iter) {
      A <- XtWX / sigma2 + diag(1e-8, p)
      b <- XtWy / sigma2
      cA <- Matrix::chol(A)
      Ainv <- Matrix::chol2inv(cA)
      mu <- Ainv %*% b
      beta <- rMVNormCovariance(1, mu = mu, Sigma = Ainv)

      resid <- y - X %*% beta
      alpha <- 2 + length(y) / 2
      rate <- 2 + as.numeric(Matrix::t(resid) %*% W %*% resid) / 2
      sigma2 <- 1 / stats::rgamma(1, shape = alpha, rate = rate)

      if (iter > burn_in) {
        beta_samples[iter - burn_in, ] <- as.numeric(beta)
        sigma2_samples[iter - burn_in] <- sigma2
      }
      utils::setTxtProgressBar(pb, iter)
    }
    return(list(beta_samples = beta_samples, sigma2_samples = sigma2_samples))
  }

  # Full branch
  gibbs_sampling(
    gdata = gdata,
    n_iter = n_iter,
    burn_in = burn_in,
    run_selection_gibbs = run_selection_gibbs,
    do.CF = do.CF,
    selection.method = selection.method,
    Z_cov_dense = Z_cov_dense,
    CF.interactions = CF.interactions
  )
}


#' Gibbs sampler for post–synthetic-control hierarchical models (internal)
#'
#' Runs the full Gibbs sampler for the second-stage Bayesian model used by
#' `scmBayesPost` when moderator and/or endogeneity components are present.
#' This function is **internal** and is called by [gibbs_postscm()].
#'
#' The sampler targets a hierarchical regression of treated-unit (or unit–event)
#' treatment-effect parameters on moderator covariates, with optional
#' endogeneity correction via (i) a selection-equation Gibbs block for binary
#' treatment assignment, and/or (ii) instrumental-variables/control-function
#' blocks for endogenous moderators.
#'
#' @param gdata A data object produced by [prepare_data_general()]. Must contain
#'   at least `Y_block`, `X_block`, and `W`. For second-stage models it must also
#'   contain `Z_block` and, if IV functionality is used, `Z.instruments`.
#' @param n_iter Integer; total number of Gibbs iterations.
#' @param burn_in Integer; number of initial iterations discarded as burn-in.
#' @param run_selection_gibbs Logical; if `TRUE`, include a Gibbs block for a
#'   binary treatment/selection equation and update the relevant column of
#'   `X_block` with the resulting residual (2SRI) or generalized residual (GR).
#'   Requires `gdata$first_stage` components (e.g., `GRX`, `GRY`) and a column in
#'   `X_block` corresponding to the residual (typically `"GR"`).
#' @param do.CF Logical; if `TRUE`, apply a control-function approach when
#'   `gdata$Z.instruments` is provided (append first-stage residuals, optionally
#'   interacted). If `FALSE`, use a 2SLS-style replacement with first-stage
#'   predicted values. Ignored when no instruments are provided.
#' @param selection.method Character; residual type used in the selection block:
#'   `"2SRI"` uses probit residual inclusion; `"GR"` uses generalized residuals.
#' @param Z_cov_dense Logical; if `TRUE`, estimate a dense covariance matrix for
#'   the first-stage IV residuals (requires `MCMCpack` for inverse-Wishart draws);
#'   if `FALSE`, assumes a diagonal covariance.
#' @param CF.interactions Logical; if `TRUE` and `do.CF = TRUE`, include
#'   interaction terms among control-function residuals.
#'
#' @return A list of posterior draws. Components include (depending on model):
#'   \describe{
#'     \item{beta_samples}{Matrix of draws for stacked second-stage coefficients.}
#'     \item{gamma_samples}{Matrix of draws for moderator effects.}
#'     \item{sigma2_samples}{Vector of draws for outcome noise variance.}
#'     \item{tau_samples}{Matrix of draws for hierarchical scale parameters.}
#'     \item{beta_fs_samples}{(Optional) draws for selection-equation coefficients.}
#'     \item{sigma2_fs_samples}{(Optional) draws for selection-equation variance.}
#'     \item{beta.Z_samples}{(Optional) list of draws for IV first-stage coefficients.}
#'     \item{sigma.Z_samples}{(Optional) draws for IV covariance (diagonal or dense).}
#'   }
#'
#' @seealso [gibbs_postscm()], [prepare_data_general()]
#'
#' @keywords internal
gibbs_sampling <- function(gdata,
                                n_iter = 1000,
                                burn_in = 500,
                                run_selection_gibbs = FALSE,
                                do.CF = TRUE,
                                selection.method = c("2SRI","GR"),
                                Z_cov_dense = FALSE,
                                CF.interactions = TRUE) {

  selection.method <- match.arg(selection.method)

  y_block <- gdata$Y_block
  X_block <- gdata$X_block
  Zbase   <- gdata$Z_block
  W_block <- gdata$W
  Zim     <- gdata$Z.instruments

  if (!do.CF) Z_cov_dense <- FALSE

  # selection equation support
  if (run_selection_gibbs) {
    if (is.null(gdata$dtaidx) || is.null(gdata$GRX) || is.null(gdata$GRY) || is.null(gdata$X_idlist)) {
      stop("run_selection_gibbs=TRUE requires gdata$dtaidx, GRX, GRY, and X_idlist.")
    }
    .dtaidx <- data.table::copy(gdata$dtaidx)
    GRX <- gdata$GRX
  }

  K  <- length(gdata$cov$Xcols)
  J0 <- gdata$cov$J0
  G  <- ncol(Zbase)

  # moderator endogeneity bookkeeping
  if (!is.null(Zim)) {
    LZ <- L <- length(Zim$dv)
    if (do.CF) {
      if (CF.interactions) L <- (2^LZ - 1)
      G <- G + L
    }
  }

  # priors
  Sigma_gamma_prior <- 10
  mu_gamma_prior <- matrix(rep(0, G), nrow = G)

  a_sigma_alpha_prior <- 2
  b_sigma_alpha_prior <- 2
  a_sigma_tau_prior <- 2
  b_sigma_tau_prior <- 2

  # init
  beta   <- matrix(0, K * J0, 1)
  gamma  <- rep(0, G)
  sigma2 <- 1
  tau    <- rep(1, K)

  if (run_selection_gibbs) {
    beta_fs <- matrix(0, ncol(GRX), 1)
    sigma2_fs <- 1
  }

  if (!is.null(Zim)) {
    if (Z_cov_dense) {
      sigma.Z <- initialize_dense_covariance(LZ)
    } else {
      sigma.Z <- diag(rep(1, LZ))
    }
  }

  # precompute
  XtW  <- Matrix::t(X_block) %*% W_block
  XtWX <- XtW %*% X_block
  XtWy <- XtW %*% y_block

  n_keep <- n_iter - burn_in
  beta_samples  <- matrix(0, n_keep, K * J0)
  gamma_samples <- matrix(0, n_keep, G)
  sigma2_samples <- numeric(n_keep)
  tau_samples <- matrix(0, n_keep, K)

  if (run_selection_gibbs) {
    beta_fs_samples <- matrix(0, n_keep, ncol(GRX))
    sigma2_fs_samples <- numeric(n_keep)
  }

  if (!is.null(Zim)) {
    beta.Z_samples <- vector("list", LZ)
    for (l in 1:LZ) beta.Z_samples[[l]] <- matrix(0, n_keep, ncol(Zim$Z.im[[l]]))
    if (Z_cov_dense) sigma.Z_samples <- array(0, dim = c(LZ, LZ, n_keep))
    else sigma.Z_samples <- matrix(0, n_keep, LZ)
  }

  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)

  for (iter in seq_len(n_iter)) {

    ## 1) selection equation block (optional)
    if (run_selection_gibbs) {

      .bd <- gibbs_binomial_probit(
        X = GRX, y = gdata$GRY,
        beta = beta_fs,
        Sigma_prior = 0.2 * diag(ncol(GRX)),
        mu_prior = rep(0, ncol(GRX)),
        sigma2 = sigma2_fs,
        alpha_0 = 4,
        beta_0 = 4
      )

      beta_fs <- matrix(.bd$beta, ncol = 1)
      sigma2_fs <- .bd$sigma2

      resid <- gdata$GRY - stats::pnorm(GRX %*% beta_fs)
      PR <- stats::pnorm(GRX %*% beta_fs)
      GR <- gdata$GRY * stats::dnorm(PR) / stats::pnorm(PR) -
        (1 - gdata$GRY) * stats::dnorm(-PR) / stats::pnorm(-PR)
      .r <- if (selection.method == "GR") GR else resid

      .dtaidx[, fs_resid := as.numeric(.r)]
      .dtaX <- gdata$X_idlist
      data.table::setkey(.dtaidx, id, wID)
      data.table::setkey(.dtaX,   id, wID)

      CR <- .dtaX[.dtaidx][["fs_resid"]]

      X_block <- replace_kth_column_block_diagonal(X_block, new_col = CR, K = K,
                                                   k = which(gdata$cov$Xcols=="GR"))

      XtW  <- Matrix::t(X_block) %*% W_block
      XtWX <- XtW %*% X_block
      XtWy <- XtW %*% y_block
    }

    ## 2) endogenous moderators first stage (optional)
    if (!is.null(Zim)) {
      Ylist <- Qzlist <- vector("list", LZ)
      for (l in 1:LZ) {
        Ylist[[l]] <- Zbase[, Zim$dv[l]]
        Qzlist[[l]] <- Zim$Z.im[[l]]
      }

      gibbsZ <- gibbs_sampler_one_draw_cov(
        y_list = Ylist,
        X_list = Qzlist,
        sigma_param = sigma.Z,
        covariance_type = if (Z_cov_dense) "dense" else "diagonal"
      )

      sigma.Z <- gibbsZ$sigma2_samples
      beta.Z  <- gibbsZ$beta_samples

      if (do.CF) {
        .lmat <- as.data.frame(gibbsZ$residuals)
        colnames(.lmat) <- paste0("lres", seq_len(LZ))
        if (CF.interactions) {
          .lmat <- stats::model.matrix(stats::formula(paste("~ 0 + .^", LZ)), data = .lmat)
        }
        Z <- Matrix::as.matrix(cbind(Zbase, .lmat))
      } else {
        .pmat <- as.matrix(gibbsZ$predicted)
        Z <- Zbase
        for (l in 1:LZ) {
          cl <- match(Zim$dv[l], colnames(Z))
          Z[, cl] <- .pmat[, l]
          colnames(Z)[cl] <- paste0("pred.", Zim$dv[l])
        }
      }
    } else {
      Z <- Zbase
    }

    ## 3) Z_star for treatment coefficient moderation
    dZ <- rep(0, K)
    dZ[gdata$cov$intX] <- 1
    Z_star <- Matrix::kronecker(Z, dZ)

    ## 4) sample beta
    Sigma_beta_prior_inv <- Matrix::kronecker(diag(J0), diag(1 / tau))
    q <- Z_star %*% gamma
    A <- XtWX / sigma2 + Sigma_beta_prior_inv
    b <- XtWy / sigma2 + Sigma_beta_prior_inv %*% q

    cA <- Matrix::chol(A)
    A_inv <- Matrix::chol2inv(cA)
    mu_beta <- A_inv %*% b
    beta <- rMVNormCovariance(1, mu = mu_beta, Sigma = A_inv)

    beta_matrix <- matrix(beta, ncol = K, byrow = TRUE)

    ## 5) sample gamma
    Sigma_gamma_prior_inv <- (1 / Sigma_gamma_prior) * diag(G)
    beta_tr <- beta_matrix[, gdata$cov$intX]
    V_gamma <- Matrix::solve(t(Z) %*% Z / tau[2]^2 + Sigma_gamma_prior_inv)
    m_gamma <- V_gamma %*% (t(Z) %*% beta_tr / tau[2]^2 + Sigma_gamma_prior_inv %*% mu_gamma_prior)
    gamma <- rMVNormCovariance(1, mu = m_gamma, Sigma = V_gamma)

    ## 6) sample sigma2
    resid_y <- y_block - X_block %*% beta
    alpha <- a_sigma_alpha_prior + length(y_block) / 2
    rate  <- b_sigma_alpha_prior + as.numeric(Matrix::t(resid_y) %*% W_block %*% resid_y) / 2
    sigma2 <- 1 / stats::rgamma(1, shape = alpha, rate = rate)

    ## 7) sample tau
    for (k in seq_len(K)) {
      mu_k <- Z_star[((seq_len(J0) - 1) * K + k), , drop = FALSE] %*% gamma
      tau[k] <- sqrt(1 / stats::rgamma(
        1,
        shape = a_sigma_tau_prior + (J0 / 2),
        rate  = b_sigma_tau_prior + (sum((beta_matrix[, k] - mu_k)^2) / 2)
      ))
    }

    ## store
    if (iter > burn_in) {
      s <- iter - burn_in
      beta_samples[s, ] <- as.numeric(beta)
      gamma_samples[s, ] <- as.numeric(gamma)
      sigma2_samples[s]  <- sigma2
      tau_samples[s, ]   <- tau

      if (run_selection_gibbs) {
        beta_fs_samples[s, ] <- as.numeric(beta_fs)
        sigma2_fs_samples[s] <- sigma2_fs
      }

      if (!is.null(Zim)) {
        for (l in seq_len(LZ)) beta.Z_samples[[l]][s, ] <- as.numeric(beta.Z[[l]])
        if (Z_cov_dense) sigma.Z_samples[, , s] <- sigma.Z
        else sigma.Z_samples[s, ] <- diag(sigma.Z)
      }
    }

    utils::setTxtProgressBar(pb, iter)
  }

  out <- list(
    beta_samples  = beta_samples,
    gamma_samples = gamma_samples,
    sigma2_samples = sigma2_samples,
    tau_samples   = tau_samples
  )
  if (run_selection_gibbs) {
    out$beta_fs_samples <- beta_fs_samples
    out$sigma2_fs_samples <- sigma2_fs_samples
  }
  if (!is.null(Zim)) {
    out$beta.Z_samples <- beta.Z_samples
    out$sigma.Z_samples <- sigma.Z_samples
  }
  out
}




#' Gibbs Sampling for a Binomial Probit Model with sigma^2
#'
#' This function performs one Gibbs sampling iteration for a binomial probit model,
#' with an additional step to update the variance parameter \eqn{\sigma^2}.
#' The model is defined as: \eqn{y_i = 1} if \eqn{z_i > 0}, and \eqn{y_i = 0} if \eqn{z_i \leq 0},
#' where \eqn{z_i = X_i \beta + \epsilon_i}, \eqn{\epsilon_i \sim N(0, \sigma^2)}.
#' A prior distribution is assumed for both \eqn{\beta} and \eqn{\sigma^2}.
#' Specifically, \eqn{\beta \sim N(\mu_{\beta}, \Sigma_{\beta})}, and \eqn{\sigma^2 \sim \text{Inv-Gamma}(\alpha_0, \beta_0)}.
#'
#' @param X A matrix of predictor variables of dimension \eqn{n \times p}, where \eqn{n} is the number of observations and \eqn{p} is the number of predictors.
#' @param y A binary vector (0/1) of response variables of length \eqn{n}, where \eqn{y_i \in \{0,1\}}.
#' @param beta A numeric vector of initial values for the regression coefficients \eqn{\beta} of length \eqn{p}.
#' @param Sigma_prior A \eqn{p \times p} matrix representing the prior covariance matrix for \eqn{\beta}.
#' @param mu_prior A numeric vector of length \eqn{p} representing the prior mean vector for \eqn{\beta}.
#' @param sigma2 A numeric value representing the initial value of \eqn{\sigma^2}.
#' @param alpha_0 A numeric value representing the shape parameter of the inverse-gamma prior for \eqn{\sigma^2}.
#' @param beta_0 A numeric value representing the scale parameter of the inverse-gamma prior for \eqn{\sigma^2}.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{beta}}{A numeric vector of updated values for the regression coefficients \eqn{\beta}.}
#'   \item{\code{sigma2}}{An updated value for the variance parameter \eqn{\sigma^2}.}
#'   \item{\code{z}}{A numeric vector of latent variables \eqn{z} sampled from the truncated normal distribution.}
#' }
#'
#' @examples
#' # Example usage:
#' set.seed(123)
#' X <- matrix(rnorm(100 * 3), 100, 3)  # 100 observations, 3 predictors
#' y <- rbinom(100, 1, 0.5)  # Binary outcome (0 or 1)
#' beta <- rep(0, 3)  # Starting beta values
#' Sigma_prior <- diag(3)  # Prior covariance matrix
#' mu_prior <- rep(0, 3)  # Prior mean for beta
#' sigma2 <- 1  # Initial value for sigma^2
#' alpha_0 <- 2  # Shape parameter for prior on sigma^2
#' beta_0 <- 2  # Scale parameter for prior on sigma^2
#'
#' result <- gibbs_binomial_probit(X, y, beta, Sigma_prior, mu_prior, sigma2, alpha_0, beta_0)
#' print(result$beta)  # Updated beta values
#' print(result$sigma2)  # Updated sigma^2 value
#'
#' @importFrom truncnorm rtruncnorm
#' @export
gibbs_binomial_probit <- function(X, y, beta, Sigma_prior,
                                  mu_prior, sigma2 = 1, alpha_0 = 2,
                                  beta_0 = 2) {

  # Number of observations and predictors
  n <- length(y)
  p <- ncol(X)

  # Step 1: Sample latent variable z from the truncated normal distribution
  z <- rep(0, n)

  for (i in 1:n) {
    mu_i <- X[i, ] %*% beta

    if (y[i] == 1) {
      # Sample z from truncated normal distribution, truncated below at 0 for y == 1
      z[i] <- rtruncnorm(1, a = 0, mean = mu_i, sd = sqrt(sigma2))
    } else {
      # Sample z from truncated normal distribution, truncated above at 0 for y == 0
      z[i] <- rtruncnorm(1, b = 0, mean = mu_i, sd = sqrt(sigma2))
    }
  }

  # Step 2: Sample beta from the normal posterior
  Sigma_posterior <- solve(t(X) %*% X / sigma2 + solve(Sigma_prior))
  mu_posterior <- Sigma_posterior %*% (t(X) %*% z / sigma2 + solve(Sigma_prior) %*% mu_prior)

  # Draw a sample for beta from the multivariate normal distribution
  beta_new <- MASS::mvrnorm(1, mu_posterior, Sigma_posterior)

  # Step 3: Sample sigma^2 from the inverse-gamma posterior
  alpha_sigma <- alpha_0 + n / 2
  beta_sigma <- beta_0 + 0.5 * sum((z - X %*% beta_new)^2)

  # Sample sigma^2 from inverse-gamma distribution
  sigma2_new <- 1 / stats::rgamma(1, shape = alpha_sigma, rate = beta_sigma)

  return(list(beta = beta_new, sigma2 = sigma2_new, z = z))
}



#' Perform one iteration of the Gibbs sampler for an OLS model
#'
#' This function performs one iteration of a Gibbs sampler for a simple OLS
#' (Ordinary Least Squares) regression model. It samples the coefficients
#' (beta) and the variance of the errors (sigma^2) given the current data and
#' parameters. Note there is currently no prior that can be passed here.
#'
#' @param y A vector of the response variable values (n x 1).
#' @param X A matrix of the explanatory variables (n x p), including the intercept.
#' @param beta A vector of current values of the regression coefficients (p x 1).
#' @param sigma2 The current value of the variance of the errors.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{beta_sample}{A vector of the sampled regression coefficients.}
#'   \item{sigma2_sample}{The sampled variance of the errors.}
#'   \item{residuals}{The residuals from the regression, i.e., y - X * beta.}
#' }
#'
#' @examples
#' # Simulated data example
#' set.seed(123)
#' n <- 100
#' p <- 2
#' X <- cbind(1, rnorm(n))  # Design matrix (including intercept term)
#' beta_true <- c(2, 3)  # True beta coefficients
#' sigma2_true <- 1  # True variance of the error
#'
#' # Generate y according to the model y = X * beta + noise
#' y <- X %*% beta_true + rnorm(n, mean = 0, sd = sqrt(sigma2_true))
#'
#' # Initial values for beta and sigma^2
#' beta_init <- rep(0, p)
#' sigma2_init <- 1
#'
#' # Run one iteration of the Gibbs sampler
#' sample_result <- gibbs_sampler_one_draw(y, X, beta_init, sigma2_init)
#'
#' # Output
#' print(sample_result$beta_sample)
#' print(sample_result$sigma2_sample)
#' print(sample_result$residuals[1:10])
#' @importFrom MASS mvrnorm
#'
#' @export
gibbs_sampler_one_draw <- function(y, X, beta, sigma2) {
  n <- nrow(X)
  p <- ncol(X)

  # Precompute some matrix operations for efficiency
  XtX_inv <- solve(t(X) %*% X)

  # 1. Sample beta | y, X, sigma^2
  beta_mean <- XtX_inv %*% t(X) %*% y
  beta_var <- sigma2 * XtX_inv
  beta_sample <- MASS::mvrnorm(1, beta_mean, beta_var)  # Draw beta from the multivariate normal distribution

  # 2. Sample sigma^2 | y, X, beta
  residuals <- y - X %*% beta_sample
  alpha <- (n / 2)
  beta_param <- sum(residuals^2) / 2
  sigma2_sample <- 1 / stats::rgamma(1, shape = alpha, rate = beta_param)  # Draw sigma^2 from the inverse-gamma

  # Return beta sample, sigma^2 sample, and residuals
  return(list(betaZ = beta_sample,
              sigmaZ2 = sigma2_sample,
              residuals = residuals))
}


#' Gibbs Sampler for Multivariate Regression with Covariance Option
#'
#' This function performs a single Gibbs sampling draw for a multivariate regression model with an option
#' to use either a dense or diagonal covariance matrix for the outcomes. Each outcome has its own set of predictors.
#'
#' The model assumes that each outcome \eqn{y_j} follows a regression model:
#' \deqn{y_j = X_j \beta_j + \epsilon_j, \quad \epsilon_j \sim \mathcal{N}(0, \Sigma)}
#' where \eqn{X_j} is the predictor matrix for outcome \eqn{y_j}, \eqn{\beta_j} is the vector of regression
#' coefficients, and \eqn{\Sigma} is the covariance matrix of the residuals. The covariance matrix can either
#' be diagonal (independent outcomes) or dense (allowing correlation between outcomes).
#'
#' @param y_list A list of outcome vectors. Each element \eqn{y_j} in the list corresponds to the outcomes for one dependent variable.
#' @param X_list A list of predictor matrices. Each element \eqn{X_j} in the list corresponds to the predictor matrix for the corresponding outcome \eqn{y_j}.
#' @param sigma_param A list of initial values for the residual variances \eqn{\sigma^2_j} for each outcome when `covariance_type = "diagonal"`, or a dense covariance matrix when `covariance_type = "dense"`.
#' @param covariance_type A string specifying the type of covariance matrix to use. Can be either `"diagonal"` (no correlation between outcomes) or `"dense"` (allows correlation between outcomes).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{beta_samples}}{A list of sampled regression coefficients \eqn{\beta_j} for each outcome.}
#'   \item{\code{sigma2_samples}}{A list of sampled residual variances \eqn{\sigma^2_j} for each outcome (diagonal case), or the sampled covariance matrix \eqn{\Sigma} (dense case).}
#'   \item{\code{residuals}}{A list of residuals for each outcome \eqn{r_j = y_j - X_j \beta_j}.}
#' }
#'
#' @examples
#' # Simulate data for two outcomes with different predictor matrices
#' set.seed(123)
#' n <- 100  # Number of observations
#' p <- 2    # Number of outcomes
#'
#' # Create predictor matrices X for each outcome
#' X1 <- cbind(1, rnorm(n))
#' X2 <- cbind(1, rnorm(n), rnorm(n))
#'
#' # True beta values
#' beta1 <- c(1, 0.5)
#' beta2 <- c(1, -0.3, 0.2)
#'
#' # Simulate outcomes
#' y1 <- X1 %*% beta1 + rnorm(n)
#' y2 <- X2 %*% beta2 + rnorm(n)
#'
#' # Prepare data for the Gibbs sampler
#' y_list <- list(y1, y2)
#' X_list <- list(X1, X2)
#' sigma_param <- list(1, 1)  # Diagonal case: Initial values for sigma^2
#'
#' # Run the Gibbs sampler with diagonal covariance
#' result_diag <-
#'     gibbs_sampler_one_draw_cov(y_list, X_list,
#'                               sigma_param,
#'                               covariance_type = "diagonal")
#' print(result_diag)
#'
#' # Simulate a dense covariance matrix
#' sigma_param_dense <- matrix(c(1, 0.5, 0.5, 1), 2, 2)  # Dense covariance matrix
#'
#' # Run the Gibbs sampler with dense covariance
#' result_dense <-
#'     gibbs_sampler_one_draw_cov(y_list, X_list,
#'                               sigma_param_dense,
#'                               covariance_type = "dense")
#' print(result_dense)
#'
#' @export
gibbs_sampler_one_draw_cov <- function(y_list, X_list, sigma_param,
                                       covariance_type = "diagonal") {
  p <- length(y_list)  # Number of outcomes
  n <- nrow(X_list[[1]])  # Number of observations (assumed to be the same for all outcomes)

  # Storage for results
  beta_samples <- vector("list", p)
  # Store residuals and predicted for all outcomes
  residuals_matrix <- predicted_matrix <- matrix(0, nrow = n, ncol = p)

  if (covariance_type == "diagonal") {
    # Diagonal case: sigma_param is a list of sigma^2 values
    sigma2_list <- sigma_param
  } else if (covariance_type == "dense") {
    # Dense case: sigma_param is the dense covariance matrix
    Sigma <- sigma_param
    Sigma_inv <- solve(Sigma)  # Inverse of the covariance matrix
  } else {
    stop("Invalid covariance type. Use 'diagonal' or 'dense'.")
  }

  # Step 1: Sample beta and sigma^2 for each outcome (or use dense covariance)
  for (j in 1:p) {
    X_j <- X_list[[j]]
    y_j <- y_list[[j]]

    # Precompute some matrix operations for efficiency
    XtX_inv_j <- solve(t(X_j) %*% X_j)

    # 1. Sample beta_j | y_j, X_j, sigma^2_j or Sigma (if dense)
    if (covariance_type == "diagonal") {
      beta_mean_j <- XtX_inv_j %*% t(X_j) %*% y_j
      beta_var_j <- sigma2_list[[j]] * XtX_inv_j
      beta_samples[[j]] <- MASS::mvrnorm(1, beta_mean_j, beta_var_j)  # Draw beta from the multivariate normal distribution
    } else if (covariance_type == "dense") {
      beta_mean_j <- XtX_inv_j %*% t(X_j) %*% y_j
      beta_var_j <- Sigma_inv[j, j] * XtX_inv_j  # Use the inverse of Sigma for variance
      beta_samples[[j]] <- MASS::mvrnorm(1, beta_mean_j, beta_var_j)  # Draw beta from the multivariate normal distribution
    }

    # 2. Compute residuals and store them
    residuals_j <- y_j - X_j %*% beta_samples[[j]]
    residuals_matrix[, j] <- residuals_j  # Store residuals for this outcome

    predicted_matrix[, j] <- X_j %*% beta_samples[[j]]

    # 3. Sample sigma^2_j | y_j, X_j, beta_j (diagonal case)
    if (covariance_type == "diagonal") {
      alpha <- n / 2
      beta_param <- sum(residuals_j^2) / 2
      sigma2_list[j] <- 1 / stats::rgamma(1, shape = alpha, rate = beta_param)  # Sample sigma^2 from inverse-gamma
    }
  }

  # Step 4: If dense covariance, sample the covariance matrix
  if (covariance_type == "dense") {
    # Degrees of freedom for inverse-Wishart
    nu <- n + p + 1
    # Scale matrix for inverse-Wishart
    S <- t(residuals_matrix) %*% residuals_matrix
    sigma_sample <- MCMCpack::riwish(nu, S)  # Sample covariance matrix from inverse-Wishart
  } else {
    sigma_sample <- sigma2_list  # For diagonal case, return the sampled sigma^2 values from diagonal
  }

  # Return results: beta samples, residuals, and covariance matrix (if dense)
  return(list(beta_samples = beta_samples,
              sigma2_samples = sigma_sample,
              residuals = as.data.frame(residuals_matrix),
              predicted = as.data.frame(predicted_matrix),
              covariance = if (covariance_type == "dense") sigma_sample else NULL))
}


