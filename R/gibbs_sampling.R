#' Gibbs sampling dispatcher for post-SCM models
#'
#' Runs an appropriate Gibbs sampler depending on the complexity level encoded in gdata$cov.
#'
#' @param gdata Output from prepare_data_general().
#' @param n_iter Number of iterations.
#' @param burn_in Burn-in iterations.
#' @param do.CF Use control function when IV blocks exist.
#' @param Z_cov_dense Dense covariance for IV first stage (requires MCMCpack if TRUE).
#' @param CF.interactions Expand CF residual interactions.
#'
#' @return Posterior draws.
#' @export
gibbs_postscm <- function(gdata,
                          n_iter = 1000,
                          burn_in = 500,
                          do.CF = TRUE,
                          Z_cov_dense = FALSE,
                          CF.interactions = TRUE) {

  lvl <- gdata$cov$second_stage

  if (is.null(lvl) || lvl == "none") {
    # Minimal: sample beta and sigma2 only (no gamma, no tau hierarchy unless you want it)
    # Here we implement a lightweight conjugate-ish Gibbs similar to your OLS sampler,
    # but respecting W and block X.
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
      # beta | sigma2
      A <- XtWX / sigma2 + diag(1e-8, p)
      b <- XtWy / sigma2
      cA <- Matrix::chol(A)
      Ainv <- Matrix::chol2inv(cA)
      mu <- Ainv %*% b
      beta <- rMVNormCovariance(1, mu = mu, Sigma = Ainv)

      # sigma2 | beta
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

  # For moderators / moderators_iv: use your full sampler (slightly refactored)
  # Here we simply call your existing gibbs_sampling() with compatible fields.
  gibbs_sampling(gdata = gdata,
                 n_iter = n_iter,
                 burn_in = burn_in,
                 run_selection_gibbs = FALSE,
                 do.CF = do.CF,
                 selection.method = "2SRI",
                 Z_cov_dense = Z_cov_dense,
                 CF.interactions = CF.interactions)
}
