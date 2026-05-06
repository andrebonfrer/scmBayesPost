# =============================================================================
#  Convergence diagnostics for gibbs_postscm() fit objects
# =============================================================================
#
#  Expected fit object structure (output of gibbs_postscm):
#    $beta_samples    [n_save x K*J0]  - unit x covariate coefficients
#    $sigma2_samples  [n_save]         - observation noise variance
#    $tau_samples     [n_save x K]     - covariate-level dispersion SDs
#    $gamma_samples   [n_save x G]     - moderator coefficients (if present)
#
#  Usage:
#    source("convergence_diagnostics.R")
#
#    # Single-chain diagnostics (fast)
#    diag <- convergence_summary(fit)
#
#    # With multi-chain R-hat (re-runs sampler n_chains times)
#    diag <- convergence_summary(fit, gdata = gdata, run_rhat = TRUE,
#                                n_iter = 2000, burn_in = 1000)
# =============================================================================


# -----------------------------------------------------------------------------
#' Extract named parameter matrix from a gibbs_postscm fit
#'
#' Combines all monitored parameters into a single draws matrix, dropping
#' beta columns by default (typically high-dimensional and not of primary
#' interest for convergence checks).
#'
#' @param fit          Output list from gibbs_postscm().
#' @param include_beta Logical. Include beta columns? Default FALSE.
#'
#' @return A numeric matrix [n_save x n_params] with column names.
#' @keywords internal
.extract_draws_matrix <- function(fit, include_beta = FALSE) {

  # Flexibly detect slot names, accepting common naming variants
  slot <- function(candidates) {
    hit <- intersect(candidates, names(fit))
    if (length(hit)) fit[[hit[1]]] else NULL
  }

  beta_draws   <- slot(c("beta_samples",   "betaSamples",   "beta"))
  gamma_draws  <- slot(c("gamma_samples",  "gammaSamples",  "gamma"))
  sigma2_draws <- slot(c("sigma2_samples", "sigma2Samples", "sigma2"))
  tau_draws    <- slot(c("tau_samples",    "tauSamples",    "tau"))

  parts <- list()

  if (include_beta && !is.null(beta_draws)) {
    b <- as.matrix(beta_draws)
    if (is.null(colnames(b)))
      colnames(b) <- paste0("beta[", seq_len(ncol(b)), "]")
    parts$beta <- b
  }

  if (!is.null(gamma_draws)) {
    g <- as.matrix(gamma_draws)
    if (is.null(colnames(g)))
      colnames(g) <- paste0("gamma[", seq_len(ncol(g)), "]")
    parts$gamma <- g
  }

  if (!is.null(sigma2_draws)) {
    s <- matrix(as.numeric(sigma2_draws), ncol = 1,
                dimnames = list(NULL, "sigma2_alpha"))
    parts$sigma2 <- s
  }

  if (!is.null(tau_draws)) {
    ta <- as.matrix(tau_draws)^2     # store variance, not SD
    if (is.null(colnames(tau_draws)))
      colnames(ta) <- paste0("tau2[", seq_len(ncol(ta)), "]")
    else
      colnames(ta) <- paste0(colnames(tau_draws), "_tau2")
    parts$tau2 <- ta
  }

  if (length(parts) == 0)
    stop(
      "No recognised parameter slots found in fit object.\n",
      "  Available names: ", paste(names(fit), collapse = ", "), "\n",
      "  Expected one or more of: beta_samples, gamma_samples, ",
      "sigma2_samples, tau_samples",
      call. = FALSE
    )

  do.call(cbind, parts)
}


# -----------------------------------------------------------------------------
#' Effective sample size for a single MCMC chain
#'
#' Uses Geyer's initial monotone sequence estimator for the autocorrelation
#' sum, which is more robust than naive truncation.
#'
#' @param x Numeric vector of post-burn-in draws.
#' @return Scalar ESS estimate, capped at length(x).
#' @keywords internal
.ess_single <- function(x) {
  n       <- length(x)
  x       <- x - mean(x)
  acf_obj <- stats::acf(x, lag.max = n - 1, plot = FALSE)$acf[-1]

  # Geyer's initial monotone sequence: sum consecutive pairs until negative
  pairs  <- acf_obj[c(TRUE, FALSE)] + acf_obj[c(FALSE, TRUE)]
  pairs  <- pairs[seq_len(floor(length(acf_obj) / 2))]
  cutoff <- which(pairs < 0)[1]
  if (is.na(cutoff)) cutoff <- length(pairs)

  tau <- 1 + 2 * sum(acf_obj[seq_len(2 * cutoff - 1)])
  min(n / max(tau, 1), n)
}


# -----------------------------------------------------------------------------
#' Gelman-Rubin R-hat for a single parameter across chains
#'
#' @param chains List of numeric vectors (one per chain), equal length.
#' @return Scalar R-hat.
#' @keywords internal
.rhat_single <- function(chains) {
  m <- length(chains)
  n <- length(chains[[1]])

  chain_means <- vapply(chains, mean, numeric(1))
  chain_vars  <- vapply(chains, stats::var, numeric(1))

  B       <- n * stats::var(chain_means)   # between-chain variance
  W       <- mean(chain_vars)              # within-chain variance

  var_hat <- (n - 1) / n * W + B / n
  sqrt(var_hat / W)
}


# =============================================================================
#' Trace plots for key parameters
#'
#' Produces trace plots for gamma, sigma2, and tau2 parameters from a
#' gibbs_postscm fit. A horizontal dashed red line marks the posterior mean.
#'
#' @param fit    Output from gibbs_postscm().
#' @param params Character vector of parameter groups to plot. One or more of
#'   "gamma", "sigma2", "tau2". Default: all three.
#' @param ask    Logical. Prompt before each page of plots when there are many
#'   parameters. Default TRUE in interactive sessions.
#'
#' @return Invisibly returns the draws matrix used for plotting.
#' @export
plot_traces <- function(fit,
                        params = c("gamma", "sigma2", "tau2"),
                        ask    = interactive()) {

  draws  <- .extract_draws_matrix(fit, include_beta = FALSE)
  params <- match.arg(params, several.ok = TRUE)

  cols_to_plot <- character(0)
  if ("gamma"  %in% params)
    cols_to_plot <- c(cols_to_plot,
                      grep("^gamma",  colnames(draws), value = TRUE))
  if ("sigma2" %in% params)
    cols_to_plot <- c(cols_to_plot,
                      grep("^sigma2", colnames(draws), value = TRUE))
  if ("tau2"   %in% params)
    cols_to_plot <- c(cols_to_plot,
                      grep("tau2$",   colnames(draws), value = TRUE))

  if (length(cols_to_plot) == 0)
    stop(
      "No columns matching params = c(",
      paste0('"', params, '"', collapse = ", "),
      ") found.\n",
      "  Available columns: ", paste(colnames(draws), collapse = ", "),
      call. = FALSE
    )

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  n_col <- min(2L, length(cols_to_plot))
  n_row <- min(ceiling(length(cols_to_plot) / n_col), 4L)
  graphics::par(mfrow = c(n_row, n_col),
                mar   = c(3, 3, 2, 1),
                mgp   = c(2, 0.6, 0),
                ask   = ask)

  iters <- seq_len(nrow(draws))

  for (nm in cols_to_plot) {
    graphics::plot(
      iters, draws[, nm],
      type = "l",
      col  = "#2166AC",
      lwd  = 0.6,
      xlab = "Iteration",
      ylab = nm,
      main = nm
    )
    graphics::abline(h   = mean(draws[, nm]),
                     col = "#D7191C",
                     lty = 2,
                     lwd = 1.2)
  }

  invisible(draws)
}


# =============================================================================
#' Effective sample sizes for all monitored parameters
#'
#' @param fit          Output from gibbs_postscm().
#' @param include_beta Logical. Include beta parameters? Default FALSE.
#' @param min_ess      ESS threshold below which parameters are flagged.
#'   Default 200.
#'
#' @return A data frame with columns: parameter, ess, adequate.
#' @export
compute_ess <- function(fit, include_beta = FALSE, min_ess = 200) {

  draws    <- .extract_draws_matrix(fit, include_beta = include_beta)
  ess_vals <- apply(draws, 2, .ess_single)

  result <- data.frame(
    parameter = names(ess_vals),
    ess       = round(ess_vals, 1),
    adequate  = ess_vals >= min_ess,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  n_low <- sum(!result$adequate)
  if (n_low > 0) {
    message(sprintf(
      "ESS < %d for %d parameter(s): %s",
      min_ess, n_low,
      paste(result$parameter[!result$adequate], collapse = ", ")
    ))
  } else {
    message(sprintf("All ESS >= %d. Adequate posterior precision.", min_ess))
  }

  result
}


# =============================================================================
#' Gelman-Rubin R-hat from multiple independent chains
#'
#' Re-runs the sampler n_chains times from different random seeds and computes
#' R-hat for all monitored parameters. Chains run in parallel via
#' parallel::mclapply on Mac/Linux, sequentially on Windows.
#'
#' @param gdata      Prepared data object passed to the sampler.
#' @param sampler_fn Function to use as the sampler. Defaults to
#'   gibbs_postscm if it can be found, otherwise must be supplied
#'   explicitly, e.g. sampler_fn = gibbs_postscm.
#' @param package    Character. Name of the package to load inside each worker
#'   so that sampler dependencies are available. E.g. package = "scmBayesPost".
#' @param n_chains   Integer. Number of independent chains. Default 4.
#' @param n_iter     Integer. Total iterations per chain.
#' @param burn_in    Integer. Burn-in iterations per chain.
#' @param control    Optional control list passed to the sampler.
#' @param seeds      Integer vector of length n_chains for reproducibility.
#'   Randomly generated if NULL.
#' @param threshold  Scalar R-hat threshold for flagging. Default 1.1.
#'
#' @return A data frame with columns: parameter, rhat, converged.
#'   The list of fitted chains is attached as attr(result, "fits") and
#'   the seeds used as attr(result, "seeds").
#' @export
compute_rhat <- function(gdata,
                         sampler_fn = NULL,
                         package    = NULL,
                         n_chains   = 4L,
                         n_iter     = 1000L,
                         burn_in    = 500L,
                         control    = NULL,
                         seeds      = NULL,
                         threshold  = 1.1) {

  # Resolve sampler: explicit fn > name lookup > informative error
  if (is.null(sampler_fn)) {
    if (exists("gibbs_postscm", mode = "function")) {
      sampler_fn <- get("gibbs_postscm")
    } else {
      stop(
        "Cannot find gibbs_postscm(). ",
        "Pass it explicitly: sampler_fn = gibbs_postscm",
        call. = FALSE
      )
    }
  }
  if (!is.function(sampler_fn))
    stop("'sampler_fn' must be a function.", call. = FALSE)

  if (is.null(seeds))
    seeds <- sample.int(1e6, n_chains)

  if (length(seeds) != n_chains)
    stop("'seeds' must have length equal to 'n_chains'.", call. = FALSE)

  n_cores <- if (.Platform$OS.type == "windows") 1L else
    max(1L, parallel::detectCores(logical = FALSE) - 1L)
  n_cores <- min(n_cores, n_chains)

  message(sprintf("Running %d chains on %d core(s) ...", n_chains, n_cores))

  fits <- parallel::mclapply(seq_len(n_chains), function(i) {
    if (!is.null(package))
      library(package, character.only = TRUE)
    set.seed(seeds[i])
    sampler_fn(
      gdata   = gdata,
      n_iter  = n_iter,
      burn_in = burn_in,
      control = control
    )
  }, mc.cores = n_cores)

  # Validate all chains returned usable results
  failed <- vapply(fits, function(f) is.null(names(f)), logical(1))
  if (any(failed))
    stop(sprintf(
      "%d chain(s) returned NULL or unrecognised output. Check gibbs_postscm().",
      sum(failed)
    ), call. = FALSE)

  draws_list  <- lapply(fits, .extract_draws_matrix, include_beta = FALSE)
  param_names <- colnames(draws_list[[1]])

  rhat_vals <- vapply(param_names, function(nm) {
    chains <- lapply(draws_list, function(d) d[, nm])
    .rhat_single(chains)
  }, numeric(1))

  result <- data.frame(
    parameter = param_names,
    rhat      = round(rhat_vals, 4),
    converged = rhat_vals < threshold,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  n_fail <- sum(!result$converged)
  if (n_fail > 0) {
    message(sprintf(
      "R-hat >= %.2f for %d parameter(s): %s",
      threshold, n_fail,
      paste(result$parameter[!result$converged], collapse = ", ")
    ))
  } else {
    message(sprintf(
      "All R-hat < %.2f. Chains have converged to a common distribution.",
      threshold
    ))
  }

  attr(result, "fits")  <- fits
  attr(result, "seeds") <- seeds
  result
}


# =============================================================================
#' Full convergence summary for a gibbs_postscm fit
#'
#' Combines trace plots, ESS, posterior summaries, and optionally Gelman-Rubin
#' R-hat into a single printed report. Mirrors the convergence assessment
#' described in the paper.
#'
#' @param fit          Output from gibbs_postscm(). Required.
#' @param gdata        Prepared data object. Required only if run_rhat = TRUE.
#' @param sampler_fn   The sampler function. Passed to compute_rhat(). Defaults
#'   to gibbs_postscm if found, otherwise supply explicitly.
#' @param run_rhat     Logical. Run multi-chain R-hat? Default FALSE.
#' @param n_chains     Passed to compute_rhat(). Default 4.
#' @param n_iter       Passed to compute_rhat(). Default 1000.
#' @param burn_in      Passed to compute_rhat(). Default 500.
#' @param control      Passed to compute_rhat().
#' @param min_ess      ESS adequacy threshold. Default 200.
#' @param rhat_thresh  R-hat convergence threshold. Default 1.1.
#' @param show_traces  Logical. Display trace plots? Default TRUE.
#'
#' @return Invisibly returns a list with elements:
#'   $ess               data frame from compute_ess()
#'   $posterior_summary data frame of means and 95% credible intervals
#'   $rhat              data frame from compute_rhat() (only if run_rhat = TRUE)
#' @export
convergence_summary <- function(fit,
                                gdata       = NULL,
                                sampler_fn  = NULL,
                                package     = NULL,
                                run_rhat    = FALSE,
                                n_chains    = 4L,
                                n_iter      = 1000L,
                                burn_in     = 500L,
                                control     = NULL,
                                min_ess     = 200,
                                rhat_thresh = 1.1,
                                show_traces = TRUE) {

  cat("=======================================================\n")
  cat("  Convergence Diagnostics: gibbs_postscm\n")
  cat("=======================================================\n\n")

  # Detect n_save robustly across slot name variants
  n_save_slot <- intersect(
    c("beta_samples", "betaSamples", "beta",
      "gamma_samples", "gammaSamples",
      "sigma2_samples", "sigma2Samples",
      "tau_samples", "tauSamples"),
    names(fit)
  )
  if (length(n_save_slot)) {
    obj    <- fit[[n_save_slot[1]]]
    n_save <- if (is.matrix(obj)) nrow(obj) else length(obj)
    cat(sprintf("  Post-burn-in draws retained : %d\n\n", n_save))
  }

  # --- Trace plots ------------------------------------------------------------
  if (show_traces) {
    cat("Generating trace plots ...\n")
    plot_traces(fit)
  }

  # --- ESS --------------------------------------------------------------------
  cat("\n--- Effective Sample Sizes (ESS) ---\n")
  ess_df <- compute_ess(fit, min_ess = min_ess)
  print(ess_df, row.names = FALSE)

  # --- Posterior summaries ----------------------------------------------------
  cat("\n--- Posterior Summaries (mean and 95% credible interval) ---\n")
  draws <- .extract_draws_matrix(fit, include_beta = FALSE)
  summ  <- t(apply(draws, 2, function(x) {
    c(mean  = mean(x),
      sd    = stats::sd(x),
      q2.5  = stats::quantile(x, 0.025),
      q97.5 = stats::quantile(x, 0.975))
  }))
  summ_df <- as.data.frame(round(summ, 4))
  summ_df <- cbind(parameter = rownames(summ_df), summ_df)
  rownames(summ_df) <- NULL
  print(summ_df, row.names = FALSE)

  out <- list(ess = ess_df, posterior_summary = summ_df)

  # --- R-hat (optional) -------------------------------------------------------
  if (run_rhat) {
    if (is.null(gdata))
      stop("'gdata' must be supplied when run_rhat = TRUE.", call. = FALSE)

    cat("\n--- Gelman-Rubin R-hat (multi-chain) ---\n")
    rhat_df <- compute_rhat(
      gdata      = gdata,
      sampler_fn = sampler_fn,
      package    = package,
      n_chains   = n_chains,
      n_iter     = n_iter,
      burn_in    = burn_in,
      control    = control,
      threshold  = rhat_thresh
    )
    print(rhat_df, row.names = FALSE)
    out$rhat <- rhat_df
  }

  cat("\n=======================================================\n")
  invisible(out)
}
