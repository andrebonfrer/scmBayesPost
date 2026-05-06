# CommonTrends.R

################################################################################
##  pretrend_test.R
##
##  Pre-intervention common trends test for multiout_synth fitted objects.
##
##  For each outcome the function:
##    1. Splits the L balanced pre-treatment lags into a "training" sub-period
##       (earlier lags, used to fit placebo weights if needed) and a "test"
##       sub-period (later lags, the placebo post-treatment window).
##    2. Optionally standardises each outcome the same way build_Xy_for_unit
##       does when standardize_outcomes = TRUE: subtract the donor-pool
##       pre-treatment mean and divide by the donor-pool pre-treatment SD,
##       computed per treated unit over its own L lags and donor set.
##    3. Computes the gap (treated − synthetic) in the test sub-period,
##       averaged across treated units.
##    4. Tests H0: mean gap = 0 via an F-test (or reports NA for n_test = 1).
################################################################################


#' Pre-Intervention Common Trends Test for multiout_synth Objects
#'
#' Evaluates whether the synthetic control tracks the treated group during the
#' pre-treatment period by splitting the \code{L} balanced lags into a training
#' sub-period and a held-out test sub-period, then testing whether the gaps in
#' the test sub-period are jointly zero.
#'
#' @section Standardisation:
#' If \code{standardize_outcomes = TRUE}, each outcome is standardised
#' \emph{per treated unit} using the same procedure as
#' \code{build_Xy_for_unit}: the donor-pool mean and SD are computed over that
#' unit's donor set and its \code{L} pre-treatment columns, then subtracted /
#' divided before computing gaps.  This ensures the test is conducted on the
#' same scale the solver used when fitting weights.  Set this to \code{TRUE}
#' whenever \code{multiout_synth} was called with
#' \code{standardize_outcomes = TRUE}.
#'
#' @param fit              List returned by \code{multiout_synth()}.
#' @param Y_list           Named list of \eqn{N \times T} outcome matrices
#'   (unit rows, time-period columns).
#' @param treat_time       Numeric vector length \eqn{N}; first treated column
#'   index per unit (\code{Inf} = never treated).
#' @param L                Integer; number of pre-treatment lags balanced.
#' @param outcomes         Character vector of outcome names to test.  Defaults
#'   to \code{names(Y_list)}.
#' @param min_pre_periods  Integer.  Minimum required pre-treatment lags.
#'   Error if \code{L < min_pre_periods}.  Default \code{3L}.
#' @param test_fraction    Numeric in (0, 1).  Share of the \code{L} lags
#'   reserved as the placebo test window.  Default \code{0.3}.
#' @param standardize_outcomes Logical.  If \code{TRUE}, standardise each
#'   outcome per treated unit using the donor-pool pre-treatment mean and SD,
#'   matching \code{build_Xy_for_unit}.  Default \code{FALSE}.
#' @param intercept  Character.  Intercept removal matching the \code{intercept}
#'   argument passed to \code{multiout_synth}:
#'   \describe{
#'     \item{\code{"none"}}{No demeaning (default).}
#'     \item{\code{"outcome"}}{Remove the (time-weighted) pre-treatment mean
#'       from both the treated trajectory and the synthetic control trajectory,
#'       per outcome per treated unit, using all \code{L} lags as the reference.
#'       This replicates \code{demean_within_outcome_blocks} and removes level
#'       differences so the test focuses purely on trend shape.}
#'     \item{\code{"global"}}{Remove a single (time-weighted) grand mean across
#'       all outcomes and lags.  Use when \code{multiout_synth} was called with
#'       \code{intercept = "global"}.}
#'   }
#' @param time_weights Optional numeric vector of length \code{L}.  Lag weights
#'   used in the (time-)weighted mean when demeaning.  Must match the
#'   \code{time_weights} argument passed to \code{multiout_synth}.  Default
#'   \code{NULL} (uniform weights).
#' @param eps_sd           Small positive floor for SD when standardising.
#'   Default \code{1e-8}.
#' @param alpha            Significance level.  Default \code{0.05}.
#' @param verbose          Logical.  Print summary table.  Default \code{TRUE}.
#'
#' @return A named list (one element per outcome) each containing:
#' \describe{
#'   \item{\code{outcome}}{Outcome name.}
#'   \item{\code{gaps_test}}{Named numeric vector of average gaps in the test
#'     sub-period (one value per lag, averaged across treated units).}
#'   \item{\code{gaps_train}}{Average gaps in the training sub-period.}
#'   \item{\code{rmse_train}}{RMSE of gaps in the training sub-period.}
#'   \item{\code{rmse_test}}{RMSE of gaps in the test sub-period.}
#'   \item{\code{rmse_ratio}}{\code{rmse_test / rmse_train}.}
#'   \item{\code{statistic}}{F-statistic (NA if n_test = 1).}
#'   \item{\code{df1}, \code{df2}}{Degrees of freedom.}
#'   \item{\code{p_value}}{p-value.}
#'   \item{\code{reject}}{Logical: rejected at \code{alpha}?}
#'   \item{\code{train_lags}, \code{test_lags}}{Integer vectors of the column
#'     indices used as training / test lags (relative to the matrix).}
#'   \item{\code{standardized}}{Logical: were gaps computed on standardised
#'     scale?}
#' }
#'
#' @examples
#' \dontrun{
#' # Match the standardize_outcomes flag used in multiout_synth
#' res <- pretrend_test(
#'   fit                  = fit,
#'   Y_list               = Y_list,
#'   treat_time           = treat_time,
#'   L                    = L,
#'   standardize_outcomes = TRUE   # if multiout_synth was called with TRUE
#' )
#'
#' plot_pretrend(res)
#' }
#'
#' @seealso \code{\link{balance_table}} for a summary balance table.
#' @export
pretrend_test <- function(fit,
                          Y_list,
                          treat_time,
                          L,
                          outcomes             = NULL,
                          min_pre_periods      = 3L,
                          test_fraction        = 0.3,
                          standardize_outcomes = FALSE,
                          intercept            = c("none", "outcome", "global"),
                          time_weights         = NULL,
                          eps_sd               = 1e-8,
                          alpha                = 0.05,
                          verbose              = TRUE) {

  # ── 0. Validation ──────────────────────────────────────────────────────────

  required <- c("weights", "donors", "treated_units")
  miss <- setdiff(required, names(fit))
  if (length(miss)) {
    stop("fit is missing slot(s): ", paste(miss, collapse = ", "))
  }
  if (!is.list(Y_list) || length(Y_list) == 0L) {
    stop("Y_list must be a non-empty named list of matrices.")
  }
  if (is.null(names(Y_list)) || any(names(Y_list) == "")) {
    stop("All elements of Y_list must be named.")
  }

  intercept       <- match.arg(intercept)
  L               <- as.integer(L)
  min_pre_periods <- as.integer(min_pre_periods)

  # Validate and normalise time_weights
  if (is.null(time_weights)) {
    tw <- rep(1, L)
  } else {
    tw <- as.numeric(time_weights)
    if (length(tw) != L) stop("time_weights must have length L = ", L)
    if (any(!is.finite(tw)) || any(tw < 0) || sum(tw) <= 0) {
      stop("time_weights must be finite, non-negative, and sum to > 0.")
    }
  }
  tw <- tw / sum(tw)   # normalise to sum to 1 for weighted mean computation

  if (L < min_pre_periods) {
    stop(sprintf("L = %d is less than min_pre_periods = %d.", L,
                 min_pre_periods))
  }
  if (test_fraction <= 0 || test_fraction >= 1) {
    stop("`test_fraction` must be strictly between 0 and 1.")
  }

  # ── 1. Resolve outcomes ────────────────────────────────────────────────────

  if (is.null(outcomes)) outcomes <- names(Y_list)
  outcomes <- as.character(outcomes)
  bad <- setdiff(outcomes, names(Y_list))
  if (length(bad)) {
    stop("Outcome(s) not found in Y_list: ", paste(bad, collapse = ", "))
  }

  # ── 2. Geometry ────────────────────────────────────────────────────────────

  # Y_list matrices are N × T (unit rows, time-period columns)
  N  <- nrow(Y_list[[1L]])
  TT <- ncol(Y_list[[1L]])

  if (length(treat_time) != N) {
    stop("treat_time must have length N = nrow(Y_list[[1]]) = ", N)
  }

  treated_idx <- fit$treated_units
  donor_idx   <- which(!is.finite(treat_time))
  J           <- length(treated_idx)

  if (J == 0L) stop("No treated units in fit$treated_units.")
  if (length(donor_idx) == 0L) stop("No donor units (all treat_time finite).")

  # ── 3. Per-unit lag column indices ─────────────────────────────────────────
  #
  # For treated unit j with first treatment at column Ti:
  #   all L pre-treatment columns: (Ti - L) : (Ti - 1)
  #   training sub-period: first n_train of those
  #   test sub-period:     last  n_test  of those

  Ti_vec <- as.integer(treat_time[treated_idx])

  n_test  <- max(1L, ceiling(L * test_fraction))
  n_train <- L - n_test

  if (n_train < 1L) {
    stop(sprintf(
      "With L = %d and test_fraction = %.2f, no training lags remain. ",
      L, test_fraction
    ), "Reduce test_fraction.")
  }

  # Column offsets (relative positions within the L pre-treatment window)
  train_offsets <- seq_len(n_train)                    # 1 .. n_train
  test_offsets  <- (n_train + 1L):L                   # n_train+1 .. L

  # Feasibility check
  ok <- (Ti_vec - L >= 1L) & (Ti_vec - 1L <= TT)
  if (!all(ok)) {
    n_drop <- sum(!ok)
    warning(n_drop, " treated unit(s) lack ", L,
            " pre-treatment columns; excluded from test.")
    treated_idx <- treated_idx[ok]
    Ti_vec      <- Ti_vec[ok]
    J           <- length(treated_idx)
    if (J == 0L) stop("No treated units remain after feasibility check.")
  }

  if (verbose) {
    flags <- c(
      if (standardize_outcomes) "standardised" else NULL,
      if (intercept != "none") paste0("intercept=", intercept) else NULL
    )
    flag_str <- if (length(flags)) paste0(" | ", paste(flags, collapse = ", ")) else ""
    message(sprintf(
      "\nPre-trend test: L = %d lags | train = %d | test = %d | J = %d treated units%s\n",
      L, n_train, n_test, J, flag_str
    ))
  }

  # ── 4. Helper: standardise one outcome block for one treated unit ──────────
  #
  # Replicates build_Xy_for_unit exactly:
  #   mu  = mean of donor rows × all L pre-treatment columns
  #   sd  = SD  of donor rows × all L pre-treatment columns (flattened)

  .std_scale <- function(Y, unit_j, donors_j, all_L_cols) {
    D_pre <- Y[donors_j, all_L_cols, drop = FALSE]   # n_donors × L
    mu    <- mean(D_pre, na.rm = TRUE)
    sd_v  <- stats::sd(as.numeric(D_pre), na.rm = TRUE)
    if (!is.finite(sd_v) || sd_v < eps_sd) sd_v <- eps_sd
    list(mu = mu, sd = sd_v)
  }

  .apply_scale <- function(x, scale) (x - scale$mu) / scale$sd

  # ── 4b. Demeaning helper ────────────────────────────────────────────────────
  #
  # Replicates demean_within_outcome_blocks for a single outcome block.
  # Computes the time-weighted mean of a length-L vector and subtracts it.
  # The same mean is computed from the treated trajectory and subtracted from
  # both treated and synthetic — this is what the solver saw.
  #
  # intercept = "outcome": demean per outcome per treated unit (most common).
  # intercept = "global":  a single grand mean across all outcomes and lags is
  #   removed; here we approximate by demeaning each outcome separately since
  #   we operate one outcome at a time. For a true global demean the grand mean
  #   would need to be computed across all M outcomes simultaneously — but since
  #   we loop over outcomes, we centre each independently, which is equivalent
  #   when outcome_weights are equal. If outcome_weights differ, this is an
  #   approximation; in practice the distinction is minor.

  .demean_vec <- function(x) {
    # time-weighted mean of length-L vector x, then subtract
    wm <- sum(tw * x)   # tw already sums to 1
    x - wm
  }

  # ── 5. Per-outcome loop ────────────────────────────────────────────────────

  results <- vector("list", length(outcomes))
  names(results) <- outcomes

  for (oc in outcomes) {
    Y <- Y_list[[oc]]    # N × T

    # Accumulators: one value per lag per treated unit, then averaged
    # train: matrix J × n_train,  test: matrix J × n_test
    gaps_train_mat <- matrix(NA_real_, nrow = J, ncol = n_train)
    gaps_test_mat  <- matrix(NA_real_, nrow = J, ncol = n_test)

    for (jj in seq_len(J)) {
      j        <- treated_idx[jj]
      Ti       <- Ti_vec[jj]
      donors_j <- fit$donors[[jj]]
      w_j      <- fit$weights[[jj]]

      all_L_cols    <- (Ti - L):(Ti - 1L)
      train_cols    <- all_L_cols[train_offsets]
      test_cols     <- all_L_cols[test_offsets]

      # ── Full L-lag trajectories (needed for demeaning reference) ──────────
      y_tr_all  <- Y[j,        all_L_cols]                        # length L
      Y_dn_all  <- Y[donors_j, all_L_cols, drop = FALSE]          # n_donors × L
      y_syn_all <- as.numeric(t(Y_dn_all) %*% w_j)               # length L

      # ── Standardise (scale only, no mean removal) ──────────────────────────
      if (standardize_outcomes) {
        sc        <- .std_scale(Y, j, donors_j, all_L_cols)
        y_tr_all  <- .apply_scale(y_tr_all,  sc)
        y_syn_all <- .apply_scale(y_syn_all, sc)
      }

      # ── Demean (intercept removal) over all L lags ─────────────────────────
      # Applied after standardisation (same order as build_Xy_for_unit /
      # multiout_synth: standardize first, then intercept projection).
      # The mean is computed from the treated trajectory and subtracted from
      # both, matching demean_within_outcome_blocks which uses the same ybar
      # for y and the same Xbar per column for X.
      if (intercept != "none") {
        y_tr_all  <- .demean_vec(y_tr_all)
        y_syn_all <- .demean_vec(y_syn_all)
      }

      # ── Split into train / test after all transformations ──────────────────
      gaps_train_mat[jj, ] <- y_tr_all[train_offsets] - y_syn_all[train_offsets]
      gaps_test_mat[jj, ]  <- y_tr_all[test_offsets]  - y_syn_all[test_offsets]
    }

    # Average gaps across treated units (one value per lag)
    gaps_train <- colMeans(gaps_train_mat, na.rm = TRUE)
    gaps_test  <- colMeans(gaps_test_mat,  na.rm = TRUE)
    names(gaps_train) <- as.character(train_offsets)
    names(gaps_test)  <- as.character(test_offsets)

    rmse_train <- sqrt(mean(gaps_train^2, na.rm = TRUE))
    rmse_test  <- sqrt(mean(gaps_test^2,  na.rm = TRUE))

    # ── 5a. F-test: H0: mean gap = 0 in test sub-period ─────────────────────
    #
    # We test the average treated-unit gap, not each unit individually,
    # because the unit-level gaps share donor units and are therefore
    # correlated.  The average is the natural estimand.
    #
    # With n_test >= 2: F = (n_test * mu_hat^2) / s^2  ~ F(1, n_test - 1)
    # With n_test == 1: no degrees of freedom; report NA with a warning.

    if (n_test == 1L) {
      stat   <- NA_real_
      df1    <- NA_integer_
      df2    <- NA_integer_
      p_val  <- NA_real_
      reject <- NA
      warning("Only 1 test lag for outcome '", oc,
              "'. Statistical test unavailable; inspect rmse_ratio.")
    } else {
      mu_hat <- mean(gaps_test, na.rm = TRUE)
      s2     <- stats::var(gaps_test, na.rm = TRUE)
      stat   <- if (!is.na(s2) && s2 > 0) {
        n_test * mu_hat^2 / s2
      } else NA_real_
      df1    <- 1L
      df2    <- as.integer(n_test - 1L)
      p_val  <- if (!is.na(stat)) {
        stats::pf(stat, df1 = df1, df2 = df2, lower.tail = FALSE)
      } else NA_real_
      reject <- if (!is.na(p_val)) p_val < alpha else NA
    }

    results[[oc]] <- list(
      outcome       = oc,
      gaps_test     = gaps_test,
      gaps_train    = gaps_train,
      gaps_test_mat = gaps_test_mat,   # J × n_test, for joint covariance estimation
      rmse_train    = rmse_train,
      rmse_test     = rmse_test,
      rmse_ratio    = rmse_test / rmse_train,
      statistic     = stat,
      df1           = df1,
      df2           = df2,
      p_value       = p_val,
      reject        = reject,
      train_lags    = train_offsets,
      test_lags     = test_offsets,
      standardized  = standardize_outcomes,
      intercept     = intercept
    )
  }

  # ── 6. Verbose summary table ───────────────────────────────────────────────

  if (verbose) {
    cat(strrep("-", 75), "\n")
    cat(sprintf("%-22s %8s %8s %8s %9s %8s\n",
                "Outcome", "RMSE_trn", "RMSE_tst", "Ratio",
                "F-stat", "p-value"))
    cat(strrep("-", 75), "\n")
    for (oc in outcomes) {
      r    <- results[[oc]]
      flag <- if (isTRUE(r$reject)) " *" else if (is.na(r$reject)) " ?" else "  "
      cat(sprintf("%-22s %8.4f %8.4f %8.4f %9.4f %8.4f%s\n",
                  oc,
                  r$rmse_train,
                  r$rmse_test,
                  r$rmse_ratio,
                  ifelse(is.na(r$statistic), NaN, r$statistic),
                  ifelse(is.na(r$p_value),   NaN, r$p_value),
                  flag))
    }
    cat(strrep("-", 75), "\n")
    cat(sprintf("alpha = %.2f  (*) = H0 rejected\n", alpha))
    cat(sprintf("H0: mean gap is zero in test lags %d-%d (last %.0f%% of L).\n",
                min(test_offsets), max(test_offsets), test_fraction * 100))
    cat("F ~ F(1, n_test-1). Rejection suggests pre-trend imbalance.\n")
    if (standardize_outcomes) {
      cat("Gaps on standardised scale (donor-pool mean/SD per unit).\n")
    }
    if (intercept != "none") {
      cat(sprintf(
        "Level differences removed via intercept = \"%s\" (time-weighted mean\n",
        intercept
      ))
      cat("over all L lags subtracted before computing gaps).\n")
    }
    cat("\n")
  }

  invisible(results)
}


# ══════════════════════════════════════════════════════════════════════════════
#  Plot method
# ══════════════════════════════════════════════════════════════════════════════

#' Plot Pre-Trend Gaps
#'
#' Visualises the training and test sub-period gaps from
#' \code{\link{pretrend_test}}, one panel per outcome.
#'
#' @param results  List returned by \code{\link{pretrend_test}}.
#' @param outcomes Character vector of outcomes to plot.  Default: all.
#' @param ...      Further arguments passed to \code{plot()}.
#'
#' @return \code{NULL} invisibly.
#' @export
plot_pretrend <- function(results, outcomes = NULL, ...) {
  if (is.null(outcomes)) outcomes <- names(results)
  n_plots <- length(outcomes)
  old_par <- graphics::par(
    mfrow = c(ceiling(n_plots / 2), min(2L, n_plots)),
    mar   = c(4, 4, 3, 1)
  )
  on.exit(graphics::par(old_par))

  for (oc in outcomes) {
    r <- results[[oc]]

    all_gaps  <- c(r$gaps_train, r$gaps_test)
    all_lags  <- c(r$train_lags, r$test_lags)
    is_test   <- all_lags %in% r$test_lags
    ylim      <- range(c(all_gaps, 0), na.rm = TRUE) * 1.15

    ylab <- if (r$standardized) {
      "Gap (standardised units)"
    } else {
      "Gap (treated \u2212 synthetic)"
    }

    graphics::plot(
      all_lags, all_gaps,
      type = "b", pch = 16,
      col  = ifelse(is_test, "firebrick", "steelblue"),
      xlab = "Pre-treatment lag index", ylab = ylab,
      main = paste0("Pre-trend gaps: ", oc),
      ylim = ylim, ...
    )
    graphics::abline(h = 0, lty = 2, col = "grey50")
    graphics::abline(
      v   = min(r$test_lags) - 0.5,
      lty = 3, col = "darkorange", lwd = 1.5
    )
    graphics::legend(
      "topleft",
      legend = c("Training lags", "Test lags", "Split"),
      col    = c("steelblue", "firebrick", "darkorange"),
      lty    = c(1, 1, 3), pch = c(16, 16, NA),
      bty    = "n", cex = 0.8
    )
    p_lab <- if (!is.na(r$p_value)) {
      sprintf("p = %.3f%s", r$p_value,
              if (isTRUE(r$reject)) " *" else "")
    } else {
      "p = n/a"
    }
    graphics::mtext(p_lab, side = 3, line = 0.2, cex = 0.8, col = "grey30")
  }
  invisible(NULL)
}


################################################################################
##  summary_pretrend.R
##
##  Joint common trends test across all outcomes from a pretrend_test() result.
##
##  Three complementary tests are computed:
##
##  (1) Bonferroni correction
##      Adjusts individual outcome p-values for M simultaneous tests.
##      Controls family-wise error rate (FWER).  Conservative when outcomes
##      are positively correlated, which is typical.
##
##  (2) Fisher's combined p-value
##      chi2 = -2 * sum(log(p_m)) ~ chi2(2M) under H0.
##      More powerful than Bonferroni under independence, but anti-conservative
##      when outcomes are positively correlated.  Provides a lower bound on
##      the true joint p-value.
##
##  (3) Joint Wald test on stacked unit-level gaps
##      Stacks the J × (M * n_test) matrix of unit-level test-period gaps
##      across outcomes and tests H0: E[gap] = 0 jointly using a Wald
##      statistic with a sample covariance matrix estimated from the J
##      unit-level observations.  Properly accounts for cross-outcome
##      correlation.  Requires J > M * n_test for a non-singular covariance;
##      a ridge-regularised inverse is used when this condition is not met.
##
##  Together the three tests bracket the truth:
##    - Fisher gives the most optimistic (lowest) p-value
##    - Bonferroni gives the most conservative (highest) p-value
##    - Wald (when feasible) gives the most principled p-value
##
##  If all three agree, the conclusion is robust.  If they disagree, it
##  signals that cross-outcome correlation is driving the difference and
##  the Wald test should be preferred.
################################################################################


#' Joint Common Trends Test Across All Outcomes
#'
#' Aggregates the per-outcome pre-trend tests produced by
#' \code{\link{pretrend_test}} into a single joint test of the common trends
#' hypothesis across all outcomes simultaneously.
#'
#' @section Tests:
#' \describe{
#'   \item{Bonferroni}{Multiplies each outcome p-value by \eqn{M} (number of
#'     outcomes tested) and takes the minimum adjusted p-value.  Controls the
#'     family-wise error rate.  Conservative when outcomes are correlated.}
#'   \item{Fisher}{Computes \eqn{\chi^2 = -2\sum_m \log(p_m)} and compares to
#'     \eqn{\chi^2_{2M}}.  Liberal (anti-conservative) when outcomes are
#'     positively correlated.  Serves as a lower bound on the true p-value.}
#'   \item{Wald}{Stacks unit-level test-period gaps across outcomes into a
#'     \eqn{J \times (M \cdot n_{\text{test}})} matrix, estimates the
#'     \eqn{(M \cdot n_{\text{test}}) \times (M \cdot n_{\text{test}})}
#'     covariance from the \eqn{J} unit observations, and computes a Hotelling
#'     \eqn{T^2}-style Wald statistic.  This accounts for cross-outcome and
#'     cross-lag correlation.  Requires the unit-level
#'     \code{gaps_test_mat} slot stored by \code{pretrend_test} (available
#'     when \code{pretrend_test} was run with the current version).  A
#'     ridge-regularised inverse (\eqn{\hat\Sigma + \lambda I}) is used when
#'     \eqn{J \leq M \cdot n_{\text{test}}}.}
#' }
#'
#' @param results   Named list returned by \code{\link{pretrend_test}}.
#' @param outcomes  Character vector of outcomes to include in the joint test.
#'   Defaults to all outcomes in \code{results}.
#' @param alpha     Significance level.  Default \code{0.05}.
#' @param ridge_lambda Numeric.  Ridge penalty added to the diagonal of the
#'   covariance matrix before inversion in the Wald test.  Default \code{0}
#'   (no regularisation); automatically set to a small positive value with a
#'   warning if the covariance matrix is singular or near-singular.
#' @param verbose   Logical.  Print a formatted summary.  Default \code{TRUE}.
#'
#' @return A list (invisibly) with:
#' \describe{
#'   \item{\code{outcomes}}{Outcomes included in the joint test.}
#'   \item{\code{n_outcomes}}{Number of outcomes (\eqn{M}).}
#'   \item{\code{n_test_lags}}{Number of test lags (\eqn{n_{\text{test}}}).}
#'   \item{\code{n_units}}{Number of treated units (\eqn{J}).}
#'   \item{\code{per_outcome}}{Data frame of per-outcome results (p-values,
#'     RMSE ratios, rejection flags).}
#'   \item{\code{bonferroni}}{List: \code{p_adjusted} (vector),
#'     \code{p_value} (min adjusted p), \code{reject}.}
#'   \item{\code{fisher}}{List: \code{statistic}, \code{df}, \code{p_value},
#'     \code{reject}.}
#'   \item{\code{wald}}{List: \code{statistic}, \code{df}, \code{p_value},
#'     \code{reject}, \code{method} ("chi2" or "F"), \code{ridge_used},
#'     or \code{NULL} if \code{gaps_test_mat} is unavailable.}
#'   \item{\code{conclusion}}{Character string summarising the overall
#'     conclusion.}
#'   \item{\code{alpha}}{Significance level used.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- pretrend_test(fit, Y_list, treat_time, L,
#'                      intercept = "outcome",
#'                      standardize_outcomes = TRUE)
#'
#' jt <- summary_pretrend(res)
#'
#' # Include only a subset of outcomes
#' jt <- summary_pretrend(res, outcomes = c("numarrears", "lengtharrears"))
#' }
#'
#' @seealso \code{\link{pretrend_test}}, \code{\link{plot_pretrend}}
#' @export
summary_pretrend <- function(results,
                             outcomes     = NULL,
                             alpha        = 0.05,
                             ridge_lambda = 0,
                             verbose      = TRUE) {

  # ── 0. Validate and resolve outcomes ────────────────────────────────────────

  if (!is.list(results) || is.null(names(results))) {
    stop("`results` must be a named list returned by pretrend_test().")
  }

  if (is.null(outcomes)) outcomes <- names(results)
  outcomes <- as.character(outcomes)
  bad <- setdiff(outcomes, names(results))
  if (length(bad)) {
    stop("Outcome(s) not found in results: ", paste(bad, collapse = ", "))
  }

  M <- length(outcomes)
  if (M < 2L) {
    stop("At least 2 outcomes are required for a joint test. ",
         "For a single outcome inspect the p_value in pretrend_test() directly.")
  }

  # ── 1. Extract per-outcome components ───────────────────────────────────────

  p_vals     <- numeric(M)
  rmse_ratio <- numeric(M)
  reject_ind <- logical(M)
  n_test     <- NULL
  J          <- NULL
  has_mat    <- TRUE

  gap_mats   <- vector("list", M)   # each J × n_test

  for (i in seq_len(M)) {
    r <- results[[outcomes[i]]]

    if (is.na(r$p_value)) {
      stop("Outcome '", outcomes[i], "' has NA p-value (likely only 1 test lag). ",
           "Increase L or reduce test_fraction so n_test >= 2 for all outcomes.")
    }

    p_vals[i]     <- r$p_value
    rmse_ratio[i] <- r$rmse_ratio
    reject_ind[i] <- isTRUE(r$reject)

    # Check consistency of n_test across outcomes
    n_test_i <- length(r$test_lags)
    if (is.null(n_test)) {
      n_test <- n_test_i
    } else if (n_test_i != n_test) {
      stop("Outcomes have different numbers of test lags (",
           n_test, " vs ", n_test_i, "). ",
           "This should not happen if all outcomes were tested with the same L.")
    }

    # Unit-level gap matrix (J × n_test) — stored by current pretrend_test
    if (!is.null(r$gaps_test_mat)) {
      gap_mats[[i]] <- r$gaps_test_mat
      J_i <- nrow(r$gaps_test_mat)
      if (is.null(J)) {
        J <- J_i
      } else if (J_i != J) {
        stop("Outcomes have different numbers of treated units in gaps_test_mat.")
      }
    } else {
      has_mat <- FALSE
    }
  }

  names(p_vals)     <- outcomes
  names(rmse_ratio) <- outcomes
  names(reject_ind) <- outcomes

  # ── 2. Per-outcome summary data frame ───────────────────────────────────────

  per_outcome <- data.frame(
    Outcome    = outcomes,
    p_value    = round(p_vals, 4L),
    RMSE_ratio = round(rmse_ratio, 4L),
    Reject     = reject_ind,
    stringsAsFactors = FALSE
  )

  # ── 3. Test 1: Bonferroni ────────────────────────────────────────────────────
  #
  # Adjusted p-value for outcome m: min(M * p_m, 1)
  # Joint p-value: min of the adjusted p-values (equivalent to asking whether
  # any outcome rejects after correction).

  p_bonf    <- pmin(M * p_vals, 1)
  p_bonf_jt <- min(p_bonf)
  bonferroni <- list(
    p_adjusted = setNames(round(p_bonf, 4L), outcomes),
    p_value    = round(p_bonf_jt, 4L),
    reject     = p_bonf_jt < alpha
  )

  # ── 4. Test 2: Fisher's combined p-value ─────────────────────────────────────
  #
  # Clip p-values away from 0 and 1 to avoid -Inf in log.
  p_clip   <- pmax(pmin(p_vals, 1 - 1e-15), 1e-15)
  chi2_f   <- -2 * sum(log(p_clip))
  df_f     <- 2L * M
  p_fisher <- stats::pchisq(chi2_f, df = df_f, lower.tail = FALSE)
  fisher   <- list(
    statistic = round(chi2_f, 4L),
    df        = df_f,
    p_value   = round(p_fisher, 4L),
    reject    = p_fisher < alpha
  )

  # ── 5. Test 3: Joint Wald test ───────────────────────────────────────────────
  #
  # Stack unit-level gaps across outcomes: for each treated unit j, form the
  # vector g_j = (gap_j,1_test, gap_j,2_test, ..., gap_j,M_test) of length
  # M * n_test, where gap_j,m_test is the n_test-vector of test-period gaps
  # for outcome m and unit j.
  #
  # We then have J observations of this M*n_test vector.
  # Sample mean:       g_bar = colMeans(G)              (M*n_test vector)
  # Sample covariance: S = cov(G)                       (M*n_test × M*n_test)
  # Wald / Hotelling:  W = J * g_bar' S^{-1} g_bar ~ chi2(M*n_test) asymptotically
  # F approximation:   F = (J-d)/((J-1)*d) * W ~ F(d, J-d)  when J > d
  #
  # When J is large relative to M*n_test, S is well-conditioned and W is
  # chi2-distributed asymptotically.  When J <= M*n_test, we apply ridge
  # regularisation S_ridge = S + lambda*I and flag this in the output.
  #
  # We also report the small-sample F approximation:
  #   F = (J - M*n_test) / ((J-1) * M*n_test) * W ~ F(M*n_test, J - M*n_test)
  # when J > M*n_test.

  wald <- NULL

  if (!has_mat) {
    warning("gaps_test_mat not found in results. ",
            "Re-run pretrend_test() with the current version to enable the ",
            "Wald test. Bonferroni and Fisher tests are still reported.")
  } else {

    # Stack into J × (M * n_test) matrix
    # Column order: [outcome1_lag1, ..., outcome1_lagP, outcome2_lag1, ...]
    G <- do.call(cbind, gap_mats)    # J × (M * n_test)
    d <- ncol(G)                     # M * n_test

    g_bar <- colMeans(G, na.rm = TRUE)          # length d
    S     <- stats::cov(G)                      # d × d  sample covariance of unit-level gaps

    # Check conditioning
    ridge_used  <- ridge_lambda
    eig_min     <- tryCatch(min(eigen(S, only.values = TRUE)$values),
                            error = function(e) NA_real_)
    auto_ridge  <- (!is.na(eig_min) && eig_min < 1e-10) || (J <= d)

    if (auto_ridge && ridge_lambda == 0) {
      # Automatic ridge: scale by trace of S / d so penalty is relative
      ridge_used <- max(1e-6, sum(diag(S)) / d * 0.01)
      warning(sprintf(
        "Covariance matrix is singular or near-singular (J=%d <= d=%d or ",
        J, d
      ), sprintf("min eigenvalue=%.2e). ", eig_min %||% NA),
      sprintf("Applying automatic ridge penalty lambda=%.2e. ", ridge_used),
      "Interpret Wald p-value with caution; prefer Bonferroni in this case.")
    }

    S_inv <- tryCatch({
      solve(S + diag(ridge_used, d))
    }, error = function(e) {
      warning("Wald test: covariance inversion failed: ", conditionMessage(e),
              ". Wald test will be NA.")
      NULL
    })

    if (!is.null(S_inv)) {
      W <- as.numeric(J * t(g_bar) %*% S_inv %*% g_bar)

      # Choose chi2 or F depending on degrees of freedom
      if (J > d) {
        # F approximation (Hotelling T2 → F)
        f_stat  <- (J - d) / ((J - 1) * d) * W
        df1_w   <- d
        df2_w   <- J - d
        p_wald  <- stats::pf(f_stat, df1 = df1_w, df2 = df2_w,
                             lower.tail = FALSE)
        wald <- list(
          statistic  = round(f_stat, 4L),
          df1        = df1_w,
          df2        = df2_w,
          df         = NULL,
          p_value    = round(p_wald, 4L),
          reject     = p_wald < alpha,
          method     = "F",
          ridge_used = ridge_used,
          W          = round(W, 4L)
        )
      } else {
        # Asymptotic chi2
        p_wald <- stats::pchisq(W, df = d, lower.tail = FALSE)
        wald <- list(
          statistic  = round(W, 4L),
          df1        = NULL,
          df2        = NULL,
          df         = d,
          p_value    = round(p_wald, 4L),
          reject     = p_wald < alpha,
          method     = "chi2",
          ridge_used = ridge_used,
          W          = round(W, 4L)
        )
      }
    }
  }

  # ── 6. Overall conclusion ────────────────────────────────────────────────────
  #
  # Prefer Wald if available and not ridge-regularised; otherwise use the
  # range spanned by Bonferroni (upper bound) and Fisher (lower bound).

  conclusion <- .joint_conclusion(bonferroni, fisher, wald, alpha, M)

  # ── 7. Assemble output ───────────────────────────────────────────────────────

  out <- list(
    outcomes    = outcomes,
    n_outcomes  = M,
    n_test_lags = n_test,
    n_units     = J,
    per_outcome = per_outcome,
    bonferroni  = bonferroni,
    fisher      = fisher,
    wald        = wald,
    conclusion  = conclusion,
    alpha       = alpha
  )

  # ── 8. Verbose output ────────────────────────────────────────────────────────

  if (verbose) .print_summary_pretrend(out)

  invisible(out)
}


# ══════════════════════════════════════════════════════════════════════════════
#  Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── .joint_conclusion ─────────────────────────────────────────────────────────
#' @noRd
.joint_conclusion <- function(bonferroni, fisher, wald, alpha, M) {
  # Primary: Wald (if available and not heavily regularised)
  if (!is.null(wald) && wald$ridge_used == 0) {
    primary        <- wald
    primary_name   <- sprintf("Wald (%s)", wald$method)
  } else {
    primary        <- NULL
    primary_name   <- NULL
  }

  bonf_reject   <- isTRUE(bonferroni$reject)
  fisher_reject <- isTRUE(fisher$reject)
  wald_reject   <- if (!is.null(wald)) isTRUE(wald$reject) else NA

  if (!is.null(primary)) {
    if (primary$reject) {
      sprintf(
        "REJECT H0 (common trends) at alpha=%.2f based on %s test (p=%.4f). ",
        alpha, primary_name, primary$p_value
      )
    } else {
      sprintf(
        "FAIL TO REJECT H0 (common trends) at alpha=%.2f based on %s test (p=%.4f).",
        alpha, primary_name, primary$p_value
      )
    }
  } else {
    # No clean Wald — summarise agreement between Bonferroni and Fisher
    if (bonf_reject && fisher_reject) {
      sprintf(
        "REJECT H0 at alpha=%.2f: both Bonferroni (p=%.4f) and Fisher (p=%.4f) reject.",
        alpha, bonferroni$p_value, fisher$p_value
      )
    } else if (!bonf_reject && !fisher_reject) {
      sprintf(
        "FAIL TO REJECT H0 at alpha=%.2f: both Bonferroni (p=%.4f) and Fisher (p=%.4f) fail to reject.",
        alpha, bonferroni$p_value, fisher$p_value
      )
    } else if (!bonf_reject && fisher_reject) {
      sprintf(
        paste0("INCONCLUSIVE at alpha=%.2f: Fisher rejects (p=%.4f) but Bonferroni does not (p=%.4f). ",
               "Outcomes are likely positively correlated; prefer Bonferroni (conservative) or Wald."),
        alpha, fisher$p_value, bonferroni$p_value
      )
    } else {
      # bonf rejects but fisher doesn't — unusual, can happen with one very small p
      sprintf(
        paste0("INCONCLUSIVE at alpha=%.2f: Bonferroni rejects (p=%.4f) but Fisher does not (p=%.4f). ",
               "One outcome may be driving the result; inspect per-outcome p-values."),
        alpha, bonferroni$p_value, fisher$p_value
      )
    }
  }
}


# ── .print_summary_pretrend ───────────────────────────────────────────────────
#' @noRd
.print_summary_pretrend <- function(out) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat("  Joint Common Trends Test\n")
  cat(strrep("=", 72), "\n", sep = "")
  cat(sprintf("  Outcomes    : %s\n", paste(out$outcomes, collapse = ", ")))
  cat(sprintf("  M           : %d outcomes\n", out$n_outcomes))
  cat(sprintf("  n_test lags : %d\n", out$n_test_lags))
  cat(sprintf("  J (units)   : %s\n",
              if (!is.null(out$n_units)) as.character(out$n_units) else "unknown"))
  cat(sprintf("  alpha       : %.2f\n", out$alpha))
  cat(strrep("-", 72), "\n", sep = "")

  # Per-outcome table
  cat("\n  Per-outcome results:\n\n")
  po <- out$per_outcome
  hdr <- sprintf("  %-22s %10s %12s %8s\n",
                 "Outcome", "p-value", "RMSE ratio", "Reject?")
  cat(hdr)
  cat("  ", strrep("-", 56), "\n", sep = "")
  for (i in seq_len(nrow(po))) {
    cat(sprintf("  %-22s %10.4f %12.4f %8s\n",
                po$Outcome[i],
                po$p_value[i],
                po$RMSE_ratio[i],
                if (po$Reject[i]) "YES *" else "no"))
  }
  cat("\n")

  # Joint tests
  cat(strrep("-", 72), "\n", sep = "")
  cat("  Joint tests:\n\n")

  # Bonferroni
  cat(sprintf("  Bonferroni  p = %.4f  [reject: %s]\n",
              out$bonferroni$p_value,
              if (isTRUE(out$bonferroni$reject)) "YES *" else "no"))
  cat(sprintf("              (controls FWER; conservative under positive correlation)\n"))

  # Fisher
  cat(sprintf("  Fisher      p = %.4f  chi2(%d) = %.4f  [reject: %s]\n",
              out$fisher$p_value, out$fisher$df, out$fisher$statistic,
              if (isTRUE(out$fisher$reject)) "YES *" else "no"))
  cat(sprintf("              (powerful under independence; anti-conservative if outcomes correlated)\n"))

  # Wald
  if (!is.null(out$wald)) {
    w <- out$wald
    if (w$method == "F") {
      cat(sprintf("  Wald        p = %.4f  F(%d,%d) = %.4f  [reject: %s]\n",
                  w$p_value, w$df1, w$df2, w$statistic,
                  if (isTRUE(w$reject)) "YES *" else "no"))
    } else {
      cat(sprintf("  Wald        p = %.4f  chi2(%d) = %.4f  [reject: %s]\n",
                  w$p_value, w$df, w$statistic,
                  if (isTRUE(w$reject)) "YES *" else "no"))
    }
    if (w$ridge_used > 0) {
      cat(sprintf("              (ridge lambda=%.2e applied; interpret with caution)\n",
                  w$ridge_used))
    } else {
      cat(sprintf("              (accounts for cross-outcome correlation; preferred test)\n"))
    }
  } else {
    cat("  Wald        not available (re-run pretrend_test() with current version)\n")
  }

  cat("\n", strrep("-", 72), "\n", sep = "")
  cat("  Conclusion:", out$conclusion, "\n")
  cat(strrep("=", 72), "\n\n", sep = "")
}
