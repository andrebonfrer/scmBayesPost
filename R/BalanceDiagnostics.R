################################################################################
##  balance_table.R
##
##  Pre-treatment balance diagnostics for multiout_synth fitted objects.
##
##  Compares, for each outcome, the pre-treatment mean of:
##    (1) the treated group   (average across all treated units)
##    (2) the synthetic control (unit-specific donor weights, averaged over
##        treated units)
##    (3) the unweighted donor pool (simple average of all donor units)
##
##  For each outcome the table reports:
##    - Treated mean
##    - Synthetic control mean
##    - Unweighted donor mean          (optional)
##    - SD of treated group            (optional)
##    - Standardised mean difference: treated vs synthetic   (SMD_syn)
##    - Standardised mean difference: treated vs donors      (SMD_donors)
##    - RMSE of pre-treatment gaps (treated − synthetic) over the L lags
##
##  Output options:
##    - Plain data.frame (always returned invisibly)
##    - Console-formatted table  (format = "console", default)
##    - LaTeX via xtable          (format = "latex")
##    - Markdown via knitr        (format = "markdown")
################################################################################

utils::globalVariables(c("SMD", "Outcome", "Comparison"))

#' Pre-Treatment Balance Diagnostics Table for multiout_synth Objects
#'
#' Computes and formats a pre-treatment balance table comparing the treated
#' group, the synthetic control, and the unweighted donor pool for each outcome.
#'
#' @section Why raw inputs are required:
#' \code{multiout_synth} does not store \code{Y_list}, \code{treat_time}, or
#' \code{L} inside the returned object, so they must be passed explicitly here
#' using the same values supplied to \code{multiout_synth}.
#'
#' @section Matrix orientation:
#' \code{Y_list} matrices must be \eqn{N \times T} (unit rows, time-period
#' columns), as produced by \code{panel_to_Ylist_int_time} (2633 units ×
#' 156 periods in your case).  \code{treat_time} is length \eqn{N}, with
#' values being 1-based column indices of the first treated period.
#'
#' @param fit        List returned by \code{multiout_synth()}.
#' @param Y_list     Named list of \eqn{T \times N} outcome matrices.
#'   Names become the row labels in the table.
#' @param treat_time Numeric vector length \eqn{N}; first treated \emph{row
#'   index} per unit (\code{Inf} = never treated).  Same vector passed to
#'   \code{multiout_synth}.
#' @param L          Integer; number of pre-treatment lags balanced.  Same
#'   value passed to \code{multiout_synth}.
#' @param outcomes   Character vector of outcome names to include.  Defaults
#'   to \code{names(Y_list)}.
#' @param standardize_outcomes Logical.  If \code{TRUE}, RMSE and SMD are
#'   computed on the standardised scale used by \code{build_Xy_for_unit}.
#'   Means are always in raw units.  Default \code{FALSE}.
#' @param intercept  Character matching the \code{intercept} argument passed
#'   to \code{multiout_synth}.  Controls whether the pre-treatment mean is
#'   removed from gaps before computing RMSE:
#'   \describe{
#'     \item{\code{"none"}}{No demeaning (default).}
#'     \item{\code{"outcome"}}{Remove the time-weighted pre-treatment mean
#'       per outcome per treated unit before computing gaps, matching
#'       \code{demean_within_outcome_blocks}.  Means and SMDs are unaffected.}
#'     \item{\code{"global"}}{Remove a single time-weighted grand mean per
#'       treated unit across all outcomes.}
#'   }
#' @param time_weights Optional numeric vector length \code{L}.  Must match
#'   the \code{time_weights} passed to \code{multiout_synth}.  Default
#'   \code{NULL} (uniform).
#' @param eps_sd     Small positive floor for SD when standardising.
#' @param std_denom  Denominator for standardised mean differences:
#'   \code{"treated"} (SD across treated-unit pre-treatment means, default),
#'   \code{"pooled"} (pooled treated + donor SD), or \code{"donor"}.
#' @param show_sd    Logical.  Add an SD (Treated) column.  Default \code{FALSE}.
#' @param show_donors_mean Logical.  Include unweighted donor mean and
#'   SMD vs donors.  Default \code{TRUE}.
#' @param digits     Integer.  Decimal places in formatted output.  Default 3.
#' @param format     One of \code{"console"} (default), \code{"latex"},
#'   \code{"markdown"}, \code{"none"}.
#' @param caption    Caption string for latex/markdown tables.
#' @param label      LaTeX label.  Default \code{"tab:balance"}.
#' @param verbose    Logical.  Print header with sample counts and L.
#'   Default \code{TRUE}.
#'
#' @return A \code{data.frame} (invisibly) with columns
#'   \code{Outcome}, \code{Mean_Treated}, \code{Mean_Synthetic},
#'   and \code{RMSE}. When \code{show_sd = TRUE}, \code{SD_Treated}
#'   is added. When \code{show_donors_mean = TRUE}, \code{Mean_Donors}
#'   and \code{SMD_Donors} are added. \code{SMD_Synthetic} is always
#'   present.
#'
#' @examples
#' \dontrun{
#' tab <- balance_table(
#'   fit        = fit,
#'   Y_list     = Y_list,
#'   treat_time = treat_time,
#'   L          = L
#' )
#'
#' # LaTeX for a paper
#' balance_table(fit, Y_list, treat_time, L,
#'               format  = "latex",
#'               caption = "Pre-treatment balance",
#'               label   = "tab:balance")
#'
#' # Love plot of SMDs (missing loveplot?)
#' tab <- balance_table(fit, Y_list, treat_time, L, show_donors_mean = TRUE)
#' }
#'
#' @export
#' @importFrom stats sd
balance_table <- function(fit,
                          Y_list,
                          treat_time,
                          L,
                          outcomes             = NULL,
                          standardize_outcomes = FALSE,
                          intercept            = c("none", "outcome", "global"),
                          time_weights         = NULL,
                          eps_sd               = 1e-8,
                          std_denom            = c("treated", "pooled", "donor"),
                          show_sd              = FALSE,
                          show_donors_mean     = TRUE,
                          digits               = 3L,
                          format               = c("console", "latex",
                                                   "markdown", "none"),
                          caption              = NULL,
                          label                = "tab:balance",
                          verbose              = TRUE) {

  # ── 0. Argument checks ──────────────────────────────────────────────────────

  std_denom <- match.arg(std_denom)
  intercept <- match.arg(intercept)
  format    <- match.arg(format)
  digits    <- as.integer(digits)
  L         <- as.integer(L)

  # Validate and normalise time_weights (same logic as pretrend_test)
  if (is.null(time_weights)) {
    tw <- rep(1 / L, L)
  } else {
    tw <- as.numeric(time_weights)
    if (length(tw) != L) stop("time_weights must have length L = ", L)
    if (any(!is.finite(tw)) || any(tw < 0) || sum(tw) <= 0) {
      stop("time_weights must be finite, non-negative, and sum to > 0.")
    }
    tw <- tw / sum(tw)
  }

  required <- c("weights", "donors", "treated_units")
  miss <- setdiff(required, names(fit))
  if (length(miss)) {
    stop("fit is missing required slot(s): ", paste(miss, collapse = ", "))
  }
  if (!is.list(Y_list) || length(Y_list) == 0L) {
    stop("Y_list must be a non-empty named list of matrices.")
  }
  if (is.null(names(Y_list)) || any(names(Y_list) == "")) {
    stop("All elements of Y_list must be named (outcome names).")
  }

  # ── 1. Resolve outcomes ─────────────────────────────────────────────────────

  if (is.null(outcomes)) outcomes <- names(Y_list)
  outcomes <- as.character(outcomes)
  bad <- setdiff(outcomes, names(Y_list))
  if (length(bad)) {
    stop("Outcome(s) not found in Y_list: ", paste(bad, collapse = ", "))
  }

  # ── 2. Geometry ─────────────────────────────────────────────────────────────

  # Y_list matrices are N × T (unit rows, time-period columns)
  N  <- nrow(Y_list[[1L]])
  TT <- ncol(Y_list[[1L]])

  if (length(treat_time) != N) {
    stop("treat_time must have length N = nrow(Y_list[[1]]) = ", N,
         ".\nGot length ", length(treat_time), ".")
  }

  treated_idx <- fit$treated_units          # integer col indices into Y matrices
  donor_idx   <- which(!is.finite(treat_time))

  n_treated <- length(treated_idx)
  n_donors  <- length(donor_idx)

  if (n_treated == 0L) stop("No treated units found in fit$treated_units.")
  if (n_donors  == 0L) stop("No donor units: all treat_time are finite.")

  J <- n_treated

  # ── 3. Per-unit pre-treatment lag row indices ────────────────────────────────

  Ti_vec   <- as.integer(treat_time[treated_idx])
  lag_rows <- lapply(Ti_vec, function(Ti) (Ti - L):(Ti - 1L))

  # Drop units that lack L pre-treatment columns
  ok <- vapply(seq_len(J), function(jj) {
    min(lag_rows[[jj]]) >= 1L && max(lag_rows[[jj]]) <= TT
  }, logical(1))

  if (!all(ok)) {
    n_drop <- sum(!ok)
    warning(n_drop, " treated unit(s) lack ", L,
            " pre-treatment rows in the matrices; excluded from balance.")
    treated_idx <- treated_idx[ok]
    Ti_vec      <- Ti_vec[ok]
    lag_rows    <- lag_rows[ok]
    J           <- length(treated_idx)
    if (J == 0L) stop("No treated units remain after feasibility check.")
  }

  # ── 4. Global donor pre-treatment window ────────────────────────────────────
  #
  # For the unweighted donor mean we use the L periods immediately before
  # the earliest treatment onset in the sample — the most conservative
  # common pre-treatment window across all cohorts.

  earliest_Ti  <- min(Ti_vec)
  global_pre   <- max(1L, earliest_Ti - L):(earliest_Ti - 1L)
  has_global   <- length(global_pre) >= 1L

  # ── 4a. Standardisation helper (mirrors build_Xy_for_unit exactly) ──────────
  #
  # mu and sd are computed from the donor-pool rows over the unit's own L
  # pre-treatment columns — the same calculation build_Xy_for_unit uses.
  # Applied only to gap computation; raw means are always in original units.

  .std_scale <- function(Y, donors_j, lag_cols_j) {
    D_pre <- Y[donors_j, lag_cols_j, drop = FALSE]
    mu    <- mean(D_pre, na.rm = TRUE)
    sd_v  <- stats::sd(as.numeric(D_pre), na.rm = TRUE)
    if (!is.finite(sd_v) || sd_v < eps_sd) sd_v <- eps_sd
    list(mu = mu, sd = sd_v)
  }

  # Replicates demean_within_outcome_blocks for a single length-L vector.
  # tw is pre-normalised to sum to 1.
  .demean_vec <- function(x) x - sum(tw * x)

  # ── 5. Per-outcome balance statistics ───────────────────────────────────────

  rows_list <- vector("list", length(outcomes))

  for (i_oc in seq_along(outcomes)) {
    oc <- outcomes[i_oc]
    Y  <- Y_list[[oc]]     # N × T

    mean_treated_j   <- numeric(J)
    mean_synthetic_j <- numeric(J)
    gap_vec          <- numeric(0)   # standardised if standardize_outcomes

    for (jj in seq_len(J)) {
      col_j   <- treated_idx[jj]
      cols    <- lag_rows[[jj]]
      donors  <- fit$donors[[jj]]
      w       <- fit$weights[[jj]]

      y_treated   <- Y[col_j,  cols]
      Y_donors    <- Y[donors, cols, drop = FALSE]   # n_donors_j × L
      y_synthetic <- as.numeric(t(Y_donors) %*% w)  # length L

      # Raw means — always in original units, never demeaned
      mean_treated_j[jj]   <- mean(y_treated,   na.rm = TRUE)
      mean_synthetic_j[jj] <- mean(y_synthetic, na.rm = TRUE)

      # Gaps for RMSE — apply same transformations solver used
      gap <- y_treated - y_synthetic

      if (standardize_outcomes) {
        sc  <- .std_scale(Y, donors, cols)
        gap <- gap / sc$sd
      }

      if (intercept != "none") {
        gap <- .demean_vec(gap)
      }

      gap_vec <- c(gap_vec, gap)
    }

    mean_treated   <- mean(mean_treated_j,   na.rm = TRUE)
    mean_synthetic <- mean(mean_synthetic_j, na.rm = TRUE)
    rmse           <- sqrt(mean(gap_vec^2,   na.rm = TRUE))

    # SMD denominator is always based on raw treated-unit means, regardless
    # of standardize_outcomes, because mean_treated_j is always in raw units.
    # This makes SMD_Synthetic directly comparable to SMD_Donors (which is
    # also in raw units).  RMSE is on the standardised scale when
    # standardize_outcomes = TRUE, so the two columns are intentionally on
    # different scales in that case — noted in the table footer.

    # Unweighted donor mean over global pre-treatment window
    if (has_global && n_donors >= 1L) {
      donor_pre_mat <- Y[donor_idx, global_pre, drop = FALSE]  # donors × window
      mean_donors   <- mean(donor_pre_mat, na.rm = TRUE)
      donor_unit_means <- rowMeans(donor_pre_mat, na.rm = TRUE)
      sd_donor      <- if (n_donors > 1L) sd(donor_unit_means, na.rm = TRUE)
      else NA_real_
    } else {
      mean_donors  <- NA_real_
      sd_donor     <- NA_real_
    }

    sd_treated <- if (J > 1L) sd(mean_treated_j, na.rm = TRUE) else NA_real_

    sd_denom_val <- switch(std_denom,
                           "treated" = sd_treated,
                           "donor"   = sd_donor,
                           "pooled"  = {
                             if (!is.na(sd_treated) && !is.na(sd_donor) && J > 1L && n_donors > 1L) {
                               sqrt(((J - 1) * sd_treated^2 + (n_donors - 1) * sd_donor^2) /
                                      (J + n_donors - 2))
                             } else {
                               sd_treated
                             }
                           }
    )

    .smd <- function(diff) {
      if (!is.na(sd_denom_val) && sd_denom_val > 0) diff / sd_denom_val
      else NA_real_
    }

    smd_synthetic <- .smd(mean_treated - mean_synthetic)
    smd_donors    <- if (show_donors_mean && !is.na(mean_donors)) {
      .smd(mean_treated - mean_donors)
    } else NA_real_

    # Assemble row
    row <- data.frame(Outcome = oc, Mean_Treated = mean_treated,
                      stringsAsFactors = FALSE)
    if (show_sd)          row$SD_Treated    <- sd_treated
    row$Mean_Synthetic    <- mean_synthetic
    if (show_donors_mean) row$Mean_Donors   <- mean_donors
    row$SMD_Synthetic     <- smd_synthetic
    if (show_donors_mean) row$SMD_Donors    <- smd_donors
    row$RMSE              <- rmse

    rows_list[[i_oc]] <- row
  }

  tab <- do.call(rbind, rows_list)
  attr(tab, "standardized") <- standardize_outcomes

  # ── 6. Round ─────────────────────────────────────────────────────────────────

  num_cols <- setdiff(names(tab), "Outcome")
  tab[num_cols] <- lapply(tab[num_cols], round, digits = digits)

  # ── 7. Verbose header ────────────────────────────────────────────────────────

  if (verbose) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("  Pre-Treatment Balance Diagnostics\n")
    cat(strrep("=", 70), "\n", sep = "")
    cat(sprintf("  Treated units  : %d\n", n_treated))
    cat(sprintf("  Donor units    : %d\n", n_donors))
    cat(sprintf("  Pre-treat lags : L = %d\n", L))
    cat(sprintf("  Outcomes       : %s\n", paste(outcomes, collapse = ", ")))
    cat(sprintf("  SMD denominator: %s SD (raw units)\n", std_denom))
    rmse_desc <- if (standardize_outcomes && intercept != "none") {
      sprintf("standardised + demeaned (intercept = \"%s\")", intercept)
    } else if (standardize_outcomes) {
      "standardised scale (donor-pool SD per unit)"
    } else if (intercept != "none") {
      sprintf("demeaned (intercept = \"%s\"), raw scale", intercept)
    } else {
      "raw units"
    }
    cat(sprintf("  RMSE           : %s\n", rmse_desc))
    cat(strrep("-", 70), "\n\n", sep = "")
  }

  # ── 8. Format and print ──────────────────────────────────────────────────────

  .format_balance_table(tab, format, caption, label, digits,
                        show_sd, show_donors_mean, std_denom)

  invisible(tab)
}


# ══════════════════════════════════════════════════════════════════════════════
#  Internal formatting helpers
# ══════════════════════════════════════════════════════════════════════════════

`%||%` <- function(a, b) if (!is.null(a)) a else b

.col_labels_bt <- function(col_names) {
  map <- c(
    Outcome        = "Outcome",
    Mean_Treated   = "Mean (Treated)",
    SD_Treated     = "SD (Treated)",
    Mean_Synthetic = "Mean (Synthetic)",
    Mean_Donors    = "Mean (Donors)",
    SMD_Synthetic  = "SMD (Syn.)",
    SMD_Donors     = "SMD (Donors)",
    RMSE           = "RMSE"
  )
  unname(map[col_names])
}

.rename_for_display_bt <- function(df) {
  names(df) <- .col_labels_bt(names(df))
  df
}

.format_balance_table <- function(tab, format, caption, label, digits,
                                  show_sd, show_donors_mean, std_denom) {
  if (format == "none") return(invisible(NULL))

  if (format == "console") {
    .print_console_table_bt(tab, digits, std_denom)
    return(invisible(NULL))
  }

  if (format == "latex") {
    if (!requireNamespace("xtable", quietly = TRUE)) {
      warning("Package 'xtable' needed for format='latex'. Falling back to console.")
      .print_console_table_bt(tab, digits, std_denom)
      return(invisible(NULL))
    }
    cap <- caption %||%
      "Pre-treatment balance: treated group vs.\\ synthetic control."
    xt <- xtable::xtable(
      .rename_for_display_bt(tab),
      caption = cap,
      label   = label,
      digits  = c(0L, rep(digits, ncol(tab) - 1L))
    )
    print(xt, include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top",
          sanitize.colnames.function = identity)
    return(invisible(NULL))
  }

  if (format == "markdown") {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      warning("Package 'knitr' needed for format='markdown'. Falling back to console.")
      .print_console_table_bt(tab, digits, std_denom)
      return(invisible(NULL))
    }
    cap <- caption %||% "Pre-treatment balance: treated vs. synthetic control."
    cat(knitr::kable(
      .rename_for_display_bt(tab),
      format  = "markdown",
      digits  = digits,
      caption = cap,
      align   = c("l", rep("r", ncol(tab) - 1L))
    ), sep = "\n")
    return(invisible(NULL))
  }
}

.print_console_table_bt <- function(tab, digits, std_denom) {
  display <- .rename_for_display_bt(tab)
  fmt     <- paste0("%.", digits, "f")

  for (nm in names(display)[-1L]) {
    display[[nm]] <- ifelse(is.na(display[[nm]]), "NA",
                            sprintf(fmt, as.numeric(display[[nm]])))
  }

  # Compute per-column widths as plain unnamed integers
  widths <- mapply(function(col_name, col_vals) {
    max(nchar(col_name), max(nchar(as.character(col_vals)), na.rm = TRUE)) + 2L
  }, names(display), display, SIMPLIFY = TRUE, USE.NAMES = FALSE)
  widths <- as.integer(widths)

  header <- paste(mapply(function(nm, w) formatC(nm, width = w, flag = "-"),
                         names(display), widths, SIMPLIFY = TRUE),
                  collapse = "")
  sep_ln <- strrep("-", nchar(header))

  cat(sep_ln, "\n")
  cat(header, "\n")
  cat(sep_ln, "\n")
  for (i in seq_len(nrow(display))) {
    cat(paste(mapply(function(v, w) formatC(as.character(v), width = w, flag = "-"),
                     v = as.list(display[i, ]), w = widths,
                     SIMPLIFY = TRUE),
              collapse = ""), "\n")
  }
  cat(sep_ln, "\n")
  cat(sprintf(
    "\nSMD = (Treated \u2212 Comparison) / SD_%s (raw units).\n",
    switch(std_denom, treated = "Treated", pooled = "Pooled", donor = "Donor")
  ))
  cat(sprintf(
    "RMSE computed over L balanced pre-treatment lags per unit (%s).\n\n",
    if (isTRUE(attr(tab, "standardized"))) "standardised scale"
    else "raw units"
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
#  Love plot
# ══════════════════════════════════════════════════════════════════════════════

#' Love Plot of Standardised Mean Differences
#'
#' Produces a Love plot (dot chart of SMDs before and after synthetic control
#' weighting) from the output of \code{\link{balance_table}}.
#'
#' @param tab       A \code{data.frame} returned by \code{\link{balance_table}}.
#'   Run with \code{show_donors_mean = TRUE} to show both before and after.
#' @param threshold Numeric.  Reference lines drawn at ±\code{threshold}.
#'   Default \code{0.1}.
#' @param title     Plot title.
#' @param xlab      x-axis label.
#'
#' @return A \code{ggplot} object (invisibly if \pkg{ggplot2} available),
#'   or a base-graphics plot.
#' @export
love_plot <- function(tab,
                      threshold = 0.1,
                      title = "Standardised Mean Differences (Love Plot)",
                      xlab  = "Standardised Mean Difference") {

  has_donors <- "SMD_Donors" %in% names(tab)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    long <- if (has_donors) {
      rbind(
        data.frame(Outcome    = tab$Outcome,
                   SMD        = tab$SMD_Synthetic,
                   Comparison = "Synthetic Control",
                   stringsAsFactors = FALSE),
        data.frame(Outcome    = tab$Outcome,
                   SMD        = tab$SMD_Donors,
                   Comparison = "Unweighted Donors",
                   stringsAsFactors = FALSE)
      )
    } else {
      data.frame(Outcome    = tab$Outcome,
                 SMD        = tab$SMD_Synthetic,
                 Comparison = "Synthetic Control",
                 stringsAsFactors = FALSE)
    }
    long$Outcome <- factor(long$Outcome, levels = rev(unique(tab$Outcome)))

    p <- ggplot2::ggplot(long,
                         ggplot2::aes(x = SMD, y = Outcome,
                                      shape = Comparison, color = Comparison)) +
      ggplot2::geom_point(size = 3.5) +
      ggplot2::geom_vline(xintercept =  threshold, linetype = "dashed",
                          color = "grey50") +
      ggplot2::geom_vline(xintercept = -threshold, linetype = "dashed",
                          color = "grey50") +
      ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
      ggplot2::scale_color_manual(
        values = c("Synthetic Control" = "#1f78b4",
                   "Unweighted Donors" = "#e31a1c")) +
      ggplot2::scale_shape_manual(
        values = c("Synthetic Control" = 16L,
                   "Unweighted Donors" = 17L)) +
      ggplot2::labs(title = title, x = xlab, y = NULL,
                    color = NULL, shape = NULL) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position  = "bottom",
                     panel.grid.minor = ggplot2::element_blank())
    print(p)
    return(invisible(p))
  }

  # Base graphics fallback
  op <- graphics::par(mar = c(4, max(nchar(tab$Outcome)) * 0.5 + 1, 3, 2))
  on.exit(graphics::par(op))
  n    <- nrow(tab)
  xlim <- range(c(tab$SMD_Synthetic,
                  if (has_donors) tab$SMD_Donors else NULL,
                  -threshold, threshold), na.rm = TRUE) * 1.15

  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = c(0.5, n + 0.5))
  graphics::abline(v = c(-threshold, 0, threshold),
                   lty = c(2, 1, 2), col = c("grey60", "black", "grey60"))
  for (i in seq_len(n)) {
    graphics::points(tab$SMD_Synthetic[i], i, pch = 16L,
                     col = "#1f78b4", cex = 1.4)
    if (has_donors)
      graphics::points(tab$SMD_Donors[i], i, pch = 17L,
                       col = "#e31a1c", cex = 1.4)
  }
  graphics::axis(1)
  graphics::axis(2, at = seq_len(n), labels = tab$Outcome, las = 1)
  graphics::title(main = title, xlab = xlab)
  if (has_donors) {
    graphics::legend("bottomright",
                     legend = c("Synthetic Control", "Unweighted Donors"),
                     pch = c(16L, 17L), col = c("#1f78b4", "#e31a1c"), bty = "n")
  }
  invisible(NULL)
}
