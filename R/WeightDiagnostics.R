################################################################################
##  weight_diagnostics.R
##
##  Donor weight distribution diagnostics for multiout_synth fitted objects.
##
##  Computes, for each treated unit and across the full sample:
##
##  Unit-level summaries:
##    - Effective number of donors  (inverse Herfindahl: 1 / sum(w^2))
##    - Maximum donor weight
##    - Number of donors with weight above a threshold
##    - Gini coefficient of the weight distribution
##    - Entropy of the weight distribution
##
##  Donor-level summaries:
##    - Number of treated units for which the donor receives non-negligible weight
##    - Total weight received across all treated units
##    - Mean weight conditional on being used
##
##  Outputs:
##    - unit_summary:  data frame, one row per treated unit
##    - donor_summary: data frame, one row per donor unit
##    - aggregate:     scalar summaries across all treated units
##    - Plots:         histogram of N_eff, max weight distribution,
##                     donor usage frequency, weight Lorenz curve
################################################################################


#' Donor Weight Distribution Diagnostics for multiout_synth Objects
#'
#' Summarises the distribution of synthetic control weights across donor units,
#' both at the level of individual treated units and aggregated across all
#' treated units.  This is a standard robustness diagnostic in the synthetic
#' control literature: highly concentrated weights (few donors receiving most
#' of the weight) indicate a fragile counterfactual.
#'
#' @section Metrics:
#' \describe{
#'   \item{N_eff}{Effective number of donors: \eqn{1 / \sum_d w_d^2}.  Equals
#'     1 when all weight is on one donor; equals the total donor count when
#'     weights are uniform.  The most interpretable concentration measure.}
#'   \item{Max weight}{Largest single donor weight.  Above ~0.3 suggests
#'     substantial reliance on a single donor.}
#'   \item{N_nonzero}{Number of donors with weight above \code{threshold}.}
#'   \item{Gini}{Gini coefficient of the weight vector.  0 = perfectly uniform;
#'     1 = all weight on one donor.}
#'   \item{Entropy}{Shannon entropy: \eqn{-\sum_d w_d \log(w_d)}.  Higher =
#'     more spread.  Normalised by \eqn{\log(N_{\text{donors}})} to lie in
#'     [0, 1].}
#' }
#'
#' @param fit         List returned by \code{multiout_synth()}.
#' @param unit_ids    Optional character/numeric vector of length \eqn{N}
#'   (all units).  Used to label donor rows in \code{donor_summary}.  If
#'   \code{NULL}, uses \code{fit$unit_ids} if present, otherwise integer
#'   indices.
#' @param threshold   Numeric.  Weight threshold for counting "non-negligible"
#'   donors.  Default \code{0.01}.
#' @param top_n_donors Integer.  Number of top donors (by total weight received)
#'   to highlight in verbose output.  Default \code{10L}.
#' @param verbose     Logical.  Print summary tables.  Default \code{TRUE}.
#'
#' @return A list (invisibly) with:
#' \describe{
#'   \item{\code{unit_summary}}{Data frame, one row per treated unit, with
#'     columns \code{unit_id}, \code{N_donors_available}, \code{N_eff},
#'     \code{max_weight}, \code{N_nonzero}, \code{gini}, \code{entropy_norm}.}
#'   \item{\code{donor_summary}}{Data frame, one row per donor unit, with
#'     columns \code{donor_id}, \code{n_used} (number of treated units for
#'     which this donor has weight > \code{threshold}), \code{total_weight}
#'     (sum across treated units), \code{mean_weight_when_used}.}
#'   \item{\code{aggregate}}{Named numeric vector of cross-unit summary
#'     statistics: mean, median, sd, p10, p90 of N_eff and max_weight;
#'     fraction of treated units with max_weight > 0.3; total unique donors
#'     used.}
#'   \item{\code{threshold}}{The threshold used.}
#' }
#'
#' @examples
#' \dontrun{
#' wd <- weight_diagnostics(fit)
#'
#' # Plots
#' plot_weight_diagnostics(wd)
#'
#' # Access unit-level table
#' head(wd$unit_summary)
#'
#' # Donors used by more than 10% of treated units
#' subset(wd$donor_summary, n_used > 0.1 * nrow(wd$unit_summary))
#' }
#'
#' @seealso \code{\link{plot_weight_diagnostics}}, \code{\link{balance_table}}
#' @export
weight_diagnostics <- function(fit,
                               unit_ids    = NULL,
                               threshold   = 0.01,
                               top_n_donors = 10L,
                               verbose     = TRUE) {

  # ── 0. Validation ──────────────────────────────────────────────────────────

  required <- c("weights", "donors", "treated_units")
  miss <- setdiff(required, names(fit))
  if (length(miss)) {
    stop("fit is missing slot(s): ", paste(miss, collapse = ", "))
  }

  J           <- length(fit$treated_units)
  treated_idx <- fit$treated_units

  if (J == 0L) stop("No treated units in fit$treated_units.")

  # Resolve unit IDs
  if (is.null(unit_ids)) {
    if (!is.null(fit$unit_ids)) {
      unit_ids <- fit$unit_ids
    } else {
      # Fall back to integer indices
      max_idx  <- max(unlist(fit$donors), treated_idx)
      unit_ids <- seq_len(max_idx)
    }
  }

  treated_unit_ids <- if (!is.null(fit$treated_unit_ids)) {
    fit$treated_unit_ids
  } else {
    unit_ids[treated_idx]
  }

  # All donor indices ever used
  all_donor_idx <- sort(unique(unlist(fit$donors)))
  donor_ids     <- unit_ids[all_donor_idx]
  D             <- length(all_donor_idx)

  # ── 1. Unit-level metrics ──────────────────────────────────────────────────

  n_donors_avail <- integer(J)
  n_eff          <- numeric(J)
  max_w          <- numeric(J)
  n_nonzero      <- integer(J)
  gini_v         <- numeric(J)
  entropy_v      <- numeric(J)

  for (jj in seq_len(J)) {
    w       <- fit$weights[[jj]]
    donors  <- fit$donors[[jj]]
    nd      <- length(donors)

    # Ensure non-negative and sum to 1 (solver should guarantee this,
    # but floating point may introduce tiny violations)
    w <- pmax(w, 0)
    ws <- sum(w)
    if (ws > 0) w <- w / ws

    n_donors_avail[jj] <- nd
    n_eff[jj]          <- .n_eff(w)
    max_w[jj]          <- max(w)
    n_nonzero[jj]      <- sum(w > threshold)
    gini_v[jj]         <- .gini(w)
    entropy_v[jj]      <- .entropy_norm(w, nd)
  }

  unit_summary <- data.frame(
    unit_id            = as.character(treated_unit_ids),
    N_donors_available = n_donors_avail,
    N_eff              = round(n_eff,     2L),
    max_weight         = round(max_w,     4L),
    N_nonzero          = n_nonzero,
    gini               = round(gini_v,    4L),
    entropy_norm       = round(entropy_v, 4L),
    stringsAsFactors   = FALSE
  )

  # ── 2. Donor-level metrics ─────────────────────────────────────────────────
  #
  # For each donor, count how many treated units use it (weight > threshold)
  # and compute total weight received (summed across treated units, where
  # weight is 0 if the donor is not in that unit's donor set).

  n_used       <- integer(D)
  total_weight <- numeric(D)

  # Index map: donor global index → position in all_donor_idx
  donor_pos <- match(all_donor_idx, all_donor_idx)   # identity, but kept for clarity
  idx_map   <- setNames(seq_len(D), as.character(all_donor_idx))

  for (jj in seq_len(J)) {
    w      <- pmax(fit$weights[[jj]], 0)
    ws     <- sum(w); if (ws > 0) w <- w / ws
    donors <- fit$donors[[jj]]

    for (kk in seq_along(donors)) {
      d_idx <- as.character(donors[kk])
      pos   <- idx_map[d_idx]
      if (!is.na(pos)) {
        total_weight[pos] <- total_weight[pos] + w[kk]
        if (w[kk] > threshold) n_used[pos] <- n_used[pos] + 1L
      }
    }
  }

  mean_w_when_used <- ifelse(n_used > 0, total_weight / n_used, NA_real_)

  donor_summary <- data.frame(
    donor_id            = as.character(donor_ids),
    n_used              = n_used,
    pct_treated         = round(100 * n_used / J, 1L),
    total_weight        = round(total_weight,      4L),
    mean_weight_used    = round(mean_w_when_used,  4L),
    stringsAsFactors    = FALSE
  )
  donor_summary <- donor_summary[order(-donor_summary$total_weight), ]

  # ── 3. Aggregate summaries ─────────────────────────────────────────────────

  .pct <- function(x, p) as.numeric(stats::quantile(x, p / 100, na.rm = TRUE))

  aggregate <- c(
    N_eff_mean       = round(mean(n_eff),    2L),
    N_eff_median     = round(stats::median(n_eff), 2L),
    N_eff_sd         = round(stats::sd(n_eff),     2L),
    N_eff_p10        = round(.pct(n_eff, 10),       2L),
    N_eff_p90        = round(.pct(n_eff, 90),       2L),
    maxw_mean        = round(mean(max_w),    4L),
    maxw_median      = round(stats::median(max_w),  4L),
    maxw_p10         = round(.pct(max_w, 10),        4L),
    maxw_p90         = round(.pct(max_w, 90),        4L),
    frac_maxw_gt30   = round(mean(max_w > 0.30),     3L),
    frac_maxw_gt50   = round(mean(max_w > 0.50),     3L),
    gini_mean        = round(mean(gini_v),   4L),
    gini_median      = round(stats::median(gini_v),  4L),
    entropy_mean     = round(mean(entropy_v),4L),
    n_unique_donors  = D,
    n_treated_units  = J
  )

  out <- list(
    unit_summary  = unit_summary,
    donor_summary = donor_summary,
    aggregate     = aggregate,
    threshold     = threshold
  )

  # ── 4. Verbose output ──────────────────────────────────────────────────────

  if (verbose) .print_weight_diagnostics(out, top_n_donors)

  invisible(out)
}


# ══════════════════════════════════════════════════════════════════════════════
#  Metric helpers
# ══════════════════════════════════════════════════════════════════════════════

# Effective number of donors (inverse Herfindahl)
.n_eff <- function(w) {
  w <- w[w > 0]
  if (length(w) == 0L) return(0)
  1 / sum(w^2)
}

# Gini coefficient
.gini <- function(w) {
  w <- sort(w[w >= 0])
  n <- length(w)
  if (n <= 1L || sum(w) == 0) return(0)
  # Standard Gini formula for sorted non-negative values
  idx <- seq_len(n)
  2 * sum(idx * w) / (n * sum(w)) - (n + 1) / n
}

# Normalised Shannon entropy (in [0,1])
.entropy_norm <- function(w, n_donors) {
  w <- w[w > 0]
  if (length(w) <= 1L || n_donors <= 1L) return(0)
  H    <- -sum(w * log(w))
  Hmax <- log(n_donors)
  if (Hmax <= 0) return(0)
  H / Hmax
}


# ══════════════════════════════════════════════════════════════════════════════
#  Verbose printer
# ══════════════════════════════════════════════════════════════════════════════

#' @noRd
.print_weight_diagnostics <- function(out, top_n) {
  ag <- out$aggregate
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("  Donor Weight Distribution Diagnostics\n")
  cat(strrep("=", 70), "\n", sep = "")
  cat(sprintf("  Treated units  : %d\n", ag["n_treated_units"]))
  cat(sprintf("  Unique donors  : %d\n", ag["n_unique_donors"]))
  cat(sprintf("  Weight threshold (non-negligible): %.3f\n\n", out$threshold))

  cat("  Effective number of donors (N_eff = 1/sum(w^2)):\n")
  cat(sprintf("    Mean   = %6.2f   Median = %6.2f   SD  = %6.2f\n",
              ag["N_eff_mean"], ag["N_eff_median"], ag["N_eff_sd"]))
  cat(sprintf("    P10    = %6.2f   P90    = %6.2f\n\n",
              ag["N_eff_p10"], ag["N_eff_p90"]))

  cat("  Maximum single-donor weight:\n")
  cat(sprintf("    Mean   = %6.4f   Median = %6.4f\n",
              ag["maxw_mean"], ag["maxw_median"]))
  cat(sprintf("    P10    = %6.4f   P90    = %6.4f\n",
              ag["maxw_p10"], ag["maxw_p90"]))
  cat(sprintf("    Fraction with max weight > 0.30 : %.1f%%\n",
              100 * ag["frac_maxw_gt30"]))
  cat(sprintf("    Fraction with max weight > 0.50 : %.1f%%\n\n",
              100 * ag["frac_maxw_gt50"]))

  cat("  Weight inequality (across donors within each treated unit):\n")
  cat(sprintf("    Gini coefficient  : mean = %.4f, median = %.4f\n",
              ag["gini_mean"], ag["gini_median"]))
  cat(sprintf("    Normalised entropy: mean = %.4f  (0=concentrated, 1=uniform)\n\n",
              ag["entropy_mean"]))

  # Top donors
  top_n  <- min(top_n, nrow(out$donor_summary))
  top_ds <- out$donor_summary[seq_len(top_n), ]
  cat(sprintf("  Top %d donors by total weight received:\n\n", top_n))
  cat(sprintf("  %-20s %8s %10s %12s %16s\n",
              "Donor ID", "N used", "% treated", "Total wt", "Mean wt (when used)"))
  cat("  ", strrep("-", 70), "\n", sep = "")
  for (i in seq_len(nrow(top_ds))) {
    cat(sprintf("  %-20s %8d %10.1f %12.4f %16.4f\n",
                top_ds$donor_id[i],
                top_ds$n_used[i],
                top_ds$pct_treated[i],
                top_ds$total_weight[i],
                ifelse(is.na(top_ds$mean_weight_used[i]),
                       NA_real_, top_ds$mean_weight_used[i])))
  }
  cat("\n", strrep("=", 70), "\n\n", sep = "")
}


# ══════════════════════════════════════════════════════════════════════════════
#  Plot function
# ══════════════════════════════════════════════════════════════════════════════

#' Plot Donor Weight Distribution Diagnostics
#'
#' Produces a four-panel diagnostic plot from the output of
#' \code{\link{weight_diagnostics}}:
#' \enumerate{
#'   \item Histogram of the effective number of donors (N_eff) across treated
#'     units.
#'   \item Histogram of the maximum single-donor weight across treated units.
#'   \item Bar chart of the top donors by total weight received, showing how
#'     many treated units each donor serves.
#'   \item Lorenz curve of the weight distribution aggregated across all
#'     treated units, illustrating overall concentration.
#' }
#'
#' @param wd         List returned by \code{\link{weight_diagnostics}}.
#' @param top_n      Integer.  Number of top donors to show in panel 3.
#'   Default \code{15L}.
#' @param bins       Integer.  Number of histogram bins.  Default \code{30L}.
#' @param threshold_lines Logical.  Draw reference lines at max_weight = 0.3
#'   and 0.5 in panel 2.  Default \code{TRUE}.
#'
#' @return \code{NULL} invisibly (plots produced as side effect).  If
#'   \pkg{ggplot2} is available, a \code{patchwork} or \code{gridExtra}
#'   layout is used; otherwise base graphics.
#' @export
plot_weight_diagnostics <- function(wd,
                                    top_n           = 15L,
                                    bins            = 30L,
                                    threshold_lines = TRUE) {

  us <- wd$unit_summary
  ds <- wd$donor_summary
  J  <- nrow(us)

  top_n  <- min(as.integer(top_n), nrow(ds))
  top_ds <- ds[seq_len(top_n), ]
  # Reverse for horizontal bar chart so highest is at top
  top_ds <- top_ds[rev(seq_len(nrow(top_ds))), ]

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .plot_gg(us, top_ds, J, bins, threshold_lines, wd$threshold)
  } else {
    .plot_base(us, top_ds, J, bins, threshold_lines, wd$threshold)
  }

  invisible(NULL)
}


# ── ggplot2 version ───────────────────────────────────────────────────────────
#' @noRd
.plot_gg <- function(us, top_ds, J, bins, threshold_lines, threshold) {

  p1 <- ggplot2::ggplot(us, ggplot2::aes(x = N_eff)) +
    ggplot2::geom_histogram(bins = bins, fill = "#1f78b4", color = "white",
                            alpha = 0.85) +
    ggplot2::geom_vline(xintercept = stats::median(us$N_eff),
                        linetype = "dashed", color = "firebrick", linewidth = 0.7) +
    ggplot2::labs(
      title    = "Effective number of donors (N_eff)",
      subtitle = sprintf("Median = %.1f  |  1 = single donor, higher = more spread",
                         stats::median(us$N_eff)),
      x = expression(N[eff] == 1 / sum(w[d]^2)),
      y = "Treated units"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  p2 <- ggplot2::ggplot(us, ggplot2::aes(x = max_weight)) +
    ggplot2::geom_histogram(bins = bins, fill = "#33a02c", color = "white",
                            alpha = 0.85) +
    ggplot2::geom_vline(xintercept = stats::median(us$max_weight),
                        linetype = "dashed", color = "firebrick", linewidth = 0.7)
  if (threshold_lines) {
    p2 <- p2 +
      ggplot2::geom_vline(xintercept = 0.30, linetype = "dotted",
                          color = "darkorange", linewidth = 0.6) +
      ggplot2::geom_vline(xintercept = 0.50, linetype = "dotted",
                          color = "red3", linewidth = 0.6) +
      ggplot2::annotate("text", x = 0.31, y = Inf, label = "0.30",
                        hjust = 0, vjust = 1.5, size = 3, color = "darkorange") +
      ggplot2::annotate("text", x = 0.51, y = Inf, label = "0.50",
                        hjust = 0, vjust = 1.5, size = 3, color = "red3")
  }
  p2 <- p2 +
    ggplot2::labs(
      title    = "Maximum single-donor weight",
      subtitle = sprintf("Median = %.3f  |  %.1f%% of units have max weight > 0.30",
                         stats::median(us$max_weight),
                         100 * mean(us$max_weight > 0.30)),
      x = "Max weight", y = "Treated units"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  top_ds$donor_id <- factor(top_ds$donor_id,
                            levels = top_ds$donor_id)  # preserve order
  p3 <- ggplot2::ggplot(top_ds,
                        ggplot2::aes(x = total_weight, y = donor_id,
                                     fill = pct_treated)) +
    ggplot2::geom_col(alpha = 0.9) +
    ggplot2::scale_fill_gradient(low = "#deebf7", high = "#08519c",
                                 name = "% treated\nunits") +
    ggplot2::labs(
      title    = sprintf("Top %d donors by total weight", nrow(top_ds)),
      subtitle = "Colour = % of treated units for which this donor has non-negligible weight",
      x        = "Total weight (sum across treated units)",
      y        = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor  = ggplot2::element_blank(),
                   legend.position   = "right")

  # Lorenz curve: aggregate weight distribution across all treated units
  lc <- .lorenz_data(us)
  p4 <- ggplot2::ggplot(lc, ggplot2::aes(x = cum_units, y = cum_weight)) +
    ggplot2::geom_line(color = "#1f78b4", linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey50") +
    ggplot2::labs(
      title    = "Lorenz curve: N_eff distribution",
      subtitle = sprintf("Gini = %.3f  (0 = uniform, 1 = fully concentrated)",
                         mean(us$gini, na.rm = TRUE)),
      x = "Cumulative share of treated units (sorted by N_eff)",
      y = "Cumulative share of total N_eff"
    ) +
    ggplot2::scale_x_continuous(labels = scales_pct) +
    ggplot2::scale_y_continuous(labels = scales_pct) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  # Arrange in 2×2 grid
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2L)
  } else {
    # Print sequentially if gridExtra not available
    print(p1); print(p2); print(p3); print(p4)
  }
}

# Simple percentage formatter without scales dependency
scales_pct <- function(x) paste0(round(x * 100), "%")


# ── Base graphics version ─────────────────────────────────────────────────────
#' @noRd
.plot_base <- function(us, top_ds, J, bins, threshold_lines, threshold) {
  old_par <- graphics::par(mfrow = c(2L, 2L), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(old_par))

  # Panel 1: N_eff histogram
  graphics::hist(us$N_eff, breaks = bins, col = "#1f78b4", border = "white",
                 main = "Effective number of donors (N_eff)",
                 xlab = expression(N[eff] == 1 / sum(w[d]^2)),
                 ylab = "Treated units")
  graphics::abline(v = stats::median(us$N_eff), lty = 2, col = "firebrick",
                   lwd = 1.5)
  graphics::legend("topright",
                   legend = sprintf("Median = %.1f", stats::median(us$N_eff)),
                   lty = 2, col = "firebrick", bty = "n", cex = 0.85)

  # Panel 2: Max weight histogram
  graphics::hist(us$max_weight, breaks = bins, col = "#33a02c", border = "white",
                 main = "Maximum single-donor weight",
                 xlab = "Max weight", ylab = "Treated units")
  graphics::abline(v = stats::median(us$max_weight), lty = 2,
                   col = "firebrick", lwd = 1.5)
  if (threshold_lines) {
    graphics::abline(v = 0.30, lty = 3, col = "darkorange", lwd = 1.5)
    graphics::abline(v = 0.50, lty = 3, col = "red3",       lwd = 1.5)
    graphics::legend("topright",
                     legend = c(sprintf("Median = %.3f", stats::median(us$max_weight)),
                                "0.30 threshold", "0.50 threshold"),
                     lty = c(2, 3, 3), col = c("firebrick", "darkorange", "red3"),
                     bty = "n", cex = 0.8)
  }

  # Panel 3: Top donors horizontal bar chart
  nd  <- nrow(top_ds)
  col_scale <- grDevices::colorRampPalette(c("#deebf7", "#08519c"))(nd)
  pct_rank  <- rank(top_ds$pct_treated)
  bar_cols  <- col_scale[pct_rank]

  graphics::barplot(top_ds$total_weight,
                    names.arg = top_ds$donor_id,
                    horiz     = TRUE,
                    col       = bar_cols,
                    border    = NA,
                    las       = 1,
                    cex.names = 0.7,
                    main      = sprintf("Top %d donors by total weight", nd),
                    xlab      = "Total weight (sum across treated units)")

  # Panel 4: Lorenz curve of N_eff
  lc <- .lorenz_data(us)
  graphics::plot(lc$cum_units, lc$cum_weight,
                 type = "l", col = "#1f78b4", lwd = 2,
                 main = "Lorenz curve: N_eff distribution",
                 xlab = "Cumulative share of treated units (sorted by N_eff)",
                 ylab = "Cumulative share of total N_eff",
                 xlim = c(0, 1), ylim = c(0, 1))
  graphics::abline(0, 1, lty = 2, col = "grey50")
  graphics::legend("topleft",
                   legend = sprintf("Mean Gini = %.3f", mean(us$gini, na.rm = TRUE)),
                   bty = "n", cex = 0.85)
}


# ── Lorenz curve data ─────────────────────────────────────────────────────────
#' @noRd
.lorenz_data <- function(us) {
  # Lorenz curve of N_eff across treated units
  x <- sort(us$N_eff)
  n <- length(x)
  data.frame(
    cum_units  = c(0, seq_len(n) / n),
    cum_weight = c(0, cumsum(x) / sum(x))
  )
}
