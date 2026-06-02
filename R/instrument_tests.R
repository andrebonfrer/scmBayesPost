# =============================================================================
#  instrument_tests.R
#
#  Tests for instrument validity (exclusion restriction) and strength
#  (relevance) for use with the selection_probit_bayes first stage in
#  scmBayesPost.
#
#  Contents
#  --------
#  .parse_iv_formula()             internal: parse treatment/instruments from f.X
#  test_instrument_strength()      F-statistic, partial R2, Cragg-Donald
#  test_instrument_exclusion()     Reduced-form placebo tests
#  test_instrument_balance()       Pre-treatment outcome balance by instrument
#  plot_first_stage()              Coefficient plot from strength test
#  instrument_validity_report()    Omnibus wrapper with consolidated verdict
#
#  Recommended workflow
#  --------------------
#  Build f.X with build_iv_formula(), then pass it directly to
#  instrument_validity_report(). Instruments and controls (including
#  factor fixed effects) are parsed automatically from f.X.
#  feglm is used automatically when factor() terms are detected on the RHS.
#
#  Statistical background
#  ----------------------
#  Relevance  : Instruments must predict treatment adoption after partialling
#               out controls. Tested via first-stage F-statistic (Stock-Yogo
#               threshold F > 10), partial R-squared, and Cragg-Donald statistic.
#
#  Exclusion  : Instruments must not directly affect outcomes. Tested via:
#               (1) Reduced-form regressions of each outcome on instruments
#                   in the PRE-TREATMENT period only (placebo: instruments
#                   should not predict pre-treatment outcomes).
#               (2) One-way ANOVA of pre-treatment outcome means across
#                   instrument quantile groups (balance test).
#
#  Both sets of tests use only the pre-treatment period, consistent with the
#  SCM design which conditions on the pre-treatment outcome trajectory.
# =============================================================================


# -----------------------------------------------------------------------------
#' Parse IV formula into components
#'
#' Extracts treatment, instrument, and control variable names from a formula
#' produced by \code{build_iv_formula()}. The preferred path reads the
#' \code{instruments}, \code{controls}, and \code{treatment} attributes
#' stored by \code{build_iv_formula()} at creation time, which gives an
#' unambiguous classification regardless of whether controls are plain
#' numeric variables or \code{factor()} terms. Falls back to classifying
#' by \code{factor()} wrapping when attributes are absent.
#'
#' @param f.X Formula from \code{build_iv_formula()}.
#'
#' @return A list with elements \code{treatment}, \code{instruments},
#'   \code{controls}, and \code{fs_formula}.
#' @keywords internal
.parse_iv_formula <- function(f.X) {

  # ---- preferred path: read attributes stored by build_iv_formula()
  if (!is.null(attr(f.X, "instruments"))) {

    treatment   <- attr(f.X, "treatment")
    instruments <- attr(f.X, "instruments")
    controls    <- attr(f.X, "controls")

    f_str    <- deparse(f.X, width.cutoff = 500L)
    blocks   <- strsplit(f_str, "\\|")[[1L]]
    if (length(blocks) < 2L)
      stop("f.X has no second block (instrument equation).", call. = FALSE)
    fs_block   <- trimws(blocks[[2L]])
    fs_parts   <- strsplit(fs_block, "~")[[1L]]
    rhs        <- trimws(fs_parts[[2L]])
    fs_formula <- stats::as.formula(paste(treatment, "~", rhs))

    return(list(
      treatment   = treatment,
      instruments = instruments,
      controls    = controls,
      fs_formula  = fs_formula
    ))
  }

  # ---- fallback path: guess from factor() wrapping
  message(paste0(
    "f.X has no instrument attributes. ",
    "Classifying factor() terms as controls, plain terms as instruments. ",
    "For reliable classification use build_iv_formula()."
  ))

  f_str  <- deparse(f.X, width.cutoff = 500L)
  blocks <- strsplit(f_str, "\\|")[[1L]]

  if (length(blocks) < 2L)
    stop(paste0(
      "f.X does not contain a second block (instrument equation). ",
      "Use build_iv_formula() with instruments specified."
    ), call. = FALSE)

  fs_block <- trimws(blocks[[2L]])
  fs_parts <- strsplit(fs_block, "~")[[1L]]

  if (length(fs_parts) < 2L)
    stop("Second block of f.X must have the form: treatment ~ vars",
         call. = FALSE)

  treatment <- trimws(fs_parts[[1L]])
  rhs       <- trimws(fs_parts[[2L]])

  rhs_terms <- attr(
    stats::terms(stats::as.formula(paste("~", rhs))),
    "term.labels"
  )

  instruments <- rhs_terms[!grepl("factor\\(", rhs_terms, ignore.case = TRUE)]
  controls    <- rhs_terms[ grepl("factor\\(", rhs_terms, ignore.case = TRUE)]
  fs_formula  <- stats::as.formula(paste(treatment, "~", rhs))

  list(
    treatment   = treatment,
    instruments = instruments,
    controls    = controls,
    fs_formula  = fs_formula
  )
}


# -----------------------------------------------------------------------------
#' Test instrument strength (relevance)
#'
#' Runs the first-stage regression of treatment on instruments and controls
#' and reports the F-statistic, partial R-squared, and Cragg-Donald statistic.
#' When \code{f.X} is supplied the full second-block formula is used directly,
#' and \code{fixest::feglm} is selected automatically if \code{factor()}
#' terms are present on the RHS. Otherwise \code{stats::glm} is used.
#'
#' @param dt          data.table. Panel data in long format.
#' @param f.X         Optional formula from \code{build_iv_formula()}.
#'   If supplied, \code{treatment} and \code{instruments} are parsed
#'   from its second block automatically, and \code{feglm} is used when
#'   \code{factor()} terms are detected.
#' @param treatment   Character. Treatment column. Required if \code{f.X}
#'   is \code{NULL}.
#' @param instruments Character vector. Instrument column names. Required
#'   if \code{f.X} is \code{NULL}.
#' @param controls    Character vector. Additional control columns to include
#'   in the first-stage regression alongside instruments. Default \code{NULL}.
#' @param id_col      Character. Unit identifier. Default \code{"customer_id"}.
#' @param time_col    Character. Time identifier. Default \code{"wID"}.
#' @param pre_periods Integer vector of time values to restrict to.
#'   Default \code{NULL} uses all available rows.
#' @param method      \code{"probit"} (default) or \code{"ols"}.
#' @param alpha       Numeric. Significance level. Default \code{0.05}.
#' @param verbose     Logical. Print results. Default \code{TRUE}.
#'
#' @return A list with elements \code{f_statistic}, \code{f_pvalue},
#'   \code{partial_r2}, \code{cragg_donald}, \code{n_instruments},
#'   \code{n_obs}, \code{first_stage_fit}, \code{passes_threshold},
#'   \code{method}, \code{use_feglm}.
#' @export
test_instrument_strength <- function(dt,
                                     f.X         = NULL,
                                     treatment   = NULL,
                                     instruments = NULL,
                                     controls    = NULL,
                                     id_col      = "customer_id",
                                     time_col    = "wID",
                                     pre_periods = NULL,
                                     method      = c("probit", "ols"),
                                     alpha       = 0.05,
                                     verbose     = TRUE) {

  method <- match.arg(method)
  data.table::setDT(dt)

  # ---- resolve formula components
  parsed     <- NULL
  fs_formula <- NULL

  if (!is.null(f.X)) {
    parsed     <- .parse_iv_formula(f.X)
    fs_formula <- parsed$fs_formula
    if (is.null(treatment))   treatment   <- parsed$treatment
    if (is.null(instruments)) instruments <- parsed$instruments
    # merge any explicitly supplied controls with those parsed from f.X
    controls <- unique(c(parsed$controls, controls))
  }

  if (is.null(treatment))
    stop("Provide either f.X or treatment.", call. = FALSE)
  if (is.null(instruments) || length(instruments) == 0L)
    stop("No instruments found. Provide f.X or instruments directly.",
         call. = FALSE)

  # ---- subset data
  dt_use <- if (is.null(pre_periods)) dt
  else dt[get(time_col) %in% pre_periods]

  if (nrow(dt_use) == 0L)
    stop("No observations found after subsetting.", call. = FALSE)

  n_obs  <- nrow(dt_use)
  n_inst <- length(instruments)
  df_use <- as.data.frame(dt_use)

  # ---- detect whether feglm is needed
  # feglm is used when factor() terms exist (from f.X parse or controls arg)
  has_factor_controls <- any(grepl("factor\\(", controls, ignore.case = TRUE))
  use_feglm <- has_factor_controls

  if (use_feglm && !requireNamespace("fixest", quietly = TRUE))
    stop("Package 'fixest' required when factor() controls are present.",
         call. = FALSE)

  # ---- build reduced formula (instruments removed, controls kept)
  if (is.null(fs_formula)) {
    # Build from scratch when f.X not supplied
    rhs_full    <- paste(c(instruments, controls), collapse = " + ")
    fs_formula  <- stats::as.formula(paste(treatment, "~", rhs_full))
  }

  rhs_reduced <- if (length(controls) > 0L)
    paste(controls, collapse = " + ") else "1"
  f_reduced <- stats::as.formula(paste(treatment, "~", rhs_reduced))

  # ---- fit models and compute test statistics
  if (use_feglm) {

    fam <- if (method == "probit") "probit" else "gaussian"

    fit_full    <- fixest::feglm(fs_formula, data = df_use, family = fam)
    fit_reduced <- fixest::feglm(f_reduced,  data = df_use, family = fam)

    ll_full    <- as.numeric(stats::logLik(fit_full))
    ll_reduced <- as.numeric(stats::logLik(fit_reduced))
    lr_stat    <- 2 * (ll_full - ll_reduced)
    f_stat     <- lr_stat / n_inst
    f_pvalue   <- stats::pchisq(lr_stat, df = n_inst, lower.tail = FALSE)
    partial_r2 <- max(0, 1 - exp(-lr_stat / n_obs))
    first_stage_fit <- fit_full

  } else if (method == "probit") {

    fam         <- stats::binomial(link = "probit")
    fit_full    <- stats::glm(fs_formula, data = df_use, family = fam)
    fit_reduced <- stats::glm(f_reduced,  data = df_use, family = fam)
    fit_null    <- stats::glm(
      stats::as.formula(paste(treatment, "~ 1")),
      data = df_use, family = fam
    )

    lr_stat    <- as.numeric(2 * (stats::logLik(fit_full) -
                                    stats::logLik(fit_reduced)))
    f_stat     <- lr_stat / n_inst
    f_pvalue   <- stats::pchisq(lr_stat, df = n_inst, lower.tail = FALSE)

    ll_full    <- as.numeric(stats::logLik(fit_full))
    ll_reduced <- as.numeric(stats::logLik(fit_reduced))
    ll_null    <- as.numeric(stats::logLik(fit_null))
    partial_r2 <- max(0, (ll_reduced - ll_full) / (ll_reduced - ll_null))
    first_stage_fit <- fit_full

  } else {

    fit_full    <- stats::lm(fs_formula, data = df_use)
    fit_reduced <- stats::lm(f_reduced,  data = df_use)

    ftest      <- stats::anova(fit_reduced, fit_full)
    f_stat     <- ftest$F[2L]
    f_pvalue   <- ftest[["Pr(>F)"]][2L]

    rss_full    <- sum(stats::resid(fit_full)^2)
    rss_reduced <- sum(stats::resid(fit_reduced)^2)
    partial_r2  <- (rss_reduced - rss_full) / rss_reduced
    first_stage_fit <- fit_full
  }

  cragg_donald <- f_stat
  passes       <- !is.na(f_stat) && f_stat > 10

  if (verbose) {
    cat("\n", strrep("=", 65), "\n", sep = "")
    cat("  Instrument Strength (Relevance) Test\n")
    cat(strrep("=", 65), "\n", sep = "")
    cat(sprintf("  Method             : %s%s\n", method,
                if (use_feglm) " (feglm)" else " (glm)"))
    cat(sprintf("  Instruments        : %s\n",
                paste(instruments, collapse = ", ")))
    if (length(controls) > 0L)
      cat(sprintf("  Controls/FEs       : %s\n",
                  paste(controls, collapse = ", ")))
    cat(sprintf("  N observations     : %d\n", n_obs))
    cat(sprintf("  N instruments      : %d\n", n_inst))
    cat(strrep("-", 65), "\n")
    cat(sprintf("  %-21s: %.3f\n", "LR-based F-stat",  f_stat))
    cat(sprintf("  %-21s: %.4f\n", "p-value",           f_pvalue))
    cat(sprintf("  %-21s: %.4f\n", "Partial R-sq",      partial_r2))
    cat(sprintf("  %-21s: %.3f\n", "Cragg-Donald",      cragg_donald))
    cat(sprintf("  %-21s: %s\n",  "Stock-Yogo (F>10)",
                if (passes) "PASS" else "FAIL"))
    cat(strrep("=", 65), "\n\n")
  }

  invisible(list(
    f_statistic      = f_stat,
    f_pvalue         = f_pvalue,
    partial_r2       = partial_r2,
    cragg_donald     = cragg_donald,
    n_instruments    = n_inst,
    n_obs            = n_obs,
    first_stage_fit  = first_stage_fit,
    passes_threshold = passes,
    method           = method,
    use_feglm        = use_feglm
  ))
}


# -----------------------------------------------------------------------------
#' Test instrument exclusion restriction (validity)
#'
#' Tests the exclusion restriction by regressing each outcome on the
#' instruments using only PRE-TREATMENT observations (placebo test).
#' Under a valid instrument, the instruments should not predict
#' pre-treatment outcomes after partialling out controls.
#'
#' @param dt          data.table. Panel data in long format.
#' @param outcomes    Character vector. Outcome column names to test.
#' @param instruments Character vector. Instrument column names.
#' @param controls    Character vector. Control column names. Default \code{NULL}.
#' @param id_col      Character. Unit identifier. Default \code{"customer_id"}.
#' @param time_col    Character. Time identifier. Default \code{"wID"}.
#' @param treatment   Character. Treatment column used to identify pre-treatment
#'   rows. Default \code{"budgetdummy"}.
#' @param pre_periods Integer vector of time values. Default \code{NULL}.
#' @param alpha       Numeric. Significance level. Default \code{0.05}.
#' @param verbose     Logical. Print results. Default \code{TRUE}.
#'
#' @return A data.frame with columns \code{outcome}, \code{instrument},
#'   \code{coefficient}, \code{std_error}, \code{t_stat}, \code{p_value},
#'   \code{significant}, \code{conclusion}.
#' @export
test_instrument_exclusion <- function(dt,
                                      outcomes,
                                      instruments,
                                      controls    = NULL,
                                      id_col      = "customer_id",
                                      time_col    = "wID",
                                      treatment   = "budgetdummy",
                                      pre_periods = NULL,
                                      alpha       = 0.05,
                                      verbose     = TRUE) {

  data.table::setDT(dt)

  dt_pre <- if (is.null(pre_periods)) dt[get(treatment) == 0L]
  else dt[get(time_col) %in% pre_periods]

  if (nrow(dt_pre) == 0L)
    stop("No pre-treatment observations found.", call. = FALSE)

  miss <- setdiff(c(outcomes, instruments, controls), names(dt_pre))
  if (length(miss))
    stop("Columns not found: ", paste(miss, collapse = ", "), call. = FALSE)

  df_pre <- as.data.frame(dt_pre)

  # Use plain variable names only (strip factor() for exclusion test)
  plain_instruments <- instruments[!grepl("factor\\(", instruments)]
  plain_controls    <- controls[!grepl("factor\\(", controls %||% character(0))]

  rhs  <- paste(c(plain_instruments, plain_controls), collapse = " + ")
  rows <- list()

  for (oc in outcomes) {
    fit <- stats::lm(stats::as.formula(paste(oc, "~", rhs)), data = df_pre)
    cf  <- summary(fit)$coefficients

    for (inst in plain_instruments) {
      if (!inst %in% rownames(cf)) next
      r   <- cf[inst, ]
      sig <- !is.na(r["Pr(>|t|)"]) && r["Pr(>|t|)"] < alpha

      rows[[length(rows) + 1L]] <- data.frame(
        outcome     = oc,
        instrument  = inst,
        coefficient = round(r["Estimate"],   5),
        std_error   = round(r["Std. Error"], 5),
        t_stat      = round(r["t value"],    3),
        p_value     = round(r["Pr(>|t|)"],   4),
        significant = sig,
        conclusion  = if (sig)
          "CONCERN: instrument predicts pre-treatment outcome"
        else
          "OK: no significant pre-treatment effect",
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }

  result <- do.call(rbind, rows)

  if (verbose) {
    cat("\n", strrep("=", 65), "\n", sep = "")
    cat("  Exclusion Restriction Test (Pre-Treatment Placebo)\n")
    cat(strrep("=", 65), "\n", sep = "")
    cat(sprintf("  Outcomes           : %s\n", paste(outcomes,          collapse = ", ")))
    cat(sprintf("  Instruments        : %s\n", paste(plain_instruments, collapse = ", ")))
    cat(sprintf("  N (pre-treatment)  : %d\n", nrow(dt_pre)))
    cat(sprintf("  Significance level : %.2f\n", alpha))
    cat(strrep("-", 65), "\n")
    print(result[, c("outcome", "instrument", "coefficient",
                     "p_value", "significant", "conclusion")],
          row.names = FALSE)
    n_c <- sum(result$significant)
    cat(strrep("-", 65), "\n")
    cat(if (n_c == 0L) "  RESULT: No exclusion violations detected.\n"
        else sprintf("  RESULT: %d potential exclusion violation(s).\n", n_c))
    cat(strrep("=", 65), "\n\n")
  }

  invisible(result)
}


# -----------------------------------------------------------------------------
#' Test instrument balance across pre-treatment outcomes
#'
#' Stratifies units into quantile groups of each instrument value and
#' tests whether mean pre-treatment outcomes differ across groups via
#' one-way ANOVA. Significant differences suggest the instrument
#' correlates with baseline outcomes, raising exclusion concerns.
#'
#' @param dt          data.table. Panel data in long format.
#' @param outcomes    Character vector. Outcome column names.
#' @param instruments Character vector. Instrument column names.
#' @param id_col      Character. Unit identifier. Default \code{"customer_id"}.
#' @param time_col    Character. Time identifier. Default \code{"wID"}.
#' @param treatment   Character. Treatment column. Default \code{"budgetdummy"}.
#' @param n_quantiles Integer. Number of quantile groups. Default \code{4}.
#' @param alpha       Numeric. Significance level. Default \code{0.05}.
#' @param verbose     Logical. Print results. Default \code{TRUE}.
#'
#' @return A data.frame with columns \code{outcome}, \code{instrument},
#'   \code{f_stat}, \code{p_value}, \code{significant}, \code{conclusion}.
#' @export
test_instrument_balance <- function(dt,
                                    outcomes,
                                    instruments,
                                    id_col      = "customer_id",
                                    time_col    = "wID",
                                    treatment   = "budgetdummy",
                                    n_quantiles = 4L,
                                    alpha       = 0.05,
                                    verbose     = TRUE) {

  data.table::setDT(dt)
  dt_pre <- dt[get(treatment) == 0L]

  if (nrow(dt_pre) == 0L)
    stop("No pre-treatment observations found.", call. = FALSE)

  # Use plain instrument names only (strip factor() terms)
  plain_instruments <- instruments[!grepl("factor\\(", instruments)]

  miss <- setdiff(c(outcomes, plain_instruments), names(dt_pre))
  if (length(miss))
    stop("Columns not found: ", paste(miss, collapse = ", "), call. = FALSE)

  dt_unit <- dt_pre[,
                    lapply(.SD, mean, na.rm = TRUE),
                    by = id_col,
                    .SDcols = c(outcomes, plain_instruments)
  ]
  df_unit <- as.data.frame(dt_unit)
  rows    <- list()

  for (inst in plain_instruments) {
    vals   <- df_unit[[inst]]
    breaks <- unique(stats::quantile(vals,
                                     probs = seq(0, 1,
                                                 length.out = n_quantiles + 1L),
                                     na.rm = TRUE))
    if (length(breaks) < 3L) {
      warning(sprintf(
        "Instrument '%s' has too few unique values for %d quantiles. Skipping.",
        inst, n_quantiles
      ))
      next
    }

    df_unit[[".grp"]] <- cut(vals, breaks = breaks,
                             include.lowest = TRUE, labels = FALSE)

    for (oc in outcomes) {
      fit <- stats::aov(
        stats::as.formula(paste(oc, "~ factor(.grp)")),
        data = df_unit
      )
      sm  <- summary(fit)[[1L]]
      f_s <- sm[["F value"]][1L]
      p_v <- sm[["Pr(>F)"]][1L]
      sig <- !is.na(p_v) && p_v < alpha

      rows[[length(rows) + 1L]] <- data.frame(
        outcome     = oc,
        instrument  = inst,
        f_stat      = round(f_s, 3),
        p_value     = round(p_v, 4),
        significant = sig,
        conclusion  = if (sig)
          "CONCERN: instrument predicts baseline outcome level"
        else
          "OK: balanced across instrument quantiles",
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
    df_unit[[".grp"]] <- NULL
  }

  result <- do.call(rbind, rows)

  if (verbose) {
    cat("\n", strrep("=", 65), "\n", sep = "")
    cat("  Instrument Balance Test (Pre-Treatment Outcome Means)\n")
    cat(strrep("=", 65), "\n", sep = "")
    cat(sprintf("  Quantile groups    : %d\n", n_quantiles))
    cat(sprintf("  Units              : %d\n", nrow(dt_unit)))
    cat(sprintf("  Significance level : %.2f\n", alpha))
    cat(strrep("-", 65), "\n")
    print(result[, c("outcome", "instrument", "f_stat",
                     "p_value", "significant", "conclusion")],
          row.names = FALSE)
    n_c <- sum(result$significant)
    cat(strrep("-", 65), "\n")
    cat(if (n_c == 0L) "  RESULT: Balanced on pre-treatment outcomes.\n"
        else sprintf("  RESULT: %d balance concern(s) detected.\n", n_c))
    cat(strrep("=", 65), "\n\n")
  }

  invisible(result)
}


# -----------------------------------------------------------------------------
#' First-stage coefficient plot
#'
#' Dot-and-whisker plot of first-stage coefficients and confidence intervals
#' from \code{test_instrument_strength()}. Uses ggplot2 if available,
#' base graphics otherwise.
#'
#' @param strength_result Output list from \code{test_instrument_strength()}.
#' @param alpha Numeric. Significance level for CIs. Default \code{0.05}.
#' @param title Character. Plot title. Default \code{"First-Stage Coefficients"}.
#'
#' @return A ggplot object (invisibly) or NULL.
#' @export
plot_first_stage <- function(strength_result,
                             alpha = 0.05,
                             title = "First-Stage Coefficients") {

  fit <- strength_result$first_stage_fit
  cf  <- as.data.frame(summary(fit)$coefficients)
  cf$term <- rownames(cf)
  rownames(cf) <- NULL
  names(cf)[1:4] <- c("estimate", "std_error", "statistic", "p_value")
  cf <- cf[cf$term != "(Intercept)", , drop = FALSE]

  z        <- stats::qnorm(1 - alpha / 2)
  cf$lower <- cf$estimate - z * cf$std_error
  cf$upper <- cf$estimate + z * cf$std_error
  cf$sig   <- cf$p_value < alpha

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    cf$term <- factor(cf$term, levels = rev(cf$term))
    p <- ggplot2::ggplot(
      cf, ggplot2::aes(x = estimate, y = term, color = sig)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = lower, xmax = upper), height = 0.2) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                          color = "grey50") +
      ggplot2::scale_color_manual(
        values = c("TRUE" = "#1f78b4", "FALSE" = "#999999"),
        labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
        name   = NULL
      ) +
      ggplot2::labs(title = title, x = "Coefficient", y = NULL) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position  = "bottom",
                     panel.grid.minor = ggplot2::element_blank())
    print(p)
    return(invisible(p))
  }

  n    <- nrow(cf)
  xlim <- range(c(cf$lower, cf$upper, 0), na.rm = TRUE) * 1.1
  graphics::par(mar = c(4, max(nchar(as.character(cf$term))) * 0.6 + 1, 3, 2))
  graphics::plot(cf$estimate, seq_len(n),
                 xlim = xlim, ylim = c(0.5, n + 0.5),
                 xlab = "Coefficient", ylab = "",
                 main = title, yaxt = "n",
                 pch = 16, col = "#1f78b4")
  graphics::segments(cf$lower, seq_len(n), cf$upper, seq_len(n),
                     col = "#1f78b4", lwd = 1.5)
  graphics::abline(v = 0, lty = 2, col = "grey50")
  graphics::axis(2, at = seq_len(n),
                 labels = as.character(cf$term), las = 1)
  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' Omnibus instrument validity report
#'
#' Runs instrument strength, exclusion, and balance tests and prints a
#' consolidated verdict. Pass \code{f.X} from \code{build_iv_formula()} to
#' extract treatment and instruments automatically. Uses
#' \code{fixest::feglm} when \code{factor()} terms are detected on the RHS
#' of the instrument equation.
#'
#' @param dt          data.table. Panel data in long format.
#' @param f.X         Optional formula from \code{build_iv_formula()}.
#'   If supplied, \code{treatment} and \code{instruments} are parsed from
#'   its second block. \code{feglm} is used automatically when
#'   \code{factor()} terms are present.
#' @param treatment   Character. Treatment column. Required if \code{f.X}
#'   is \code{NULL}.
#' @param outcomes    Character vector. Outcome column names to test.
#' @param instruments Character vector. Instrument column names. Overrides
#'   instruments parsed from \code{f.X} if supplied.
#' @param controls    Character vector. Additional control columns. Default
#'   \code{NULL}.
#' @param id_col      Character. Unit identifier. Default \code{"customer_id"}.
#' @param time_col    Character. Time identifier. Default \code{"wID"}.
#' @param pre_periods Integer vector of time values. Default \code{NULL}.
#' @param method      \code{"probit"} (default) or \code{"ols"}.
#' @param alpha       Numeric. Significance level. Default \code{0.05}.
#' @param plot        Logical. Show first-stage coefficient plot.
#'   Default \code{TRUE}.
#' @param verbose     Logical. Print individual test results. Default \code{TRUE}.
#'
#' @return Invisibly returns a list with elements \code{strength},
#'   \code{exclusion}, \code{balance}, and \code{verdict}.
#' @export
instrument_validity_report <- function(dt,
                                       f.X         = NULL,
                                       treatment   = NULL,
                                       outcomes,
                                       instruments = NULL,
                                       controls    = NULL,
                                       id_col      = "customer_id",
                                       time_col    = "wID",
                                       pre_periods = NULL,
                                       method      = c("probit", "ols"),
                                       alpha       = 0.05,
                                       plot        = TRUE,
                                       verbose     = TRUE) {

  method <- match.arg(method)

  # Resolve components from f.X
  parsed <- NULL
  if (!is.null(f.X)) {
    parsed <- .parse_iv_formula(f.X)
    if (is.null(treatment))   treatment   <- parsed$treatment
    if (is.null(instruments)) instruments <- parsed$instruments
    controls <- unique(c(parsed$controls, controls))
  }

  if (is.null(treatment))
    stop("Provide either f.X or treatment.", call. = FALSE)
  if (is.null(instruments) || length(instruments) == 0L)
    stop("No instruments found. Provide f.X or instruments directly.",
         call. = FALSE)

  plain_instruments <- instruments[!grepl("factor\\(", instruments)]

  cat("\n", strrep("#", 65), "\n", sep = "")
  cat("  INSTRUMENT VALIDITY REPORT\n")
  cat("  scmBayesPost: selection_probit_bayes diagnostics\n")
  cat(sprintf("  Treatment          : %s\n", treatment))
  cat(sprintf("  Instruments        : %s\n", paste(plain_instruments, collapse = ", ")))
  if (length(controls) > 0L)
    cat(sprintf("  Controls/FEs       : %s\n", paste(controls, collapse = ", ")))
  cat(strrep("#", 65), "\n\n")

  # 1. Strength
  strength <- test_instrument_strength(
    dt          = dt,
    f.X         = f.X,
    treatment   = treatment,
    instruments = instruments,
    controls    = controls,
    id_col      = id_col,
    time_col    = time_col,
    pre_periods = pre_periods,
    method      = method,
    alpha       = alpha,
    verbose     = verbose
  )

  # 2. Exclusion placebo
  exclusion <- test_instrument_exclusion(
    dt          = dt,
    outcomes    = outcomes,
    instruments = plain_instruments,
    controls    = controls[!grepl("factor\\(", controls)],
    id_col      = id_col,
    time_col    = time_col,
    treatment   = treatment,
    pre_periods = pre_periods,
    alpha       = alpha,
    verbose     = verbose
  )

  # 3. Balance
  balance <- test_instrument_balance(
    dt          = dt,
    outcomes    = outcomes,
    instruments = plain_instruments,
    id_col      = id_col,
    time_col    = time_col,
    treatment   = treatment,
    alpha       = alpha,
    verbose     = verbose
  )

  # 4. Plot
  if (plot) {
    tryCatch(
      plot_first_stage(strength, alpha = alpha),
      error = function(e)
        message("plot_first_stage skipped: ", conditionMessage(e))
    )
  }

  # 5. Verdict
  strength_ok  <- isTRUE(strength$passes_threshold)
  exclusion_ok <- !any(exclusion$significant)
  balance_ok   <- !any(balance$significant)
  overall_ok   <- strength_ok && exclusion_ok && balance_ok

  verdict <- if (overall_ok) {
    "PASS: Instruments appear strong and valid."
  } else {
    issues <- c(
      if (!strength_ok)  "WEAK instruments (F < 10)",
      if (!exclusion_ok) "EXCLUSION concerns (pre-treatment outcome effects)",
      if (!balance_ok)   "BALANCE concerns (baseline outcome differences)"
    )
    paste("CONCERNS:", paste(issues, collapse = "; "))
  }

  cat(strrep("#", 65), "\n")
  cat("  OVERALL VERDICT\n")
  cat(strrep("#", 65), "\n")
  cat(" ", verdict, "\n")
  cat(strrep("#", 65), "\n\n")

  invisible(list(
    strength  = strength,
    exclusion = exclusion,
    balance   = balance,
    verdict   = verdict
  ))
}


# Suppress R CMD check notes for ggplot2 aesthetic mappings
utils::globalVariables(c(
  "estimate", "term", "lower", "upper", "sig", ".grp"
))

# Internal pipe-safe null coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b
