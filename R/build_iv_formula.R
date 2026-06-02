#' Build IV formula for use with prepare_data_general
#'
#' Constructs a two-block formula of the form
#' \code{outcome ~ 0 + treatment | treatment ~ instruments + controls}
#' for use with \code{prepare_data_general()} and
#' \code{instrument_validity_report()}. Instrument and control variable
#' names are stored as formula attributes for reliable downstream parsing
#' by \code{.parse_iv_formula()}.
#'
#' @param outcome    Character. Outcome variable name.
#' @param treatment  Character. Treatment variable name.
#'   Default \code{"budgetdummy"}.
#' @param instruments Character vector. Instrument variable names.
#'   Default \code{NULL}.
#' @param controls   Character vector. Control variable names to include
#'   in the instrument equation only (not the outcome equation).
#'   Default \code{NULL}.
#' @param intercept  Logical. Include intercept in outcome equation.
#'   Default \code{FALSE} (recommended for SCM — weights absorb the level).
#'
#' @return A formula with attributes \code{treatment}, \code{instruments},
#'   and \code{controls}.
#' @export
build_iv_formula <- function(outcome,
                             treatment   = "budgetdummy",
                             instruments = NULL,
                             controls    = NULL,
                             intercept   = FALSE) {

  # Outcome RHS: treatment only, no controls
  outcome_rhs <- if (!intercept) paste("0 +", treatment) else treatment

  # Instrument equation RHS: instruments + controls
  instr_vars <- c(instruments, controls)
  instr_part <- if (length(instr_vars) > 0L) {
    paste(treatment, "~", paste(instr_vars, collapse = " + "))
  } else {
    NULL
  }

  # Assemble formula
  f <- if (!is.null(instr_part)) {
    stats::as.formula(
      paste0(outcome, " ~ ", outcome_rhs, " | ", instr_part)
    )
  } else {
    stats::as.formula(paste0(outcome, " ~ ", outcome_rhs))
  }

  # Store attributes for reliable parsing by .parse_iv_formula()
  attr(f, "treatment")   <- treatment
  attr(f, "instruments") <- instruments
  attr(f, "controls")    <- controls

  f
}
