#' Parse complex formula specifications
#'
#' Supports formulas with '|' separators to denote additional model blocks.
#' For example:
#' - y ~ x1 + x2
#' - y ~ x1 + x2 | tvg.dummy ~ q1 + q2
#' - tvg.dummy ~ 1 + z1 + z2 | z1 ~ q1 + q2 | z2 ~ q3 + q4
#'
#' @param f A formula or formula-like object.
#' @return List of formula blocks.
#' @export
parse_complex_formula <- function(f) {
  if (is.null(f)) stop("Formula is NULL.")
  f <- stats::as.formula(f)

  # split on '|'
  rhs_txt <- paste(deparse(f), collapse = "")
  parts <- strsplit(rhs_txt, "\\|", fixed = FALSE)[[1]]
  parts <- trimws(parts)

  # First part must be a standard formula text
  blocks <- lapply(parts, function(s) stats::as.formula(s))
  blocks
}
