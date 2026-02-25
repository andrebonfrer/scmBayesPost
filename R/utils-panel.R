#' Validate required columns in panel data
#'
#' @param dta data.table
#' @param id_col unit id column name
#' @param time_col time column name
#' @param tr_col treatment column name
#' @export
check_panel_cols <- function(dta, id_col, time_col, tr_col) {
  if (!data.table::is.data.table(dta)) stop("dta must be a data.table.")
  need <- c(id_col, time_col, tr_col)
  miss <- setdiff(need, names(dta))
  if (length(miss)) stop(sprintf("Missing columns in dta: %s", paste(miss, collapse = ", ")))
}

#' Compute treated unit ids from treatment column
#'
#' @param dta data.table
#' @param id_col name of unit id column
#' @param tr_col name of treatment column
#' @param treat_type "binary" or "continuous"
#' @param treat_threshold threshold to define treated when continuous (default >0)
#' @return vector of treated ids
#' @export
treated_ids_from_panel <- function(dta, id_col, tr_col,
                                   treat_type = c("binary", "continuous"),
                                   treat_threshold = 0) {
  treat_type <- match.arg(treat_type)
  if (treat_type == "binary") {
    tlist <- dta[, max(get(tr_col), na.rm = TRUE), by = ..id_col][V1 > 0, get(id_col)]
  } else {
    tlist <- dta[, max(get(tr_col), na.rm = TRUE), by = ..id_col][V1 > treat_threshold, get(id_col)]
  }
  as.vector(tlist)
}
