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

  id_col <- as.character(id_col)[1]
  tr_col <- as.character(tr_col)[1]

  if (!id_col %in% names(dta)) stop(sprintf("id_col '%s' not found.", id_col))
  if (!tr_col %in% names(dta)) stop(sprintf("tr_col '%s' not found.", tr_col))

  ids <- dta[[id_col]]
  tr  <- dta[[tr_col]]

  # max by id, NA-safe
  max_tr <- tapply(tr, ids, function(v) {
    if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  })

  # filter treated
  if (treat_type == "binary") {
    names(max_tr)[!is.na(max_tr) & max_tr > 0]
  } else {
    names(max_tr)[!is.na(max_tr) & max_tr > treat_threshold]
  }
}
