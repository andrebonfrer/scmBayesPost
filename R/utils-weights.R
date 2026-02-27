#' Validate and standardize a weights matrix (canonical SCM weights)
#'
#' Canonical format:
#' - rows: all units in the panel (unit ids) in a known order
#' - cols: treated units (treated ids)
#' - entries: donor weights for each treated unit (0 allowed)
#'
#' @param W Numeric matrix of weights.
#' @param require_dimnames Logical; require row/col names.
#' @param nonneg Logical; require nonnegative weights.
#' @param allow_self_weight Logical; if TRUE, allows treated unit row to be nonzero.
#'
#' @return A validated numeric matrix.
#' @export
validate_weights_matrix <- function(W,
                                    require_dimnames = TRUE,
                                    nonneg = TRUE,
                                    allow_self_weight = TRUE) {
  if (is.null(W) || !is.matrix(W)) stop("W must be a numeric matrix.")
  storage.mode(W) <- "double"

  if (require_dimnames) {
    if (is.null(rownames(W)) || is.null(colnames(W))) {
      stop("W must have rownames (unit ids) and colnames (treated unit ids).")
    }
    if (anyNA(rownames(W)) || anyNA(colnames(W))) stop("W dimnames cannot contain NA.")
  }

  if (nonneg && any(W < -1e-12, na.rm = TRUE)) {
    stop("W contains negative weights; set nonneg=FALSE if this is intended.")
  }

  if (!allow_self_weight) {
    if (!is.null(rownames(W)) && !is.null(colnames(W))) {
      common <- intersect(rownames(W), colnames(W))
      if (length(common)) {
        diag_vals <- vapply(common, function(id) W[id, id], numeric(1))
        if (any(abs(diag_vals) > 1e-12)) stop("Self-weights not allowed but W has diagonal mass.")
      }
    }
  }

  W
}

#' Align weights matrix to a desired unit id universe and ordering
#'
#' @param W Canonical weights matrix (rows=units, cols=treated).
#' @param id_universe Vector of unit ids (desired row order).
#' @param treated_ids Optional vector of treated ids (desired column order).
#' @param fill_value Value to fill for missing rows/cols (default 0).
#'
#' @return Matrix W aligned to id_universe (and treated_ids if provided).
#' @export
align_weights_matrix <- function(W,
                                 id_universe,
                                 treated_ids = NULL,
                                 fill_value = 0) {
  W <- validate_weights_matrix(W, require_dimnames = TRUE)

  id_universe <- as.character(id_universe)
  if (anyDuplicated(id_universe)) stop("id_universe contains duplicates.")

  if (is.null(treated_ids)) treated_ids <- colnames(W)
  treated_ids <- as.character(treated_ids)
  if (anyDuplicated(treated_ids)) stop("treated_ids contains duplicates.")

  # Build aligned matrix
  W2 <- matrix(fill_value,
               nrow = length(id_universe),
               ncol = length(treated_ids),
               dimnames = list(id_universe, treated_ids))

  # Copy overlap
  r <- intersect(rownames(W), id_universe)
  c <- intersect(colnames(W), treated_ids)
  if (length(r) && length(c)) {
    W2[r, c] <- W[r, c, drop = FALSE]
  }

  W2
}

#' Infer treated ids from a weights matrix
#'
#' @param W Canonical weights matrix.
#' @return Character vector of treated unit ids (colnames).
#' @export
treated_ids_from_W <- function(W) {
  W <- validate_weights_matrix(W, require_dimnames = TRUE)
  colnames(W)
}

#' Convert augMultiSynth-style donor lists into a dense weights matrix
#'
#' Expects a fit object containing:
#' - treated_unit_ids: vector length J0
#' - donor_ids: list length J0 (unit ids for donors per treated)
#' - weights: list length J0 (numeric donor weights aligned to donor_ids)
#'
#' Produces canonical W with rownames = id_universe and colnames = treated_unit_ids.
#'
#' @param fit augMultiSynth-like object.
#' @param id_universe Vector of all unit ids (row universe + order).
#' @param self_weight Numeric; value set on the treated unit row for each treated column.
#'
#' @return Canonical weights matrix.
#' @export
build_W_from_augMultiSynth <- function(fit, id_universe, self_weight = 1) {
  if (is.null(fit$treated_unit_ids) || is.null(fit$donor_ids) || is.null(fit$weights)) {
    stop("fit must contain treated_unit_ids, donor_ids, and weights.")
  }

  id_universe <- as.character(id_universe)
  tlist <- as.character(fit$treated_unit_ids)
  J0 <- length(tlist)

  W <- matrix(0, nrow = length(id_universe), ncol = J0,
              dimnames = list(id_universe, tlist))

  for (j0 in seq_len(J0)) {
    tr <- tlist[j0]
    donors <- as.character(fit$donor_ids[[j0]])
    wj <- as.numeric(fit$weights[[j0]])
    if (length(donors) != length(wj)) stop("donor_ids/weights length mismatch.")

    ridx <- match(donors, id_universe)
    if (anyNA(ridx)) {
      miss <- donors[is.na(ridx)]
      stop(sprintf("Some donor_ids not found in id_universe. Example: %s",
                   paste(utils::head(miss, 5), collapse = ", ")))
    }
    W[ridx, j0] <- wj

    tr_idx <- match(tr, id_universe)
    if (is.na(tr_idx)) stop(sprintf("treated id %s not in id_universe.", tr))
    W[tr_idx, j0] <- self_weight
  }

  validate_weights_matrix(W, require_dimnames = TRUE)
}

#' Coerce augsynth::multisynth object into canonical weights matrix
#'
#' multisynth documents that each column corresponds to a treated unit.
#' In implementation, rownames are set to unit IDs.
#'
#' @param msynth multisynth object.
#' @return Canonical weights matrix.
#' @export
weights_from_multisynth <- function(msynth) {
  if (is.null(msynth$weights)) stop("msynth must contain $weights.")
  W <- msynth$weights
  if (!is.matrix(W)) W <- as.matrix(W)
  validate_weights_matrix(W, require_dimnames = TRUE)
}
