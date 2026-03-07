#' Prepare data blocks for post-SCM estimation (generalized mscPost)
#'
#' Build stacked outcome/design matrices and weights for downstream estimation.
#' Accepts a canonical weights matrix with rownames=unit ids and colnames=treated ids.
#'
#' @param dta data.table panel containing id/time/treatment and covariates.
#' @param W weights matrix (rownames unit ids, colnames treated ids) OR NULL if res provided.
#' @param res optional fit object (e.g., multisynth, augMultiSynth-like) used to build W.
#' @param y_name outcome column in dta.
#' @param f.X formula for the per-observation design matrix (must include treatment column if desired).
#' @param f.Z optional moderators formula (treated-unit level) for second stage.
#' @param id_col unit id column name.
#' @param time_col time column name.
#' @param tr_col treatment column name (can be binary or continuous).
#' @param treat_type "binary" or "continuous".
#' @param second_stage "none","moderators","moderators_iv".
#' @param first_stage "none","selection_probit","treat_iv".
#' @param treat_threshold threshold for continuous treated definition when inferring treated set.
#' @param verbose logical.
#'
#' @return list with Y_block, X_block, W, Z_block (optional), instruments (optional), and metadata.
#' @export
prepare_data_general <- function(dta,
                                 W = NULL,
                                 res = NULL,
                                 y_name,
                                 f.X,
                                 f.Z = NULL,
                                 id_col = "id",
                                 time_col = "wID",
                                 tr_col = "tvg.dummy",
                                 treat_type = c("binary", "continuous"),
                                 second_stage = c("none", "moderators", "moderators_iv"),
                                 first_stage = c("none", "selection_probit", "treat_iv"),
                                 treat_threshold = 0,
                                 verbose = TRUE) {

  treat_type <- match.arg(treat_type)
  second_stage <- match.arg(second_stage)
  first_stage <- match.arg(first_stage)

  check_panel_cols(dta, id_col, time_col, tr_col)
  if (!y_name %in% names(dta)) stop("y_name not found in dta.")
  if (is.null(f.X)) stop("Provide f.X.")

  # ---- canonical unit universe + ordering
  id_universe <- as.character(unique(dta[[id_col]]))
  Tn <- data.table::uniqueN(dta[[time_col]])

  # ---- build / validate W
  if (is.null(W)) {
    if (is.null(res)) stop("Provide either W or res.")
    # try known res types
    if (!is.null(res$weights) && is.matrix(res$weights)) {
      W <- res$weights
    } else if (!is.null(res$treated_unit_ids) && !is.null(res$donor_ids) && !is.null(res$weights)) {
      W <- build_W_from_augMultiSynth(res, id_universe, self_weight = 1)
    } else {
      stop("Unrecognized res structure. Provide W directly or add a coercer.")
    }
  }
  W <- validate_weights_matrix(W, require_dimnames = TRUE)
  W <- align_weights_matrix(W, id_universe = id_universe)

  treated_ids <- colnames(W)

  # ---- derive treated ids from panel if desired and cross-check
  tlist_panel <- treated_ids_from_panel(dta, id_col, tr_col, treat_type, treat_threshold)
  if (verbose) {
    message(sprintf("prepare_data_general: |units|=%d, |treated(W)|=%d, |treated(panel)|=%d",
                    length(id_universe), length(treated_ids), length(tlist_panel)))
  }

  # optional consistency check (soft)
  if (length(intersect(treated_ids, as.character(tlist_panel))) == 0) {
    warning("No overlap between treated ids implied by W columns and treated ids implied by panel treatment column.")
  }

  # ---- parse formulas
  fX_blocks <- parse_complex_formula(f.X)
  X_mm <- stats::model.matrix(fX_blocks[[1]], data = dta)

  # index of treatment regressor inside X
  intX <- match(tr_col, colnames(X_mm))

  if (is.na(intX)) {
    stop(sprintf(
      "Treatment column '%s' not found in X design matrix. X columns are: %s",
      tr_col, paste(colnames(X_mm), collapse = ", ")
    ))
  }

  # outcome vector
  y_vec <- dta[[y_name]]

  # ---- First-stage functionality
  first_stage_obj <- NULL
  if (first_stage == "selection_probit") {
    if (treat_type != "binary") stop("selection_probit first_stage requires treat_type='binary'.")
    if (length(fX_blocks) < 2) stop("selection_probit requires a second formula block: | tr ~ ...")

    if (verbose) message("Running first-stage probit selection equation.")
    fs_fit <- stats::glm(fX_blocks[[2]], family = stats::binomial(link = "probit"), data = dta)
    PR <- stats::predict(fs_fit, type = "response")

    # generalized residual (GR) for binary treatment
    GR <- dta[[tr_col]] * stats::dnorm(PR) / stats::pnorm(PR) -
      (1 - dta[[tr_col]]) * stats::dnorm(-PR) / stats::pnorm(-PR)

    dta[, GR := GR]
    first_stage_obj <- list(fit = fs_fit,
                            GRX = stats::model.matrix(fX_blocks[[2]], data = dta),
                            GRY = dta[[tr_col]])

    # append GR to X_mm if GR column exists in model.matrix of primary equation
    # (you can also force it in f.X)
  } else if (first_stage == "treat_iv") {
    # Placeholder hook:
    # - if binary: could do 2SRI linear probability / probit with IV
    # - if continuous: classic first-stage regression with instruments
    # Implementation depends on how you want to specify instruments (formula blocks).
    if (length(fX_blocks) < 2) stop("treat_iv requires a second formula block: | tr ~ instruments")
    if (verbose) message("first_stage='treat_iv' is a hook: implement your preferred IV/CF first-stage here.")
    first_stage_obj <- list(hook = TRUE, formula = fX_blocks[[2]])
  }

  # ---- Build treated-specific pseudo-panels (treated + all controls)
  # mscPost-style: for each treated column j0, stack rows of treated + donors
  # Here we use: (treated unit) + (all non-treated units) — same as your existing approach.
  # Weight vector for each treated column defines the weights for those rows.
  id_controls <- setdiff(id_universe, treated_ids)  # pool
  J0 <- length(treated_ids)
  J <- length(id_universe)

  # fixed pseudo-panel order per treated unit: treated first, then controls
  pseudo_ids <- lapply(treated_ids, function(tr) c(tr, id_controls))
  n_per <- length(pseudo_ids[[1]]) * Tn

  # Precompute indexing table once
  dta_idx <- data.table::data.table(
    id  = as.character(dta[[id_col]]),
    wID = dta[[time_col]]
  )

  # Build X_list, y_list, w_list per treated unit
  X_list <- vector("list", J0)
  y_list <- vector("list", J0)
  w_list <- vector("list", J0)
  X_idlist <- NULL

  if (verbose) message("Creating treated-specific pseudo-panels...")

  # Defensive: ensure types and alignment
  data.table::setDT(dta_idx)
  stopifnot(nrow(dta_idx) == nrow(X_mm))
  stopifnot(length(y_vec) == nrow(dta_idx))
  stopifnot(all(c("id", "wID") %in% names(dta_idx)))

  # Ensure id is character (matching W rownames)
  dta_idx[, id := as.character(id)]
  wn <- rownames(W)
  if (is.null(wn)) stop("W must have rownames = unit ids.")
  if (!is.character(wn)) wn <- as.character(wn)
  rownames(W) <- wn

  # Build a fast lookup from id -> row weight for each j0 (avoid names indexing in a loop)
  # Also precompute row indices for each pseudo panel
  row_idx_list <- vector("list", J0)

  id_vec <- dta_idx[["id"]]  # plain vector
  for (j0 in seq_len(J0)) {
    ids_j <- pseudo_ids[[j0]]
    # integer row positions (fast + unambiguous)
    row_idx_list[[j0]] <- which(id_vec %chin% ids_j)
    if (length(row_idx_list[[j0]]) == 0L) {
      stop(sprintf("Pseudo-panel %d has zero rows. Check pseudo_ids and dta_idx$id.", j0))
    }
  }

  # Collect id/wID rows to rbind once (faster and avoids repeated rbind coercion)
  xid_chunks <- vector("list", J0)

  for (j0 in seq_len(J0)) {
    tr   <- as.character(treated_ids[j0])
    ridx <- row_idx_list[[j0]]

    # Subset X and y using integer rows
    Xj <- X_mm[ridx, , drop = FALSE]
    yj <- y_vec[ridx]

    # weights for rows in this pseudo panel
    # W[, j0] must be named by rownames(W) (unit ids)
    wcol <- W[, j0]
    # map row ids -> weights
    wj <- as.numeric(wcol[id_vec[ridx]])

    # convention: treated unit weight = 1
    wj[id_vec[ridx] == tr] <- 1

    # store
    X_list[[j0]] <- Xj
    y_list[[j0]] <- yj
    w_list[[j0]] <- wj

    # store id mapping rows
    xid_chunks[[j0]] <- dta_idx[ridx, .(id, wID)]
  }

  # One safe bind at end
  X_idlist <- data.table::rbindlist(xid_chunks, use.names = TRUE, fill = TRUE)

  # block objects
  X_block <- Matrix::bdiag(X_list)
  y_long <- unlist(y_list, use.names = FALSE)
  w_long <- unlist(w_list, use.names = FALSE)

  Y_block <- Matrix::Matrix(y_long, ncol = 1)
  W_block <- Matrix::Diagonal(x = w_long)

  # ---- second-stage moderators
  Z_block <- NULL
  Z_instruments <- NULL

  cov_meta <- list(
    Xcols = colnames(X_mm),
    intX = intX,
    J0 = J0,
    J = J,
    T = Tn,
    treated_ids = treated_ids,
    id_universe = id_universe,
    second_stage = second_stage,
    first_stage = first_stage,
    treat_type = treat_type,
    tr_col = tr_col
  )


  if (second_stage != "none") {
    if (is.null(f.Z)) stop("second_stage requires f.Z.")

    fZ_blocks <- parse_complex_formula(f.Z)

    # treated-unit level slice: take last time period for each treated id
    # (you can replace with your own aggregator)
    # treated-unit level slice
    last_t <- max(dta[[time_col]], na.rm = TRUE)

    Zdt <- dta[get(time_col) == last_t & as.character(get(id_col)) %in% treated_ids]

    # keep one row per treated unit
    Zdt <- Zdt[!duplicated(as.character(Zdt[[id_col]]))]

    # force exact treated-unit order
    Zdt <- Zdt[match(treated_ids, as.character(Zdt[[id_col]]))]

    # build moderator matrix from RHS only
    fZ_blocks <- parse_complex_formula(f.Z)
    tt_Z <- stats::terms(fZ_blocks[[1]])
    tt_Z_rhs <- stats::delete.response(tt_Z)

    Z_block <- stats::model.matrix(tt_Z_rhs, data = Zdt)


    if (second_stage == "moderators_iv") {
      if (length(fZ_blocks) < 2) {
        stop("moderators_iv requires IV blocks: f.Z like '... | Z1 ~ Q1 + Q2 | ...'")
      }
      dv <- character(0)
      Zim <- list()
      for (b in 2:length(fZ_blocks)) {
        dv <- c(dv, as.character(fZ_blocks[[b]])[2])
        Zim[[b - 1]] <- stats::model.matrix(fZ_blocks[[b]], data = Zdt)
      }
      Z_instruments <- list(dv = dv, Z.im = Zim)
    }
  }

  list(
    Y_block = Y_block,
    X_block = X_block,
    W = W_block,
    Z_block = Z_block,
    Z.instruments = Z_instruments,
    X_idlist = X_idlist,
    dtaidx = dta_idx,
    first_stage = first_stage_obj,
    cov = cov_meta
  )
}
