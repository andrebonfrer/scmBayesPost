# =============================================================================
#  prepare_data_general — updated first-stage handling
# =============================================================================
#
#  Changes vs the original:
#
#  1.  New first_stage option: "selection_probit_bayes"
#      Runs the frequentist probit to get starting values for delta, but
#      crucially populates gdata$first_stage with:
#        $X_fs   — first-stage design matrix [n_obs x p_fs]
#        $d      — treatment indicator vector [n_obs]
#        $delta0 — MLE starting values for delta [p_fs]
#        $fit    — the glm fit object (for diagnostics)
#
#  2.  The existing "selection_probit" path is unchanged (frequentist GR
#      correction; user appends GR manually via f.X).
#
#  3.  gdata$cov$first_stage now stores the matched first_stage string so
#      the Gibbs dispatcher can detect it.
#
#  Only the first-stage section and cov_meta are modified. Everything else
#  is identical to the original.
# =============================================================================


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
#' @param first_stage "none","selection_probit","selection_probit_bayes","treat_iv".
#'   \describe{
#'     \item{"none"}{No first stage.}
#'     \item{"selection_probit"}{Frequentist probit. Computes GR; user includes
#'       it in f.X manually. No Bayesian first-stage Gibbs step.}
#'     \item{"selection_probit_bayes"}{Fully Bayesian probit first stage via
#'       Albert-Chib augmentation. Populates gdata$first_stage with X_fs, d,
#'       and MLE starting values. The Gibbs sampler samples delta and z* jointly
#'       with the outcome equation parameters. Requires a second formula block
#'       in f.X: \code{outcome ~ covars | treatment ~ selection_covars}.}
#'     \item{"treat_iv"}{Placeholder hook for IV/CF first stage.}
#'   }
#' @param treat_threshold threshold for continuous treated definition.
#' @param verbose logical.
#'
#' @return list with Y_block, X_block, W, Z_block (optional), instruments
#'   (optional), and metadata. When first_stage = "selection_probit_bayes",
#'   gdata$first_stage contains X_fs, d, delta0, and fit.
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
                                 treat_type   = c("binary", "continuous"),
                                 second_stage = c("none", "moderators",
                                                  "moderators_iv"),
                                 first_stage  = c("none", "selection_probit",
                                                  "selection_probit_bayes",
                                                  "treat_iv"),
                                 treat_threshold = 0,
                                 verbose = TRUE) {

  treat_type   <- match.arg(treat_type)
  second_stage <- match.arg(second_stage)
  first_stage  <- match.arg(first_stage)

  check_panel_cols(dta, id_col, time_col, tr_col)
  if (!y_name %in% names(dta)) stop("y_name not found in dta.")
  if (is.null(f.X)) stop("Provide f.X.")

  # ---- validate treatment column
  tr_vec <- dta[[tr_col]]
  if (!is.numeric(tr_vec))
    stop(sprintf("Treatment column '%s' must be numeric.", tr_col),
         call. = FALSE)

  if (treat_type == "binary") {
    tr_vals <- unique(stats::na.omit(tr_vec))
    if (!all(tr_vals %in% c(0, 1)))
      stop(sprintf(
        "Treatment column '%s' must contain only 0/1 values when treat_type = 'binary'.",
        tr_col), call. = FALSE)
  }

  # ---- canonical unit universe + ordering
  id_universe <- as.character(unique(dta[[id_col]]))
  Tn <- data.table::uniqueN(dta[[time_col]])

  # ---- build / validate W
  if (is.null(W)) {
    if (is.null(res)) stop("Provide either W or res.")
    if (!is.null(res$weights) && is.matrix(res$weights)) {
      W <- res$weights
    } else if (!is.null(res$treated_unit_ids) && !is.null(res$donor_ids) &&
               !is.null(res$weights)) {
      W <- build_W_from_augMultiSynth(res, id_universe, self_weight = 1)
    } else {
      stop("Unrecognized res structure. Provide W directly or add a coercer.")
    }
  }
  W <- validate_weights_matrix(W, require_dimnames = TRUE)
  W <- align_weights_matrix(W, id_universe = id_universe)

  treated_ids <- colnames(W)

  # ---- derive treated ids from panel + cross-check
  tlist_panel <- treated_ids_from_panel(dta, id_col, tr_col,
                                        treat_type, treat_threshold)
  if (verbose) {
    message(sprintf(
      "prepare_data_general: |units|=%d, |treated(W)|=%d, |treated(panel)|=%d",
      length(id_universe), length(treated_ids), length(tlist_panel)
    ))
  }
  if (length(intersect(treated_ids, as.character(tlist_panel))) == 0) {
    warning("No overlap between treated ids implied by W columns and treated ids implied by panel treatment column.")
  }

  # ---- parse formulas
  fX_blocks <- parse_complex_formula(f.X)
  X_mm <- stats::model.matrix(fX_blocks[[1]], data = dta)

  intX <- match(tr_col, colnames(X_mm))
  if (is.na(intX))
    stop(sprintf(
      "Treatment column '%s' not found in X design matrix. X columns are: %s",
      tr_col, paste(colnames(X_mm), collapse = ", ")
    ))

  y_vec <- dta[[y_name]]

  # =========================================================================
  # ---- First-stage handling
  # =========================================================================
  first_stage_obj <- NULL

  if (first_stage == "selection_probit") {
    # ------------------------------------------------------------------
    # Frequentist probit + GR.  GR is stored but NOT appended here.
    # User should include GR in f.X if they want the control function.
    # ------------------------------------------------------------------
    if (treat_type != "binary")
      stop("selection_probit first_stage requires treat_type='binary'.")
    if (length(fX_blocks) < 2)
      stop("selection_probit requires a second formula block: f.X = outcome ~ ... | treatment ~ ...")

    if (verbose) message("Running first-stage probit (frequentist).")
    fs_fit <- stats::glm(fX_blocks[[2]],
                         family = stats::binomial(link = "probit"),
                         data   = dta)
    PR <- stats::predict(fs_fit, type = "response")
    GR <- dta[[tr_col]] * stats::dnorm(PR) / stats::pnorm(PR) -
      (1 - dta[[tr_col]]) * stats::dnorm(-PR) / stats::pnorm(-PR)

    first_stage_obj <- list(
      fit  = fs_fit,
      GR   = GR,
      GRX  = stats::model.matrix(fX_blocks[[2]], data = dta),
      GRY  = dta[[tr_col]]
    )

  } else if (first_stage == "selection_probit_bayes") {
    # ------------------------------------------------------------------
    # Fully Bayesian probit via Albert-Chib augmentation.
    # Populates first_stage_obj with everything the Gibbs sampler needs.
    # ------------------------------------------------------------------
    if (treat_type != "binary")
      stop("selection_probit_bayes requires treat_type='binary'.")
    if (length(fX_blocks) < 2)
      stop(paste0(
        "selection_probit_bayes requires a second formula block in f.X:\n",
        "  f.X = outcome ~ covars | treatment ~ selection_covars"
      ))

    if (verbose)
      message("Building first-stage design matrix for Bayesian probit.")

    # First-stage design matrix (RHS of second block)
    fs_formula <- fX_blocks[[2]]
    X_fs <- stats::model.matrix(fs_formula, data = dta)
    d    <- as.numeric(dta[[tr_col]])

    # Run frequentist probit to get MLE starting values for delta
    if (verbose)
      message("Fitting frequentist probit for starting values (delta0).")
    fs_fit <- stats::glm(fs_formula,
                         family = stats::binomial(link = "probit"),
                         data   = dta)
    delta0 <- as.numeric(stats::coef(fs_fit))

    if (verbose) {
      message(sprintf(
        "  First-stage probit: %d obs, %d covariates, pseudo-R2 = %.3f",
        length(d), ncol(X_fs),
        1 - fs_fit$deviance / fs_fit$null.deviance
      ))
    }

    first_stage_obj <- list(
      X_fs   = X_fs,      # [n_obs x p_fs] first-stage design matrix
      d      = d,         # [n_obs]        treatment indicator
      delta0 = delta0,    # [p_fs]         MLE starting values
      fit    = fs_fit     # glm object for diagnostics
    )

  } else if (first_stage == "treat_iv") {
    if (length(fX_blocks) < 2)
      stop("treat_iv requires a second formula block: | tr ~ instruments")
    if (verbose)
      message("first_stage='treat_iv' is a hook: implement IV/CF first-stage here.")
    first_stage_obj <- list(hook = TRUE, formula = fX_blocks[[2]])
  }

  # =========================================================================
  # ---- Build treated-specific pseudo-panels
  # =========================================================================
  id_controls <- setdiff(id_universe, treated_ids)
  J0 <- length(treated_ids)
  J  <- length(id_universe)

  pseudo_ids <- lapply(treated_ids, function(tr) c(tr, id_controls))

  dta_idx <- data.table::data.table(
    id  = as.character(dta[[id_col]]),
    wID = dta[[time_col]]
  )

  X_list <- vector("list", J0)
  y_list <- vector("list", J0)
  w_list <- vector("list", J0)

  if (verbose) message("Creating treated-specific pseudo-panels...")

  data.table::setDT(dta_idx)
  stopifnot(nrow(dta_idx) == nrow(X_mm))
  stopifnot(length(y_vec) == nrow(dta_idx))
  stopifnot(all(c("id", "wID") %in% names(dta_idx)))

  dta_idx[, id := as.character(id)]
  wn <- rownames(W)
  if (is.null(wn)) stop("W must have rownames = unit ids.")
  if (!is.character(wn)) wn <- as.character(wn)
  rownames(W) <- wn

  id_vec <- dta_idx[["id"]]

  row_idx_list <- vector("list", J0)
  for (j0 in seq_len(J0)) {
    ids_j <- pseudo_ids[[j0]]
    row_idx_list[[j0]] <- which(id_vec %chin% ids_j)
    if (length(row_idx_list[[j0]]) == 0L)
      stop(sprintf("Pseudo-panel %d has zero rows.", j0))
  }

  xid_chunks <- vector("list", J0)

  for (j0 in seq_len(J0)) {
    tr   <- as.character(treated_ids[j0])
    ridx <- row_idx_list[[j0]]

    Xj <- X_mm[ridx, , drop = FALSE]
    yj <- y_vec[ridx]

    wcol <- W[, j0]
    wj   <- as.numeric(wcol[id_vec[ridx]])
    wj[id_vec[ridx] == tr] <- 1

    X_list[[j0]] <- Xj
    y_list[[j0]] <- yj
    w_list[[j0]] <- wj

    xid_chunks[[j0]] <- dta_idx[ridx, .(id, wID)]
  }

  X_idlist <- data.table::rbindlist(xid_chunks, use.names = TRUE, fill = TRUE)

  # Do NOT materialise X_block / Y_block / W_block.
  # With J0 ~ 1879 and n_obs ~ 83k the stacked matrix has ~157M rows
  # making every %*% operation catastrophically slow.
  # All Gibbs samplers operate on X_list / y_list / w_list directly.

  # =========================================================================
  # ---- Second-stage moderators
  # =========================================================================
  Z_block      <- NULL
  Z_instruments <- NULL

  # Store first_stage string in cov_meta so the Gibbs dispatcher can detect it
  cov_meta <- list(
    Xcols        = colnames(X_mm),
    intX         = intX,
    J0           = J0,
    J            = J,
    T            = Tn,
    treated_ids  = treated_ids,
    id_universe  = id_universe,
    second_stage = second_stage,
    first_stage  = first_stage,   # <-- now stored as string
    treat_type   = treat_type,
    tr_col       = tr_col
  )

  if (second_stage != "none") {
    if (is.null(f.Z)) stop("second_stage requires f.Z.")

    fZ_blocks <- parse_complex_formula(f.Z)
    last_t    <- max(dta[[time_col]], na.rm = TRUE)

    Zdt <- dta[get(time_col) == last_t &
                 as.character(get(id_col)) %in% treated_ids]
    Zdt <- Zdt[!duplicated(as.character(Zdt[[id_col]]))]
    Zdt <- Zdt[match(treated_ids, as.character(Zdt[[id_col]]))]

    tt_Z     <- stats::terms(fZ_blocks[[1]])
    tt_Z_rhs <- stats::delete.response(tt_Z)
    Z_block  <- stats::model.matrix(tt_Z_rhs, data = Zdt)

    if (second_stage == "moderators_iv") {
      if (length(fZ_blocks) < 2)
        stop("moderators_iv requires IV blocks in f.Z.")
      dv  <- character(0)
      Zim <- list()
      for (b in 2:length(fZ_blocks)) {
        dv       <- c(dv, as.character(fZ_blocks[[b]])[2])
        Zim[[b - 1]] <- stats::model.matrix(fZ_blocks[[b]], data = Zdt)
      }
      Z_instruments <- list(dv = dv, Z.im = Zim)
    }
  }

  list(
    X_list         = X_list,          # list of J0 matrices [n_j x K]
    y_list         = y_list,          # list of J0 vectors  [n_j]
    w_list         = w_list,          # list of J0 vectors  [n_j]
    row_idx_list   = row_idx_list,    # list of J0 integer vectors (original dta rows)
    Z_block        = Z_block,
    Z.instruments  = Z_instruments,
    X_idlist       = X_idlist,
    dtaidx         = dta_idx,
    first_stage    = first_stage_obj,
    cov            = cov_meta
  )
}
