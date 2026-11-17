# Conditional (gap-series) likelihood for prepare_data_general()
#
# Synthetic-control model for treated unit j with donor weights w (sum 1):
#   y_j(t) = mu_j + sum_d w_d y_d(t) + x~_j(t)' beta + e_j(t),
# with donor paths and weights conditioned on. The pseudo-panel likelihood
# factorises as p(y_j | y_D; beta, sigma^2) p(y_D), and only the first factor
# involves the parameters. It is obtained by applying the within-period
# contrast c = (1, -w) to the outcome and to every column of the design:
#   * sum(c) = 0, so period effects drop out;
#   * the intercept column contrasts to 0 and is replaced by 1, the unit-effect
#     contrast mu_j = alpha_j - sum_d w_d alpha_d;
#   * the treatment column contrasts to D_j(t) when donors are untreated in the
#     window (checked below);
#   * other covariates enter as treated-minus-synthetic differences.
# Each treated unit contributes L + length(post) rows with weight 1.

#' Build gap-series (conditional-likelihood) panels
#'
#' Internal worker for `prepare_data_general(likelihood = "conditional")`.
#'
#' @param id_vec,t_vec,d_vec Unit id (character), time, and treatment vectors
#'   aligned with the rows of `X_mm`.
#' @param X_mm Design matrix from `f.X`.
#' @param y_vec Outcome vector.
#' @param W Weights matrix (rownames unit ids, colnames treated ids).
#' @param treated_ids Column names of `W`.
#' @param L Number of pre-treatment periods (event times -L..-1).
#' @param post Post event times (>= 0) defining the estimand.
#' @param treat_threshold Values of `d_vec` above this count as treated.
#' @param mass_tol Tolerance on `sum(w) == 1` for each treated unit's donors.
#' @return list with X_list, y_list, w_list, row_idx_list, X_idlist.
#' @keywords internal
.build_conditional_panels <- function(id_vec, t_vec, d_vec, X_mm, y_vec, W,
                                      treated_ids, L, post,
                                      treat_threshold = 0, mass_tol = 1e-6) {
  key <- paste(id_vec, t_vec, sep = "\r")
  if (anyDuplicated(key)) stop("dta has duplicated id-time rows; the conditional likelihood needs one row per unit and period.")

  ev   <- c(-L:-1, post)
  n_ev <- length(ev)
  int_col <- match("(Intercept)", colnames(X_mm))
  if (is.na(int_col)) {
    warning("f.X has no intercept: the conditional likelihood then omits the unit-effect contrast mu_j, ",
            "so tau is not the de-meaned gap estimator.", call. = FALSE)
  }

  J0 <- length(treated_ids)
  X_list <- y_list <- w_list <- row_idx_list <- xid <- vector("list", J0)
  bad_mass <- character(0)
  wn <- rownames(W)

  for (j0 in seq_len(J0)) {
    tr <- as.character(treated_ids[j0])
    rows_tr <- which(id_vec == tr)
    on <- rows_tr[!is.na(d_vec[rows_tr]) & d_vec[rows_tr] > treat_threshold]
    if (!length(on)) stop(sprintf("Treated unit '%s' has no treated periods in dta.", tr))
    T0  <- min(t_vec[on])
    per <- T0 + ev

    r_tr <- match(paste(tr, per, sep = "\r"), key)
    if (anyNA(r_tr)) {
      stop(sprintf("Event window %d..%d around adoption (time %s) is not fully observed for treated unit '%s'.",
                   min(ev), max(ev), format(T0), tr))
    }

    wcol <- W[, j0]
    keep <- is.finite(wcol) & wcol != 0 & wn != tr          # negative weights allowed (augmentation)
    dn <- wn[keep]; wd <- as.numeric(wcol[keep]); nd <- length(dn)
    if (!nd) stop(sprintf("Treated unit '%s' has no nonzero-weight donors in W.", tr))
    if (abs(sum(wd) - 1) > mass_tol) bad_mass <- c(bad_mass, tr)

    R_d <- match(paste(rep(dn, times = n_ev), rep(per, each = nd), sep = "\r"), key)
    if (anyNA(R_d)) {
      miss <- unique(rep(dn, times = n_ev)[is.na(R_d)])
      stop(sprintf("Donor(s) %s of treated unit '%s' are not observed over its event window.",
                   paste(utils::head(miss, 5), collapse = ", "), tr))
    }
    dd <- d_vec[R_d]
    if (any(!is.na(dd) & dd > treat_threshold)) {
      stop(sprintf(paste0("Donor(s) of treated unit '%s' are treated inside its event window; the contrast would not ",
                          "isolate tau. Restrict W to donors untreated through the last post period."), tr))
    }

    # contrast: treated row minus weighted donor rows, period by period
    Yd <- matrix(y_vec[R_d], nrow = nd)                       # nd x n_ev
    yj <- y_vec[r_tr] - as.numeric(crossprod(wd, Yd))
    Xj <- X_mm[r_tr, , drop = FALSE]
    Xd <- X_mm[R_d, , drop = FALSE]
    for (e in seq_len(n_ev)) {
      idx <- ((e - 1L) * nd + 1L):(e * nd)
      Xj[e, ] <- Xj[e, ] - as.numeric(crossprod(wd, Xd[idx, , drop = FALSE]))
    }
    if (!is.na(int_col)) Xj[, int_col] <- 1
    if (anyNA(yj) || anyNA(Xj)) stop(sprintf("Missing outcome or covariate values in the event window of '%s'.", tr))

    X_list[[j0]] <- Xj
    y_list[[j0]] <- yj
    w_list[[j0]] <- rep(1, n_ev)
    row_idx_list[[j0]] <- r_tr
    xid[[j0]] <- data.table::data.table(id = tr, wID = per, event_time = ev)
  }

  if (length(bad_mass)) {
    warning(sprintf(paste0("Donor weights do not sum to 1 for %d treated unit(s) (e.g. %s); the contrast then does ",
                           "not remove period effects."), length(bad_mass),
                    paste(utils::head(bad_mass, 3), collapse = ", ")), call. = FALSE)
  }

  list(X_list = X_list, y_list = y_list, w_list = w_list,
       row_idx_list = row_idx_list, X_idlist = data.table::rbindlist(xid))
}

.scmBayesPost_env <- new.env(parent = emptyenv())

.warn_pseudo_panel_once <- function() {
  if (isFALSE(getOption("scmBayesPost.warn_pseudo_panel", TRUE))) return(invisible())
  if (isTRUE(.scmBayesPost_env$pseudo_panel_warned)) return(invisible())
  .scmBayesPost_env$pseudo_panel_warned <- TRUE
  warning(paste0(
    "prepare_data_general(likelihood = \"pseudo_panel\") treats donor rows as observations and, unless f.X ",
    "includes unit and period effects, lets tau absorb level differences between each treated unit and its ",
    "synthetic control. Posterior intervals are then not those of the synthetic-control model. ",
    "likelihood = \"conditional\" uses the gap-series likelihood, which conditions on donor paths and weights. ",
    "(Shown once per session; silence with options(scmBayesPost.warn_pseudo_panel = FALSE).)"), call. = FALSE)
}
