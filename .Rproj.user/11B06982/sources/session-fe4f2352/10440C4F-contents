# R/utils-math.R

#' Draw from multivariate normal with (possibly sparse) covariance
#'
#' Robust MVN sampler that accepts covariance matrices in either base R or
#' Matrix classes. Uses a Cholesky factorization with jitter if needed.
#'
#' @param n Integer number of draws.
#' @param mu Mean vector (numeric).
#' @param Sigma Covariance matrix (base matrix or Matrix).
#' @param jitter Nonnegative numeric; diagonal jitter added if chol fails.
#' @param max_tries Integer; number of jitter escalations.
#'
#' @return If n=1, a numeric vector. If n>1, a matrix with n rows.
#' @export
rMVNormCovariance <- function(n = 1L, mu, Sigma, jitter = 1e-8, max_tries = 6L) {
  n <- as.integer(n)
  if (n < 1L) stop("n must be >= 1.")
  mu <- as.numeric(mu)

  p <- length(mu)
  if (p == 0L) stop("mu is empty.")
  if (is.null(Sigma)) stop("Sigma is NULL.")

  # ensure matrix-like
  if (inherits(Sigma, "Matrix")) {
    # keep as Matrix for chol where possible
  } else {
    Sigma <- as.matrix(Sigma)
  }

  if (nrow(Sigma) != p || ncol(Sigma) != p) stop("Sigma dimension mismatch with mu.")

  L <- safe_chol(Sigma, jitter = jitter, max_tries = max_tries)

  Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  X <- sweep(Z %*% t(L), 2, mu, FUN = "+") # L is upper-tri from chol; using t(L) gives lower
  if (n == 1L) return(as.numeric(X))
  X
}

#' Safe Cholesky factorization with diagonal jitter
#'
#' @param A Symmetric PSD matrix (base or Matrix).
#' @param jitter Starting jitter.
#' @param max_tries Number of escalations.
#' @return Upper-triangular Cholesky factor.
#' @keywords internal
safe_chol <- function(A, jitter = 1e-8, max_tries = 6L) {
  if (inherits(A, "Matrix")) {
    # prefer Matrix::chol (gives "Cholesky" object) but base chol often fine after as.matrix
    # We'll use base chol on dense if needed.
  }

  j <- 0
  while (TRUE) {
    out <- tryCatch({
      if (inherits(A, "Matrix")) {
        # ensure symmetry
        # Matrix::chol can return "Cholesky" object; we want numeric upper-tri matrix
        ch <- Matrix::chol(A + (10^j) * jitter * Matrix::Diagonal(n = nrow(A)))
        as.matrix(ch)
      } else {
        chol(A + (10^j) * jitter * diag(nrow(A)))
      }
    }, error = function(e) e)

    if (!inherits(out, "error")) return(out)

    j <- j + 1
    if (j >= max_tries) {
      stop("Cholesky failed even after jitter. Matrix may be indefinite or badly scaled.")
    }
  }
}

#' Initialize a dense covariance matrix for LZ endogenous variables
#'
#' @param LZ Integer number of endogenous variables.
#' @param diag_value Positive numeric diagonal initial value.
#' @return Dense covariance matrix (LZ x LZ).
#' @export
initialize_dense_covariance <- function(LZ, diag_value = 1) {
  LZ <- as.integer(LZ)
  if (LZ < 1L) stop("LZ must be >= 1.")
  if (!is.finite(diag_value) || diag_value <= 0) stop("diag_value must be > 0.")
  diag(diag_value, LZ)
}

#' Safe inverse-Wishart draw wrapper
#'
#' Uses MCMCpack::riwish if available; otherwise errors with a clear message.
#'
#' @param nu Degrees of freedom.
#' @param S Scale matrix.
#' @return Sample from inverse-Wishart.
#' @export
riwish_safe <- function(nu, S) {
  if (!requireNamespace("MCMCpack", quietly = TRUE)) {
    stop("Dense covariance draws require MCMCpack (suggested dependency). Install MCMCpack or set Z_cov_dense=FALSE.")
  }
  MCMCpack::riwish(nu, S)
}

#' Replace a single column in a block-diagonal matrix efficiently
#'
#' Replaces the k-th column within each block of a block-diagonal matrix.
#' This is used to update X_block when a within-iteration regressor changes
#' (e.g., selection residuals).
#'
#' @param X_block A block-diagonal Matrix (typically from Matrix::bdiag).
#' @param new_col Numeric vector containing the replacement column values for
#'   the full stacked matrix (same nrow as X_block).
#' @param K Integer number of columns per block (assumed constant across blocks).
#' @param k Integer which column (1..K) within each block to replace.
#'
#' @return Updated X_block (same class as input where possible).
#' @export
replace_kth_column_block_diagonal <- function(X_block, new_col, K, k) {
  if (!inherits(X_block, "Matrix")) X_block <- Matrix::Matrix(X_block, sparse = TRUE)
  K <- as.integer(K); k <- as.integer(k)
  if (k < 1L || k > K) stop("k must be between 1 and K.")
  n <- nrow(X_block)
  if (length(new_col) != n) stop("new_col length must equal nrow(X_block).")

  # Determine block structure from X_block@Dim and assuming equal K columns per block
  p <- ncol(X_block)
  if (p %% K != 0L) stop("ncol(X_block) must be a multiple of K.")
  J0 <- p %/% K

  # We need to replace columns: (block 1 col k), (block 2 col k), ..., (block J0 col k)
  # global column indices:
  col_idx <- (0:(J0 - 1L)) * K + k

  # We also need row ranges per block. For bdiag of blocks with same number of rows nb:
  # We can infer rows per block if blocks are identical size; otherwise you should pass row_index list.
  # Here assume equal rows per block:
  if (n %% J0 != 0L) stop("nrow(X_block) must be divisible by J0 under equal-block-row assumption.")
  nb <- n %/% J0

  # Convert to a mutable sparse triplet
  Xs <- Matrix::as(X_block, "dgCMatrix")

  # For each block, replace entries in that column for that block's row range
  for (j in seq_len(J0)) {
    r0 <- (j - 1L) * nb + 1L
    r1 <- j * nb
    rows <- r0:r1
    cj <- col_idx[j]

    # Zero out existing column entries in those rows
    # easiest: assign vector in a dense way to that slice; Matrix supports this for dgCMatrix
    Xs[rows, cj] <- new_col[rows]
  }
  Xs
}
