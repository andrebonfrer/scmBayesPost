# Self-contained: a small simulated staggered-adoption panel, a weights matrix
# with many exact zeros, and the Gibbs sampler run twice with the same seed.

make_panel <- function(n_tr = 6, n_dn = 40, Tn = 60, seed = 1) {
  set.seed(seed)
  ids <- seq_len(n_tr + n_dn); treat_t <- c(sample(30:40, n_tr, TRUE), rep(NA, n_dn))
  d <- data.table::CJ(customer_id = ids, wID = seq_len(Tn))
  d[, unit_mean := rep(stats::rnorm(length(ids), 10, 2), each = Tn)]
  d[, tt := treat_t[customer_id]]
  d[, budgetdummy := as.integer(!is.na(tt) & wID >= tt)]
  d[, y := unit_mean + 0.1 * wID + stats::rnorm(.N) - 0.5 * budgetdummy]
  d[, z1 := rep(stats::rnorm(length(ids)), each = Tn)]
  d[, tt := NULL]
  d
}

make_W <- function(d, n_tr, zero_share = 0.6, seed = 2) {
  set.seed(seed)
  ids <- as.character(unique(d$customer_id)); tr <- ids[seq_len(n_tr)]; dn <- setdiff(ids, tr)
  W <- matrix(0, length(ids), n_tr, dimnames = list(ids, tr))
  for (j in seq_len(n_tr)) {
    w <- stats::runif(length(dn)); w[sample(length(dn), floor(zero_share * length(dn)))] <- 0
    W[dn, j] <- w / sum(w)
  }
  W
}

run_sampler <- function(gdata, seed = 11) {
  set.seed(seed)
  gibbs_postscm(gdata, n_iter = 80, burn_in = 20,
                control = list(Sigma_gamma_prior = 1000, a_sigma_tau_prior = 5, b_sigma_tau_prior = 1))
}

prep <- function(d, W, ...) prepare_data_general(
  dta = d, W = W, y_name = "y", f.X = y ~ 1 + budgetdummy, f.Z = budgetdummy ~ z1,
  id_col = "customer_id", time_col = "wID", tr_col = "budgetdummy", treat_type = "binary",
  second_stage = "moderators", first_stage = "none", verbose = FALSE, ...)

test_that("pseudo-panels contain only positive-weight donors and give identical draws to the all-controls build", {
  n_tr <- 6; d <- make_panel(n_tr = n_tr); W <- make_W(d, n_tr)
  g <- prep(d, W)
  # rows per unit = (1 + positive-weight donors) x T
  pos <- colSums(W > 0)
  expect_equal(unname(sapply(g$X_list, nrow)), unname((1 + pos) * 60))
  # rebuild the old all-controls pseudo-panels by hand and compare draws
  ids <- as.character(d$customer_id); X_mm <- model.matrix(~ 1 + budgetdummy, d)
  ctrl <- setdiff(unique(ids), colnames(W))
  g0 <- g
  for (j in seq_len(n_tr)) {
    tr <- colnames(W)[j]; ridx <- which(ids %in% c(tr, ctrl))
    g0$X_list[[j]] <- X_mm[ridx, , drop = FALSE]; g0$y_list[[j]] <- d$y[ridx]
    w <- as.numeric(W[ids[ridx], j]); w[ids[ridx] == tr] <- 1; g0$w_list[[j]] <- w
  }
  p1 <- run_sampler(g); p0 <- run_sampler(g0)
  expect_equal(p1$beta_samples, p0$beta_samples, tolerance = 1.0e6)
  expect_equal(p1$gamma_samples, p0$gamma_samples, tolerance = 1.0e6)
})

test_that("w_min drops small donors, renormalises to the same mass, and leaves w_min = 0 untouched", {
  n_tr <- 6; d <- make_panel(n_tr = n_tr); W <- make_W(d, n_tr, zero_share = 0)
  g0 <- prep(d, W); g1 <- prep(d, W, w_min = 0.03)
  expect_true(all(sapply(g1$X_list, nrow) < sapply(g0$X_list, nrow)))
  for (j in seq_len(n_tr)) {
    w1 <- g1$w_list[[j]]; w0 <- g0$w_list[[j]]
    expect_equal(sum(w1[w1 < 1]), sum(w0[w0 < 1]), tolerance = 1e-10)   # donor mass preserved (treated unit has w = 1)
    expect_true(all(w1[w1 < 1] >= 0.03 * sum(w0[w0 < 1]) / sum(w0[w0 >= 0.03 & w0 < 1]) - 1e-12))
  }
  expect_error(prep(d, W, w_min = 0), NA)
})
