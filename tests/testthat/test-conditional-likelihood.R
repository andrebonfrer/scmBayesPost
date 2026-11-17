# Conditional (gap-series) likelihood: prepare_data_general(likelihood = "conditional")

cl_panel <- function(seed = 3, n_tr = 5, n_late = 1, n_never = 30, Tn = 60) {
  set.seed(seed)
  N <- n_tr + n_late + n_never
  tt <- c(sample(20:30, n_tr, TRUE), rep(50L, n_late), rep(NA, n_never))
  f <- cumsum(stats::rnorm(Tn))                         # common factor
  lam <- stats::runif(N, 0.5, 1.5); a <- stats::rnorm(N, 0, 2)
  d <- data.table::CJ(id = seq_len(N), t = seq_len(Tn))
  d[, D := as.integer(!is.na(tt[id]) & t >= tt[id])]
  d[, z := stats::rnorm(.N)]
  d[, y := a[id] + lam[id] * f[t] + 0.3 * z + stats::rnorm(.N, sd = 0.5) + 1.5 * D]
  d[, id := as.character(id)]
  list(d = d, tt = tt, n_tr = n_tr, n_late = n_late, N = N)
}

cl_W <- function(p, seed = 4) {
  set.seed(seed)
  ids <- as.character(seq_len(p$N)); cols <- ids[seq_len(p$n_tr + p$n_late)]
  never <- ids[(p$n_tr + p$n_late + 1):p$N]; late <- ids[p$n_tr + seq_len(p$n_late)]
  W <- matrix(0, p$N, length(cols), dimnames = list(ids, cols))
  for (j in seq_along(cols)) {
    pool <- if (cols[j] %in% late) never else c(late, never)
    w <- stats::runif(length(pool)); w[sample(length(pool), floor(length(pool) / 2))] <- 0
    if (!(cols[j] %in% late)) w[1] <- 0.2             # make sure the late adopter is a donor
    W[pool, j] <- w / sum(w)
  }
  W
}

cl_prep <- function(d, W, f.X = y ~ 1 + D, ...) prepare_data_general(
  dta = d, W = W, y_name = "y", f.X = f.X, id_col = "id", time_col = "t", tr_col = "D",
  treat_type = "binary", second_stage = "none", first_stage = "none", verbose = FALSE,
  likelihood = "conditional", ...)

hand_gap <- function(p, W, L, post) {
  Y <- matrix(p$d$y, nrow = p$N, byrow = TRUE)           # d is sorted by id, t (CJ)
  ids <- as.character(seq_len(p$N)); ev <- c(-L:-1, post)
  data.table::rbindlist(lapply(colnames(W), function(tr) {
    j <- as.integer(tr); per <- p$tt[j] + ev
    data.table::data.table(id = tr, wID = seq_along(ev),
                           y = Y[j, per] - as.numeric(crossprod(W[, tr], Y[, per])),
                           D = as.integer(ev >= 0))
  }))
}

test_that("conditional panels are the gap series with an intercept and the unit's own treatment", {
  p <- cl_panel(); W <- cl_W(p); L <- 10; post <- 1:5
  g <- cl_prep(p$d, W, L = L, post = post)
  hg <- hand_gap(p, W, L, post)
  expect_equal(unlist(g$y_list), hg$y, tolerance = 1e-12)
  X <- do.call(rbind, g$X_list)
  expect_equal(unname(X[, "(Intercept)"]), rep(1, nrow(X)))
  expect_equal(unname(X[, "D"]), hg$D)
  expect_equal(unlist(g$w_list), rep(1, nrow(X)))
  expect_identical(g$cov$likelihood, "conditional")

  # the late adopter (a column of W) is used as a donor for the early units
  W0 <- W; W0[as.character(p$n_tr + 1), seq_len(p$n_tr)] <- 0
  W0[, seq_len(p$n_tr)] <- sweep(W0[, seq_len(p$n_tr)], 2, colSums(W0[, seq_len(p$n_tr)]), "/")
  g0 <- cl_prep(p$d, W0, L = L, post = post)
  expect_false(isTRUE(all.equal(g0$y_list[[1]], g$y_list[[1]])))
})

test_that("OLS on each unit's rows is the post-minus-pre mean gap and equals weighted two-way FE", {
  p <- cl_panel(); W <- cl_W(p); L <- 10; post <- 0:5
  g <- cl_prep(p$d, W, L = L, post = post)
  for (j in seq_len(p$n_tr)) {
    b <- stats::lm.fit(g$X_list[[j]], g$y_list[[j]])$coefficients[["D"]]
    yg <- g$y_list[[j]]; pre <- seq_len(L)
    expect_equal(b, mean(yg[-pre]) - mean(yg[pre]), tolerance = 1e-12)

    # weighted two-way FE on the pseudo-panel (treated + donors) gives the same point estimate
    tr <- colnames(W)[j]; dn <- rownames(W)[W[, j] > 0]
    per <- p$tt[as.integer(tr)] + c(-L:-1, post)
    pp <- p$d[id %in% c(tr, dn) & t %in% per]
    pp[, wt := ifelse(id == tr, 1, W[id, j])]
    pp[, Dj := as.integer(id == tr & t >= p$tt[as.integer(tr)])]
    b_fe <- stats::coef(stats::lm(y ~ factor(id) + factor(t) + Dj, data = pp, weights = wt))[["Dj"]]
    expect_equal(b_fe, b, tolerance = 1e-8)
  }
})

test_that("posterior draws equal those from the gap series passed directly", {
  p <- cl_panel(); W <- cl_W(p); L <- 10; post <- 1:5
  g <- cl_prep(p$d, W, L = L, post = post)
  hg <- hand_gap(p, W, L, post)
  Wi <- diag(ncol(W)); dimnames(Wi) <- list(colnames(W), colnames(W))
  op <- options(scmBayesPost.warn_pseudo_panel = FALSE); on.exit(options(op))
  gh <- prepare_data_general(dta = hg, W = Wi, y_name = "y", f.X = y ~ 1 + D, id_col = "id",
                             time_col = "wID", tr_col = "D", treat_type = "binary",
                             second_stage = "none", first_stage = "none", verbose = FALSE)
  run <- function(gd) { set.seed(7); utils::capture.output(r <- gibbs_postscm(gd, n_iter = 60, burn_in = 10)); r }
  r1 <- run(g); r2 <- run(gh)
  expect_equal(r1$beta_samples, r2$beta_samples, tolerance = 1e-10)
})

test_that("covariates enter as treated-minus-synthetic differences", {
  p <- cl_panel(); W <- cl_W(p); L <- 8; post <- 1:4
  g <- cl_prep(p$d, W, f.X = y ~ 1 + D + z, L = L, post = post)
  Z <- matrix(p$d$z, nrow = p$N, byrow = TRUE)
  j <- 2L; tr <- colnames(W)[j]; per <- p$tt[as.integer(tr)] + c(-L:-1, post)
  expect_equal(unname(g$X_list[[j]][, "z"]),
               Z[as.integer(tr), per] - as.numeric(crossprod(W[, tr], Z[, per])), tolerance = 1e-12)
})

test_that("conditional likelihood checks its inputs", {
  p <- cl_panel(); W <- cl_W(p)
  expect_error(cl_prep(p$d, W, L = 10), "needs L")
  expect_error(prepare_data_general(dta = p$d, W = W, y_name = "y", f.X = y ~ 1 + D | D ~ z,
                                    id_col = "id", time_col = "t", tr_col = "D",
                                    first_stage = "selection_probit", likelihood = "conditional",
                                    L = 10, post = 1:5, verbose = FALSE), "first_stage")
  # window running past the panel
  expect_error(cl_prep(p$d, W, L = 10, post = 1:40), "not fully observed")
  # a donor treated inside the window: put the late adopter (adopts at 50) on a window reaching 50
  expect_error(cl_prep(p$d, W, L = 10, post = 1:29), "treated inside")
  # weights not summing to one
  W2 <- W; W2[, 1] <- W2[, 1] * 0.9
  expect_warning(cl_prep(p$d, W2, L = 10, post = 1:5), "do not sum to 1")
})

test_that("the pseudo-panel likelihood warns once per session", {
  p <- cl_panel(); W <- cl_W(p)
  op <- options(scmBayesPost.warn_pseudo_panel = TRUE); on.exit(options(op))
  env <- get(".scmBayesPost_env", envir = asNamespace("scmBayesPost"))
  env$pseudo_panel_warned <- FALSE
  pp <- function() prepare_data_general(dta = p$d, W = W, y_name = "y", f.X = y ~ 1 + D,
                                        id_col = "id", time_col = "t", tr_col = "D", verbose = FALSE)
  expect_warning(pp(), "conditional")
  expect_silent(pp())
  expect_silent(cl_prep(p$d, W, L = 10, post = 1:5))
})

test_that("negative (augmented) weights work in the conditional likelihood only", {
  p <- cl_panel(); W <- cl_W(p); L <- 10; post <- 1:5
  # augmented-style weights: add a sum-zero correction with negative entries
  set.seed(11)
  never <- rownames(W)[(p$n_tr + p$n_late + 1):p$N]
  for (j in seq_len(p$n_tr)) {
    corr <- stats::rnorm(length(never), sd = 0.05); corr <- corr - mean(corr)
    W[never, j] <- W[never, j] + corr
  }
  expect_true(any(W < 0))
  g <- cl_prep(p$d, W, L = L, post = post)
  hg <- hand_gap(p, W, L, post)
  expect_equal(unlist(g$y_list), hg$y, tolerance = 1e-12)
  expect_error(prepare_data_general(dta = p$d, W = W, y_name = "y", f.X = y ~ 1 + D, id_col = "id",
                                    time_col = "t", tr_col = "D", verbose = FALSE),
               "negative weights")
  fit <- list(treated_unit_ids = colnames(W)[1:2],
              donor_ids = list(never, never), weights = list(W[never, 1], W[never, 2]),
              augment = list(method = "ridge"))
  Wb <- build_W_from_augMultiSynth(fit, id_universe = rownames(W))
  expect_true(any(Wb < 0))
  expect_error(build_W_from_augMultiSynth(fit, id_universe = rownames(W), allow_negative = FALSE),
               "negative")
})
