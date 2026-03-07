library(augMultiSynth)
library(scmBayesPost)

simulate_gamma_recovery_data <- function(
    J0 = 100,          # treated units
    Jc = 200,          # control/donor units
    Tpre = 20,
    Tpost = 10,
    gamma_true = c(0.5, 1.0, -0.75),   # intercept, z1, z2
    sigma_tau = 0.3,
    sigma_y = 1,
    seed = 123
) {
  set.seed(seed)

  Ttot <- Tpre + Tpost
  ids_t <- paste0("T", seq_len(J0))
  ids_c <- paste0("C", seq_len(Jc))
  ids <- c(ids_t, ids_c)

  # unit-level moderators for treated units
  Z_t <- data.frame(
    id = ids_t,
    z1 = rnorm(J0),
    z2 = rnorm(J0)
  )

  # true heterogeneous treatment effects
  tau_true <- gamma_true[1] +
    gamma_true[2] * Z_t$z1 +
    gamma_true[3] * Z_t$z2 +
    rnorm(J0, sd = sigma_tau)

  # baseline intercepts for all units
  alpha <- rnorm(J0 + Jc, mean = 0, sd = 1)
  names(alpha) <- ids

  # build panel
  dt <- expand.grid(
    id = ids,
    wID = seq_len(Ttot)
  )
  dt <- data.table::as.data.table(dt)

  # treatment indicator
  dt[, treated := as.integer(id %in% ids_t)]
  dt[, tvg.dummy := as.integer(treated == 1 & wID > Tpre)]

  # attach moderators to treated units only
  dt <- merge(dt, Z_t, by = "id", all.x = TRUE)
  dt[is.na(z1), `:=`(z1 = 0, z2 = 0)]   # controls can have zero/unused Z

  # generate outcome
  dt[, tau_unit := 0]
  dt[id %in% ids_t, tau_unit := tau_true[match(id[id %in% ids_t], ids_t)]]

  dt[, y := alpha[id] + tau_unit * tvg.dummy + rnorm(.N, sd = sigma_y)]

  list(
    dt = dt,
    treated_ids = ids_t,
    donor_ids = ids_c,
    tau_true = stats::setNames(tau_true, ids_t),
    gamma_true = gamma_true
  )
}


build_test_weights <- function(dt, treated_ids) {
  id_universe <- unique(as.character(dt$id))
  donor_ids <- setdiff(id_universe, treated_ids)

  W <- matrix(
    0,
    nrow = length(id_universe),
    ncol = length(treated_ids),
    dimnames = list(id_universe, treated_ids)
  )

  w_donor <- rep(1 / length(donor_ids), length(donor_ids))

  for (j in seq_along(treated_ids)) {
    W[donor_ids, j] <- w_donor
    W[treated_ids[j], j] <- 1
  }

  W
}


run_gamma_recovery_example <- function(
    J0 = 100,
    Jc = 200,
    Tpre = 20,
    Tpost = 10,
    gamma_true = c(0.5, 1.0, -0.75),
    n_iter = 2000,
    burn_in = 1000,
    seed = 123
) {
  sim <- simulate_gamma_recovery_data(
    J0 = J0,
    Jc = Jc,
    Tpre = Tpre,
    Tpost = Tpost,
    gamma_true = gamma_true,
    seed = seed
  )

  dt <- sim$dt
  W <- build_test_weights(dt, sim$treated_ids)

  gdata <- prepare_data_general(
    dta = dt,
    W = W,
    y_name = "y",
    f.X = y ~ 1 + tvg.dummy,
    f.Z = ~ z1 + z2,
    id_col = "id",
    time_col = "wID",
    tr_col = "tvg.dummy",
    treat_type = "binary",
    second_stage = "moderators",
    first_stage = "none",
    verbose = TRUE
  )

  post <- gibbs_postscm(
    gdata = gdata,
    n_iter = n_iter,
    burn_in = burn_in
  )

  list(
    sim = sim,
    W = W,
    gdata = gdata,
    post = post
  )
}


summarise_gamma_recovery <- function(post, gamma_true, probs = c(0.025, 0.975)) {
  Gs <- post$gamma_samples
  if (is.null(Gs)) stop("gamma_samples not found in posterior output.")

  out <- data.frame(
    term = colnames(Gs),
    truth = gamma_true,
    mean = colMeans(Gs),
    median = apply(Gs, 2, stats::median),
    lo = apply(Gs, 2, stats::quantile, probs = probs[1]),
    hi = apply(Gs, 2, stats::quantile, probs = probs[2]),
    row.names = NULL
  )

  out$covered <- out$truth >= out$lo & out$truth <= out$hi
  out$bias <- out$mean - out$truth
  out
}



res <- run_gamma_recovery_example()

gamma_summary <- summarise_gamma_recovery(
  post = res$post,
  gamma_true = c("(Intercept)" = 0.5, "z1" = 1.0, "z2" = -0.75)
)

gamma_summary


plot_gamma_recovery <- function(post, gamma_true) {
  Gs <- post$gamma_samples
  gs <- data.frame(
    term = colnames(Gs),
    mean = colMeans(Gs),
    lo = apply(Gs, 2, stats::quantile, probs = 0.025),
    hi = apply(Gs, 2, stats::quantile, probs = 0.975),
    truth = gamma_true
  )

  y <- seq_len(nrow(gs))
  op <- par(mar = c(5, 10, 4, 2) + 0.1)
  on.exit(par(op), add = TRUE)

  plot(gs$mean, y,
       xlim = range(c(gs$lo, gs$hi, gs$truth)),
       yaxt = "n",
       ylab = "",
       xlab = "Gamma",
       pch = 16)
  axis(2, at = y, labels = gs$term, las = 1)
  segments(gs$lo, y, gs$hi, y)
  points(gs$truth, y, pch = 4, cex = 1.3, lwd = 2)
  abline(v = 0, lty = 2)
}


plot_gamma_recovery(
  res$post,
  gamma_true = c("(Intercept)" = 0.5, "z1" = 1.0, "z2" = -0.75)
)
