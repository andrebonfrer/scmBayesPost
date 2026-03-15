library(scmBayesPost)
library(data.table)

simulate_continuous_demo_data <- function(
    J0 = 100,
    Jc = 200,
    Tpre = 20,
    Tpost = 10,
    gamma_true = c(5, -2, 3),
    sigma_tau = 1,
    sigma_y = 2,
    seed = 123
) {

  set.seed(seed)

  Ttot <- Tpre + Tpost

  ids_t <- paste0("T", seq_len(J0))
  ids_c <- paste0("C", seq_len(Jc))
  ids   <- c(ids_t, ids_c)

  # moderators for treated units
  Z_t <- data.frame(
    id = ids_t,
    z1 = rnorm(J0),
    z2 = rnorm(J0)
  )

  tau_true <- gamma_true[1] +
    gamma_true[2]*Z_t$z1 +
    gamma_true[3]*Z_t$z2 +
    rnorm(J0, sd = sigma_tau)

  alpha <- rnorm(J0 + Jc, 100, 10)
  names(alpha) <- ids

  dt <- expand.grid(
    id  = ids,
    wID = seq_len(Ttot)
  )

  dt <- as.data.table(dt)

  dt[, treated := as.integer(id %in% ids_t)]

  # continuous intervention intensity
  dt[, budget_size := 0]

  post_rows <- dt$treated == 1 & dt$wID > Tpre

  dt[post_rows, budget_size := rgamma(sum(post_rows), 2, 1)]

  dt <- merge(dt, Z_t, by = "id", all.x = TRUE)

  dt[is.na(z1), `:=`(z1 = 0, z2 = 0)]

  dt[, tau_unit := 0]

  dt[id %in% ids_t, tau_unit :=
       tau_true[match(id[id %in% ids_t], ids_t)]]

  dt[, y := alpha[id] + tau_unit * budget_size + rnorm(.N, 0, sigma_y)]

  list(
    dt = dt,
    treated_ids = ids_t,
    tau_true = setNames(tau_true, ids_t),
    gamma_true = gamma_true
  )
}


build_test_weights_continuous <- function(dt, treated_ids) {

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


run_continuous_demo <- function() {

  sim <- simulate_continuous_demo_data()

  dt <- sim$dt

  W <- build_test_weights_continuous(dt, sim$treated_ids)

  gdata <- prepare_data_general(
    dta = dt,
    W   = W,
    y_name = "y",
    f.X = y ~ 1 + budget_size,
    f.Z = ~ z1 + z2,
    id_col = "id",
    time_col = "wID",
    tr_col = "budget_size",
    treat_type = "continuous",
    second_stage = "moderators",
    first_stage = "none",
    verbose = TRUE
  )

  post <- gibbs_postscm(
    gdata,
    n_iter = 2000,
    burn_in = 1000
  )

  list(
    sim = sim,
    gdata = gdata,
    post = post,
    W = W
  )
}


extract_tau_draws <- function(post, gdata, tr_name = NULL) {

  B <- post$beta_samples

  K  <- length(gdata$cov$Xcols)
  J0 <- gdata$cov$J0

  if (is.null(tr_name)) {
    k_tr <- gdata$cov$intX
  } else {
    k_tr <- match(tr_name, gdata$cov$Xcols)
  }

  if (is.na(k_tr)) stop("Treatment variable not found in Xcols")

  n_draws <- nrow(B)

  tau_draws <- matrix(NA_real_, n_draws, J0)

  for (d in seq_len(n_draws)) {

    beta_mat <- matrix(B[d, ], ncol = K, byrow = TRUE)

    tau_draws[d, ] <- beta_mat[, k_tr]

  }

  colnames(tau_draws) <- gdata$cov$treated_ids

  tau_draws
}

