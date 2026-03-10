# examples to run the model

library(scmBayesPost)
library(augMultiSynth)
library(data.table)

rm(list=ls())

res <- run_demo(N=1000, T=100, M=3, treated_eval=300, L=40, K=10,
                max_donors=500, pooled_adjustment=TRUE, nu=1)

fit <- res$fit
sim <- res$sim
unit_ids <- res$unit_ids

m <- 1  # outcome index to analyze

# Build long panel
Ymat <- sim$Yobs[[m]]  # N x T
N <- nrow(Ymat); TT <- ncol(Ymat)

dt <- data.table(
  id  = rep(unit_ids, times = TT),
  wID = rep(seq_len(TT), each = N),
  y   = as.numeric(Ymat)
)

# Treatment column:
# If you have binary in sim, use it. Commonly you'd construct it from treat_time.
# Here: tvg.dummy = 1 if t >= treat_time for that unit, else 0
tt <- sim$treat_time
dt[, tvg.dummy := {
  i <- match(id, unit_ids)
  as.integer(is.finite(tt[i]) & wID >= tt[i])
}]

W <- build_W_from_augMultiSynth(fit, id_universe = unit_ids, self_weight = 1)
# W now has rownames = unit_ids, colnames = treated_unit_ids

gdata <- prepare_data_general(
  dta = dt,
  W   = W,
  y_name = "y",
  f.X = y ~ 1 + tvg.dummy,     # minimal: intercept + treatment
  # no f.Z
  id_col = "id",
  time_col = "wID",
  tr_col = "tvg.dummy",
  treat_type = "binary",
  second_stage = "none",
  first_stage = "none",
  verbose = TRUE
)

gdata$Z_block <- matrix(1, nrow = gdata$cov$J0, ncol = 1)
colnames(gdata$Z_block) <- "(Intercept)"
gdata$Z.instruments <- NULL

# ensure cov fields exist
gdata$cov$intX <- match("tvg.dummy", gdata$cov$Xcols)  # index of treatment coefficient in X
if (is.na(gdata$cov$intX)) stop("tvg.dummy not found in Xcols")

set.seed(1)
post <- gibbs_postscm(
  gdata,
  n_iter = 200,
  burn_in = 100
)

# run diagnostics and comparisons

extract_unit_effect_draws <- function(post, gdata) {
  stopifnot(!is.null(post$beta_samples))
  B <- post$beta_samples

  K  <- length(gdata$cov$Xcols)
  J0 <- gdata$cov$J0
  k_tr <- gdata$cov$intX

  stopifnot(ncol(B) == K * J0)

  # column indices for treatment coefficient in each treated-unit block
  idx <- (0:(J0 - 1)) * K + k_tr

  tau_draws <- B[, idx, drop = FALSE]   # draws x J0
  colnames(tau_draws) <- as.character(gdata$tlist)  # if you stored treated ids here
  tau_draws
}


posterior_unit_summaries <- function(tau_draws, probs = c(0.025, 0.5, 0.975),
                                     treated_ids = NULL) {

  stopifnot(is.matrix(tau_draws) || is.data.frame(tau_draws))
  tau_draws <- as.matrix(tau_draws)

  J0 <- ncol(tau_draws)
  if (J0 < 1) stop("tau_draws has 0 columns.")

  # choose treated_ids in priority order
  if (is.null(treated_ids)) treated_ids <- colnames(tau_draws)
  if (is.null(treated_ids) || length(treated_ids) != J0) {
    treated_ids <- paste0("tr", seq_len(J0))
  }

  q <- t(apply(tau_draws, 2, stats::quantile, probs = probs, na.rm = TRUE))

  out <- data.frame(
    treated_id = treated_ids,
    mean = colMeans(tau_draws, na.rm = TRUE),
    sd   = apply(tau_draws, 2, stats::sd, na.rm = TRUE),
    q,
    row.names = NULL,
    check.names = FALSE
  )
  out
}




compare_to_truth <- function(summ, res, m = 1) {
  # treated IDs used in gdata
  tr_ids <- summ$treated_id

  # map treated IDs to sim indices
  idx <- match(tr_ids, res$unit_ids)
  if (anyNA(idx)) stop("Some treated_id not found in res$unit_ids. Check ID plumbing.")

  tau_true <- res$sim$tau_unit[idx, m]
  df <- summ
  df$tau_true <- tau_true
  df$error_mean <- df$mean - df$tau_true
  df$abs_error_mean <- abs(df$error_mean)

  df
}

tau_draws <- extract_unit_effect_draws(post, gdata)
colnames(tau_draws) <- colnames(W)   # or gdata$cov$treated_ids if you stored it
summ <- posterior_unit_summaries(tau_draws, probs=c(0.025,0.5,0.975))
head(summ)

cmp <- compare_to_truth(summ, res, m = 1)
with(cmp, cor(tau_true, mean))
with(cmp, sqrt(mean((mean - tau_true)^2)))
head(cmp[order(cmp$abs_error_mean, decreasing = TRUE), ], 10)

plot(cmp$tau_true, cmp$mean,
     xlab = "True unit effect",
     ylab = "Posterior mean unit effect")
abline(0, 1)

# columns named "2.5%" and "97.5%" if you used probs=c(0.025,0.5,0.975)
low  <- cmp[["2.5%"]]
high <- cmp[["97.5%"]]
mean(low <= cmp$tau_true & cmp$tau_true <= high)
