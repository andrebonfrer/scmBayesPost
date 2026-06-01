# scmBayesPost 0.4.0

## Breaking changes

- `prepare_data_general()` no longer returns `X_block`, `Y_block`, or
  `W` (the stacked block-diagonal matrices). These have been replaced
  by `X_list`, `y_list`, `w_list`, and `row_idx_list` — lists of J0
  unit-level matrices and vectors. Code that accesses `gdata$X_block`
  directly will need to be updated.

## Major architectural change — list-based samplers

All Gibbs samplers have been rewritten to operate on unit-level lists
directly, without ever constructing the stacked block-diagonal design
matrix `X_block`. On the motivating dataset (J0 = 1,879 treated units,
N_obs = 292,199 observations) the stacked matrix had 157 million rows,
making every `%*%` operation inside the Gibbs loop the dominant cost.

The new implementation eliminates:

- `Matrix::bdiag()` call in `prepare_data_general()` — the 157M-row
  sparse matrix is never built
- `Matrix::t(X_block)` inside the Gibbs loop
- `X_block %*% beta_bd` over 157M rows — replaced by J0 scalar
  multiplies (K=1 case)
- `t(nu_col) %*% W_block %*% nu_col` — replaced by unit-level
  `sum(wj * nuj^2)`
- All sparse matrix construction and manipulation inside the loop

## New internal helpers

- `.precompute_unit_stats()` — computes XtWX and XtWy per unit once
  before the loop. For K=1 these are scalars.
- `.compute_XtWy_units()` — recomputes XtWy per unit when y_tilde
  changes each iteration.
- `.sample_beta_units()` — for K=1: scalar arithmetic only (one
  division, one rnorm per unit). For K>1: Cholesky of [K x K] per unit.
- `.compute_wss()` — weighted sum of squared residuals, unit by unit.
- `.sample_z_star()` — fast inverse-CDF truncated normal sampler.

## API unchanged

`prepare_data_general()` and `gibbs_postscm()` function signatures and
the Gibbs fit output structure are unchanged. Only the intermediate
`gdata` object structure changes.
