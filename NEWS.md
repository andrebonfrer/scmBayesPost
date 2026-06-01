# scmBayesPost 0.2.9

## Performance improvements

- `gibbs_sampling_selection()` and `gibbs_sampling_selection_moderators()`
  are now substantially faster when the first-stage probit is used with
  large datasets (tested on N = 292,199 observations):

  - **Truncated normal sampling**: replaced `truncnorm::rtruncnorm()` with
    a fast inverse-CDF sampler (`.sample_z_star()`). Treated and untreated
    index vectors are precomputed once before the loop, eliminating repeated
    `ifelse` dispatch across 290k+ observations per iteration. Approx 2-4x
    speedup on the z* block.

  - **Outcome cross-products**: `X_block'WX_block` and `X_block'W` are
    now precomputed once before the Gibbs loop since `X_block` is fixed
    across iterations. Only `X_block'W y_tilde` is recomputed per
    iteration. `Z_star` (the Kronecker moderator selector) is also
    precomputed in `gibbs_sampling_selection_moderators()`.

  - **nu_hat alignment**: replaced the per-iteration `data.table` merge
    in `.build_nu_hat_stacked()` with a one-time integer index map
    (`.build_nu_hat_index()`). Inside the loop, stacking nu_hat is now
    a single vector lookup (`nu_hat[nu_row_idx]`) rather than a 290k-row
    join. Approx 10-50x speedup on that operation.

## New internal functions

- `.sample_z_star()`: fast inverse-CDF truncated normal sampler.
- `.build_nu_hat_index()`: precomputes the integer alignment map from
  `X_block` rows to original `dta` rows, used in place of the
  per-iteration `data.table` merge.
