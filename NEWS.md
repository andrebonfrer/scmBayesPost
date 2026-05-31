# scmBayesPost 0.2.6

## New features
- Bayesian probit first stage (`first_stage = "selection_probit_bayes"`) 
  via Albert-Chib augmentation in `gibbs_sampling_selection()`
- `prepare_data_general()` now populates `gdata$first_stage` with `X_fs`, 
  `d`, and MLE starting values when `first_stage = "selection_probit_bayes"`

## Bug fixes
- Fixed broken `@return` Rd link in `BalanceDiagnostics.R`
- `resolve_sampler_control()` now accepts `mu_delta_prior` and 
  `Sigma_delta_prior`
- Progress bar now closed properly after Gibbs iterations
