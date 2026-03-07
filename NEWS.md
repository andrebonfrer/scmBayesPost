# scmBayesPost 0.2.0

Major internal refactor.

* Replaced legacy `gibbs_sampling()` with modular samplers
  - `gibbs_sampling_simple()`
  - `gibbs_sampling_moderators()`

* Generalized `prepare_data_general()` to support moderators and future multi-outcome models.

* Clean separation between observation model (`f.X`) and heterogeneity model (`f.Z`).

* Simplified Gibbs dispatch logic in `gibbs_postscm()`.

* sim_example.R in scripts/ replaced to run a model with second stage included
