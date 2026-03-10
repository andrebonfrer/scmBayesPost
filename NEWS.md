# scmBayesPost 0.2.3.9000

Development version.

* Added sampler `control` list support in `gibbs_postscm()`.
* Priors and tuning parameters can now be passed from the top-level API to internal samplers.
* `gibbs_sampling_simple()` and `gibbs_sampling_moderators()` now use resolved 
sampler controls rather than hard-coded prior values.
