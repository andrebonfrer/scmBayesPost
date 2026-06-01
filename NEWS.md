## scmBayesPost v0.2.8

### What was broken

When using both `first_stage = "selection_probit_bayes"` and
`second_stage = "moderators"` together, `gibbs_postscm()` threw a
non-conformable matrix error and an erroneous warning that moderators
were not supported alongside the Bayesian probit first stage. Both
issues stemmed from the dispatcher routing the combined model to the
wrong internal sampler.

### What is fixed

A new internal sampler `gibbs_sampling_selection_moderators()` now
handles the full joint model correctly. The dispatcher routes
explicitly on all four combinations of first stage and moderator
presence. The combined specification is now fully supported with no
warnings.

### The joint model

**First stage (Albert-Chib probit):**  
`z*_i = x_fs_i' delta + nu_i`,  `d_i = 1[z*_i > 0]`

**Moderator equation (treatment heterogeneity):**  
`beta_{j, k_tr} = z_j' gamma + eta_j`

**Outcome equation:**  
`y_i = x_i' beta + rho * nu_hat_i + eps_i`

The `f.Z` moderators affect only the unit-specific treatment
coefficients and are entirely independent of the selection equation.

### Installation

```r
remotes::install_github("andrebonfrer/scmBayesPost@v0.2.8")
```
