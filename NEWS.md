## scmBayesPost 0.4.2

### New features

- `test_instrument_strength()` gains a `fixed_effects` argument that
  switches from `stats::glm` to `fixest::feglm` when fixed effects are
  specified. This resolves convergence failures with large datasets
  (N = 292k observations) where standard `glm` is unable to handle
  high-dimensional fixed effects. Example usage:

```r
  test_instrument_strength(
    dt            = bdt_horizon,
    treatment     = "budgetdummy",
    instruments   = "peer_adoption_rate",
    fixed_effects = "event_time",
    method        = "probit"
  )
```

  When `fixed_effects` is `NULL` (default), behaviour is unchanged and
  `stats::glm` is used as before.

- `fixest` added to `Suggests`. It is only required when
  `fixed_effects` is non-NULL in `test_instrument_strength()`. A
  clear error is raised if `fixest` is not installed and fixed effects
  are requested.

### Notes

- Unit (`customer_id`) fixed effects should not be included in the
  first-stage strength test. The incidental parameters problem causes
  severe bias in nonlinear probit models with many unit dummies, and
  unit-level heterogeneity is already absorbed by the SCM weights.
  Only time fixed effects (e.g. `fixed_effects = "event_time"`) are
  appropriate here.
  
