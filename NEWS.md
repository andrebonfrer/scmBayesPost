# scmBayesPost 0.4.8

## New features

- `prepare_data_general(likelihood = "conditional", L, post)`: the synthetic-control
  (gap-series) likelihood, conditional on donor paths and weights. The pseudo-panel
  likelihood factorises as p(y_j | y_D; beta, sigma^2) p(y_D), and only the first factor
  involves the parameters. For each treated unit the within-period contrast (1, -w) is
  applied to the outcome and to every column of the design over event times `-L..-1` and
  `post`: period effects drop out, the intercept column becomes the unit-effect contrast
  (set to 1), the treatment column becomes the unit's own indicator, and covariates enter
  as treated-minus-synthetic differences. With `f.X = y ~ 1 + treatment`, OLS on each
  unit's rows is the post-minus-pre mean of the gap series, i.e. augMultiSynth's
  `mean(fit$tau[j, m, post + 1])`, and `gibbs_postscm()` draws are identical to those from
  the gap series passed in directly. All positive-weight donors are used (including
  treated units that adopt after the window); a donor treated inside the window is an
  error. Requires `first_stage = "none"`.

## Behaviour changes

- `likelihood = "pseudo_panel"` (the default, unchanged) now warns once per session that
  donor rows enter as observations and that, without unit and period effects in `f.X`, tau
  absorbs the pre-period level gap between a treated unit and its synthetic control. In an
  augMultiSynth simulation (8 outcomes, 300 treated units) this spec gave RMSE 0.38-0.54
  against the true unit effects and 4-9% coverage of 95% intervals, against 0.11-0.12 and
  91-95% under the conditional likelihood. Silence with
  `options(scmBayesPost.warn_pseudo_panel = FALSE)`.
- Known issue in the pseudo-panel spec, left unchanged for reproducibility: donors that
  are themselves columns of `W` are dropped without renormalising the remaining weights.

# scmBayesPost 0.4.6

## Bug fixes

- Fixed `build_iv_formula()` losing formula attributes (`treatment`,
  `instruments`, `controls`) after construction. `as.formula()` strips
  custom attributes, causing `.parse_iv_formula()` to fall back to the
  unreliable `factor()` heuristic and misclassify plain numeric controls
  as instruments. Fixed by using `structure()` to preserve attributes
  on the formula object after `as.formula()` is called.
