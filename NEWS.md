# scmBayesPost 0.4.3

## New features

- `instrument_validity_report()` and `test_instrument_strength()` now
  accept `f.X` directly from `build_iv_formula()`. Treatment variable
  and instruments are parsed automatically from the second block of
  `f.X` (the instrument equation `treatment ~ instruments + controls`),
  eliminating the need to re-specify them separately:

```r
  f.X <- build_iv_formula(
    outcome     = outcomevar,
    instruments = "peer_adoption_rate",
    controls    = "factor(event_time)"
  )

  validity <- instrument_validity_report(
    dt       = bdt_horizon,
    f.X      = f.X,
    outcomes = c(outcomevar),
    id_col   = "customer_id",
    time_col = "event_time",
    method   = "probit"
  )
```

- `fixest::feglm` is now selected automatically in
  `test_instrument_strength()` when `factor()` terms are detected on
  the RHS of the instrument equation. The `fixed_effects` argument
  introduced in 0.4.2 is no longer needed and has been removed.
  Passing `controls = "factor(event_time)"` in `build_iv_formula()`
  is sufficient.

- New internal helper `.parse_iv_formula()` classifies RHS terms from
  the instrument equation: plain variable names are treated as
  instruments, `factor()`-wrapped terms as controls or fixed effects.
  The full parsed formula is passed directly to `feglm` or `glm`
  without reconstruction.

- `test_instrument_exclusion()` and `test_instrument_balance()` now
  strip `factor()` terms from the instrument list automatically before
  running their tests, since these tests operate on plain numeric
  instrument values only.

## Breaking changes

- The `fixed_effects` argument added to `test_instrument_strength()`
  and `instrument_validity_report()` in 0.4.2 has been removed.
  Include fixed effects via `factor()` in the `controls` argument of
  `build_iv_formula()` instead.
