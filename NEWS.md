# scmBayesPost 0.4.1

## New features

- Added `instrument_tests.R` with five exported functions for a priori
  instrument validity and strength testing, intended to be run before
  `prepare_data_general()` and the Gibbs sampler:

  - `test_instrument_strength()`: tests instrument relevance via
    first-stage F-statistic (or likelihood-ratio statistic for probit),
    partial R^2, and Cragg-Donald statistic. Flags instruments that fail
    the Stock-Yogo threshold of F > 10.

  - `test_instrument_exclusion()`: pre-treatment placebo test. Regresses
    each outcome on the instruments using only pre-treatment observations
    (`budgetdummy == 0`). A significant coefficient indicates the
    instrument may have a direct effect on outcomes, violating the
    exclusion restriction.

  - `test_instrument_balance()`: stratifies units into quantile groups
    of the instrument and tests whether pre-treatment outcome means
    differ across groups via one-way ANOVA. Significant differences
    suggest the instrument is correlated with baseline outcome levels.

  - `plot_first_stage()`: coefficient plot with confidence intervals
    from a `test_instrument_strength()` result, showing which instruments
    drive first-stage relevance.

  - `instrument_validity_report()`: omnibus wrapper that runs all three
    tests in sequence and prints a consolidated PASS/FAIL verdict.
    Intended as the primary pre-estimation diagnostic when using
    `first_stage = "selection_probit_bayes"`.

  The recommended workflow is to run `instrument_validity_report()` on
  the raw panel data before calling `prepare_data_general()`, and only
  proceed to full Bayesian estimation if the instrument passes strength,
  exclusion, and balance checks.


