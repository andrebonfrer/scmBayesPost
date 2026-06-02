# scmBayesPost 0.4.6

## Bug fixes

- Fixed `build_iv_formula()` losing formula attributes (`treatment`,
  `instruments`, `controls`) after construction. `as.formula()` strips
  custom attributes, causing `.parse_iv_formula()` to fall back to the
  unreliable `factor()` heuristic and misclassify plain numeric controls
  as instruments. Fixed by using `structure()` to preserve attributes
  on the formula object after `as.formula()` is called.
