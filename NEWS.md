# scmBayesPost 0.4.4

## Bug fixes and improvements

- `build_iv_formula()` now stores `treatment`, `instruments`, and
  `controls` as attributes on the returned formula object. This allows
  `.parse_iv_formula()` to reliably identify which variables are
  instruments versus controls without guessing from `factor()` wrapping,
  which was fragile when plain numeric controls were present alongside
  instruments.

- `.parse_iv_formula()` updated to read formula attributes as the
  preferred path. Falls back to `factor()` classification heuristic
  for hand-built formulas that lack attributes, with a message
  informing the user.
