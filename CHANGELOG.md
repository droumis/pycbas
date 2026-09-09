# Changelog

## 0.2.0

Two new capabilities, and one fix that can change results you have already reported.
Read the Fixed section before comparing new output against old.

### Added

- **Higher-order criteria.** `CBASParams.criterion_order` sets what `criterion`
  counts: `0` trials, as before; `1` rewarded trials; `k` runs of `k` consecutive
  rewarded choices. Each subject's stream is truncated where it reaches that level of
  performance rather than at a fixed length, which matches subjects on performance
  instead of on exposure.

  A subject who never reaches the criterion is **not** truncated and contributes
  everything it has, so with a performance-based criterion the weakest subjects can
  contribute the most data, unevenly across groups. `subject_criteria` reports who fell
  short, and the app shows the same thing. This is a property of the design, not a bug,
  but it is worth measuring before drawing conclusions.

- **Several contingencies at once.** `run_cbas_multicontingency` counts each task
  contingency as its own set of hypotheses and corrects across all of them jointly, so
  the same arm sequence under two contingencies is two hypotheses and the
  false-positive budget is shared rather than spent separately in each.
  `result.sequences` is keyed by `(block, sequence)`.

  New: `load_cohort_with_contingencies`, `load_subject_data_with_contingencies`,
  `shared_contingency_blocks`, `build_multicontingency_count_matrix`,
  `criterion_trial`, `subject_criteria`, `reached_criterion`, `record_criteria`.

  Counting several contingencies multiplies the hypothesis space, and the step-down's
  memory is linear in it, so check `estimate_resources` first.

- **GUI support for both.** The criterion order is selectable, with the criterion's
  label, help text and step size following it. Multi-contingency folders are detected
  on load and get a contingency block selector and, where the info table offers one, a
  subject filter; results are labelled with their block.

  The criterion shortfall report works for multi-contingency data too, via
  `record_criteria`. It counts subject-contingency pairs and says how many subjects are
  affected in at least one contingency, which is the group most exposed to the
  falls-short hazard and previously the only one with no way to check it.

  The app's output is verified identical to the equivalent API call, both in the
  hypothesis set and bitwise in the adjusted p-values, on a single-contingency and a
  multi-contingency run.

- `pycbas --version`, and `python -m pycbas` as an equivalent to the console script.
- `pycbas.__version__`, and the version now has one source, `pycbas/__init__.py`, which
  `pyproject.toml` reads. Reading it the other way round reports whatever was last
  installed, which is misleading when a release can change numerical output.

- `CHANGELOG.md`, which did not exist before.

### Fixed

- **The observed and bootstrap statistics now agree bitwise, which changes results at
  ties.** The step-down asks whether a bootstrap statistic reaches the observed one
  using `>=`, and for a sequence seen in only a handful of subjects of one group that is
  an equality test: the statistic is an exact small rational and every resample drawing
  the same pattern reproduces the identical value. The two were computed by different
  routes, `std(ddof=1) / sqrt(n)` squared against `var / (n * (n - 1))`, which are
  algebraically identical and round two units in the last place apart. An entire block
  of resamples then failed the comparison at once.

  The effect was one-directional: it could only ever **add** significant sequences. In
  one analysis it moved a tied block's adjusted p-value by 0.105 and reported 21
  sequences as significant that should not have been.

  Both paths now build the statistic from raw sums through one function,
  `pycbas/_moments.py`. **If you re-run an analysis from 0.1.0 you may get fewer
  significant sequences.** That is the fix working. Results are unchanged wherever no
  tie sat at the alpha boundary, including the published rat validation, whose
  significance counts and Igor agreement are identical.

  The guarantee holds for **integer count matrices**. A rate-normalised matrix has no
  such guarantee; `compute_test_stats` now warns once when it sees one, and the
  step-down functions accept a `tie_rtol` for that case.

- `find_k_fwer_k1` raised `NameError` on every call.
- `python -m pycbas.cli` did nothing: the module defined `main()` and never called it,
  so it exited silently with no output and no error.
- The GUI launcher passed panel's `--autoreload`, a development flag, so the first thing
  a new user saw was a `FutureWarning` about the watchfiles package.
- The correlative statistic's tau numerator used `sum(X_dev**2 * Y_dev**2)` where the
  bootstrap used `sum((x_dev * y_dev)**2)`; these are not bitwise equal. Aligned. Full
  exactness is not reachable for a real-valued covariate, so the correlative path
  remains more exposed to ties than the comparative one.
- `estimate_resources` overstated chunked memory by about 1.8x: it doubled the column
  count and ignored the direction array. The real allocation is
  `resample_number x hypotheses x 9` bytes.
- The GUI resource estimate ignored `criterion_order` and `block_aware`, so switching to
  a performance-based criterion changed the real hypothesis count while the estimate did
  not move.
- An empty hypothesis space rendered in the GUI as "0 significant" after dividing by
  zero out of sight. It is a configuration error, usually a contingency filter matching
  no trials, and now says so.
- `load_cohort_with_contingencies` and `load_subject_data_with_contingencies` raise on a
  contingency change partway through a session, rather than silently disagreeing with
  the Igor reference, which assigns the whole session to the new block.
- Documentation: the parameter table listed `block_aware` as a `CBASParams` field when it
  is a pipeline argument; `docs/algorithm.md` gave the superseded formula for the
  standard error; the criterion's block-aware semantics were described as a cap on
  positions counted rather than a maximum start position; `find_k_fwer_chunked` was said
  to match `find_k_fwer` unconditionally when it always applies direction-conditional
  removal; the chunked memory saving was quoted as three different figures.

### Changed

- Prose throughout refers to artefacts rather than to people, and no longer explains
  itself in terms of one specific dataset's composition. `scripts/compare_with_david.py`
  is now `scripts/compare_with_reference.py`.

## 0.1.0

First release: comparative and correlative CBAS, the Romano-Wolf step-down with k-FWER
iteration, the chunked pipeline, the GUI, and validation against the Igor reference on
flies, humans and rats.
