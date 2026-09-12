# Changelog

## 0.2.0

Two additions: criteria that count performance instead of trials, and analysis across
several task contingencies at once, both usable from the GUI. One fix changes numerical
output. Re-running a 0.1.0 analysis can report fewer significant sequences, because the
step-down mishandled ties at the alpha boundary and resolved them towards significance.

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
- `NonIntegerCountWarning`, the category `compute_test_stats` warns with, so it can be
  filtered on its own.
- Tests for the GUI, which had none. They run on a synthetic cohort, so they need no
  unpublished data, and each pins one of the display and estimate faults listed below.
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
- **`estimate_resources` overstated the unchunked figure by the ratio of the enumerable
  sequence space to the observed hypothesis count**, which is two orders of magnitude on
  a sparse space. Both matrices are sized by the count matrix's columns, so the unchunked
  peak is about twice the chunked one, as the documentation says. `n_observed` now
  overrides the worst case for both figures, as it always claimed to.
- `estimate_resources` left `resample_number` out of its time estimate, reporting the
  same seconds for a run five times the size.
- The time estimate was a single measured point fitted through the origin, so it
  under-reported every small run: the GUI offered "~0s" for an analysis it then took
  several seconds to finish, on the same screen. It now carries a fixed cost as well as
  a per-hypothesis one and reproduces both published validation runs.
- **The GUI's plots took over the mouse wheel.** Bokeh makes wheel-zoom the active
  scroll tool, so scrolling the page over a plot silently rescaled its axes instead of
  moving the page, and a reader could not tell a wrecked view from a real one. This is
  how the committed overview screenshot came to show a sequence axis running to ten
  thousand and a negative `-log10(g-value)`, neither of which the code can produce. Zoom
  remains in the toolbar, deliberately, with reset alongside it.
- The Manhattan plot's axes were left to auto-range, spending most of the plot area on
  regions that cannot hold a point: a g-value is at most 1, so `-log10(g)` is never
  negative, and default padding on a log axis extended the sequence axis well past the
  last sequence. Both axes are now fitted to the data.
- The GUI resource estimate ignored `criterion_order` and `block_aware`, so switching to
  a performance-based criterion changed the real hypothesis count while the estimate did
  not move.
- **The GUI never counted hypotheses for a multi-contingency run at all.** The estimate
  was gated on the single-contingency subject list, which that path leaves empty, so the
  pane reported the enumerable sequence space: a figure that depends only on the alphabet
  and the sequence length, and that moves for neither the contingency blocks nor the
  subjects. It now counts what the run will test, and follows both controls.
- **The GUI crashed when given a small performance criterion.** Choosing an order above
  trials lowers the criterion field's minimum to 1, because a performance level is
  usually a small number, but the underlying parameter kept a lower bound of 10, so
  entering one raised `ValueError` inside the widget's watcher.
- **Multi-contingency results reported a sequence length of 2 for every hypothesis**, in
  the Manhattan plot's length bands, the significant-sequences table and the exported
  CSV's `length` column. Those keys are `(block, sequence)` pairs, so the length of the
  key is not the length of the sequence. All four now go through one helper.
- The GUI resource estimate and criterion shortfall report were not recomputed when the
  subject filter or the contingency block selection changed, and were never populated by
  a multi-contingency load at all, so they stayed blank or showed a previous selection's
  numbers at the moment they were most worth reading.
- After loading multi-contingency data the GUI's arm count, reward encoding and
  block-aware settings were applied to the analysis but not to the widgets, so the panel
  displayed something other than what the run would use, and editing one of those widgets
  silently reverted the detected value.
- A failed analysis reported nothing in the GUI, leaving the button on "Running...". The
  failure callback closed over the `except ... as` variable, which Python unbinds at the
  end of the block, and the callback is deferred onto the server's event loop.
- **`pycbas gui` launched panel through whatever `panel` was first on `PATH`**, which is
  not always the environment pycbas is installed in: a conda base env ahead of an
  activated venv is enough. The server started, the page loaded, and then every callback
  failed with "No module named 'pycbas'". It now runs panel through the current
  interpreter.
- **The GUI reported more sequences to test than were possible.** A multi-contingency run
  counts the sequence space once per contingency block, but the ceiling shown alongside
  the observed count was the space for a single block, so the pane read like
  "3,542 (of 1,884 possible)". `estimate_resources` takes `n_hypothesis_sets` for this.
- The multi-contingency loader set its detected arm count, reward encoding, block-aware
  flag and subject filter with parameter events discarded. Panel pushes a widget's value
  to the browser through those events, so the server was updated and the page was not:
  the arm count read 2, the checkboxes were clear and the subject filter rendered empty
  while the analysis used the detected values. The expensive recomputation is now
  deferred instead of the events being suppressed.
- A problem with multi-contingency data was reported as "info file is empty or
  malformed". The loader's diagnosis, such as a contingency changing partway through a
  session, was discarded by a blanket `except` and the fall-through then blamed the info
  file. Format detection is now structural, so a recognised file that cannot be loaded
  reports why.
- `pycbas._moments.tie_rtol_for` derived its tolerance from the number of finite matrix
  cells rather than the number of subjects whenever any entry was non-finite, inflating
  it by roughly the column count, and returned zero for a matrix of all NaN. A non-finite
  entry now defeats the exactness claim rather than supporting it.
- The non-integer count matrix warning was gated on a module-level flag that could not be
  reset, so it fired once per process and was unobservable to any test that did not run
  first. Deduplication is the warnings module's job now, and the warning has its own
  category, `NonIntegerCountWarning`, so a caller working in rates can silence exactly it.
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
