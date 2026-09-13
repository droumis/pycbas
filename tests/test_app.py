"""Regression tests for the GUI, on a synthetic cohort so they run anywhere.

Every failure pinned here reached a release because the app had no tests at all.
They are the kind that a reader of the code does not see: the displayed number is
plausible, the export is well-formed, the widgets look right, and the answer is
wrong anyway.

The app module is a singleton holding global widget state, which is why these tests
load data explicitly rather than sharing a fixture instance. They exercise the same
functions the callbacks do.
"""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("panel", reason="the GUI is an optional extra")

import pycbas._app as app  # noqa: E402

# A synthetic multi-contingency cohort. Exploration is block 0, then three
# alternation contingencies, each spanning whole sessions so the loader's
# mid-session guard is satisfied. Two genotype values, two lesion groups.
ARMS = {0: (-1, -1), 1: (2, 1), 2: (3, 2), 3: (4, 3)}
SESSION_BLOCKS = [0, 1, 1, 2, 2, 3, 3]


def _write_cohort(directory, n_subjects=12, trials=40, mid_session_change=False,
                  header=True, two_header_info=True):
    rng = np.random.default_rng(0)
    info = ["anInfo", "experiment,sex,genotype,lesion"] if two_header_info else \
           ["experiment,sex,genotype,lesion"]
    for subject in range(1, n_subjects + 1):
        rows = []
        for session, block in enumerate(SESSION_BLOCKS):
            centre, left = ARMS[block]
            for trial in range(trials):
                c, l = centre, left
                if mid_session_change and session == 1 and trial >= trials // 2:
                    c, l = ARMS[2]
                rows.append(f"{session},{rng.integers(0, 6)},{rng.integers(0, 2)},"
                            f"{c},{l},0,0,0,0")
        body = "\n".join(rows)
        prefix = "session,choice,reward,centre,left,t1,t2,f1,f2\n" if header else ""
        (directory / f"an{subject}.txt").write_text(prefix + body + "\n")
        genotype = "WT" if subject % 2 else "Mut"
        # Deliberately generic, and deliberately not an exact match for one of the
        # loader's keywords: the group mapper falls back to a substring test, and
        # that is the branch a real info table's wordier labels take.
        lesion = "Control" if subject <= n_subjects // 2 else "Lesion Group"
        info.append(f"0,Male,{genotype},{lesion}")
    (directory / "anInfo.txt").write_text("\n".join(info) + "\n")
    return directory


@pytest.fixture
def cohort_dir(tmp_path):
    return _write_cohort(tmp_path / "cohort", )


@pytest.fixture
def cohort_dir_factory(tmp_path):
    def make(name="cohort", **kwargs):
        target = tmp_path / name
        target.mkdir()
        return _write_cohort(target, **kwargs)
    return make


@pytest.fixture
def loaded(cohort_dir_factory):
    """Load a synthetic multi-contingency cohort through the folder callback."""
    directory = cohort_dir_factory()
    app.folder_selector.value = [str(directory)]
    app.load_from_folder(None)
    assert app.data_status.alert_type == "success", app.data_status.object
    return directory


# ---------------------------------------------------------------------------
# Sequence keys: (block, symbols) vs a bare sequence
# ---------------------------------------------------------------------------

class TestSequenceKeys:
    """`len(entry)` is 2 for every multi-contingency hypothesis, whatever its length.

    That is how the Manhattan plot's colour bands, the significance table's Length
    column and the exported CSV's `length` column all came to report 2 for every row
    of a multi-contingency run. `_seq_len` exists so no call site can get it wrong.
    """

    @pytest.mark.parametrize("entry,expected", [
        ((1, (3, 4, 5)), 3),
        ((2, (7,)), 1),
        ((3, (1, 2)), 2),
        ((3, 4, 5), 3),
        ((7,), 1),
        ((1, 2), 2),
    ])
    def test_length_counts_symbols_not_tuple_elements(self, entry, expected):
        assert app._seq_len(entry) == expected

    def test_label_carries_the_block_only_when_there_is_one(self):
        assert app._seq_label((1, (3, 4))) == "c1: 3-4"
        assert app._seq_label((3, 4)) == "3-4"

    def test_length_agrees_with_the_label_for_both_key_shapes(self):
        for entry in ((1, (3, 4, 5)), (3, 4, 5), (2, (9,)), (9,)):
            symbols = app._seq_label(entry).split(": ")[-1].split("-")
            assert len(symbols) == app._seq_len(entry)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

class TestMultiContingencyLoad:
    def test_widgets_match_the_state_the_run_will_use(self, loaded):
        """There is no reverse binding from state to widget, so the loader must set both.

        Setting only `app_state` left the panel showing something else, and the next
        edit of one of these widgets silently reverted the detected value.
        """
        assert app.num_arms_widget.value == app.app_state.num_arms
        assert app.encode_reward_widget.value == app.app_state.encode_reward is True
        assert app.block_aware_widget.value == app.app_state.block_aware is True

    def test_the_sync_emits_events_so_the_browser_hears_about_it(self, cohort_dir_factory):
        """Setting a widget inside `discard_events` never reaches the page.

        Panel pushes a value to the browser through the parameter event, so
        suppressing the event updates the server and leaves the rendered control
        showing the old value: a desynchronised panel that looks synchronised from
        Python. Asserting the event fires is the testable half of "the user sees it".
        """
        # The widgets are module globals shared across tests, so put them somewhere
        # other than the values the loader will set, or setting them is a no-op and
        # emits nothing whether the events are suppressed or not.
        app.num_arms_widget.value = 2
        app.encode_reward_widget.value = False
        app.block_aware_widget.value = False

        seen = []
        watchers = [
            widget.param.watch(lambda event, w=widget: seen.append(w), "value")
            for widget in (app.num_arms_widget, app.encode_reward_widget,
                           app.block_aware_widget)
        ]
        try:
            directory = cohort_dir_factory(name="events")
            app.folder_selector.value = [str(directory)]
            app.load_from_folder(None)
        finally:
            for widget, watcher in zip((app.num_arms_widget, app.encode_reward_widget,
                                        app.block_aware_widget), watchers):
                widget.param.unwatch(watcher)

        assert app.encode_reward_widget in seen
        assert app.block_aware_widget in seen

    def test_the_subject_filter_is_populated_visibly(self, cohort_dir_factory):
        """Its options and value must be set with events live, for the same reason.

        Asserting the values alone would not catch this: `discard_events` still
        performs the assignment, it only stops anyone hearing about it, so the server
        state looks right while the rendered selector stays empty.
        """
        app.subject_filter.options = []
        app.subject_filter.value = []

        seen = []
        watcher = app.subject_filter.param.watch(lambda event: seen.append(event), "value")
        try:
            directory = cohort_dir_factory(name="filtervisible")
            app.folder_selector.value = [str(directory)]
            app.load_from_folder(None)
        finally:
            app.subject_filter.param.unwatch(watcher)

        assert list(app.subject_filter.options) == ["Mut", "WT"]
        assert list(app.subject_filter.value) == ["Mut", "WT"]
        assert seen, "the selector's value was set without an event, so the page never updates"

    def test_the_filter_keeps_every_subject_with_its_own_group(self, loaded):
        """Filtering must not shift the correspondence between subject and group.

        The fixture makes the two facts independent: genotype alternates by subject
        number and group splits the cohort in half, so a filter on genotype keeps an
        interleaved subset and any off-by-one in the pairing shows up. Asserting only
        that the subject count fell, as the resource-estimate test does, passes just
        as happily when every surviving label belongs to the wrong animal.
        """
        truth = {s.id: s.meta["lesion"] for s in app.app_state._all_cohort}
        wanted = [v for v in app.subject_filter.options if v == "WT"]
        assert wanted, "expected a WT genotype in the fixture"

        app.subject_filter.value = wanted
        cohort = app.app_state.cohort
        assert 0 < len(cohort) < len(truth), "expected a proper subset"

        labels = app.app_state.group_labels
        assert len(labels) == len(cohort)
        for subject, label in zip(cohort, labels):
            expected = 0 if "control" in truth[subject.id].lower() else 1
            assert label == expected, (
                f"{subject.id} came from a {truth[subject.id]!r} row but is "
                f"labelled {label}")
            assert subject.meta["genotype"] == "WT"

    def test_block_and_filter_controls_appear(self, loaded):
        assert app.block_row.visible
        assert app.subject_filter_row.visible
        assert list(app.block_selector.value) == [1, 2, 3]

    def test_a_broken_cohort_reports_its_real_cause(self, cohort_dir_factory):
        """A mid-session contingency change must not read as a malformed info file.

        The loader raises a precise ValueError for this. It used to be swallowed by a
        blanket except, and the fall-through to the single-header parser then blamed
        `anInfo.txt`, which is the wrong file and the wrong diagnosis.
        """
        directory = cohort_dir_factory(name="broken", mid_session_change=True)
        app.folder_selector.value = [str(directory)]
        app.load_from_folder(None)
        message = str(app.data_status.object)
        assert app.data_status.alert_type == "danger"
        assert "partway through session" in message
        assert "empty or malformed" not in message

    def test_a_single_header_table_is_not_claimed_by_this_loader(self, cohort_dir_factory):
        """Detection must not swallow folders belonging to the other loaders."""
        directory = cohort_dir_factory(name="single", two_header_info=False)
        assert app._has_two_header_info_table(directory / "anInfo.txt") is False

    def test_the_two_header_table_is_recognised(self, loaded):
        assert app._has_two_header_info_table(loaded / "anInfo.txt") is True

    def test_headerless_subject_files_still_load(self, cohort_dir_factory):
        directory = cohort_dir_factory(name="noheader", header=False)
        app.folder_selector.value = [str(directory)]
        app.load_from_folder(None)
        assert app.data_status.alert_type == "success", app.data_status.object


# ---------------------------------------------------------------------------
# The resource estimate
# ---------------------------------------------------------------------------

class TestResourceEstimate:
    def test_it_counts_observed_hypotheses_not_the_enumerable_space(self, loaded):
        """The estimate gated on `subjects_data`, which this path sets to [].

        So no multi-contingency run ever counted its hypotheses, and the pane
        reported the enumerable sequence space instead: a number that depends only on
        the alphabet and the length, and moves for neither the block selection nor
        the subjects.
        """
        est = app.app_state.get_resource_estimate()
        assert est["observed_sequences"] is not None
        assert est["observed_sequences"] < est["total_sequences"]
        assert app.resource_estimate_pane.visible

    def test_it_follows_the_block_selection(self, loaded):
        """Each block is its own hypothesis set, so this is the biggest control there is."""
        app.block_selector.value = [1, 2, 3]
        with_all = app.app_state.get_resource_estimate()["observed_sequences"]
        app.block_selector.value = [1]
        with_one = app.app_state.get_resource_estimate()["observed_sequences"]
        assert with_one < with_all, (with_one, with_all)

    def test_it_follows_the_subject_filter(self, loaded):
        app.subject_filter.value = list(app.subject_filter.options)
        everyone = app.app_state.n_subjects
        app.subject_filter.value = list(app.subject_filter.options)[:1]
        assert app.app_state.n_subjects < everyone

    def test_the_reported_ceiling_is_not_below_the_observed_count(self, loaded):
        """Each block counts the sequence space again, so the ceiling must scale.

        Without that the pane reported more sequences to test than were possible.
        """
        for blocks in ([1], [1, 2], [1, 2, 3]):
            app.block_selector.value = blocks
            est = app.app_state.get_resource_estimate()
            assert est["observed_sequences"] <= est["total_sequences"], (blocks, est)

    def test_the_pane_is_populated_by_the_load_itself(self, loaded):
        """Not left blank until the user happens to touch a parameter widget."""
        assert app.resource_estimate_pane.visible
        assert "Sequences to test" in str(app.resource_estimate_pane.object)


# ---------------------------------------------------------------------------
# The criterion
# ---------------------------------------------------------------------------

class TestCriterion:
    def test_a_small_performance_criterion_does_not_raise(self, loaded):
        """The widget's minimum drops to 1 for a performance criterion.

        `app_state.criterion` kept a lower bound of 10, so entering the small count
        that a performance criterion is for raised ValueError inside the widget's
        watcher and took the app down.
        """
        app.criterion_order_widget.value = "4 rewarded in a row"
        assert app.criterion_widget.start == 1
        app.criterion_widget.value = 3
        assert app.app_state.criterion == 3

    def test_trial_order_restores_the_trial_floor(self, loaded):
        app.criterion_order_widget.value = "4 rewarded in a row"
        app.criterion_order_widget.value = "Trials (standard)"
        assert app.criterion_widget.start == 10
        assert app.criterion_widget.step == 10

    def test_shortfall_pane_reports_pairs_and_subjects(self, loaded):
        app.criterion_order_widget.value = "4 rewarded in a row"
        app.criterion_widget.value = 3
        assert app.criterion_shortfall_pane.visible
        text = str(app.criterion_shortfall_pane.object)
        assert "subject-contingency pairs" in text or "reach this criterion" in text


# ---------------------------------------------------------------------------
# Results and export
# ---------------------------------------------------------------------------

class TestResults:
    @pytest.fixture
    def result(self, loaded):
        app.seq_len_max_widget.value = 3
        app.resample_widget.value = 200
        app.block_selector.value = [1, 2, 3]
        app.subject_filter.value = list(app.subject_filter.options)
        app.app_state.run_analysis()
        return app.app_state.result

    def test_keys_are_block_sequence_pairs(self, result):
        block, symbols = app._split_sequence(result.sequences[0])
        assert block in (1, 2, 3)
        assert isinstance(symbols, tuple)

    def test_reported_lengths_span_the_real_range(self, result):
        lengths = {app._seq_len(s) for s in result.sequences}
        assert lengths == {1, 2, 3}, lengths

    def test_exported_csv_length_matches_its_own_label(self, result):
        section = app.make_download_section(result)
        downloads = list(section.select(app.pn.widgets.FileDownload))
        assert downloads, "expected a download button"
        frame = pd.read_csv(downloads[0].callback())
        assert sorted(frame["length"].unique()) == [1, 2, 3]
        symbols = frame["sequence"].str.split(": ").str[-1].str.split("-")
        assert (symbols.str.len() == frame["length"]).all()

    def test_result_panes_build(self, result):
        assert app.make_sig_table(result) is not None
        assert app.build_results_tabs() is not None


class TestManhattanPlot:
    """The figure in README.md, so its faults are the first thing anyone sees."""

    @pytest.fixture
    def rendered(self, loaded):
        import holoviews as hv
        app.seq_len_max_widget.value = 3
        app.resample_widget.value = 200
        app.block_selector.value = [1, 2, 3]
        app.app_state.run_analysis()
        column = app.make_manhattan_plot(app.app_state.result)
        return hv.render(column[1].object), app.app_state.result

    def test_the_wheel_is_not_bound_to_zoom(self, rendered):
        """Otherwise scrolling the page over the plot rescales its axes instead.

        This is not a cosmetic preference. Bokeh takes the scroll slot by default, so
        a reader scrolling past the plot silently zooms it, and cannot tell the
        wrecked view from the real one. It is how the committed overview screenshot
        came to show a sequence axis running to ten thousand and a negative
        -log10(g-value), neither of which this code can produce.
        """
        figure, _ = rendered
        assert figure.toolbar.active_scroll is None
        # still offered deliberately, with a reset to undo it
        names = [type(t).__name__ for t in figure.toolbar.tools]
        assert "WheelZoomTool" in names
        assert "ResetTool" in names

    def test_the_axes_are_fitted_to_the_data(self, rendered):
        """Auto-ranging spent most of the plot area on regions with no points in them."""
        figure, result = rendered
        n_seq = len(result.sequences)

        # a g-value is at most 1, so -log10(g) is never negative
        assert figure.y_range.start > -0.1
        assert figure.y_range.start <= 0, "room for the many sequences sitting at 0"

        # the sequence axis stops near the last sequence rather than decades past it
        assert figure.x_range.end < n_seq * 1.5, (figure.x_range.end, n_seq)
        assert figure.x_range.end >= n_seq
        assert 0 < figure.x_range.start <= 1, "the first sequence must be visible"


# ---------------------------------------------------------------------------
# Failure reporting
# ---------------------------------------------------------------------------

def test_a_failed_run_is_reported_even_when_the_callback_is_deferred(loaded, monkeypatch):
    """`except ... as e` unbinds `e` at the end of the block.

    On a served app the completion callback is deferred onto the event loop, by which
    point a lambda closing over `e` raises NameError instead of reporting anything:
    the run appeared to hang with the button stuck on "Running...". Deferral is
    simulated here by collecting the callbacks and invoking them after the except
    block has returned, which is what the real scheduler does.
    """
    deferred = []
    monkeypatch.setattr(app.pn.state, "execute",
                        lambda callback, *a, **kw: deferred.append(callback))

    def explode():
        raise RuntimeError("synthetic failure")

    monkeypatch.setattr(app.app_state, "run_analysis", explode)

    class Inline:
        """Run the thread body synchronously so the test can observe the outcome."""

        def __init__(self, target=None, daemon=None):
            self._target = target

        def start(self):
            self._target()

    monkeypatch.setattr(app.threading, "Thread", Inline)
    monkeypatch.setattr(app, "BROWSER_MODE", False)

    app.on_run_click(None)

    assert deferred, "the failure path scheduled no callback"
    for callback in deferred:
        callback()  # must not raise NameError

    assert "synthetic failure" in str(app.run_status.object)
