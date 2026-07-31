import pandas as pd
import hvplot.pandas  # noqa: F401
import holoviews as hv

hv.extension("bokeh")


def ranked_gvalues_html(csv_path, output_path, alpha=0.5, title=None):
    """Create interactive ranked g-values plot from a sequence_pvalues CSV.

    Parameters
    ----------
    csv_path : path-like
        Path to sequence_pvalues.csv (columns: rank, sequence_encoded,
        sequence_decoded, length, g_value, significant, direction).
    output_path : path-like
        Where to write the .html file.
    alpha : float
        Significance threshold line.
    title : str or None
        Plot title.
    """
    df = pd.read_csv(csv_path)

    if title is None:
        n_sig = int(df["significant"].sum())
        title = f"Ranked g-values — {n_sig}/{len(df)} significant (α={alpha})"

    df["length_str"] = "len=" + df["length"].astype(str)

    scatter = df.hvplot.scatter(
        x="rank",
        y="g_value",
        by="length_str",
        hover_cols=["sequence_decoded", "direction", "g_value", "length", "rank"],
        title=title,
        xlabel="Rank",
        ylabel="g-value (ζ)",
        width=1000,
        height=500,
        size=20,
        alpha=0.6,
    )

    threshold = hv.HLine(alpha).opts(
        color="red", line_dash="dashed", line_width=1.5
    )

    plot = scatter * threshold

    hvplot.save(plot, str(output_path))
    return str(output_path)
