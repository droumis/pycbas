import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def direction_counts(n_positive, n_negative, pos_label, neg_label, ax=None):
    """Bar chart of significant sequences by direction.

    Parameters
    ----------
    n_positive : int
        Count of significant sequences in the positive direction.
    n_negative : int
        Count of significant sequences in the negative direction.
    pos_label : str
        Label for the positive direction bar.
    neg_label : str
        Label for the negative direction bar.
    ax : matplotlib Axes or None
        If None, creates a new figure.

    Returns
    -------
    fig, ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 3))
    else:
        fig = ax.get_figure()

    bars = ax.bar([pos_label, neg_label], [n_positive, n_negative],
                  color=["#0066cc", "#cc6600"])
    ax.set_ylabel("# significant sequences")
    ax.set_title("Significant Sequences by Direction")

    for bar, val in zip(bars, [n_positive, n_negative]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
                str(val), ha="center", fontsize=10)

    fig.tight_layout()
    return fig, ax


def paper_comparison(our_fig_path, paper_fig_path, our_title, paper_title, ax=None):
    """Side-by-side image comparison of our figure vs a published paper figure.

    Parameters
    ----------
    our_fig_path : str or Path
        Path to our generated figure image.
    paper_fig_path : str or Path
        Path to the paper's figure image.
    our_title : str
        Title for our figure panel.
    paper_title : str
        Title for the paper figure panel.
    ax : matplotlib Axes or None
        If None, creates a new figure with two subplots. If provided, it is
        replaced by a 1x2 subplot grid on the same figure.

    Returns
    -------
    fig, ax : fig and array of two axes
    """
    if ax is None:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    else:
        fig = ax.get_figure()
        ax.remove()
        axes = fig.subplots(1, 2)

    our_img = mpimg.imread(str(our_fig_path))
    paper_img = mpimg.imread(str(paper_fig_path))

    axes[0].imshow(our_img)
    axes[0].set_title(our_title)
    axes[0].axis("off")

    axes[1].imshow(paper_img)
    axes[1].set_title(paper_title)
    axes[1].axis("off")

    fig.tight_layout()
    return fig, axes
