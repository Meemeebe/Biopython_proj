
import matplotlib.pyplot as plt


def create_gc_plot(sequence, window_size=100):
    """Create a sliding-window GC content plot.

    Returns a Matplotlib figure, or None if no valid windows exist.
    """
    sequence = str(sequence).upper()

    if window_size < 1:
        raise ValueError("Window size must be at least 1.")

    if len(sequence) < window_size:
        return None

    positions = []
    gc_values = []

    for start in range(len(sequence) - window_size + 1):
        window = sequence[start:start + window_size]

        # Calculate GC only for windows containing A, T, G and C.
        if not all(base in "ATGC" for base in window):
            continue

        gc = (
            (window.count("G") + window.count("C"))
            / window_size
        ) * 100

        positions.append(start + 1)
        gc_values.append(gc)

    if not gc_values:
        return None

    fig, ax = plt.subplots(figsize=(10, 4))

    ax.plot(positions, gc_values)
    ax.set_xlabel("Window start position (bp)")
    ax.set_ylabel("GC content (%)")
    ax.set_title("Sliding-Window GC Content")
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    return fig
