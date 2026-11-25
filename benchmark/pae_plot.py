"""Plot the AlphaFold predicted aligned error (PAE) for PGK1."""

import argparse
from pathlib import Path

import matplotlib
import numpy as np

from .snapshot import ROOT, verify_snapshot
from .structure import C_CORE, N_CORE, compare, region_map

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "pgk1_pae.png")
    args = parser.parse_args()
    verify_snapshot()
    comparison = compare()
    matrix = comparison.pae
    masks = region_map(np.arange(1, len(comparison.canonical) + 1))
    n_positions = np.flatnonzero(masks[N_CORE])
    c_positions = np.flatnonzero(masks[C_CORE])
    nc = float(matrix[np.ix_(n_positions, c_positions)].mean())
    cn = float(matrix[np.ix_(c_positions, n_positions)].mean())

    figure, axis = plt.subplots(figsize=(7.2, 6.4), layout="constrained")
    image = axis.imshow(matrix, origin="lower", cmap="viridis_r", vmin=0,
                        extent=(0.5, len(matrix) + 0.5, 0.5, len(matrix) + 0.5))
    for boundary in (185.5, 200.5, 401.5):
        axis.axvline(boundary, color="white", alpha=0.7, linewidth=0.8, linestyle="--")
        axis.axhline(boundary, color="white", alpha=0.7, linewidth=0.8, linestyle="--")
    axis.set(title="AlphaFold PAE · PGK1", xlabel="UniProt residue", ylabel="UniProt residue")
    figure.colorbar(image, ax=axis, label="Predicted aligned error (Å)")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180)
    plt.close(figure)

    print(f"N-core to C-core mean PAE: {nc:.2f} angstroms")
    print(f"C-core to N-core mean PAE: {cn:.2f} angstroms")
    print(f"Saved PAE plot: {args.output}")


if __name__ == "__main__":
    main()
