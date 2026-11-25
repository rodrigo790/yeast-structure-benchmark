"""Plot the C-alpha distances from the global PGK1 superposition."""

import argparse
import csv
from pathlib import Path

import matplotlib
import numpy as np

from .snapshot import ROOT
from .structure import REGIONS, rmsd
from .superposition import DEFAULT_OUTPUT as DISTANCES_CSV

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_OUTPUT = ROOT / "results" / "pgk1_residue_distances.png"
REGION_COLORS = ("#d9eaf3", "#f3e5c6", "#dceee4", "#eee2ef")


def load_distances(path: Path) -> tuple[np.ndarray, dict[int, str]]:
    distances = np.full(REGIONS[-1][2], np.nan)
    labels = {}
    with path.open(encoding="utf-8", newline="") as source:
        for row in csv.DictReader(source):
            position = int(row["uniprot_position"])
            distance = float(row["ca_distance_angstrom"])
            if not 1 <= position <= len(distances) or position in labels or not np.isfinite(distance) or distance < 0:
                raise ValueError(f"Invalid distance entry for UniProt position {position}")
            distances[position - 1] = distance
            labels[position] = row["amino_acid"]
    if not labels:
        raise ValueError("The residue table is empty")
    return distances, labels


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Destination for the PNG")
    args = parser.parse_args()

    if not DISTANCES_CSV.is_file():
        parser.error("Run `python -m benchmark.superposition` first to create the residue table.")
    distances, labels = load_distances(DISTANCES_CSV)
    overall = rmsd(distances[np.isfinite(distances)])
    maximum = int(np.nanargmax(distances)) + 1
    highest_displacement = float(distances[maximum - 1])

    figure, axis = plt.subplots(figsize=(10.4, 4.8), layout="constrained")
    for (name, first, last), color in zip(REGIONS[1:], REGION_COLORS, strict=True):
        axis.axvspan(first - 0.5, last + 0.5, color=color, linewidth=0, label=name)

    axis.plot(np.arange(1, len(distances) + 1), distances, color="#234a60", linewidth=1.2, zorder=3)
    axis.axhline(overall, color="#95552f", linestyle="--", linewidth=1.2,
                 label=f"Global RMSD ({overall:.3f} Å)")
    axis.scatter(maximum, highest_displacement, color="#234a60", s=24, zorder=4)
    label_offset = -34 if maximum > len(distances) / 2 else 30
    axis.annotate(f"{labels[maximum]}{maximum} · {highest_displacement:.2f} Å",
                  xy=(maximum, highest_displacement), xytext=(label_offset, 12), textcoords="offset points",
                  ha="right" if label_offset < 0 else "left", color="#234a60",
                  arrowprops={"arrowstyle": "-", "color": "#234a60"})

    axis.set(title="PGK1 Cα distances after global superposition",
             xlabel="UniProt P00560 residue position", ylabel="AlphaFold vs. 3PGK distance (Å)",
             xlim=(0.5, len(distances) + 0.5), ylim=(0, highest_displacement * 1.23))
    axis.grid(axis="y", color="#526975", alpha=0.16, linewidth=0.7)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncols=3, frameon=False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180)
    plt.close(figure)
    print(f"Plotted C-alpha pairs: {len(labels)}")
    print(f"Largest distance: {labels[maximum]}{maximum} at {highest_displacement:.3f} angstroms")
    print(f"Saved plot: {args.output}")


if __name__ == "__main__":
    main()
