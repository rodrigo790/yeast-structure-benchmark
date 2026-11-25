"""Compare AlphaFold local confidence with agreement against experimental 3PGK."""

import argparse
import csv
from pathlib import Path

import matplotlib
import numpy as np

from .snapshot import ROOT, verify_snapshot
from .structure import Comparison, compare, deviations

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PREDICTED_COLOR = "#126276"
OBSERVED_COLOR = "#9a4b20"


def save_plot(comparison: Comparison, distances: np.ndarray, peak: int, destination: Path) -> None:
    positions = np.arange(1, len(comparison.canonical) + 1)
    plddt, lddt, displacement = np.full((3, len(positions)), np.nan)
    indices = comparison.positions - 1
    plddt[indices], lddt[indices], displacement[indices] = comparison.plddt, comparison.lddt, distances
    focus = int(comparison.positions[peak])

    figure, (confidence_axis, distance_axis) = plt.subplots(
        2, 1, figsize=(10.2, 6.3), sharex=True, height_ratios=(1.35, 1), layout="constrained",
    )
    figure.suptitle("PGK1: AlphaFold confidence and agreement with 3PGK", x=0.08, ha="left", fontsize=15)
    confidence_axis.plot(positions, plddt, color=PREDICTED_COLOR, label="AlphaFold pLDDT", linewidth=1.6)
    confidence_axis.plot(positions, lddt, color=OBSERVED_COLOR, label="lDDT-Cα vs. 3PGK", linewidth=1.35)
    confidence_axis.scatter([focus, focus], [comparison.plddt[peak], comparison.lddt[peak]],
                            c=[PREDICTED_COLOR, OBSERVED_COLOR], s=36, zorder=4)
    confidence_axis.set(ylabel="Local score (0–100)", ylim=(0, 105))
    confidence_axis.legend(loc="lower left", frameon=False, ncols=2)

    distance_axis.plot(positions, displacement, color="#234a60", linewidth=1.1)
    distance_axis.scatter(focus, distances[peak], color="#234a60", s=30, zorder=4)
    distance_axis.annotate(f"{comparison.canonical[focus - 1]}{focus} · {distances[peak]:.2f} Å",
                           (focus, distances[peak]), xytext=(12, 6), textcoords="offset points", color="#234a60")
    distance_axis.set(xlabel="UniProt P00560 residue position", ylabel="Distance after global fit (Å)",
                      xlim=(0.5, len(positions) + 0.5), ylim=(0, distances[peak] * 1.22))
    for axis in (confidence_axis, distance_axis):
        axis.axvline(focus, linestyle=":", color="#6a7880", linewidth=1)
        axis.grid(axis="y", color="#526975", alpha=0.16, linewidth=0.7)
        axis.spines[["right", "top"]].set_visible(False)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results", help="Directory for the PNG and CSV")
    args = parser.parse_args()

    verify_snapshot()
    comparison = compare()
    distances = deviations(comparison.global_fit, comparison.experimental)
    peak = int(np.argmax(distances))
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    with (out / "pgk1_confidence.csv").open("w", encoding="utf-8", newline="") as sink:
        writer = csv.writer(sink, lineterminator="\n")
        writer.writerow(("uniprot_position", "amino_acid", "alphafold_plddt",
                         "lddt_ca_vs_3pgk", "global_fit_ca_distance_angstrom"))
        for position, confidence, agreement, distance in zip(
            comparison.positions, comparison.plddt, comparison.lddt, distances, strict=True,
        ):
            writer.writerow((position, comparison.canonical[position - 1], f"{confidence:.2f}",
                             f"{agreement:.6f}" if np.isfinite(agreement) else "", f"{distance:.6f}"))
    destination = out / "pgk1_confidence.png"
    save_plot(comparison, distances, peak, destination)

    position = int(comparison.positions[peak])
    print(f"Residues compared: {len(comparison.positions)}")
    print(f"Largest distance after global fit: {comparison.canonical[position - 1]}{position}, {distances[peak]:.2f} angstroms")
    print(f"At that residue: pLDDT {comparison.plddt[peak]:.2f}/100 | lDDT-Ca vs. 3PGK {comparison.lddt[peak]:.2f}/100")
    print(f"Saved plot and CSV in: {out}")


if __name__ == "__main__":
    main()
