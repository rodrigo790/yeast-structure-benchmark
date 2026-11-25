"""Superpose AlphaFold PGK1 onto 3PGK using matched C-alpha coordinates."""

import argparse
import csv
from pathlib import Path

from .snapshot import ROOT, verify_snapshot
from .structure import compare, deviations, rmsd

DEFAULT_OUTPUT = ROOT / "results" / "pgk1_superposition.csv"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Destination for the residue table")
    args = parser.parse_args()

    verify_snapshot()
    comparison = compare()
    distances = deviations(comparison.global_fit, comparison.experimental)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow((
            "uniprot_position", "amino_acid", "pdb_residue",
            "alphafold_x_aligned_angstrom", "alphafold_y_aligned_angstrom", "alphafold_z_aligned_angstrom",
            "experimental_x_angstrom", "experimental_y_angstrom", "experimental_z_angstrom",
            "ca_distance_angstrom",
        ))
        for position, pdb_residue, aligned, experimental, distance in zip(
            comparison.positions, comparison.experimental_numbers, comparison.global_fit,
            comparison.experimental, distances, strict=True,
        ):
            writer.writerow((position, comparison.canonical[position - 1], pdb_residue,
                             *(f"{value:.6f}" for value in (*aligned, *experimental, distance))))

    print(f"Matched C-alpha pairs: {len(distances)}")
    print(f"Global C-alpha RMSD: {rmsd(distances):.3f} angstroms")
    print(f"Saved per-residue coordinates and distances: {args.output}")


if __name__ == "__main__":
    main()
