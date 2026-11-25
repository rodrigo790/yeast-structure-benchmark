"""Compare the PGK1 domain cores after global and independent C-alpha fits."""

import argparse
import csv
from pathlib import Path

from .snapshot import ROOT, verify_snapshot
from .structure import C_CORE, N_CORE, compare, deviations, region_map, rmsd

DEFAULT_OUTPUT = ROOT / "results" / "pgk1_domain_rmsd.csv"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Destination for the domain table")
    args = parser.parse_args()

    verify_snapshot()
    comparison = compare()
    masks = region_map(comparison.positions)
    globally_fitted_distances = deviations(comparison.global_fit, comparison.experimental)

    rows = []
    for name, independent_fit in ((N_CORE, comparison.n_core_fit), (C_CORE, comparison.c_core_fit)):
        mask = masks[name]
        global_fit_rmsd = rmsd(globally_fitted_distances[mask])
        independent_fit_rmsd = rmsd(deviations(independent_fit[mask], comparison.experimental[mask]))
        rows.append((name, int(mask.sum()), global_fit_rmsd, independent_fit_rmsd))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(("region", "matched_ca_pairs", "global_fit_rmsd_angstrom", "independent_fit_rmsd_angstrom"))
        for name, pairs, global_value, independent_value in rows:
            writer.writerow((name, pairs, f"{global_value:.6f}", f"{independent_value:.6f}"))

    print("RMSD in angstroms for the same C-alpha pairs under both fits")
    print(f"{'Region':<17} {'Pairs':>5} {'Global fit':>11} {'Own fit':>10}")
    for name, pairs, global_value, independent_value in rows:
        print(f"{name:<17} {pairs:>5} {global_value:>11.3f} {independent_value:>10.3f}")
    print(f"Saved domain RMSD: {args.output}")


if __name__ == "__main__":
    main()
