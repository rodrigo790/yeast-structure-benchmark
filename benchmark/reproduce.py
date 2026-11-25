"""Reproduce all outputs from the frozen input snapshot."""

import os
import subprocess
import sys

from .snapshot import ROOT


def main() -> None:
    environment = os.environ.copy()
    environment["MPLBACKEND"] = "Agg"
    for module in ("snapshot", "structure", "superposition", "domains", "residue_plot", "confidence", "pae_plot", "view_3d"):
        print(f"\nRunning benchmark.{module}", flush=True)
        subprocess.run([sys.executable, "-m", f"benchmark.{module}"], cwd=ROOT, env=environment, check=True)
    print("\nReproduction completed successfully.")


if __name__ == "__main__":
    main()
