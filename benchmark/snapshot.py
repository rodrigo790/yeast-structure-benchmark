"""Download the PGK1 source files once and verify the snapshot for offline runs."""

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RAW = DATA / "raw"
MANIFEST = DATA / "manifest.json"
SOURCES = {
    "P00560.fasta": "https://rest.uniprot.org/uniprotkb/P00560.fasta",
    "AF-P00560-F1-model_v6.cif": "https://alphafold.ebi.ac.uk/files/AF-P00560-F1-model_v6.cif",
    "AF-P00560-F1-predicted_aligned_error_v6.json": "https://alphafold.ebi.ac.uk/files/AF-P00560-F1-predicted_aligned_error_v6.json",
    "3PGK.cif": "https://files.rcsb.org/download/3PGK.cif",
}


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def create_snapshot() -> None:
    if MANIFEST.exists():
        verify_snapshot()
        return
    downloads = {}
    for name, url in SOURCES.items():
        with urlopen(url, timeout=35) as response:
            downloads[name] = response.read()
        expected_start = {".fasta": b">", ".cif": b"data_", ".json": b"["}[Path(name).suffix]
        if not downloads[name].lstrip().startswith(expected_start):
            raise ValueError(f"Wrong format returned for {name}")

    RAW.mkdir(parents=True, exist_ok=True)
    for name, raw in downloads.items():
        (RAW / name).write_bytes(raw)
    metadata = {
        "protein": "PGK1",
        "uniprot": "P00560",
        "organism": "Saccharomyces cerevisiae",
        "predicted_structure": {"entry_id": "AF-P00560-F1", "file_version": 6, "chain": "A"},
        "experiment": {"pdb_id": "3PGK", "chain": "A", "method": "X-ray diffraction", "resolution_angstrom": 2.5},
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "files": {name: {"url": SOURCES[name], "sha256": digest(raw)} for name, raw in downloads.items()},
    }
    MANIFEST.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8", newline="\n")


def verify_snapshot() -> dict:
    metadata = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if set(metadata["files"]) != set(SOURCES):
        raise ValueError("Snapshot does not list the expected source files")
    for name, info in metadata["files"].items():
        if digest((RAW / name).read_bytes()) != info["sha256"]:
            raise ValueError(f"Checksum mismatch for {name}")
    return metadata


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--download", action="store_true", help="Create the snapshot if it does not exist")
    if cli.parse_args().download:
        create_snapshot()
    print(f"Verified {len(verify_snapshot()['files'])} source files for PGK1.")
