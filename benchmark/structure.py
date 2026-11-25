"""Map C-alpha observations to the UniProt sequence and compare their geometry."""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1

from .snapshot import RAW, verify_snapshot

# Sherman et al. (1997), doi:10.1002/pro.5560060415, numbers from mature Ser1.
# UniProt numbering includes Met1. The cores exclude the linker and the tail
# that folds back against the N-terminal domain.
REGIONS = (
    ("Initiator Met", 1, 1),
    ("N-terminal core", 2, 185),
    ("Interdomain helix", 186, 200),
    ("C-terminal core", 201, 401),
    ("C-terminal tail", 402, 416),
)
N_CORE, C_CORE = "N-terminal core", "C-terminal core"


@dataclass(frozen=True)
class Residue:
    amino_acid: str
    number: str
    ca: np.ndarray | None
    b_factor: float | None


def read_chain(path: Path, chain_id: str) -> list[Residue]:
    structure = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True).get_structure(path.stem, path)
    chain = structure[0][chain_id]
    residues = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        amino_acid = seq1(residue.resname)
        atom = residue["CA"] if "CA" in residue else None
        number = f"{residue.id[1]}{residue.id[2].strip()}"
        residues.append(Residue(amino_acid, number,
                                np.asarray(atom.coord, dtype=float) if atom is not None else None,
                                float(atom.bfactor) if atom is not None else None))
    if not residues:
        raise ValueError(f"No standard amino-acid residues found in chain {chain_id} of {path.name}")
    return residues


def map_to_reference(canonical: str, residues: list[Residue]) -> dict[int, Residue]:
    """Key matching residues by 1-based UniProt position, independently of PDB numbering."""
    observed = "".join(residue.amino_acid for residue in residues)
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2
    aligner.mismatch_score = -2
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(canonical, observed)[0]
    return {
        int(i + 1): residues[j] for i, j in alignment.indices.T
        if i >= 0 and j >= 0 and canonical[i] == observed[j] and canonical[i] != "X"
    }


def superpose(moving: np.ndarray, fixed: np.ndarray, fit_mask: np.ndarray | slice = slice(None)) -> np.ndarray:
    """Apply a least-squares rigid-body fit with a proper rotation."""
    source, target = moving[fit_mask], fixed[fit_mask]
    if len(source) < 3:
        raise ValueError("A rigid-body fit requires at least three matched Cα coordinates")
    origin, destination = source.mean(axis=0), target.mean(axis=0)
    u, _, vh = np.linalg.svd((source - origin).T @ (target - destination))
    if np.linalg.det(u @ vh) < 0:
        u[:, -1] *= -1
    return (moving - origin) @ (u @ vh) + destination


def deviations(fitted: np.ndarray, fixed: np.ndarray) -> np.ndarray:
    return np.linalg.norm(fitted - fixed, axis=1)


def rmsd(distances: np.ndarray) -> float:
    return float(np.sqrt(np.mean(distances**2)))


def local_agreement(moving: np.ndarray, fixed: np.ndarray) -> np.ndarray:
    """Return per-residue lDDT-Cα on a 0–100 scale for matched Cα pairs."""
    target = np.linalg.norm(fixed[:, None] - fixed[None, :], axis=2)
    source = np.linalg.norm(moving[:, None] - moving[None, :], axis=2)
    neighbors = (target < 15.0) & ~np.eye(len(fixed), dtype=bool)
    errors = np.abs(target - source)
    retained = ((errors[..., None] < np.array([0.5, 1.0, 2.0, 4.0])).mean(axis=2) * neighbors).sum(axis=1)
    counts = neighbors.sum(axis=1)
    return np.divide(100 * retained, counts, out=np.full(len(fixed), np.nan), where=counts != 0)


def region_map(positions: np.ndarray) -> dict[str, np.ndarray]:
    return {name: (positions >= first) & (positions <= last) for name, first, last in REGIONS}


@dataclass(frozen=True)
class Comparison:
    positions: np.ndarray
    canonical: str
    experimental_numbers: list[str]
    predicted: np.ndarray
    experimental: np.ndarray
    global_fit: np.ndarray
    n_core_fit: np.ndarray
    c_core_fit: np.ndarray
    plddt: np.ndarray
    lddt: np.ndarray
    pae: np.ndarray


def map_pgk1() -> tuple[str, dict[int, Residue], dict[int, Residue]]:
    canonical = str(SeqIO.read(RAW / "P00560.fasta", "fasta").seq)
    if len(canonical) != 416 or not canonical.startswith("MS"):
        raise ValueError("The snapshot does not match the expected full-length PGK1 reference")
    model = map_to_reference(canonical, read_chain(RAW / "AF-P00560-F1-model_v6.cif", "A"))
    crystal = map_to_reference(canonical, read_chain(RAW / "3PGK.cif", "A"))
    return canonical, model, crystal


def matched_ca_positions(model: dict[int, Residue], crystal: dict[int, Residue]) -> np.ndarray:
    return np.array(sorted(i for i in model.keys() & crystal.keys()
                           if model[i].ca is not None and crystal[i].ca is not None), dtype=int)


def compare() -> Comparison:
    canonical, model, crystal = map_pgk1()
    positions = matched_ca_positions(model, crystal)
    if len(positions) < 0.95 * len(canonical):
        raise ValueError("Fewer than 95% of canonical residues have matched Cα coordinates")
    if len(model) != len(canonical) or any(model[i].number != str(i) for i in model):
        raise ValueError("AlphaFold residue numbering must match UniProt positions for the PAE plot")
    predicted = np.array([model[i].ca for i in positions])
    experimental = np.array([crystal[i].ca for i in positions])
    confidence = np.array([model[i].b_factor for i in positions], dtype=float)
    pae_file = RAW / "AF-P00560-F1-predicted_aligned_error_v6.json"
    pae = np.array(json.loads(pae_file.read_text(encoding="utf-8"))[0]["predicted_aligned_error"], dtype=float)
    if pae.shape != (len(canonical), len(canonical)) or not np.isfinite(pae).all() or (pae < 0).any():
        raise ValueError("PAE matrix is invalid for the canonical sequence")
    if not (np.isfinite(confidence) & (confidence >= 0) & (confidence <= 100)).all():
        raise ValueError("Invalid model pLDDT values")
    masks = region_map(positions)
    return Comparison(
        positions=positions, canonical=canonical, experimental_numbers=[crystal[i].number for i in positions],
        predicted=predicted, experimental=experimental, global_fit=superpose(predicted, experimental),
        n_core_fit=superpose(predicted, experimental, masks[N_CORE]),
        c_core_fit=superpose(predicted, experimental, masks[C_CORE]),
        plddt=confidence, lddt=local_agreement(predicted, experimental), pae=pae,
    )


if __name__ == "__main__":
    verify_snapshot()
    canonical, model, crystal = map_pgk1()
    positions = matched_ca_positions(model, crystal)
    print(f"UniProt reference: {len(canonical)} residues")
    print(f"AlphaFold chain A: {len(model)} residues mapped to UniProt")
    print(f"3PGK chain A: {len(crystal)} residues mapped to UniProt")
    print(f"Comparable C-alpha pairs: {len(positions)}")
    print(f"\n{'UniProt':>8}  {'AA':>2}  {'AlphaFold':>9}  {'3PGK':>6}")
    for i in positions[:6]:
        print(f"{i:>8}  {canonical[i - 1]:>2}  {model[i].number:>9}  {crystal[i].number:>6}")
