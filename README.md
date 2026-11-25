# Yeast Structure Benchmark

An educational structural bioinformatics project comparing a published AlphaFold prediction of yeast phosphoglycerate kinase (PGK1) with the experimentally determined structure 3PGK.

[![Reproducibility](https://github.com/rodrigo790/yeast-structure-benchmark/actions/workflows/reproducibility.yml/badge.svg)](https://github.com/rodrigo790/yeast-structure-benchmark/actions/workflows/reproducibility.yml)

## Research question

How closely does the predicted PGK1 structure match the experimental structure, and how do AlphaFold's confidence scores compare with the observed structural agreement?

## Data

- Reference sequence: [UniProt P00560](https://www.uniprot.org/uniprotkb/P00560/entry), *Saccharomyces cerevisiae*.
- Predicted model: [AlphaFold AF-P00560-F1](https://alphafold.ebi.ac.uk/entry/P00560), file version 6, chain A.
- Experimental structure: [RCSB PDB 3PGK](https://www.rcsb.org/structure/3PGK), chain A, X-ray diffraction at 2.5 Å resolution.
- Predicted aligned error (PAE): the AlphaFold version 6 PAE matrix.

Files and their SHA-256 checksums are recorded in [`data/manifest.json`](data/manifest.json). The project uses the saved snapshot when reproducing the analysis.

## Methods

1. Align each structure's sequence to the UniProt reference using Biopython `PairwiseAligner` and map matching residues with available Cα atoms.
2. Superpose the predicted coordinates on the experimental structure with a least-squares rigid-body fit. Calculate Cα RMSD and per-residue distances.
3. Compare the N-terminal and C-terminal cores after a global fit and after independent fits of each core. Domain boundaries are based on the mature protein numbering described by [Sherman et al. (1997)](https://doi.org/10.1002/pro.5560060415), shifted by one residue to account for the initiator methionine in UniProt. The interdomain helix and C-terminal tail are reported separately from the two cores.
4. Compare AlphaFold per-residue pLDDT with a Cα-only local distance agreement score computed against the matched 3PGK residues. The local score considers reference Cα neighbors within 15 Å and distance thresholds of 0.5, 1, 2 and 4 Å. It is an implementation for this comparison, not a complete all-atom lDDT validation.
5. Visualize the PAE matrix and the globally aligned matched Cα coordinates.

## Results from the saved snapshot

| Measure | Result |
| --- | ---: |
| UniProt reference length | 416 residues |
| Matched Cα pairs | 408 |
| Global Cα RMSD | 3.242 Å |
| N-terminal core RMSD, global fit / independent fit | 3.212 / 2.855 Å |
| C-terminal core RMSD, global fit / independent fit | 3.410 / 3.170 Å |
| Mean pLDDT for matched residues | 96.75 / 100 |
| Mean per-residue Cα local agreement against 3PGK | 64.79 / 100 |

S254 has the largest Cα displacement after the global fit at approximately 10.52 Å. Its predicted pLDDT is 98.56 and its calculated local Cα agreement with 3PGK is 40.54. These measurements describe this particular comparison and do not by themselves identify the cause of the difference.

PAE is a model confidence estimate and is distinct from the observed displacement after fitting to 3PGK. The core-to-core PAE means in the two matrix directions are approximately 5.28 and 5.48 Å.

## Reproduce the analysis

On Windows, start Docker Desktop and double-click `Run benchmark.bat`. The initial build requires internet access to obtain the Python image and locked dependencies. Subsequent analysis uses the included source snapshot. Results are written to `reproduced-results/`.

Alternatively, run `docker compose run --rm --build benchmark` in this folder. To verify the scientific checks without network access after building, run `docker run --rm --network none yeast-structure-benchmark python -m unittest discover -s tests -v` if you built the image with `docker build -t yeast-structure-benchmark .`.

Requires Python 3.12. In PowerShell from this folder:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.lock
.\.venv\Scripts\python.exe -m benchmark.snapshot
.\.venv\Scripts\python.exe -m benchmark.structure
.\.venv\Scripts\python.exe -m benchmark.superposition
.\.venv\Scripts\python.exe -m benchmark.domains
.\.venv\Scripts\python.exe -m benchmark.residue_plot
.\.venv\Scripts\python.exe -m benchmark.confidence
.\.venv\Scripts\python.exe -m benchmark.pae_plot
.\.venv\Scripts\python.exe -m benchmark.view_3d
```

Outputs are saved to `results/`. Open `results/pgk1_superposition_3d.html` in a browser to inspect the Cα view.

## Interpretation

This is a comparison of one published prediction with one experimental structure, not an independent benchmark of AlphaFold's general performance. Different structures, sequence differences, conformations, crystallographic conditions and modeling choices can affect the comparison. pLDDT and PAE are confidence estimates, and the local Cα agreement score is calculated against the selected experimental structure.

Reference paper: [Watson et al. (1982), Sequence and structure of yeast phosphoglycerate kinase](https://doi.org/10.1002/j.1460-2075.1982.tb01366.x).
