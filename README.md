# PoseBusters-Benchmark
This repo contains the structure prediction benchmarking dataset from "PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences". (https://doi.org/10.48550/arXiv.2308.05777)

## The Dataset

### FASTAs
The full set of FASTA files for the 308-structure benchmark set is under `sequences/all/`. To test on a set that only contains single-chain structures, a 134-structure subset is copied into `sequences/monomer/`. Finally, a merged FASTA file of all monomers is found in `sequences/monomer_combined/`

In each directory, there is a `no_ligand/` and a `ligand` subdirectory. The `no_ligand/` directory contains just the protein sequences as is typical. To also store ligand information, the `ligand/` subdirectory contains FASTA-like files which contain the ligand SMILES as an additional chain. This is the input format for Chai-1.

### Structures
Similarly, PDB files for the 308-structure benchmark are stored in `structures/all/`, with the 134-structure monomer subset under `structures/monomer/`. These files have been stripped of all atoms except for the protein and ligand specified in the original PoseBusters list.

The raw structure files from the protein databank are stored in `structures/raw_from_pdb/`.

### ID List
The original list of PDB/ligand IDs is stored in `posebuster_benchmark.txt`.
Note: 7D6O was replaced with 8J79, as it was a duplicate and marked as obsolete.

## Reproducing

To reproduce this dataset, run `python scripts/generate_posebusters_benchmark.py` in and environment with `rdkit` and `biopython` installed.
