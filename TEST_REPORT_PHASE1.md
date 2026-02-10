# AutoScan Phase 1: Structural Blind-Docking Benchmark
**Date:** [Insert Date]
**Executor:** AutoScan Validation Suite

1. Scientific Objective
To quantify the structural accuracy of AutoScan's docking engine by performing a "Blind Redocking" experiment on 10 diverse protein-ligand systems.

2. The Dataset (N=10)
We selected 10 targets representing different pharmacological classes to ensure robustness:

| Class | PDB ID | Target Name | Ligand | Challenge |
|---|---|---|---|---|
| Antibiotic | 2XCT | S. aureus Gyrase | CPF | Baseline validation |
| Antiviral | 1HXB | HIV Protease | ROC | Large peptidomimetic |
| Cancer | 1IEP | Abl Kinase | STI | Induced-fit (Gleevec) |
| Hormone | 3ERT | Estrogen Receptor | OHT | Hydrophobic pocket |
| Cancer | 1M17 | EGFR Kinase | AQ4 | ATP-competitive |
| Steroid | 1SQN | Progesterone Rec. | PRG | Rigid core |
| Bio-tool | 1OQA | Streptavidin | STR | Ultra-high affinity |
| Cancer | 4DJU | B-Raf Kinase | 032 | V600E Mutant |
| Pain | 3LN1 | COX-2 | CEL | NSAID (Celecoxib) |
| Flu | 2HU4 | Neuraminidase | G39 | Polar/Charged |

3. Protocol
Retrieval: Download experimentally solved crystal structures (PDB).

Sanitization: Split protein and ligand; clean waters/ions.

Blind Prep: Convert ligand to PDBQT. Critical: The starting conformation is randomized/unbiased.

Docking: Run AutoScan (Vina) with a dynamically sized grid box (Ligand Size + 10A Buffer).

Analysis: Compute RMSD (Root Mean Square Deviation) between the Docked Pose and the Crystal Reference.

4. Success Metrics
Individual Pass: RMSD < 2.5 A AND Affinity < -7.0 kcal/mol.

Overall Success: > 80% of targets pass.

## Phase 1 Results
| Class | PDB ID | Ligand | Docked Energy (kcal/mol) | RMSD (A) | Status |
|---|---|---|---:|---:|---|
| Antibiotic | 2XCT | CPF | -11.080 | 0.001 | PASS |
| Antiviral | 1HXB | ROC | NA | NA | ERROR (Vina docking failed: Atom type Nd is not a valid AutoDock type) |
| Cancer | 1IEP | STI | -18.030 | 0.001 | PASS |
| Hormone | 3ERT | OHT | -14.100 | 0.001 | PASS |
| Cancer | 1M17 | AQ4 | -10.680 | 0.001 | PASS |
| Steroid | 1SQN | PRG | NA | NA | ERROR (Ligand PRG not found in workspace\phase1\1SQN\1SQN.pdb) |
| Bio-tool | 1OQA | STR | NA | NA | ERROR (Ligand STR not found in workspace\phase1\1OQA\1OQA.pdb) |
| Cancer | 4DJU | 032 | NA | NA | ERROR (Ligand 032 not found in workspace\phase1\4DJU\4DJU.pdb) |
| Pain | 3LN1 | CEL | -15.250 | 0.001 | PASS |
| Flu | 2HU4 | G39 | -9.224 | 0.001 | PASS |

Mean RMSD: 0.001 A
Success rate: 60.0%

## Round 2 Results (Phase 1.1 Fixed)
| Class | PDB ID | Ligand | Docked Energy (kcal/mol) | RMSD (A) | Status |
|---|---|---|---:|---:|---|
| Antibiotic | 2XCT | CPF | NA | NA | ERROR (Could not parse binding affinity from Vina output) |
| Antiviral | 1HSV | MK1 | NA | NA | ERROR (HTTP Error 404: Not Found) |
| Cancer | 1IEP | STI | NA | NA | ERROR (Could not parse binding affinity from Vina output) |
| Hormone | 3ERT | OHT | -14.080 | 0.001 | PASS |
| Cancer | 1M17 | AQ4 | -10.680 | 0.001 | PASS |
| Steroid | 1SQN | PRG | NA | NA | ERROR (Ligand PRG not found in workspace\phase1_fixed\1SQN\1SQN.pdb) |
| Bio-tool | 1OQA | BTN | NA | NA | ERROR (Ligand BTN not found in workspace\phase1_fixed\1OQA\1OQA.pdb) |
| Cancer | 4DJU | 032 | NA | NA | ERROR (Ligand 032 not found in workspace\phase1_fixed\4DJU\4DJU.pdb) |
| Pain | 3LN1 | CEL | NA | NA | ERROR (Could not parse binding affinity from Vina output) |
| Flu | 2HU4 | G39 | NA | NA | ERROR (Could not parse binding affinity from Vina output) |

Mean RMSD: 0.001 A
Success rate: 20.0%
