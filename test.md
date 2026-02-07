# AutoScan Test Report

## Test 1: Structural Biology Validation (Redocking Benchmark)

### Scientific Question
Can AutoScan reproduce the crystallographic pose of ciprofloxacin in 2XCT?

### Protocol Summary
1. **Fetch** PDB 2XCT (S. aureus gyrase + CPF).
2. **Split** protein chains A/B and ligand CPF into separate PDBs.
3. **Shuffle** ligand pose with random rotation and translation.
4. **Dock** the shuffled ligand back into the pocket using Vina.
5. **Compare** docked ligand vs crystal ligand using heavy-atom RMSD (Kabsch).
6. **Pass/Fail** threshold: RMSD < 2.5 A.

### Implementation Notes
- Script: tests/validation/validate_structure.py
- Ligand name: CPF
- Pocket: GyrA_pocket (from config/pockets.yaml)
- RMSD computed with heavy-atom matching and Kabsch alignment.

### Results
- **RMSD**: 0.910 A
- **Status**: PASS

### Problems Encountered and Fixes

1. **OpenBabel missing on Windows**
   - **Symptom**: Pipeline failed with "obabel not found".
   - **Fix**: Created a conda environment and installed OpenBabel:
     - `conda create -n autoscan-env -c conda-forge python=3.11 openbabel`

2. **AutoDock Vina not available via conda on Windows**
   - **Symptom**: `autodock-vina` not found in conda channels.
   - **Fix**: Downloaded official Windows Vina binary and placed it in the conda env:
     - `vina_1.2.7_win.exe` -> `CONDA_PREFIX\Scripts\vina.exe`
     - Added `autoscan-vina.exe` wrapper by copying `vina.exe`.

3. **Ligand not found in 2XCT.pdb**
   - **Symptom**: Script raised "Ligand CPF not found".
   - **Cause**: The downloaded 2XCT.pdb lacked HETATM records in the trimmed file.
   - **Fix**: Prefer `data/receptors/2XCT_orig.pdb` when present.

4. **Atom count mismatch during RMSD**
   - **Symptom**: RMSD calculation failed (24 vs 28 atoms).
   - **Cause**: Hydrogen atoms present in docked PDBQT but not in the reference PDB.
   - **Fix**: RMSD now matches heavy atoms by name with per-atom indexing.

### Conclusion
The redocking benchmark passed with RMSD 0.910 A, validating AutoScan's structural accuracy on the 2XCT reference system.

## Test 2: Thermodynamic/Kinetic Validation (Specificity Benchmark)

### Scientific Question
Can AutoScan distinguish a true binder (ciprofloxacin) from a non-specific decoy (benzene)?

### Protocol Summary
1. **Fetch** PDB 2XCT (S. aureus gyrase + CPF).
2. **Split** protein chains A/B and ligand CPF into separate PDBs.
3. **Prepare** benzene decoy (RDKit, 3D embed + optimize).
4. **Dock** ciprofloxacin and benzene into the same pocket.
5. **Compare** affinity scores and compute Delta Delta G.

### Implementation Notes
- Script: tests/validation/validate_thermodynamics.py
- Ligand A: Ciprofloxacin (CPF from 2XCT)
- Ligand B: Benzene (RDKit SMILES: c1ccccc1)
- Pocket: GyrA_pocket (from config/pockets.yaml)
- PASS criteria:
   - Cipro affinity < -8.0 kcal/mol
   - Benzene affinity > -5.5 kcal/mol
   - Delta Delta G >= 2.5 kcal/mol (Cipro stronger)

### Results
- **Ciprofloxacin affinity**: -3.662 kcal/mol
- **Benzene affinity**: -1.702 kcal/mol
- **Delta Delta G**: 1.960 kcal/mol
- **Status**: FAIL

### Problems Encountered and Fixes
1. **Toolchain prerequisites**
    - **Symptom**: Docking requires OpenBabel and Vina binaries.
    - **Fix**: Added a preflight check in the script and documented Windows conda setup in README.

2. **Binding window mismatch**
   - **Symptom**: Ciprofloxacin did not reach the expected < -8.0 kcal/mol threshold.
   - **Hypothesis**: The small grid (6x10x15 A) may be too tight for CPF to find the native pose, leading to weak scores.
   - **Next step**: Enlarge the grid or use the ligand center extracted from 2XCT to guide docking.

### Conclusion
The specificity benchmark ran successfully but did not meet the affinity thresholds, so it currently fails the PASS criteria.
