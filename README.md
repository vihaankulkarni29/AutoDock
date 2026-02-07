# AutoScan: Molecular Docking for Binding Affinity Prediction

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Vision: The "Lock and Key" Philosophy

AutoScan simulates the fundamental biochemical principle of molecular recognition:

- **The Lock**: Receptor Protein (e.g., DNA Gyrase, GyrA)
- **The Key**: Ligand / Small Molecule (e.g., Ciprofloxacin, an antibiotic)
- **The Goal**: Calculate if the "Key" still fits the "Lock" after a mutation

This tool predicts **binding affinity (Delta G)** to detect how mutations affect drug efficacy.

## Features

- **Automated Virtual Screening**: Rapidly screen small molecule libraries against bacterial targets.
- **Physics-Based Scoring**: Uses Vina to calculate Delta G binding affinities.
- **Resistance Prediction (Beta)**: Currently detects steric clashes for large structural changes. Note: single-point mutation resistance prediction is under development (integration with Molecular Dynamics planned).

## Architectural Design

Built as a **standalone tool** with a **seamless API** for future integration into the **MutationScan** pipeline as **Module 7**.

### Directory Structure

```
AutoScan/
├── .github/workflows/          # CI/CD pipelines (linting, Docker builds)
├── config/                     # Pocket definitions (receptor binding sites)
│   └── pockets.yaml            # Grid box coordinates for Vina
├── data/                       # Local storage (gitignored)
│   ├── receptors/              # PDB files (Locks)
│   └── ligands/                # PDB/PDBQT files (Keys)
├── docker/                     # Docker build context
├── src/autoscan/               # Main package
│   ├── core/
│   │   ├── fetcher.py          # RCSB PDB downloader
│   │   └── prep.py             # PDB → PDBQT converter (OpenBabel wrapper)
│   ├── engine/
│   │   ├── vina.py             # Vina wrapper (physics engine)
│   │   └── grid.py             # Grid box calculator
│   ├── utils/
│   │   └── logger.py           # Logging utility
│   └── main.py                 # CLI entrypoint
├── tests/                      # Unit tests
├── Dockerfile                  # Production container (openbabel + vina)
├── pyproject.toml              # PEP 621 packaging
└── README.md                   # This file
```

## Installation

### Local Development

```bash
# Clone repository
git clone https://github.com/vihaankulkarni29/AutoScan.git
cd AutoScan

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -e .

# Install system dependencies (Linux/macOS)
# Ubuntu/Debian:
sudo apt-get install openbabel autoscan-vina

# macOS (via Homebrew):
brew install open-babel autoscan-vina
```

### Docker (Recommended for Production)

```bash
# Build image
docker build -t autoscan:latest .

# Run container
docker run --rm autoscan:latest \
  dock --receptor-pdb 3NUU --ligand-name ciprofloxacin --mutation A:87:D:G
```

## Usage

### Command-Line Interface

#### Basic Docking (Wild-Type)

```bash
autoscan dock \
  --receptor-pdb 3NUU \
  --ligand-name ciprofloxacin \
  --pocket GyrA_pocket
```

#### Docking with Mutation (Mutant)

```bash
autoscan dock \
  --receptor-pdb 3NUU \
  --ligand-name ciprofloxacin \
  --mutation A:87:D:G \
  --pocket GyrA_D87_pocket \
  --output-format json
```

**Mutation Format**: `CHAIN:RESIDUE_NUMBER:NEW_AMINO_ACID`
- `A` = Chain A
- `87` = Residue number 87
- `G` = New amino acid (Glycine)

#### Output

**JSON (for pipeline integration)**:
```json
{
  "binding_affinity_kcal_mol": -9.4,
  "rmsd_lb": 0.0,
  "rmsd_ub": 0.5,
  "receptor_pdbqt": "/app/data/receptors/3NUU_mutant.pdbqt",
  "ligand_pdbqt": "/app/data/ligands/ciprofloxacin.pdbqt"
}
```

**Text**:
```
Binding Affinity: -9.4 kcal/mol
RMSD LB: 0.0
RMSD UB: 0.5
```

### Python API

```python
from autoscan.core import PDBFetcher, PrepareVina
from autoscan.engine import GridCalculator, VinaWrapper

# Fetch receptor
fetcher = PDBFetcher()
receptor_pdb = fetcher.fetch("3NUU")

# Mutate residue
prep = PrepareVina()
mutant_pdb = prep.mutate_residue(receptor_pdb, "A", 87, "G")

# Convert to PDBQT
receptor_pdbqt = prep.pdb_to_pdbqt(mutant_pdb)
ligand_pdbqt = prep.pdb_to_pdbqt("ligand.pdb")

# Execute docking
grid_calc = GridCalculator("config/pockets.yaml")
grid_box = grid_calc.get_grid("GyrA_D87_pocket")

vina = VinaWrapper()
result = vina.dock(receptor_pdbqt, ligand_pdbqt, grid_box.to_vina_args())

print(f"Binding Affinity: {result.binding_affinity} kcal/mol")
```

## Configuration

### Pocket Definitions (`config/pockets.yaml`)

Define receptor binding pockets as grid boxes for Vina:

```yaml
pockets:
  GyrA_pocket:
    center_x: 10.5
    center_y: 20.3
    center_z: 15.8
    size_x: 25.0
    size_y: 25.0
    size_z: 25.0

  GyrA_D87_pocket:
    center_x: 10.8
    center_y: 21.1
    center_z: 16.2
    size_x: 20.0
    size_y: 20.0
    size_z: 20.0
```

## Pipeline Integration (MutationScan Module 7)

AutoScan is designed for seamless integration into MutationScan:

```python
from autoscan.core import PDBFetcher, PrepareVina
from autoscan.engine import VinaWrapper, GridCalculator

class MutationScanModule7:
  """AutoScan integration into MutationScan pipeline."""

    def run(self, pdb_id: str, ligand: str, mutations: List[str]) -> Dict:
        results = {}
        
        for mut in mutations:
            # AutoScan handles mutation, prep, and docking
            result = docking_pipeline(pdb_id, ligand, mut)
            results[mut] = result.binding_affinity
        
        return results
```

## System Requirements

- **Python**: 3.9+
- **System Tools** (CRITICAL):
  - `openbabel` (PDB ↔ PDBQT conversion)
  - `autoscan-vina` (Vina docking engine)
  - `curl` (optional, for automated PDB download)

## Anti-Hallucination Constraints

✓ **No invented libraries**: Uses subprocess wrappers around proven binary tools (`vina`, `obabel`)  
✓ **Standard Python packaging**: PEP 621 compliant `pyproject.toml`  
✓ **PDB-first**: Input/output is PDB format (PDBQT is internal)  
✓ **JSON API**: Output designed for downstream pipeline parsing  
✓ **Modular architecture**: Each component (fetch, prep, dock) is independent

## Development

### Running Tests

```bash
pytest tests/ -v --cov=src/autoscan
```

### Code Quality

```bash
black src/ tests/
pylint src/
mypy src/
```

### Build Distribution

```bash
pip install build
python -m build
```

## References

- **PDB**: https://www.rcsb.org/
- **Vina**: https://vina.scripps.edu/
- **OpenBabel**: https://openbabel.org/
- **BioPython**: https://biopython.org/

## License

MIT License. See `LICENSE` file.

## Contact

For questions or contributions, contact the Bioinformatics Team.

---

**Version**: 0.1.0 (Alpha)  
**Status**: Under active development


