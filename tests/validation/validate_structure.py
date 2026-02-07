"""
Redocking benchmark for AutoScan using PDB 2XCT (S. aureus gyrase + ciprofloxacin).
"""

from __future__ import annotations

import math
import shutil
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, PDBList, PDBParser, Select

from autoscan.core import PrepareVina
from autoscan.engine import GridCalculator, VinaWrapper


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
RECEPTOR_DIR = DATA_DIR / "receptors"
LIGAND_DIR = DATA_DIR / "ligands"
OUTPUT_DIR = Path(__file__).resolve().parent / "output"

PDB_ID = "2XCT"
LIGAND_RESNAME = "CPF"
CHAINS = {"A", "B"}
POCKET_NAME = "GyrA_pocket"


class ProteinSelect(Select):
    def accept_chain(self, chain):
        return chain.id in CHAINS

    def accept_residue(self, residue):
        return residue.id[0] == " "


class LigandSelect(Select):
    def __init__(self, chain_id: str, residue_id):
        self.chain_id = chain_id
        self.residue_id = residue_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id == self.residue_id


def fetch_pdb(pdb_id: str, output_path: Path) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    RECEPTOR_DIR.mkdir(parents=True, exist_ok=True)
    pdbl = PDBList()
    downloaded = pdbl.retrieve_pdb_file(pdb_id, pdir=str(OUTPUT_DIR), file_format="pdb")
    downloaded_path = Path(downloaded)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(downloaded_path, output_path)
    return output_path


def split_structure(pdb_path: Path, receptor_out: Path, ligand_out: Path) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(PDB_ID, str(pdb_path))

    ligand_chain_id = None
    ligand_residue_id = None
    for chain in structure.get_chains():
        for residue in chain:
            if residue.id[0] != " " and residue.resname == LIGAND_RESNAME:
                ligand_chain_id = chain.id
                ligand_residue_id = residue.id
                break
        if ligand_chain_id is not None:
            break

    if ligand_chain_id is None:
        raise RuntimeError(f"Ligand {LIGAND_RESNAME} not found in {pdb_path}")

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(receptor_out), ProteinSelect())
    io.save(str(ligand_out), LigandSelect(ligand_chain_id, ligand_residue_id))


def random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(3, 3))
    q, r = np.linalg.qr(matrix)
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def shuffle_ligand(ligand_in: Path, ligand_out: Path, seed: int = 7) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("ligand", str(ligand_in))
    atoms = list(structure.get_atoms())
    coords = np.array([atom.get_coord() for atom in atoms], dtype=float)

    rng = np.random.default_rng(seed)
    rotation = random_rotation_matrix(rng)
    translation = rng.uniform(-10.0, 10.0, size=(3,))

    shuffled = (coords @ rotation.T) + translation
    for atom, coord in zip(atoms, shuffled):
        atom.set_coord(coord)

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(ligand_out))


def parse_pdb_coords(pdb_path: Path) -> np.ndarray:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("ref", str(pdb_path))
    atoms = list(structure.get_atoms())
    if not atoms:
        raise RuntimeError(f"No atoms found in {pdb_path}")
    return np.array([atom.get_coord() for atom in atoms], dtype=float)


def parse_pdbqt_coords(pdbqt_path: Path) -> np.ndarray:
    coords = []
    with open(pdbqt_path, "r") as handle:
        for line in handle:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                coords.append([x, y, z])
    if not coords:
        raise RuntimeError(f"No atom coordinates found in {pdbqt_path}")
    return np.array(coords, dtype=float)


def kabsch_rmsd(reference: np.ndarray, mobile: np.ndarray) -> float:
    if reference.shape != mobile.shape:
        raise RuntimeError(
            f"Atom count mismatch: {reference.shape[0]} vs {mobile.shape[0]}"
        )

    ref_center = reference.mean(axis=0)
    mob_center = mobile.mean(axis=0)
    ref = reference - ref_center
    mob = mobile - mob_center

    covariance = mob.T @ ref
    u, _, vt = np.linalg.svd(covariance)
    rotation = u @ vt

    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = u @ vt

    aligned = mob @ rotation
    diff = aligned - ref
    rmsd = math.sqrt((diff * diff).sum() / ref.shape[0])
    return rmsd


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    missing = []
    if shutil.which("obabel") is None:
        missing.append("obabel (OpenBabel)")
    if shutil.which("autoscan-vina") is None:
        missing.append("autoscan-vina")
    if missing:
        print("Missing required system tools:")
        for tool in missing:
            print(f"- {tool}")
        print("Install them and re-run this validation script.")
        sys.exit(1)

    receptor_ref = OUTPUT_DIR / "receptor_ref.pdb"
    ligand_ref = OUTPUT_DIR / "ligand_ref.pdb"
    ligand_shuffled = OUTPUT_DIR / "ligand_shuffled.pdb"
    receptor_pdbqt = OUTPUT_DIR / "receptor_ref.pdbqt"
    ligand_pdbqt = OUTPUT_DIR / "ligand_shuffled.pdbqt"
    docked_pdbqt = OUTPUT_DIR / "ligand_shuffled_docked.pdbqt"

    pdb_path = RECEPTOR_DIR / f"{PDB_ID}.pdb"
    if not pdb_path.exists():
        pdb_path = fetch_pdb(PDB_ID, pdb_path)

    split_structure(pdb_path, receptor_ref, ligand_ref)
    shuffle_ligand(ligand_ref, ligand_shuffled)

    prep = PrepareVina(use_meeko=True, ph=7.4)
    prep.pdb_to_pdbqt(receptor_ref, receptor_pdbqt, molecule_type="receptor")
    prep.pdb_to_pdbqt(ligand_shuffled, ligand_pdbqt, molecule_type="ligand")

    grid_calc = GridCalculator(str(ROOT / "config" / "pockets.yaml"))
    grid_box = grid_calc.get_grid(POCKET_NAME)

    vina = VinaWrapper()
    vina.dock(
        receptor_pdbqt,
        ligand_pdbqt,
        grid_box.to_vina_args(),
        output_pdbqt=docked_pdbqt,
        num_modes=1,
    )

    reference_coords = parse_pdb_coords(ligand_ref)
    docked_coords = parse_pdbqt_coords(docked_pdbqt)

    rmsd = kabsch_rmsd(reference_coords, docked_coords)
    status = "PASS" if rmsd < 2.5 else "FAIL"

    print(f"RMSD: {rmsd:.3f} A")
    print(status)


if __name__ == "__main__":
    main()
