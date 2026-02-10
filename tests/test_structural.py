import os
import shutil
from pathlib import Path
from typing import Iterable

import numpy as np
import pytest
from Bio.PDB import PDBIO, PDBList, PDBParser, Select

from autoscan.docking.utils import calculate_grid_box
from autoscan.docking.vina import VinaEngine


WORKSPACE = Path("workspace/test_structural")
PDB_ID = "2XCT"
LIGAND_CODE = "CPF"
CHAINS = {"A", "B"}
HEAVY_ELEMENTS = {"C", "N", "O", "F", "CL"}


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


class _AtomCoord:
    def __init__(self, coords: Iterable[float]):
        self.coords = list(coords)


class _LigandMol:
    def __init__(self, coords: Iterable[Iterable[float]]):
        self.atoms = [_AtomCoord(coord) for coord in coords]


def split_structure(pdb_path: Path, receptor_out: Path, ligand_out: Path) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(PDB_ID, str(pdb_path))

    ligand_chain_id = None
    ligand_residue_id = None
    for chain in structure.get_chains():
        for residue in chain:
            if residue.id[0] != " " and residue.resname == LIGAND_CODE:
                ligand_chain_id = chain.id
                ligand_residue_id = residue.id
                break
        if ligand_chain_id is not None:
            break

    if ligand_chain_id is None:
        raise RuntimeError(f"Ligand {LIGAND_CODE} not found in {pdb_path}")

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(receptor_out), ProteinSelect())
    io.save(str(ligand_out), LigandSelect(ligand_chain_id, ligand_residue_id))


def random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(3, 3))
    q, _ = np.linalg.qr(matrix)
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


def run_obabel(input_path: Path, output_path: Path, obabel_exe: str) -> None:
    cmd = [
        obabel_exe,
        str(input_path),
        "-O",
        str(output_path),
        "-h",
        "-xr",
    ]
    result = shutil.which(obabel_exe)
    if result is None and not Path(obabel_exe).exists():
        raise RuntimeError("obabel not found; set OBABEL_EXE or add to PATH")

    completed = os.spawnv(os.P_WAIT, obabel_exe, cmd)
    if completed != 0:
        raise RuntimeError(f"OpenBabel failed for {input_path}")


def ensure_pdbqt_root(pdbqt_path: Path) -> None:
    lines = pdbqt_path.read_text(encoding="utf-8").splitlines()
    if any(line.startswith("ROOT") for line in lines):
        return

    remarks = [line for line in lines if line.startswith("REMARK")]
    atoms = [line for line in lines if line.startswith("ATOM") or line.startswith("HETATM")]
    if not atoms:
        raise RuntimeError(f"No ATOM records found in {pdbqt_path}")

    wrapped = []
    wrapped.extend(remarks)
    wrapped.append("ROOT")
    wrapped.extend(atoms)
    wrapped.append("ENDROOT")
    wrapped.append("TORSDOF 0")
    pdbqt_path.write_text("\n".join(wrapped) + "\n", encoding="utf-8")


def parse_pdb_coords(pdb_path: Path) -> list[tuple[str, np.ndarray]]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("ref", str(pdb_path))
    atoms = []
    for atom in structure.get_atoms():
        element = atom.element.strip().upper() if atom.element else atom.get_name()[:2].upper()
        if element in HEAVY_ELEMENTS:
            atoms.append(atom)
    if not atoms:
        raise RuntimeError(f"No heavy atoms found in {pdb_path}")

    labels = []
    counts: dict[str, int] = {}
    for atom in atoms:
        name = atom.get_name().strip()
        counts[name] = counts.get(name, 0) + 1
        labels.append((f"{name}:{counts[name]}", atom.get_coord()))
    return [(label, np.array(coord, dtype=float)) for label, coord in labels]


def parse_pdbqt_coords(pdbqt_path: Path) -> list[tuple[str, np.ndarray]]:
    coords = []
    counts: dict[str, int] = {}
    with open(pdbqt_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                element = line[76:78].strip().upper()
                if not element:
                    element = atom_name[:2].strip().upper()
                if element not in HEAVY_ELEMENTS:
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                counts[atom_name] = counts.get(atom_name, 0) + 1
                label = f"{atom_name}:{counts[atom_name]}"
                coords.append((label, np.array([x, y, z], dtype=float)))
    if not coords:
        raise RuntimeError(f"No heavy atom coordinates found in {pdbqt_path}")
    return coords


def match_coords(
    reference: list[tuple[str, np.ndarray]],
    docked: list[tuple[str, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray]:
    docked_map = {label: coord for label, coord in docked}
    ref_coords = []
    docked_coords = []
    for label, coord in reference:
        if label in docked_map:
            ref_coords.append(coord)
            docked_coords.append(docked_map[label])
    if len(ref_coords) < 3:
        raise RuntimeError("Insufficient matched atoms for RMSD calculation")
    return np.array(ref_coords), np.array(docked_coords)


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
    return float(np.sqrt((diff * diff).sum() / ref.shape[0]))


@pytest.fixture(scope="module")
def setup_workspace():
    if WORKSPACE.exists():
        shutil.rmtree(WORKSPACE)
    WORKSPACE.mkdir(parents=True, exist_ok=True)
    yield WORKSPACE
    if os.environ.get("KEEP_WORKSPACE") != "1":
        shutil.rmtree(WORKSPACE)


def test_redocking_accuracy(setup_workspace):
    print(f"\n[INFO] Starting Redocking Test in {setup_workspace}...")

    obabel_exe = os.environ.get("OBABEL_EXE", "obabel")
    vina_exe = os.environ.get("VINA_EXE", "tools/vina.exe")

    pdbl = PDBList()
    pdb_file = pdbl.retrieve_pdb_file(PDB_ID, pdir=str(setup_workspace), file_format="pdb")
    pdb_path = Path(pdb_file)

    receptor_pdb = setup_workspace / "receptor.pdb"
    ligand_crystal = setup_workspace / "ligand_crystal.pdb"
    ligand_input = setup_workspace / "ligand_input.pdb"
    receptor_pdbqt = setup_workspace / "receptor.pdbqt"
    ligand_pdbqt = setup_workspace / "ligand_input.pdbqt"
    docked_pdbqt = setup_workspace / "ligand_docked.pdbqt"

    split_structure(pdb_path, receptor_pdb, ligand_crystal)
    shuffle_ligand(ligand_crystal, ligand_input)

    run_obabel(receptor_pdb, receptor_pdbqt, obabel_exe)
    run_obabel(ligand_input, ligand_pdbqt, obabel_exe)
    ensure_pdbqt_root(ligand_pdbqt)

    reference_atoms = parse_pdb_coords(ligand_crystal)
    ref_coords = np.array([coord for _, coord in reference_atoms])
    center = ref_coords.mean(axis=0).tolist()

    ligand_mol = _LigandMol(ref_coords)
    calculate_grid_box(center, ligand_mol=ligand_mol)

    engine = VinaEngine(str(receptor_pdbqt), str(ligand_pdbqt), vina_executable=vina_exe)
    score = engine.run(center=center, ligand_mol=ligand_mol, output_pdbqt=str(docked_pdbqt))

    docked_atoms = parse_pdbqt_coords(docked_pdbqt)
    reference_coords, docked_coords = match_coords(reference_atoms, docked_atoms)
    rmsd = kabsch_rmsd(reference_coords, docked_coords)

    assert score < 0.0
    assert rmsd < 2.5
    print("[PASS] Structural Validation Complete.")
