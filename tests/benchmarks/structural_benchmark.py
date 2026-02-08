"""Multi-target redocking benchmark for AutoScan (structural validation)."""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
from Bio.PDB import PDBIO, PDBList, PDBParser, Select

from autoscan.docking.utils import calculate_grid_box
from autoscan.docking.vina import VinaEngine


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(__file__).resolve().parent / "data"
RESULTS_DIR = Path(__file__).resolve().parent / "results"

TARGETS: dict[str, dict[str, object]] = {
    "2XCT": {"ligand": "CPF", "chains": {"A", "B"}},
    "1HXB": {"ligand": "ROC", "chains": None},
    "1IEP": {"ligand": "STI", "chains": None},
    "3ERT": {"ligand": "OHT", "chains": None},
    "5UH6": {"ligand": "RFP", "chains": None},
}

HEAVY_ELEMENTS = {"C", "N", "O", "F", "CL"}


class ProteinSelect(Select):
    def __init__(self, chains: set[str] | None):
        self.chains = chains

    def accept_chain(self, chain):
        if self.chains is None:
            return True
        return chain.id in self.chains

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


def fetch_pdb(pdb_id: str, output_path: Path) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    pdbl = PDBList()
    try:
        downloaded = pdbl.retrieve_pdb_file(pdb_id, pdir=str(DATA_DIR), file_format="pdb")
        downloaded_path = Path(downloaded)
        shutil.copyfile(downloaded_path, output_path)
        return output_path
    except Exception:
        import urllib.request

        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        urllib.request.urlretrieve(url, output_path)
        return output_path


def split_structure(
    pdb_path: Path,
    receptor_out: Path,
    ligand_out: Path,
    pdb_id: str,
    ligand_code: str,
    chains: set[str] | None,
) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(pdb_path))

    ligand_chain_id = None
    ligand_residue_id = None
    for chain in structure.get_chains():
        for residue in chain:
            if residue.id[0] != " " and residue.resname == ligand_code:
                ligand_chain_id = chain.id
                ligand_residue_id = residue.id
                break
        if ligand_chain_id is not None:
            break

    if ligand_chain_id is None:
        raise RuntimeError(f"Ligand {ligand_code} not found in {pdb_path}")

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(receptor_out), ProteinSelect(chains))
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
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "OpenBabel failed: "
            f"{result.stdout.strip()} {result.stderr.strip()}"
        )


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


def _element_from_atom(atom_name: str, element_field: str) -> str:
    if element_field:
        return element_field.upper()
    return atom_name[:2].strip().upper()


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
    with open(pdbqt_path, "r") as handle:
        for line in handle:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                element_field = line[76:78].strip()
                element = _element_from_atom(atom_name, element_field)
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
    rmsd = math.sqrt((diff * diff).sum() / ref.shape[0])
    return rmsd


def _write_csv(rows: list[dict[str, object]], output_path: Path) -> None:
    import csv

    header = ["Target", "Crystal_Energy", "Docked_Energy", "RMSD", "Status"]
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    poses_dir = RESULTS_DIR / "poses"
    poses_dir.mkdir(parents=True, exist_ok=True)

    obabel_exe = os.environ.get("OBABEL_EXE", "obabel")
    vina_exe = os.environ.get("VINA_EXE", "vina")

    for tool_name, tool_cmd in ("obabel", obabel_exe), ("vina", vina_exe):
        if shutil.which(tool_cmd) is None and not Path(tool_cmd).exists():
            print(f"Missing required system tool: {tool_name} ({tool_cmd})")
            print("Set OBABEL_EXE or VINA_EXE if the binaries are not on PATH.")
            sys.exit(1)

    rows: list[dict[str, object]] = []
    rmsd_values: list[float] = []
    hard_targets: list[str] = []
    failures: list[str] = []

    for pdb_id, meta in TARGETS.items():
        ligand_code = str(meta["ligand"])
        chains = meta["chains"]
        target_dir = poses_dir / pdb_id
        target_dir.mkdir(parents=True, exist_ok=True)

        pdb_path = DATA_DIR / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            fetch_pdb(pdb_id, pdb_path)

        receptor_pdb = target_dir / "receptor.pdb"
        ligand_crystal = target_dir / "ligand_crystal.pdb"
        ligand_input = target_dir / "ligand_input.pdb"
        receptor_pdbqt = target_dir / "receptor.pdbqt"
        ligand_pdbqt = target_dir / "ligand_input.pdbqt"
        docked_pdbqt = target_dir / "ligand_docked.pdbqt"

        try:
            split_structure(pdb_path, receptor_pdb, ligand_crystal, pdb_id, ligand_code, chains)
            shuffle_ligand(ligand_crystal, ligand_input)

            run_obabel(receptor_pdb, receptor_pdbqt, obabel_exe)
            run_obabel(ligand_input, ligand_pdbqt, obabel_exe)
            ensure_pdbqt_root(ligand_pdbqt)

            reference_atoms = parse_pdb_coords(ligand_crystal)
            ref_coords = np.array([coord for _, coord in reference_atoms])
            center = ref_coords.mean(axis=0).tolist()

            ligand_mol = _LigandMol(ref_coords)
            calculate_grid_box(center, ligand_mol=ligand_mol)

            engine = VinaEngine(
                str(receptor_pdbqt),
                str(ligand_pdbqt),
                vina_executable=vina_exe,
            )
            affinity = engine.run(
                center=center,
                ligand_mol=ligand_mol,
                output_pdbqt=str(docked_pdbqt),
            )

            docked_atoms = parse_pdbqt_coords(docked_pdbqt)
            reference_coords, docked_coords = match_coords(reference_atoms, docked_atoms)
            rmsd = kabsch_rmsd(reference_coords, docked_coords)

            status = "PASS" if rmsd < 2.5 else "FAIL"
            if status == "PASS":
                rmsd_values.append(rmsd)
            else:
                hard_targets.append(pdb_id)

            rows.append(
                {
                    "Target": f"{pdb_id}:{ligand_code}",
                    "Crystal_Energy": "NA",
                    "Docked_Energy": f"{affinity:.3f}",
                    "RMSD": f"{rmsd:.3f}",
                    "Status": status,
                }
            )
        except Exception as exc:
            failures.append(f"{pdb_id}:{ligand_code} -> {exc}")
            rows.append(
                {
                    "Target": f"{pdb_id}:{ligand_code}",
                    "Crystal_Energy": "NA",
                    "Docked_Energy": "NA",
                    "RMSD": "NA",
                    "Status": "ERROR",
                }
            )

    csv_path = RESULTS_DIR / "benchmark_summary.csv"
    _write_csv(rows, csv_path)

    mean_rmsd = float(np.mean(rmsd_values)) if rmsd_values else float("nan")
    success_rate = (len([r for r in rows if r["Status"] == "PASS"]) / len(rows)) * 100.0

    print("\nSummary")
    print("| Target | Crystal_Energy | Docked_Energy | RMSD (A) | Status |")
    print("|---|---:|---:|---:|---|")
    for row in rows:
        print(
            f"| {row['Target']} | {row['Crystal_Energy']} | {row['Docked_Energy']} | "
            f"{row['RMSD']} | {row['Status']} |"
        )

    print("\nAggregate")
    print(f"Mean RMSD: {mean_rmsd:.3f} A" if rmsd_values else "Mean RMSD: NA")
    print(f"Success rate: {success_rate:.1f}%")
    if hard_targets:
        print("Hard targets (RMSD > 2.5 A): " + ", ".join(hard_targets))
    if failures:
        print("Failures:")
        for failure in failures:
            print(f"- {failure}")

    suite_status = "✅ PASS" if rmsd_values and mean_rmsd < 2.5 else "❌ FAIL"
    print(f"\nSuite verdict: {suite_status}")


if __name__ == "__main__":
    main()
