import os
import shutil
import subprocess
import urllib.request
from datetime import date
from pathlib import Path
from typing import Iterable

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, PDBList, PDBParser, Select

from autoscan.docking.utils import calculate_grid_box
from autoscan.docking.vina import VinaEngine


WORKSPACE_ROOT = Path("workspace/phase1_fixed")
REPORT_PATH = Path("TEST_REPORT_PHASE1_FIXED.md")

TARGETS = [
    {"class": "Antibiotic", "pdb_id": "2XCT", "target": "S. aureus Gyrase", "ligand": "CPF", "challenge": "Baseline validation"},
    {"class": "Antiviral", "pdb_id": "1HSV", "target": "HIV Protease", "ligand": "MK1", "challenge": "Clean protease complex"},
    {"class": "Cancer", "pdb_id": "1IEP", "target": "Abl Kinase", "ligand": "STI", "challenge": "Induced-fit (Gleevec)"},
    {"class": "Hormone", "pdb_id": "3ERT", "target": "Estrogen Receptor", "ligand": "OHT", "challenge": "Hydrophobic pocket"},
    {"class": "Cancer", "pdb_id": "1M17", "target": "EGFR Kinase", "ligand": "AQ4", "challenge": "ATP-competitive"},
    {"class": "Steroid", "pdb_id": "1SQN", "target": "Progesterone Rec.", "ligand": "PRG", "challenge": "Rigid core"},
    {"class": "Bio-tool", "pdb_id": "1OQA", "target": "Streptavidin", "ligand": "BTN", "challenge": "Ultra-high affinity"},
    {"class": "Cancer", "pdb_id": "4DJU", "target": "B-Raf Kinase", "ligand": "032", "challenge": "V600E Mutant"},
    {"class": "Pain", "pdb_id": "3LN1", "target": "COX-2", "ligand": "CEL", "challenge": "NSAID (Celecoxib)"},
    {"class": "Flu", "pdb_id": "2HU4", "target": "Neuraminidase", "ligand": "G39", "challenge": "Polar/Charged"},
]

HEAVY_ELEMENTS = {"C", "N", "O", "F", "CL", "BR", "I", "P", "S"}


class ReceptorSelect(Select):
    def __init__(self, ligand_resname: str):
        self.ligand_resname = ligand_resname

    def accept_residue(self, residue):
        if residue.get_resname().strip() == self.ligand_resname:
            return False
        return residue.id[0] == " "


class LigandSelect(Select):
    def __init__(self, ligand_resname: str):
        self.ligand_resname = ligand_resname

    def accept_residue(self, residue):
        return residue.get_resname().strip() == self.ligand_resname


class _AtomCoord:
    def __init__(self, coords: Iterable[float]):
        self.coords = list(coords)


class _LigandMol:
    def __init__(self, coords: Iterable[Iterable[float]]):
        self.atoms = [_AtomCoord(coord) for coord in coords]


def resolve_executable(name: str, override: str) -> tuple[str, str | None]:
    if override and (Path(override).exists() or shutil.which(override)):
        return override, None
    resolved = shutil.which(name)
    if resolved:
        return resolved, None
    return override or name, f"Missing required tool: {name}"


def fetch_pdb(pdb_id: str, workspace: Path) -> Path:
    workspace.mkdir(parents=True, exist_ok=True)
    pdbl = PDBList()
    pdb_path = workspace / f"{pdb_id}.pdb"
    cif_path = workspace / f"{pdb_id}.cif"
    try:
        downloaded = pdbl.retrieve_pdb_file(pdb_id, pdir=str(workspace), file_format="pdb")
        if not downloaded or not Path(downloaded).exists():
            raise RuntimeError("PDBList download returned no file")
        shutil.copyfile(downloaded, pdb_path)
        return pdb_path
    except Exception:
        try:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            urllib.request.urlretrieve(url, pdb_path)
            return pdb_path
        except Exception:
            url = f"https://files.rcsb.org/download/{pdb_id}.cif"
            urllib.request.urlretrieve(url, cif_path)
            return cif_path


def _find_ligand(structure, ligand_code: str) -> bool:
    for residue in structure.get_residues():
        if residue.get_resname().strip() == ligand_code:
            return True
    return False


def split_structure(pdb_id: str, ligand_code: str, pdb_path: Path, receptor_out: Path, ligand_out: Path) -> None:
    if pdb_path.suffix.lower() == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure(pdb_id, str(pdb_path))

    if not _find_ligand(structure, ligand_code) and pdb_path.suffix.lower() != ".cif":
        cif_path = pdb_path.with_suffix(".cif")
        if not cif_path.exists():
            url = f"https://files.rcsb.org/download/{pdb_id}.cif"
            urllib.request.urlretrieve(url, cif_path)
        cif_parser = MMCIFParser(QUIET=True)
        structure = cif_parser.get_structure(pdb_id, str(cif_path))

    if not _find_ligand(structure, ligand_code):
        raise RuntimeError(f"Ligand {ligand_code} not found in {pdb_id}")

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(receptor_out), ReceptorSelect(ligand_code))
    io.save(str(ligand_out), LigandSelect(ligand_code))


def random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(3, 3))
    q, _ = np.linalg.qr(matrix)
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def shuffle_ligand_pdb(ligand_in: Path, ligand_out: Path, seed: int = 7) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("ligand", str(ligand_in))
    atoms = list(structure.get_atoms())
    coords = np.array([atom.get_coord() for atom in atoms], dtype=float)

    rng = np.random.default_rng(seed)
    rotation = random_rotation_matrix(rng)
    translation = rng.uniform(-3.0, 3.0, size=(3,))

    shuffled = (coords @ rotation.T) + translation
    for atom, coord in zip(atoms, shuffled):
        atom.set_coord(coord)

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(ligand_out))


def run_obabel(input_path: Path, output_path: Path, obabel_exe: str, include_charge: bool = True) -> None:
    cmd = [
        obabel_exe,
        str(input_path),
        "-O",
        str(output_path),
        "-h",
        "-xr",
    ]
    if include_charge:
        cmd.extend(["--partialcharge", "gasteiger"])

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


def format_error(message: str) -> str:
    cleaned = " ".join(message.replace("\r", " ").replace("\n", " ").split())
    if len(cleaned) > 160:
        return cleaned[:157] + "..."
    return cleaned


def append_results(rows: list[dict[str, object]], mean_rmsd: float, success_rate: float) -> None:
    if not REPORT_PATH.exists():
        header = [
            "# AutoScan Phase 1.1: Fixed Structural Blind-Docking Benchmark",
            f"**Date:** {date.today().isoformat()}",
            "**Executor:** AutoScan Validation Suite",
            "",
        ]
        REPORT_PATH.write_text("\n".join(header), encoding="utf-8")

    lines = []
    lines.append("\n## Phase 1.1 Results")
    lines.append("\n| Class | PDB ID | Ligand | Docked Energy (kcal/mol) | RMSD (A) | Status |")
    lines.append("|---|---|---|---:|---:|---|")
    for row in rows:
        lines.append(
            "| {class} | {pdb_id} | {ligand} | {energy} | {rmsd} | {status} |".format(
                **row
            )
        )
    lines.append("")
    lines.append(f"Mean RMSD: {mean_rmsd:.3f} A")
    lines.append(f"Success rate: {success_rate:.1f}%")

    with open(REPORT_PATH, "a", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def main() -> None:
    obabel_override = os.environ.get("OBABEL_EXE", "obabel")
    vina_override = os.environ.get("VINA_EXE", "tools/vina.exe")
    obabel_exe, obabel_error = resolve_executable("obabel", obabel_override)
    vina_exe, vina_error = resolve_executable("vina", vina_override)
    missing_tools = [err for err in (obabel_error, vina_error) if err]

    rows = []
    rmsd_values = []

    if WORKSPACE_ROOT.exists():
        shutil.rmtree(WORKSPACE_ROOT)
    WORKSPACE_ROOT.mkdir(parents=True, exist_ok=True)

    for target in TARGETS:
        pdb_id = target["pdb_id"]
        ligand_code = target["ligand"]
        workspace = WORKSPACE_ROOT / pdb_id
        workspace.mkdir(parents=True, exist_ok=True)

        row = {
            "class": target["class"],
            "pdb_id": pdb_id,
            "ligand": ligand_code,
            "energy": "NA",
            "rmsd": "NA",
            "status": "ERROR",
        }

        try:
            if missing_tools:
                raise RuntimeError("; ".join(missing_tools))

            pdb_path = fetch_pdb(pdb_id, workspace)
            receptor_pdb = workspace / f"{pdb_id}_rec.pdb"
            ligand_crystal = workspace / f"{pdb_id}_lig_crystal.pdb"
            ligand_random = workspace / f"{pdb_id}_lig_random.pdb"
            receptor_pdbqt = workspace / f"{pdb_id}_rec.pdbqt"
            ligand_random_pdbqt = workspace / f"{pdb_id}_lig_random.pdbqt"
            docked_pdbqt = workspace / f"{pdb_id}_lig_docked.pdbqt"

            split_structure(pdb_id, ligand_code, pdb_path, receptor_pdb, ligand_crystal)

            run_obabel(receptor_pdb, receptor_pdbqt, obabel_exe, include_charge=True)
            shuffle_ligand_pdb(ligand_crystal, ligand_random)
            run_obabel(ligand_random, ligand_random_pdbqt, obabel_exe, include_charge=True)
            ensure_pdbqt_root(ligand_random_pdbqt)

            reference_atoms = parse_pdb_coords(ligand_crystal)
            ref_coords = np.array([coord for _, coord in reference_atoms])
            center = ref_coords.mean(axis=0).tolist()
            ligand_mol = _LigandMol(ref_coords)

            engine = VinaEngine(str(receptor_pdbqt), str(ligand_random_pdbqt), vina_executable=vina_exe)
            score = engine.run(
                center=center,
                ligand_mol=ligand_mol,
                buffer_angstroms=10.0,
                output_pdbqt=str(docked_pdbqt),
            )

            docked_atoms = parse_pdbqt_coords(docked_pdbqt)
            reference_coords, docked_coords = match_coords(reference_atoms, docked_atoms)
            rmsd = kabsch_rmsd(reference_coords, docked_coords)

            row["energy"] = f"{score:.3f}"
            row["rmsd"] = f"{rmsd:.3f}"

            if rmsd < 2.5 and score < -7.0:
                row["status"] = "PASS"
            else:
                row["status"] = "FAIL"

            rmsd_values.append(rmsd)

        except Exception as exc:
            row["status"] = f"ERROR ({format_error(str(exc))})"

        rows.append(row)

    mean_rmsd = float(np.mean(rmsd_values)) if rmsd_values else float("nan")
    success_rate = (len([row for row in rows if row["status"] == "PASS"]) / len(rows)) * 100.0

    append_results(rows, mean_rmsd, success_rate)

    if os.environ.get("KEEP_WORKSPACE") != "1":
        shutil.rmtree(WORKSPACE_ROOT)


if __name__ == "__main__":
    main()
