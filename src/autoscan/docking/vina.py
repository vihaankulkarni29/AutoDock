"""Vina engine wrapper for CLI usage."""

from pathlib import Path
from typing import Optional

from autoscan.docking.utils import calculate_grid_box
from autoscan.engine.vina import VinaWrapper


class VinaEngine:
    """Minimal Vina engine wrapper for the CLI."""

    def __init__(self, receptor_pdbqt: str, ligand_pdbqt: str, vina_executable: str = "vina"):
        self.receptor_pdbqt = Path(receptor_pdbqt)
        self.ligand_pdbqt = Path(ligand_pdbqt)
        self.vina = VinaWrapper(vina_executable=vina_executable)

    def run(
        self,
        center: list,
        ligand_mol=None,
        buffer_angstroms: float = 6.0,
        cpu: int = 4,
        num_modes: int = 9,
        output_pdbqt: Optional[str] = None,
    ) -> float:
        grid = calculate_grid_box(center, ligand_mol=ligand_mol, buffer_angstroms=buffer_angstroms)
        grid_args = [
            "--center_x",
            str(grid["center_x"]),
            "--center_y",
            str(grid["center_y"]),
            "--center_z",
            str(grid["center_z"]),
            "--size_x",
            str(grid["size_x"]),
            "--size_y",
            str(grid["size_y"]),
            "--size_z",
            str(grid["size_z"]),
        ]

        result = self.vina.dock(
            self.receptor_pdbqt,
            self.ligand_pdbqt,
            grid_args,
            output_pdbqt=output_pdbqt,
            cpu=cpu,
            num_modes=num_modes,
        )
        return result.binding_affinity
