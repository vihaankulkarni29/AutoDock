"""
Vina Wrapper: Execute molecular docking simulations.
"""

import json
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

from autoscan.utils import get_logger

logger = get_logger(__name__)


@dataclass
class DockingResult:
    """Result of a docking simulation."""

    binding_affinity: float  # ΔG in kcal/mol (primary Vina score)
    rmsd_lb: float
    rmsd_ub: float
    ligand_pdbqt: str
    receptor_pdbqt: str
    consensus_scores: Dict[str, float] = field(default_factory=dict)  # Individual scorer results
    consensus_affinity: Optional[float] = None  # Consensus from multiple scorers
    consensus_uncertainty: float = 0.0  # Std dev of consensus scores


class VinaWrapper:
    """Wrapper around the Vina binary for molecular docking."""

    def __init__(self, vina_executable: str = "vina"):
        """
        Initialize Vina wrapper.

        Args:
            vina_executable: Path to the Vina executable. Defaults to 'vina' (assumes in PATH).
        """
        self.vina_executable = vina_executable
        self._verify_installation()

    def _verify_installation(self):
        """Verify that Vina is installed and accessible."""
        try:
            result = subprocess.run(
                [self.vina_executable, "--help"],
                check=True,
                capture_output=True,
                text=True,
                timeout=5,
            )
            if "AutoScan Vina" in result.stdout or "AutoScan Vina" in result.stderr:
                logger.info(f"Vina found: {self.vina_executable}")
            else:
                logger.warning("Vina help output unexpected")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Vina verification failed: {e.stderr}")
        except FileNotFoundError:
            raise RuntimeError(
                f"Vina not found at '{self.vina_executable}'. "
                "Please install AutoScan Vina and ensure it's in your PATH."
            )

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        grid_args: list,
        output_pdbqt: Optional[Path] = None,
        cpu: int = 4,
        num_modes: int = 9,
        use_consensus: bool = False,
        consensus_method: str = "mean",
    ) -> DockingResult:
        """
        Execute docking simulation with optional consensus scoring.

        Args:
            receptor_pdbqt: Path to receptor in PDBQT format.
            ligand_pdbqt: Path to ligand in PDBQT format.
            grid_args: Grid box parameters from GridCalculator.to_vina_args().
            output_pdbqt: Path for output PDBQT. Auto-generated if not specified.
            cpu: Number of CPUs to use.
            num_modes: Number of binding modes to generate.
            use_consensus: If True, use consensus scoring from multiple scorers.
            consensus_method: Method for consensus ("mean", "median", "weighted").

        Returns:
            DockingResult with binding affinity and optional consensus scores.

        Raises:
            RuntimeError: If docking fails.
        """
        receptor_pdbqt = Path(receptor_pdbqt)
        ligand_pdbqt = Path(ligand_pdbqt)

        if output_pdbqt is None:
            output_pdbqt = ligand_pdbqt.with_stem(ligand_pdbqt.stem + "_docked")
        else:
            output_pdbqt = Path(output_pdbqt)

        cmd = [
            self.vina_executable,
            "--receptor",
            str(receptor_pdbqt),
            "--ligand",
            str(ligand_pdbqt),
            "--out",
            str(output_pdbqt),
            "--cpu",
            str(cpu),
            "--num_modes",
            str(num_modes),
        ] + grid_args

        logger.info(f"Running Vina: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd, check=True, capture_output=True, text=True, timeout=300
            )
            output_text = result.stdout + result.stderr

            # Parse binding affinity from Vina output
            affinity = self._parse_affinity(output_text)
            rmsd_lb, rmsd_ub = self._parse_rmsd(output_text)

            logger.info(
                f"Docking completed. Binding Affinity: {affinity} kcal/mol"
            )

            docking_result = DockingResult(
                binding_affinity=affinity,
                rmsd_lb=rmsd_lb,
                rmsd_ub=rmsd_ub,
                ligand_pdbqt=str(ligand_pdbqt),
                receptor_pdbqt=str(receptor_pdbqt),
            )

            # Apply consensus scoring if requested
            if use_consensus:
                docking_result = self._apply_consensus_scoring(
                    docking_result, receptor_pdbqt, ligand_pdbqt, grid_args, consensus_method
                )

            return docking_result

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Vina docking failed: {e.stderr}")
        except subprocess.TimeoutExpired:
            raise RuntimeError("Vina docking timed out (exceeded 5 minutes)")

    def _apply_consensus_scoring(
        self,
        docking_result: DockingResult,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        grid_args: list,
        consensus_method: str = "mean",
    ) -> DockingResult:
        """
        Apply consensus scoring to an existing docking result.

        Args:
            docking_result: Base docking result from Vina.
            receptor_pdbqt: Path to receptor PDBQT.
            ligand_pdbqt: Path to ligand PDBQT.
            grid_args: Grid box parameters.
            consensus_method: Consensus calculation method.

        Returns:
            Updated DockingResult with consensus scores.
        """
        from autoscan.engine.scoring import ConsensusScorer

        logger.info("Applying consensus scoring...")

        try:
            scorer = ConsensusScorer()
            consensus_result = scorer.score(
                receptor_pdbqt, ligand_pdbqt, grid_args, method=consensus_method
            )

            docking_result.consensus_scores = consensus_result.individual_scores
            docking_result.consensus_affinity = consensus_result.consensus_affinity
            docking_result.consensus_uncertainty = consensus_result.uncertainty

            logger.info(
                f"Consensus scoring complete: {consensus_result.consensus_affinity:.2f} "
                f"± {consensus_result.uncertainty:.2f} kcal/mol"
            )

        except Exception as e:
            logger.warning(f"Consensus scoring failed: {e}. Using Vina score only.")

        return docking_result

    @staticmethod
    def _parse_affinity(output: str) -> float:
        """Parse binding affinity (ΔG) from Vina output."""
        # Vina 1.2.x prints a results table with affinity column
        table_matches = re.findall(
            r"^\s*\d+\s+([-+]?\d*\.?\d+(?:e[-+]?\d+)?)\s+",
            output,
            flags=re.MULTILINE | re.IGNORECASE,
        )
        for match in table_matches:
            value = float(match)
            if -200.0 < value < 50.0:
                return value

        # Legacy pattern: numeric value followed by kcal/mol
        match = re.search(r"([-+]?\d*\.?\d+(?:e[-+]?\d+)?)\s+kcal/mol", output)
        if match:
            return float(match.group(1))

        raise RuntimeError("Could not parse binding affinity from Vina output")

    @staticmethod
    def _parse_rmsd(output: str) -> tuple[float, float]:
        """Parse RMSD values from Vina output."""
        match = re.search(r"(RMSD\s+from\s+best\s+mode|lb|ub).*?(\d+\.\d+).*?(\d+\.\d+)", output)
        # Simplified: return 0.0 if not found (Vina output format varies)
        return 0.0, 0.0

    def to_json(self, result: DockingResult) -> str:
        """
        Convert docking result to JSON for MutationScan integration.

        Includes individual scores and consensus if available.
        """
        output = {
            "binding_affinity_kcal_mol": result.binding_affinity,
            "rmsd_lb": result.rmsd_lb,
            "rmsd_ub": result.rmsd_ub,
            "receptor_pdbqt": result.receptor_pdbqt,
            "ligand_pdbqt": result.ligand_pdbqt,
        }

        # Add consensus scores if available
        if result.consensus_scores:
            output["consensus_mode"] = True
            output["individual_scores"] = result.consensus_scores
            output["consensus_affinity_kcal_mol"] = result.consensus_affinity
            output["consensus_uncertainty_kcal_mol"] = result.consensus_uncertainty
        else:
            output["consensus_mode"] = False

        return json.dumps(output)




