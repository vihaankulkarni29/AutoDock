"""
Consensus Scoring Framework for Multiple Docking Scoring Functions.

This module provides infrastructure for combining multiple scoring functions
to improve binding affinity prediction accuracy. Currently implements Vina scoring
with extensible architecture for additional scorers (GNINA, RF-Score, etc.).
"""

import subprocess
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from autoscan.utils import get_logger

logger = get_logger(__name__)


@dataclass
class ScoringResult:
    """Individual scoring function result."""

    scorer_name: str
    binding_affinity: float  # ΔG in kcal/mol
    available: bool = True
    error: Optional[str] = None


@dataclass
class ConsensusScoringResult:
    """Consensus result from multiple scorers."""

    individual_scores: Dict[str, float]  # {scorer_name: affinity}
    consensus_affinity: float
    consensus_method: str  # "mean", "median", "weighted"
    uncertainty: float  # Standard deviation or confidence interval
    all_available: bool


class Scorer(ABC):
    """Abstract base class for scoring functions."""

    def __init__(self, executable: str):
        """
        Initialize scorer.

        Args:
            executable: Path or name of the scoring executable.
        """
        self.executable = executable
        self.available = self._check_availability()

    @abstractmethod
    def score(self, receptor_pdbqt: Path, ligand_pdbqt: Path) -> float:
        """
        Score a docked complex.

        Args:
            receptor_pdbqt: Path to receptor PDBQT file.
            ligand_pdbqt: Path to ligand PDBQT file.

        Returns:
            Binding affinity in kcal/mol.

        Raises:
            RuntimeError: If scoring fails.
        """
        pass

    def _check_availability(self) -> bool:
        """Check if scorer executable is available."""
        try:
            subprocess.run(
                [self.executable, "--help"],
                capture_output=True,
                timeout=5,
                check=False,
            )
            return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False


class VinaScorer(Scorer):
    """Vina scoring function."""

    def __init__(self, executable: str = "vina"):
        super().__init__(executable)

    def score(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        grid_args: Optional[List[str]] = None,
    ) -> float:
        """
        Score using Vina.

        Args:
            receptor_pdbqt: Path to receptor PDBQT.
            ligand_pdbqt: Path to ligand PDBQT.
            grid_args: Grid box parameters (required for Vina).

        Returns:
            Binding affinity in kcal/mol.
        """
        if grid_args is None:
            raise ValueError("Vina scoring requires grid_args")

        import re

        cmd = [
            self.executable,
            "--receptor",
            str(receptor_pdbqt),
            "--ligand",
            str(ligand_pdbqt),
            "--center_x",
            grid_args[1],
            "--center_y",
            grid_args[3],
            "--center_z",
            grid_args[5],
            "--size_x",
            grid_args[7],
            "--size_y",
            grid_args[9],
            "--size_z",
            grid_args[11],
            "--score_only",
        ]

        try:
            result = subprocess.run(
                cmd, check=True, capture_output=True, text=True, timeout=60
            )
            output = result.stdout + result.stderr

            # Parse binding affinity from stdout
            match = re.search(r"([-+]?\d+\.\d+)\s+kcal/mol", output)
            if match:
                affinity = float(match.group(1))
                logger.debug(f"Vina score: {affinity} kcal/mol")
                return affinity
            else:
                raise RuntimeError("Could not parse Vina output")

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Vina scoring failed: {e.stderr}")
        except subprocess.TimeoutExpired:
            raise RuntimeError("Vina scoring timed out")


class GNINAScorer(Scorer):
    """
    GNINA (GPU Accelerated Molecular Docking with Deep Learning) scoring.

    Note: GNINA must be installed separately:
    https://github.com/gnina/gnina
    """

    def __init__(self, executable: str = "gnina"):
        super().__init__(executable)
        if self.available:
            logger.info("GNINA scorer available")
        else:
            logger.warning(
                "GNINA not found. Install from: https://github.com/gnina/gnina"
            )

    def score(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        grid_args: Optional[List[str]] = None,
    ) -> float:
        """
        Score using GNINA (CNN scoring function).

        Args:
            receptor_pdbqt: Path to receptor PDBQT.
            ligand_pdbqt: Path to ligand PDBQT.
            grid_args: Grid box parameters (optional for GNINA).

        Returns:
            Binding affinity in kcal/mol.
        """
        if not self.available:
            raise RuntimeError("GNINA not available")

        import re

        cmd = [
            self.executable,
            "-r",
            str(receptor_pdbqt),
            "-l",
            str(ligand_pdbqt),
            "--score_only",
        ]

        try:
            result = subprocess.run(
                cmd, check=True, capture_output=True, text=True, timeout=120
            )
            output = result.stdout + result.stderr

            # Parse CNN affinity score
            match = re.search(r"CNNaffinity\s+:\s+([-+]?\d+\.\d+)", output)
            if match:
                affinity = float(match.group(1))
                logger.debug(f"GNINA CNN score: {affinity} kcal/mol")
                return affinity
            else:
                raise RuntimeError("Could not parse GNINA output")

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"GNINA scoring failed: {e.stderr}")
        except subprocess.TimeoutExpired:
            raise RuntimeError("GNINA scoring timed out")


class RFScoreScorer(Scorer):
    """
    RF-Score (Random Forest Scoring Function).

    Note: RF-Score must be installed separately.
    """

    def __init__(self, executable: str = "rf_score"):
        super().__init__(executable)
        if self.available:
            logger.info("RF-Score available")
        else:
            logger.warning("RF-Score not found. Consider installing for better accuracy.")

    def score(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        grid_args: Optional[List[str]] = None,
    ) -> float:
        """
        Score using RF-Score.

        Args:
            receptor_pdbqt: Path to receptor PDBQT.
            ligand_pdbqt: Path to ligand PDBQT.
            grid_args: Unused for RF-Score.

        Returns:
            Binding affinity in kcal/mol.
        """
        if not self.available:
            raise RuntimeError("RF-Score not available")

        import re

        cmd = [self.executable, str(receptor_pdbqt), str(ligand_pdbqt)]

        try:
            result = subprocess.run(
                cmd, check=True, capture_output=True, text=True, timeout=120
            )
            output = result.stdout + result.stderr

            # Parse RF-Score prediction (format varies)
            match = re.search(r"([-+]?\d+\.\d+)", output)
            if match:
                affinity = float(match.group(1))
                logger.debug(f"RF-Score: {affinity} kcal/mol")
                return affinity
            else:
                raise RuntimeError("Could not parse RF-Score output")

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"RF-Score failed: {e.stderr}")
        except subprocess.TimeoutExpired:
            raise RuntimeError("RF-Score timed out")


class ConsensusScorer:
    """Manager for consensus scoring across multiple functions."""

    def __init__(self):
        """Initialize all available scorers."""
        self.scorers: Dict[str, Scorer] = {}
        self._initialize_scorers()

    def _initialize_scorers(self):
        """Initialize scorer instances."""
        # Always initialize Vina (must be available)
        self.scorers["vina"] = VinaScorer()

        # Initialize optional scorers (may not be available)
        gnina = GNINAScorer()
        if gnina.available:
            self.scorers["gnina"] = gnina

        rf_score = RFScoreScorer()
        if rf_score.available:
            self.scorers["rf_score"] = rf_score

        logger.info(
            f"Consensus scorer initialized with: {', '.join(self.scorers.keys())}"
        )

    def score(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        grid_args: List[str],
        method: str = "mean",
        use_all: bool = False,
    ) -> ConsensusScoringResult:
        """
        Run consensus scoring.

        Args:
            receptor_pdbqt: Path to receptor PDBQT.
            ligand_pdbqt: Path to ligand PDBQT.
            grid_args: Grid box parameters for Vina.
            method: Consensus method ("mean", "median", "weighted").
            use_all: If True, try all scorers. If False, use only available ones.

        Returns:
            ConsensusScoringResult with individual and consensus scores.
        """
        individual_scores: Dict[str, float] = {}
        failed_scorers = []

        logger.info(f"Running consensus scoring (method={method})")

        for scorer_name, scorer in self.scorers.items():
            try:
                if not scorer.available:
                    logger.warning(f"Scorer {scorer_name} not available, skipping")
                    continue

                if scorer_name == "vina":
                    score = scorer.score(receptor_pdbqt, ligand_pdbqt, grid_args)
                else:
                    score = scorer.score(receptor_pdbqt, ligand_pdbqt)

                individual_scores[scorer_name] = score
                logger.info(f"{scorer_name}: {score:.2f} kcal/mol")

            except Exception as e:
                logger.warning(f"Scorer {scorer_name} failed: {e}")
                failed_scorers.append(scorer_name)
                continue

        if not individual_scores:
            raise RuntimeError("No scorers produced valid results")

        # Calculate consensus
        consensus_affinity = self._calculate_consensus(
            list(individual_scores.values()), method
        )

        # Calculate uncertainty (standard deviation)
        uncertainty = self._calculate_uncertainty(
            list(individual_scores.values())
        )

        all_available = len(failed_scorers) == 0

        result = ConsensusScoringResult(
            individual_scores=individual_scores,
            consensus_affinity=consensus_affinity,
            consensus_method=method,
            uncertainty=uncertainty,
            all_available=all_available,
        )

        logger.info(
            f"Consensus: {consensus_affinity:.2f} ± {uncertainty:.2f} kcal/mol"
        )

        return result

    @staticmethod
    def _calculate_consensus(scores: List[float], method: str = "mean") -> float:
        """Calculate consensus from multiple scores."""
        if method == "mean":
            return sum(scores) / len(scores)
        elif method == "median":
            import statistics

            return statistics.median(scores)
        elif method == "weighted":
            # Weight Vina higher (primary scorer)
            if len(scores) == 1:
                return scores[0]
            # Weighted average: Vina=0.5, others=0.5/n
            weights = [0.5] + [0.5 / (len(scores) - 1)] * (len(scores) - 1)
            return sum(s * w for s, w in zip(scores, weights))
        else:
            raise ValueError(f"Unknown consensus method: {method}")

    @staticmethod
    def _calculate_uncertainty(scores: List[float]) -> float:
        """Calculate uncertainty (std dev) from multiple scores."""
        if len(scores) < 2:
            return 0.0
        import statistics

        return statistics.stdev(scores)




