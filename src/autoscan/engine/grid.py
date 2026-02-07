"""
Grid Box Calculator: Generate Vina search space from pocket definitions.
"""

from dataclasses import dataclass
from typing import Dict, List

import yaml

from autoscan.utils import get_logger

logger = get_logger(__name__)


@dataclass
class GridBox:
    """Represents a Vina grid box."""

    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float

    def to_vina_args(self) -> List[str]:
        """
        Convert grid box to Vina command-line arguments.

        Returns:
            List of command-line arguments for Vina.
        """
        return [
            "--center_x",
            str(self.center_x),
            "--center_y",
            str(self.center_y),
            "--center_z",
            str(self.center_z),
            "--size_x",
            str(self.size_x),
            "--size_y",
            str(self.size_y),
            "--size_z",
            str(self.size_z),
        ]


class GridCalculator:
    """Load and manage grid boxes for receptors."""

    def __init__(self, config_file: str):
        """
        Initialize grid calculator from YAML config.

        Args:
            config_file: Path to config/pockets.yaml.
        """
        self.config_file = config_file
        self.pockets: Dict[str, GridBox] = {}
        self._load_pockets()

    def _load_pockets(self):
        """Load pocket definitions from YAML."""
        logger.info(f"Loading pocket definitions from {self.config_file}")
        try:
            with open(self.config_file, "r") as f:
                data = yaml.safe_load(f)
            for pocket_name, coords in data.get("pockets", {}).items():
                self.pockets[pocket_name] = GridBox(
                    center_x=coords["center_x"],
                    center_y=coords["center_y"],
                    center_z=coords["center_z"],
                    size_x=coords["size_x"],
                    size_y=coords["size_y"],
                    size_z=coords["size_z"],
                )
            logger.info(f"Loaded {len(self.pockets)} pocket definitions")
        except Exception as e:
            raise RuntimeError(f"Failed to load pockets config: {e}")

    def get_grid(self, pocket_name: str) -> GridBox:
        """
        Retrieve grid box for a pocket.

        Args:
            pocket_name: Name of the pocket (e.g., "GyrA_pocket").

        Returns:
            GridBox instance.

        Raises:
            KeyError: If pocket not found.
        """
        if pocket_name not in self.pockets:
            raise KeyError(
                f"Pocket '{pocket_name}' not found. Available: {list(self.pockets.keys())}"
            )
        return self.pockets[pocket_name]




