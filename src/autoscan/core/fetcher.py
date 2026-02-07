"""
RCSB PDB Fetcher: Download receptor proteins from the Protein Data Bank.
"""

import os
import subprocess
from pathlib import Path
from typing import Optional

from autoscan.utils import get_logger

logger = get_logger(__name__)


class PDBFetcher:
    """Fetch PDB structures from RCSB PDB."""

    def __init__(self, output_dir: Optional[Path] = None):
        """
        Initialize PDB fetcher.

        Args:
            output_dir: Directory to store downloaded PDB files.
                       Defaults to data/receptors.
        """
        if output_dir is None:
            output_dir = Path(__file__).parent.parent.parent.parent / "data" / "receptors"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.base_url = "https://files.rcsb.org/download"

    def fetch(self, pdb_id: str, force: bool = False) -> Path:
        """
        Fetch a PDB structure by ID.

        Args:
            pdb_id: PDB ID (e.g., "3NUU").
            force: If True, re-download even if file exists.

        Returns:
            Path to the downloaded PDB file.

        Raises:
            RuntimeError: If download fails.
        """
        pdb_id = pdb_id.upper()
        pdb_file = self.output_dir / f"{pdb_id}.pdb"

        if pdb_file.exists() and not force:
            logger.info(f"PDB file already exists: {pdb_file}")
            return pdb_file

        url = f"{self.base_url}/{pdb_id}.pdb"
        logger.info(f"Downloading {pdb_id} from {url}")

        try:
            result = subprocess.run(
                ["curl", "-o", str(pdb_file), url],
                check=True,
                capture_output=True,
                text=True,
            )
            logger.info(f"Successfully downloaded: {pdb_file}")
            return pdb_file
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to download {pdb_id}: {e.stderr}")
        except FileNotFoundError:
            raise RuntimeError(
                "curl not found. Please install curl or use a PDB file directly."
            )




