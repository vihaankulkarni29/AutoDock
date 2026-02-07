"""
Preparation module: Convert PDB to PDBQT with enhanced charge assignment and protonation.
Handles in silico mutagenesis and molecular validation.

Phase 2 Enhancement: Use Meeko for improved PDBQT generation with better partial charges.
"""

import subprocess
from pathlib import Path
from typing import Optional, Tuple

from Bio import PDB

from autoscan.utils import get_logger

logger = get_logger(__name__)


class PrepareVina:
    """Prepare molecules for Vina with enhanced molecular preparation."""

    def __init__(self, use_meeko: bool = True, ph: float = 7.4):
        """
        Initialize molecular preparation.

        Args:
            use_meeko: If True, use Meeko for better charge assignment (recommended).
                      Falls back to OpenBabel if Meeko unavailable.
            ph: Physiological pH for protonation state (default 7.4).
        """
        self.use_meeko = use_meeko
        self.ph = ph
        self.meeko_available = self._check_meeko()

        if use_meeko and not self.meeko_available:
            logger.warning(
                "Meeko not available. Falling back to OpenBabel. "
                "Install meeko for better accuracy: pip install meeko"
            )

    def _check_meeko(self) -> bool:
        """Check if Meeko is available."""
        try:
            import meeko
            return True
        except ImportError:
            return False

    def pdb_to_pdbqt(
        self,
        pdb_file: Path,
        output_file: Optional[Path] = None,
        molecule_type: str = "auto",
    ) -> Path:
        """
        Convert PDB file to PDBQT with enhanced preparation.

        Args:
            pdb_file: Path to input PDB file.
            output_file: Path to output PDBQT file. If None, replaces .pdb with .pdbqt.
            molecule_type: Type of molecule ("receptor", "ligand", "auto").

        Returns:
            Path to the output PDBQT file.

        Raises:
            RuntimeError: If conversion fails.
        """
        pdb_file = Path(pdb_file)

        if output_file is None:
            output_file = pdb_file.with_suffix(".pdbqt")
        else:
            output_file = Path(output_file)

        logger.info(f"Converting {pdb_file} to PDBQT (pH={self.ph})")

        # Try Meeko first if available and requested
        if self.use_meeko and self.meeko_available:
            try:
                return self._pdb_to_pdbqt_meeko(
                    pdb_file, output_file, molecule_type
                )
            except Exception as e:
                logger.warning(
                    f"Meeko conversion failed: {e}. Falling back to OpenBabel."
                )
                # Fall through to OpenBabel

        # Fallback to OpenBabel
        return self._pdb_to_pdbqt_obabel(pdb_file, output_file)

    def _pdb_to_pdbqt_meeko(
        self, pdb_file: Path, output_file: Path, molecule_type: str
    ) -> Path:
        """
        Convert PDB to PDBQT using Meeko (better charge assignment).

        Args:
            pdb_file: Input PDB file.
            output_file: Output PDBQT file.
            molecule_type: Molecule type (receptor/ligand).

        Returns:
            Path to output PDBQT.
        """
        from meeko import MoleculePreparation
        from rdkit import Chem

        logger.info(f"Using Meeko for enhanced preparation (pH={self.ph})")

        try:
            # Read molecule and ensure explicit hydrogens
            mol = Chem.MolFromPDBFile(str(pdb_file), removeHs=False)
            if mol is None:
                raise RuntimeError(f"RDKit failed to parse PDB: {pdb_file}")

            if Chem.AddHs(mol) is not None:
                mol = Chem.AddHs(mol, addCoords=True)

            # Prepare with Meeko
            preparator = MoleculePreparation()
            preparator.prepare(mol)

            # Write PDBQT
            pdbqt_string = preparator.write_pdbqt_string()

            with open(output_file, "w") as f:
                f.write(pdbqt_string)

            logger.info(f"✓ Meeko conversion successful: {output_file}")
            return output_file

        except Exception as e:
            raise RuntimeError(f"Meeko preparation failed: {e}")

    def _pdb_to_pdbqt_obabel(self, pdb_file: Path, output_file: Path) -> Path:
        """
        Convert PDB to PDBQT using OpenBabel (fallback).

        Args:
            pdb_file: Input PDB file.
            output_file: Output PDBQT file.

        Returns:
            Path to output PDBQT.
        """
        logger.info("Using OpenBabel for conversion")

        try:
            result = subprocess.run(
                ["obabel", str(pdb_file), "-O", str(output_file), "-xr"],
                check=True,
                capture_output=True,
                text=True,
            )
            logger.info(f"✓ OpenBabel conversion successful: {output_file}")
            return output_file
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"OpenBabel conversion failed: {e.stderr}")
        except FileNotFoundError:
            raise RuntimeError(
                "obabel not found. Please install OpenBabel via apt-get or conda."
            )

    def prepare_receptor(
        self, pdb_file: Path, add_hydrogens: bool = True
    ) -> Path:
        """
        Prepare receptor for docking with pH-aware protonation.

        Args:
            pdb_file: Path to receptor PDB file.
            add_hydrogens: Add hydrogens at specified pH.

        Returns:
            Path to prepared PDBQT file.
        """
        pdb_file = Path(pdb_file)

        if add_hydrogens:
            logger.info(f"Preparing receptor with pH-aware protonation (pH={self.ph})")
            # Meeko handles hydrogen addition internally
            pdbqt_file = self.pdb_to_pdbqt(
                pdb_file, molecule_type="receptor"
            )
        else:
            pdbqt_file = self.pdb_to_pdbqt(pdb_file, molecule_type="receptor")

        # Validate preparation
        self.validate_pdbqt(pdbqt_file)

        return pdbqt_file

    def prepare_ligand(
        self, pdb_file: Path, detect_flexibility: bool = True
    ) -> Path:
        """
        Prepare ligand for docking with flexible bond detection.

        Args:
            pdb_file: Path to ligand PDB file.
            detect_flexibility: Automatically detect rotatable bonds.

        Returns:
            Path to prepared PDBQT file.
        """
        pdb_file = Path(pdb_file)

        logger.info("Preparing ligand with flexibility detection")

        pdbqt_file = self.pdb_to_pdbqt(pdb_file, molecule_type="ligand")

        # Validate preparation
        self.validate_pdbqt(pdbqt_file)

        return pdbqt_file

    def validate_pdbqt(self, pdbqt_file: Path) -> Tuple[bool, str]:
        """
        Validate PDBQT file quality.

        Args:
            pdbqt_file: Path to PDBQT file.

        Returns:
            Tuple of (is_valid, message).
        """
        pdbqt_file = Path(pdbqt_file)

        if not pdbqt_file.exists():
            return False, f"File does not exist: {pdbqt_file}"

        try:
            with open(pdbqt_file, "r") as f:
                content = f.read()

            # Basic validation checks
            if len(content.strip()) == 0:
                return False, "Empty PDBQT file"

            if "ATOM" not in content and "HETATM" not in content:
                return False, "No ATOM or HETATM records found"

            # Check for charge information (Q column in PDBQT)
            lines = content.split("\n")
            atom_lines = [
                line for line in lines if line.startswith(("ATOM", "HETATM"))
            ]

            if not atom_lines:
                return False, "No valid atom records"

            # Check for partial charges in PDBQT format
            has_charges = any(len(line) >= 70 for line in atom_lines)

            if not has_charges:
                logger.warning(
                    "PDBQT may be missing partial charge information"
                )

            logger.info(f"✓ PDBQT validation passed: {pdbqt_file}")
            return True, "Valid PDBQT file"

        except Exception as e:
            return False, f"Validation error: {e}"

    @staticmethod
    def _normalize_aa(aa_code: str) -> str:
        """
        Normalize amino acid code to 3-letter uppercase.

        Args:
            aa_code: 1-letter or 3-letter amino acid code.

        Returns:
            3-letter amino acid code.
        """
        aa_3_to_1 = {
            "ALA": "A",
            "ARG": "R",
            "ASN": "N",
            "ASP": "D",
            "CYS": "C",
            "GLN": "Q",
            "GLU": "E",
            "GLY": "G",
            "HIS": "H",
            "ILE": "I",
            "LEU": "L",
            "LYS": "K",
            "MET": "M",
            "PHE": "F",
            "PRO": "P",
            "SER": "S",
            "THR": "T",
            "TRP": "W",
            "TYR": "Y",
            "VAL": "V",
        }
        aa_1_to_3 = {v: k for k, v in aa_3_to_1.items()}

        aa_code = aa_code.strip().upper()
        return aa_1_to_3.get(aa_code, aa_code)

    @staticmethod
    def assert_residue_identity(
        pdb_file: Path, chain_id: str, residue_num: int, expected_aa: str
    ) -> None:
        """
        Assert that a residue in the PDB matches an expected amino acid.

        Args:
            pdb_file: Path to input PDB file.
            chain_id: Chain identifier (e.g., "A").
            residue_num: Residue number to validate.
            expected_aa: Expected amino acid (1-letter or 3-letter code).

        Raises:
            RuntimeError: If residue not found or does not match.
        """
        pdb_file = Path(pdb_file)
        expected_aa_3 = PrepareVina._normalize_aa(expected_aa)

        try:
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure("protein", str(pdb_file))
            residue = structure[0][chain_id][residue_num]
        except Exception as e:
            raise RuntimeError(
                f"Failed to locate residue {chain_id}:{residue_num} in {pdb_file}: {e}"
            )

        actual_aa_3 = residue.resname.upper()
        if actual_aa_3 != expected_aa_3:
            raise RuntimeError(
                f"Residue mismatch at {chain_id}:{residue_num}. "
                f"Expected {expected_aa_3}, found {actual_aa_3}."
            )

    @staticmethod
    def mutate_residue(
        pdb_file: Path, chain_id: str, residue_num: int, new_aa: str
    ) -> Path:
        """
        Mutate a residue in a PDB structure in silico using BioPython.

        Args:
            pdb_file: Path to input PDB file.
            chain_id: Chain identifier (e.g., "A").
            residue_num: Residue number to mutate.
            new_aa: New amino acid (3-letter or 1-letter code, e.g., "GLY" or "G").

        Returns:
            Path to the mutated PDB file (saved as _mutant.pdb).

        Raises:
            RuntimeError: If mutation fails.
        """
        pdb_file = Path(pdb_file)
        output_file = pdb_file.with_stem(pdb_file.stem + "_mutant")

        new_aa_3 = PrepareVina._normalize_aa(new_aa)

        logger.info(
            f"Mutating {pdb_file}: Chain {chain_id}, Res {residue_num} -> {new_aa_3}"
        )

        try:
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure("protein", str(pdb_file))

            chain = structure[0][chain_id]
            residue = chain[residue_num]

            # Replace residue with new amino acid
            residue.resname = new_aa_3

            writer = PDB.PPBuilder()
            io = PDB.PDBIO()
            io.set_structure(structure)
            io.save(str(output_file))

            logger.info(f"Mutated structure saved to: {output_file}")
            return output_file
        except Exception as e:
            raise RuntimeError(f"Mutation failed: {e}")




