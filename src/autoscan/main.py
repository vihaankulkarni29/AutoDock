"""
AutoScan CLI: Main entry point for the docking pipeline.

Usage:
    autoscan dock --receptor-pdb 3NUU --ligand-name ciprofloxacin --mutation D87G
"""

import json
import sys
from pathlib import Path
from typing import Optional

import typer

from autoscan.core import PDBFetcher, PrepareVina
from autoscan.engine import GridCalculator, VinaWrapper
from autoscan.utils import get_logger

app = typer.Typer(help="AutoScan: Molecular docking for binding affinity prediction")
logger = get_logger(__name__)

# Determine config directory
CONFIG_DIR = Path(__file__).parent.parent.parent / "config"


@app.command()
def dock(
    receptor_pdb: str = typer.Option(
        ..., help="PDB ID of the receptor (e.g., 3NUU)"
    ),
    ligand_name: str = typer.Option(
        ..., help="Name of the ligand (e.g., ciprofloxacin)"
    ),
    mutation: Optional[str] = typer.Option(
        None,
        help=(
            "Mutation format: RES:NEW, CHAIN:RES:NEW, or CHAIN:RES:OLDAA:NEWAA "
            "(e.g., 87:G, A:87:G, A:87:D:G)"
        ),
    ),
    pocket: str = typer.Option(
        "default", help="Pocket name from config/pockets.yaml (default: 'default')"
    ),
    output_format: str = typer.Option(
        "json", help="Output format: json or text"
    ),
    use_consensus: bool = typer.Option(
        False, help="Enable consensus scoring (uses multiple scorers if available)"
    ),
    consensus_method: str = typer.Option(
        "mean", help="Consensus method: mean, median, or weighted"
    ),
    ph: float = typer.Option(
        7.4, help="Physiological pH for protonation state (default: 7.4)"
    ),
    use_meeko: bool = typer.Option(
        True, help="Use Meeko for enhanced preparation (better charges)"
    ),
    verbose: bool = typer.Option(False, help="Enable verbose logging"),
):
    """
    Execute the docking pipeline ("Lock and Key" principle).

    Step 1: Fetch receptor PDB
    Step 2: Optionally mutate residue
    Step 3: Prepare molecules (PDB -> PDBQT) with pH-aware protonation
    Step 4: Execute docking (with optional consensus scoring)
    Step 5: Return binding affinity (ΔG)
    
    Phase 2 Enhancement: Uses Meeko for better charge assignment and pH-aware prep.
    """
    if verbose:
        logger.setLevel("DEBUG")

    try:
        logger.info("=== AutoScan Pipeline Started ===")

        # Step 1: Fetch receptor (The Lock)
        logger.info(f"Step 1: Fetching receptor {receptor_pdb}")
        fetcher = PDBFetcher(output_dir=Path("data/receptors"))
        receptor_pdb_file = fetcher.fetch(receptor_pdb)

        # Step 2: Apply mutation if specified
        if mutation:
            logger.info(f"Step 2: Applying mutation {mutation}")
            parts = mutation.split(":")
            if len(parts) == 2:
                res_num, new_aa = int(parts[0]), parts[1]
                chain_id = "A"
                old_aa = None
            elif len(parts) == 3:
                chain_id, res_num, new_aa = parts[0], int(parts[1]), parts[2]
                old_aa = None
            elif len(parts) == 4:
                chain_id, res_num, old_aa, new_aa = (
                    parts[0],
                    int(parts[1]),
                    parts[2],
                    parts[3],
                )
            else:
                raise ValueError(
                    "Mutation format: RES:NEW, CHAIN:RES:NEW, or CHAIN:RES:OLDAA:NEWAA"
                )

            if old_aa:
                PrepareVina.assert_residue_identity(
                    receptor_pdb_file, chain_id, res_num, old_aa
                )

            receptor_pdb_file = PrepareVina.mutate_residue(
                receptor_pdb_file, chain_id, res_num, new_aa
            )

        # Step 3: Prepare molecules (Phase 2 enhanced)
        logger.info(f"Step 3: Preparing molecules (pH={ph}, meeko={use_meeko})")
        prep = PrepareVina(use_meeko=use_meeko, ph=ph)
        
        logger.info("Preparing receptor with pH-aware protonation...")
        receptor_pdbqt = prep.prepare_receptor(receptor_pdb_file)

        # For ligand, allow pre-generated PDBQT to avoid re-prep
        ligand_dir = Path("data/ligands")
        ligand_pdbqt_path = ligand_dir / f"{ligand_name}.pdbqt"
        ligand_pdb = ligand_dir / f"{ligand_name}.pdb"

        if ligand_pdbqt_path.exists():
            logger.info("Using existing ligand PDBQT file")
            ligand_pdbqt = ligand_pdbqt_path
        else:
            if not ligand_pdb.exists():
                raise FileNotFoundError(
                    f"Ligand {ligand_name} not found in {ligand_dir}. "
                    "Please provide a PDB or PDBQT file for the ligand."
                )

            logger.info("Preparing ligand with flexibility detection...")
            ligand_pdbqt = prep.prepare_ligand(ligand_pdb)

        # Step 4: Execute docking (The Fit)
        logger.info("Step 4: Executing docking simulation")
        if use_consensus:
            logger.info(f"Consensus scoring ENABLED (method: {consensus_method})")
        
        config_file = CONFIG_DIR / "pockets.yaml"
        grid_calc = GridCalculator(str(config_file))
        grid_box = grid_calc.get_grid(pocket)

        vina = VinaWrapper()
        result = vina.dock(
            receptor_pdbqt,
            ligand_pdbqt,
            grid_box.to_vina_args(),
            use_consensus=use_consensus,
            consensus_method=consensus_method,
        )

        # Step 5: Return verdict (The Verdict)
        logger.info("Step 5: Returning binding affinity")

        if output_format == "json":
            output = vina.to_json(result)
            print(output)
        else:
            if use_consensus and result.consensus_affinity:
                print(
                    f"Vina Binding Affinity: {result.binding_affinity} kcal/mol\n"
                    f"Consensus Binding Affinity: {result.consensus_affinity:.2f} ± {result.consensus_uncertainty:.2f} kcal/mol\n"
                    f"Individual Scores: {result.consensus_scores}\n"
                    f"RMSD LB: {result.rmsd_lb}\n"
                    f"RMSD UB: {result.rmsd_ub}"
                )
            else:
                print(
                    f"Binding Affinity: {result.binding_affinity} kcal/mol\n"
                    f"RMSD LB: {result.rmsd_lb}\n"
                    f"RMSD UB: {result.rmsd_ub}"
                )

        logger.info("=== AutoScan Pipeline Completed Successfully ===")

    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=verbose)
        sys.exit(1)


@app.command()
def prepare(
    pdb_file: str = typer.Option(..., help="Path to PDB file"),
    output_format: str = typer.Option(
        "pdbqt", help="Output format: pdbqt or mutant"
    ),
    mutation: Optional[str] = typer.Option(
        None, help="Mutation in format CHAIN:RESIDUE:NEWAA"
    ),
):
    """Prepare a PDB file (mutation and PDBQT conversion)."""
    try:
        pdb_path = Path(pdb_file)
        if not pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")

        prep = PrepareVina()

        if mutation:
            logger.info(f"Applying mutation: {mutation}")
            parts = mutation.split(":")
            chain_id, res_num, new_aa = parts[0], int(parts[1]), parts[2]
            pdb_path = prep.mutate_residue(pdb_path, chain_id, res_num, new_aa)

        if output_format == "pdbqt":
            output_path = prep.pdb_to_pdbqt(pdb_path)
            print(f"Prepared file: {output_path}")
        else:
            print(f"Prepared file: {pdb_path}")

    except Exception as e:
        logger.error(f"Preparation failed: {e}")
        sys.exit(1)


def main():
    """Entry point for CLI."""
    app()


if __name__ == "__main__":
    main()




