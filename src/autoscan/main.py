import typer
from autoscan.docking.vina import VinaEngine
from autoscan.utils.error_handler import ErrorHandler

app = typer.Typer()
console = ErrorHandler() # The Robustness Layer

@app.command()
def dock(
    receptor: str = typer.Option(..., help="Path to Receptor PDBQT"),
    ligand: str = typer.Option(..., help="Path to Ligand PDBQT"),
    center_x: float = typer.Option(..., help="Pocket Center X"),
    center_y: float = typer.Option(..., help="Pocket Center Y"),
    center_z: float = typer.Option(..., help="Pocket Center Z")
):
    """
    Run the AutoScan Docking Protocol.
    """
    try:
        console.log("Initializing AutoScan Docking Module...")

        # 1. Physics Check (Grid Box)
        # Note: In a real run, you'd load the ligand to calculate size
        # For CLI, we pass explicit center, but logic handles size internally

        engine = VinaEngine(receptor, ligand)

        # 2. Run Validation Logic
        score = engine.run(center=[center_x, center_y, center_z])

        console.success(f"Docking Complete! Affinity: {score} kcal/mol")

    except Exception as e:
        console.error(f"Docking Failed: {str(e)}")
        raise typer.Exit(code=1)

if __name__ == "__main__":
    app()




