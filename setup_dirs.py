"""Setup script to create the required directory structure."""
import os
from pathlib import Path

ROOT = Path(__file__).parent

# Create all required directories
directories = [
    ROOT / "data" / "receptors",
    ROOT / "data" / "ligands",
    ROOT / "data" / "results",
    ROOT / "src" / "autoscan" / "dynamics",
    ROOT / "tests" / "unit",
    ROOT / "tests" / "benchmarks",
]

for directory in directories:
    directory.mkdir(parents=True, exist_ok=True)
    print(f"✓ Created: {directory.relative_to(ROOT)}")

# Create __init__.py files where needed
init_files = [
    ROOT / "src" / "autoscan" / "dynamics" / "__init__.py",
    ROOT / "tests" / "unit" / "__init__.py",
    ROOT / "tests" / "benchmarks" / "__init__.py",
]

for init_file in init_files:
    if not init_file.exists():
        init_file.write_text("\"\"\"Module initialization.\"\"\"\n")
        print(f"✓ Created: {init_file.relative_to(ROOT)}")

print("\n✅ Directory structure setup complete!")
