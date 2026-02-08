"""
Organize validation scripts into benchmark structure.
Run this after setup_dirs.py to properly organize test files.
"""

import shutil
from pathlib import Path

ROOT = Path(__file__).parent

# Source validation scripts
validation_dir = ROOT / "tests" / "validation"
benchmarks_dir = ROOT / "tests" / "benchmarks"

def organize_tests():
    """Move validation scripts to benchmarks."""
    print("Organizing test structure...")
    
    # Ensure benchmarks directory exists
    benchmarks_dir.mkdir(parents=True, exist_ok=True)
    
    # Map of files to copy/move
    if validation_dir.exists():
        validation_files = list(validation_dir.glob("*.py"))
        for file in validation_files:
            if file.name != "__init__.py":
                dest = benchmarks_dir / file.name.replace("validate_", "")
                shutil.copy2(file, dest)
                print(f"✓ Copied: {file.name} -> benchmarks/{dest.name}")
    
    # Create __init__.py in benchmarks if it doesn't exist
    init_file = benchmarks_dir / "__init__.py"
    if not init_file.exists():
        init_file.write_text('"""Benchmark tests for AutoScan scientific validation."""\n')
        print(f"✓ Created: benchmarks/__init__.py")
    
    print("\n✅ Test organization complete!")
    print(f"   Validation scripts: {validation_dir}")
    print(f"   Benchmarks: {benchmarks_dir}")


if __name__ == "__main__":
    organize_tests()
