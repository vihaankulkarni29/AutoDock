"""
Pre-flight check script for AutoScan.
Verifies all dependencies and critical paths before running tests.
"""

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent

def check_dependencies():
    """Verify system dependencies are installed."""
    print("=" * 60)
    print("1. DEPENDENCY CHECK")
    print("=" * 60)
    
    dependencies = {
        "obabel": ["obabel", "--version"],
        "vina": ["vina", "--version"],
    }
    
    missing = []
    for name, cmd in dependencies.items():
        if shutil.which(cmd[0]) is None:
            print(f"❌ {name}: NOT FOUND")
            missing.append(name)
        else:
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
                version = result.stdout.split('\n')[0] if result.stdout else result.stderr.split('\n')[0]
                print(f"✓ {name}: {version}")
            except Exception as e:
                print(f"⚠ {name}: Found but version check failed - {e}")
    
    # Check Python packages
    try:
        import openmm
        platform = openmm.Platform.getPlatformByName('CPU')
        print(f"✓ openmm: {platform.getName()} platform available")
    except ImportError:
        print("❌ openmm: NOT INSTALLED")
        missing.append("openmm")
    except Exception as e:
        print(f"⚠ openmm: Installed but platform check failed - {e}")
    
    try:
        from Bio import PDB
        print(f"✓ biopython: Available")
    except ImportError:
        print("❌ biopython: NOT INSTALLED")
        missing.append("biopython")
    
    try:
        import yaml
        print(f"✓ pyyaml: Available")
    except ImportError:
        print("❌ pyyaml: NOT INSTALLED")
        missing.append("pyyaml")
    
    if missing:
        print(f"\n❌ FAILED: Missing dependencies: {', '.join(missing)}")
        return False
    else:
        print("\n✅ All dependencies found!")
        return True


def check_paths():
    """Verify critical file paths exist."""
    print("\n" + "=" * 60)
    print("2. PATH CHECK")
    print("=" * 60)
    
    critical_paths = {
        "Config file": ROOT / "config" / "pockets.yaml",
        "Receptor dir": ROOT / "data" / "receptors",
        "Ligand dir": ROOT / "data" / "ligands",
        "Results dir": ROOT / "data" / "results",
    }
    
    missing = []
    for name, path in critical_paths.items():
        if path.exists():
            print(f"✓ {name}: {path}")
        else:
            print(f"❌ {name}: NOT FOUND - {path}")
            missing.append(str(path))
    
    # Check for 2XCT specifically
    receptor_2xct = ROOT / "data" / "receptors" / "2XCT.pdb"
    if receptor_2xct.exists():
        print(f"✓ Validation receptor: {receptor_2xct}")
    else:
        print(f"⚠ Validation receptor NOT FOUND: {receptor_2xct}")
        print(f"  (Will be downloaded automatically during validation)")
    
    if missing:
        print(f"\n❌ FAILED: Missing paths")
        print("Run: python setup_dirs.py to create directory structure")
        return False
    else:
        print("\n✅ All critical paths exist!")
        return True


def check_config():
    """Verify pocket configuration."""
    print("\n" + "=" * 60)
    print("3. CONFIG CHECK")
    print("=" * 60)
    
    config_file = ROOT / "config" / "pockets.yaml"
    if not config_file.exists():
        print(f"❌ Config file not found: {config_file}")
        return False
    
    try:
        import yaml
        with open(config_file, 'r') as f:
            data = yaml.safe_load(f)
        
        pockets = data.get('pockets', {})
        if 'GyrA_pocket' in pockets:
            pocket = pockets['GyrA_pocket']
            print(f"✓ GyrA_pocket found:")
            print(f"  Center: ({pocket['center_x']}, {pocket['center_y']}, {pocket['center_z']})")
            print(f"  Size: ({pocket['size_x']}, {pocket['size_y']}, {pocket['size_z']})")
        else:
            print("❌ GyrA_pocket NOT FOUND in config")
            return False
        
        print(f"\n✅ Configuration valid! ({len(pockets)} pockets defined)")
        return True
    except Exception as e:
        print(f"❌ Config parsing failed: {e}")
        return False


def main():
    """Run all pre-flight checks."""
    print("\n" + "╔" + "=" * 58 + "╗")
    print("║" + " " * 15 + "AUTOSCAN PRE-FLIGHT CHECK" + " " * 18 + "║")
    print("╚" + "=" * 58 + "╝\n")
    
    results = {
        "Dependencies": check_dependencies(),
        "Paths": check_paths(),
        "Config": check_config(),
    }
    
    print("\n" + "=" * 60)
    print("FINAL STATUS")
    print("=" * 60)
    
    for name, status in results.items():
        emoji = "✅" if status else "❌"
        print(f"{emoji} {name}: {'PASS' if status else 'FAIL'}")
    
    if all(results.values()):
        print("\n🚀 GREEN FOR LAUNCH! All systems go!")
        return 0
    else:
        print("\n🛑 NOT READY. Fix the issues above before proceeding.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
