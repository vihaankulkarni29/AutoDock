import json
import os
import re
import shutil
import subprocess
import sys
from importlib import metadata
from pathlib import Path


_REQUIRED_PYTHON = {
    "openmm": {"op": "==", "version": "8.0.0"},
    "numpy": {"op": "==", "version": "1.24.3"},
    "biopython": {"op": ">=", "version": "1.81"},
    "scipy": {"op": ">=", "version": "1.10"},
    "pandas": {"op": ">=", "version": "2.0"},
    "typer": {"op": ">=", "version": "0.9.0"},
    "mdtraj": {"op": ">=", "version": "1.9.9"},
}

_REQUIRED_BINARIES = {
    "vina": {
        "label": "AutoDock Vina",
        "op": ">=",
        "version": "1.2.3",
    },
    "obabel": {
        "label": "OpenBabel",
        "op": ">=",
        "version": "3.1.1",
    },
    "pdbfixer": {
        "label": "PDBFixer",
        "op": ">=",
        "version": "1.8",
    },
}

_CONDA_ENV_NAME = "autoscan-tools"
_CONDA_PYTHON_VERSION = "3.10"

_PINNED_PYTHON = {
    "openmm": "8.0.0",
    "numpy": "1.24.3",
}


def find_repo_root(start: Path) -> Path:
    for candidate in [start] + list(start.parents):
        if (candidate / "pyproject.toml").exists():
            return candidate
    return start


def _read_dependencies(pyproject_path: Path) -> list[str]:
    text = pyproject_path.read_text(encoding="utf-8")
    deps: list[str] = []
    in_block = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("dependencies") and "[" in stripped:
            in_block = True
            continue
        if in_block:
            if stripped.startswith("]"):
                break
            if stripped.startswith(("\"", "'")):
                deps.append(stripped.strip().strip(",").strip("\"'"))
    return deps


def _parse_requirement(requirement: str) -> tuple[str, str | None, str | None]:
    pattern = r"^([A-Za-z0-9_.-]+)(\[[^\]]+\])?\s*([<>=!~]{1,2})\s*([0-9][^;]*)"
    match = re.match(pattern, requirement.strip())
    if match:
        name = match.group(1)
        op = match.group(3)
        version = match.group(4).strip()
    else:
        name = re.split(r"[<=>]", requirement, 1)[0]
        op = None
        version = None
    name = name.strip().lower().replace("_", "-")
    return name, op, version


def _version_tuple(version: str) -> tuple[int, ...]:
    parts = re.split(r"[^0-9]+", version)
    return tuple(int(part) for part in parts if part)


def _compare_versions(installed: str, required: str, op: str) -> bool:
    try:
        from packaging.version import Version

        left = Version(installed)
        right = Version(required)
    except Exception:
        left = _version_tuple(installed)
        right = _version_tuple(required)

    if op == "==":
        return left == right
    if op == ">=":
        return left >= right
    return False


def _check_python_dependencies() -> list[str]:
    issues: list[str] = []
    for name, spec in _REQUIRED_PYTHON.items():
        op = spec["op"]
        version = spec["version"]
        try:
            installed = metadata.version(name)
        except metadata.PackageNotFoundError:
            issues.append(
                f"Missing Python dependency: {name} (required {op}{version})"
            )
            continue
        if not _compare_versions(installed, version, op):
            issues.append(
                f"Python dependency mismatch: {name} {installed} does not satisfy {op}{version}"
            )
    return issues


def _check_manifest_in_pyproject(pyproject_path: Path) -> list[str]:
    if not pyproject_path.exists():
        return ["pyproject.toml not found; cannot verify dependency manifest."]

    deps = _read_dependencies(pyproject_path)
    dep_map: dict[str, tuple[str | None, str | None]] = {}
    for requirement in deps:
        name, op, version = _parse_requirement(requirement)
        dep_map[name] = (op, version)

    issues: list[str] = []
    for name, spec in _REQUIRED_PYTHON.items():
        required_op = spec["op"]
        required_version = spec["version"]
        if name not in dep_map:
            issues.append(
                f"pyproject.toml missing required dependency: {name} ({required_op}{required_version})"
            )
            continue
        present_op, present_version = dep_map[name]
        if present_op is None or present_version is None:
            issues.append(
                f"pyproject.toml missing version spec for {name} (required {required_op}{required_version})"
            )
            continue
        if required_op == "==":
            if not (present_op == "==" and present_version == required_version):
                issues.append(
                    f"pyproject.toml spec mismatch for {name}: found {present_op}{present_version}, required {required_op}{required_version}"
                )
        else:
            if present_op == "==":
                if not _compare_versions(present_version, required_version, ">="):
                    issues.append(
                        f"pyproject.toml spec too strict for {name}: found {present_op}{present_version}, required {required_op}{required_version}"
                    )
            elif present_op == ">=":
                if not _compare_versions(present_version, required_version, ">="):
                    issues.append(
                        f"pyproject.toml spec too low for {name}: found {present_op}{present_version}, required {required_op}{required_version}"
                    )
            else:
                issues.append(
                    f"pyproject.toml spec unsupported for {name}: found {present_op}{present_version}, required {required_op}{required_version}"
                )
    return issues


def _pip_executable() -> list[str]:
    return [sys.executable, "-m", "pip"]


def _python_install_commands() -> list[list[str]]:
    commands: list[list[str]] = []
    for name, spec in _REQUIRED_PYTHON.items():
        if spec["op"] == "==":
            requirement = f"{name}=={spec['version']}"
        else:
            requirement = f"{name}>={spec['version']}"
        commands.append(_pip_executable() + ["install", requirement])
    return commands


def _binary_install_hints(repo_root: Path) -> list[str]:
    hints: list[str] = []
    vina_path = _vina_binary_path(repo_root)
    if not vina_path.exists():
        hints.append("Run: python setup_env.py (downloads AutoDock Vina into tools/)")
    hints.append("Install Miniconda/Anaconda and ensure `conda` is on PATH.")
    hints.append(
        f"Run: conda create -y -n {_CONDA_ENV_NAME} python={_CONDA_PYTHON_VERSION} -c conda-forge openbabel pdbfixer"
    )
    hints.append("Set OBABEL_EXE/PDBFIXER_EXE or use `conda run -n autoscan-tools`.")
    return hints


def _conda_executable() -> str | None:
    conda_exe = os.environ.get("CONDA_EXE")
    if conda_exe and (Path(conda_exe).exists() or shutil.which(conda_exe)):
        return conda_exe
    resolved = shutil.which("conda")
    return resolved


def _conda_run(command: list[str]) -> subprocess.CompletedProcess[str]:
    conda_exe = _conda_executable()
    if not conda_exe:
        raise FileNotFoundError("Conda not available")
    cmd = [conda_exe, "run", "-n", _CONDA_ENV_NAME, *command]
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def _conda_package_version(package: str) -> str | None:
    conda_exe = _conda_executable()
    if not conda_exe:
        return None
    cmd = [conda_exe, "list", "-n", _CONDA_ENV_NAME, package, "--json"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        return None
    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError:
        return None
    for entry in data:
        if entry.get("name") == package and entry.get("version"):
            return entry["version"]
    return None


def _check_conda() -> list[str]:
    issues: list[str] = []
    if not _conda_executable():
        issues.append("Conda not found. Install Miniconda/Anaconda and add conda to PATH.")
    return issues


def _install_with_conda(packages: list[str]) -> None:
    conda_exe = _conda_executable()
    if not conda_exe:
        return
    create_cmd = [
        conda_exe,
        "create",
        "-y",
        "-n",
        _CONDA_ENV_NAME,
        f"python={_CONDA_PYTHON_VERSION}",
        "-c",
        "conda-forge",
        *packages,
    ]
    subprocess.run(create_cmd, check=False)


def _vina_binary_path(repo_root: Path) -> Path:
    if sys.platform.startswith("win"):
        return repo_root / "tools" / "vina.exe"
    return repo_root / "tools" / "vina"


def _detect_version_from_output(output: str) -> str | None:
    match = re.search(r"(\d+\.\d+\.\d+)", output)
    if match:
        return match.group(1)
    return None


def _detect_command_version(command: str) -> str | None:
    for args in ([command, "--version"], [command, "--help"]):
        try:
            result = subprocess.run(
                args,
                check=True,
                capture_output=True,
                text=True,
                timeout=5,
            )
        except Exception:
            continue
        output = (result.stdout or "") + (result.stderr or "")
        version = _detect_version_from_output(output)
        if version:
            return version
    return None


def _check_vina(repo_root: Path) -> list[str]:
    issues: list[str] = []
    vina_path = _vina_binary_path(repo_root)
    if not vina_path.exists():
        issues.append(
            f"Missing Vina binary in repo: {vina_path}. Run setup_env.py to install."
        )
        return issues

    version = _detect_command_version(str(vina_path))
    required = _REQUIRED_BINARIES["vina"]
    if version is None:
        issues.append("Could not detect Vina version from tools binary.")
    elif not _compare_versions(version, required["version"], required["op"]):
        issues.append(
            f"Vina version mismatch: found {version}, required {required['op']}{required['version']}."
        )
    return issues


def _check_obabel() -> list[str]:
    issues: list[str] = []
    obabel_exe = os.environ.get("OBABEL_EXE", "obabel")
    if not (Path(obabel_exe).exists() or shutil.which(obabel_exe)):
        try:
            version = _conda_package_version("openbabel")
            if not version:
                result = _conda_run(["obabel", "-V"])
                output = (result.stdout or "") + (result.stderr or "")
                version = _detect_version_from_output(output)
            if not version:
                result = _conda_run(["obabel", "--version"])
                output = (result.stdout or "") + (result.stderr or "")
                version = _detect_version_from_output(output)
        except FileNotFoundError:
            version = None
        if version is None:
            issues.append("OpenBabel executable not found. Set OBABEL_EXE or add to PATH.")
            return issues
    else:
        version = _detect_command_version(obabel_exe)

    required = _REQUIRED_BINARIES["obabel"]
    if version is None:
        issues.append("Could not detect OpenBabel version.")
    elif not _compare_versions(version, required["version"], required["op"]):
        issues.append(
            f"OpenBabel version mismatch: found {version}, required {required['op']}{required['version']}."
        )
    return issues


def _check_pdbfixer() -> list[str]:
    issues: list[str] = []
    pdbfixer_exe = os.environ.get("PDBFIXER_EXE", "pdbfixer")
    if not (Path(pdbfixer_exe).exists() or shutil.which(pdbfixer_exe)):
        try:
            version = _conda_package_version("pdbfixer")
            if not version:
                result = _conda_run(
                    ["python", "-c", "import pdbfixer; print(pdbfixer.__version__)"]
                )
                output = (result.stdout or "") + (result.stderr or "")
                version = _detect_version_from_output(output)
        except FileNotFoundError:
            version = None
        if version is None:
            issues.append("PDBFixer executable not found. Install pdbfixer or add to PATH.")
            return issues
    else:
        version = None
        try:
            version = metadata.version("pdbfixer")
        except metadata.PackageNotFoundError:
            version = _detect_command_version(pdbfixer_exe)

    required = _REQUIRED_BINARIES["pdbfixer"]
    if version is None:
        issues.append("Could not detect PDBFixer version.")
    elif not _compare_versions(version, required["version"], required["op"]):
        issues.append(
            f"PDBFixer version mismatch: found {version}, required {required['op']}{required['version']}."
        )
    return issues


def ensure_dependencies() -> None:
    repo_root = find_repo_root(Path(__file__).resolve())
    pyproject_path = repo_root / "pyproject.toml"

    issues: list[str] = []
    issues.extend(_check_conda())
    issues.extend(_check_manifest_in_pyproject(pyproject_path))
    issues.extend(_check_python_dependencies())

    issues.extend(_check_vina(repo_root))
    issues.extend(_check_obabel())
    issues.extend(_check_pdbfixer())

    if issues:
        details = "\n".join(f"- {issue}" for issue in issues)
        raise RuntimeError(f"Dependency check failed:\n{details}")


def build_dependencies(*, install_python: bool = True, quiet: bool = False) -> None:
    repo_root = find_repo_root(Path(__file__).resolve())
    pyproject_path = repo_root / "pyproject.toml"

    issues: list[str] = []
    issues.extend(_check_conda())
    issues.extend(_check_manifest_in_pyproject(pyproject_path))
    issues.extend(_check_python_dependencies())
    issues.extend(_check_vina(repo_root))
    issues.extend(_check_obabel())
    issues.extend(_check_pdbfixer())

    if not issues:
        if not quiet:
            print("Dependency check: PASS")
        return

    if install_python:
        for cmd in _python_install_commands():
            subprocess.run(cmd, check=False)

    _install_with_conda(["openbabel", "pdbfixer"])

    hints = _binary_install_hints(repo_root)
    details = "\n".join(f"- {issue}" for issue in issues)
    hint_text = "\n".join(f"- {hint}" for hint in hints)
    raise RuntimeError(
        "Dependency check failed after build attempt:\n"
        f"{details}\n"
        "Remediation hints:\n"
        f"{hint_text}"
    )
