"""Docking utilities for defining search space geometry."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
from Bio.PDB.Residue import Residue

from autoscan.utils import get_logger

logger = get_logger(__name__)


@dataclass(frozen=True)
class GridBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float


def _iter_coords(residues: Iterable[Residue]) -> list[list[float]]:
    coords: list[list[float]] = []
    for residue in residues:
        for atom in residue.get_atoms():
            coords.append(list(atom.get_coord()))
    return coords


def calculate_grid_box(
    residues: Sequence[Residue],
    buffer_angstrom: float = 6.0,
    max_volume: float = 30000.0,
) -> GridBox:
    """
    Calculate a docking grid box from active-site residues.

    Args:
        residues: List of Bio.PDB residues defining the active site.
        buffer_angstrom: Padding added to each dimension (in A).
        max_volume: Upper bound for grid volume (in A^3).

    Returns:
        GridBox with center and size values.

    Raises:
        ValueError: If residues are empty, buffer is negative, or volume is invalid.
    """
    if buffer_angstrom < 0:
        raise ValueError("buffer_angstrom must be non-negative")
    if not residues:
        raise ValueError("residues must be a non-empty sequence")

    coords = _iter_coords(residues)
    if not coords:
        raise ValueError("No atom coordinates found for the provided residues")

    coord_array = np.array(coords, dtype=float)
    min_xyz = coord_array.min(axis=0)
    max_xyz = coord_array.max(axis=0)
    center = coord_array.mean(axis=0)

    sizes = (max_xyz - min_xyz) + (2.0 * buffer_angstrom)
    if np.any(sizes <= 0):
        raise ValueError("Computed grid box has non-positive dimension")

    volume = float(sizes[0] * sizes[1] * sizes[2])
    logger.info(
        "Grid box center=(%.3f, %.3f, %.3f) size=(%.3f, %.3f, %.3f) volume=%.1f A^3",
        center[0],
        center[1],
        center[2],
        sizes[0],
        sizes[1],
        sizes[2],
        volume,
    )

    if volume > max_volume:
        raise ValueError(
            "Grid box volume too large: %.1f A^3 exceeds limit %.1f A^3"
            % (volume, max_volume)
        )

    return GridBox(
        center_x=float(center[0]),
        center_y=float(center[1]),
        center_z=float(center[2]),
        size_x=float(sizes[0]),
        size_y=float(sizes[1]),
        size_z=float(sizes[2]),
    )
