import numpy as np


def calculate_grid_box(center_coords: list, ligand_mol=None, buffer_angstroms=6.0):
    """
    Mathematically robust grid generation.

    Args:
        center_coords (list): [x, y, z] center of the pocket.
        ligand_mol (OBMol, optional): The ligand object to size the box.
        buffer_angstroms (float): Physics padding (default 6.0 A).

    Returns:
        dict: {'center_x': ..., 'size_x': ..., ...}
    """
    # 1. Default Size (Small molecule avg) if no ligand provided
    size = [20.0, 20.0, 20.0]

    # 2. Dynamic Sizing (Physics-based)
    if ligand_mol:
        # Get min/max bounds of the ligand
        coords = [atom.coords for atom in ligand_mol.atoms]
        min_c = np.min(coords, axis=0)
        max_c = np.max(coords, axis=0)

        # Box size = (Max - Min) + Buffer
        size = (max_c - min_c) + buffer_angstroms

        # Sanity Check: Don't let box get massive (Entropy penalty)
        size = np.clip(size, 10.0, 35.0)

    return {
        'center_x': center_coords[0],
        'center_y': center_coords[1],
        'center_z': center_coords[2],
        'size_x': float(size[0]),
        'size_y': float(size[1]),
        'size_z': float(size[2])
    }
