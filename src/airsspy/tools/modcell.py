"""
Tools for modifying CASTEP cell files.

Provides functions to replace blocks in cell files and update lattice
and position blocks using ASE Atoms objects.
"""

import re
from pathlib import Path

import ase


def replace_block(
    lines: list[str],
    block_name: str,
    block_pattern: str,
    new_value: list[str],
    check: bool = True,
) -> list[str]:
    """
    Replace a block in a list of cell file lines with new values.

    Args:
        lines: Original lines of the cell file.
        block_name: Name for the replacement block header.
        block_pattern: Regex pattern to match the block name in the file.
        new_value: Lines to insert as the block content.
        check: If True, raise RuntimeError when the block is not found.

    Returns:
        A new list of lines with the block replaced.
    """
    in_block = False
    out_from_block = False
    matches_found = 0
    new_lines: list[str] = []
    for line in lines:
        if re.search(r"%BLOCK " + f"{block_pattern.upper()}", line.upper()):
            in_block = True
            matches_found += 1
            if matches_found > 1:
                raise RuntimeError(
                    f"Found multiple blocks matching '{block_pattern}'. "
                    f"Cell file may have duplicate lattice or positions blocks."
                )
            new_lines.append("%BLOCK " + block_name.lower())
            new_lines.extend(new_value)
            new_lines.append("%ENDBLOCK " + block_name.lower())
            continue
        if (
            re.search(r"%ENDBLOCK " + f"{block_pattern.upper()}", line.upper())
            and in_block
        ):
            out_from_block = True
            in_block = False
            continue
        if not in_block:
            new_lines.append(line)
    if matches_found == 0 and check:
        raise RuntimeError(f"Did not find start of the block {block_name}")
    if out_from_block is False and in_block is True:
        raise RuntimeError(f"Did not find end of the block {block_name}")
    return new_lines


def modify_cell(base_cell: str, atoms: ase.Atoms) -> list[str]:
    """
    Modify a cell file using the structure given by an ASE Atoms object.

    Replaces the ``LATTICE_CART`` and ``POSITIONS_ABS`` blocks in the
    base cell file with data from ``atoms``, leaving all other content
    unchanged.

    Args:
        base_cell: Path to the base .cell file.
        atoms: ASE Atoms object with the new structure.

    Returns:
        A list of lines for the modified cell file.
    """
    cell_lines: list[str] = []
    cell = atoms.cell
    for i in range(3):
        cell_lines.append(f"{cell[i, 0]:.10f} {cell[i, 1]:.10f} {cell[i, 2]:.10f}")

    pos = atoms.positions
    pos_lines: list[str] = []
    for symbol, i in zip(atoms.get_chemical_symbols(), range(pos.shape[0])):
        pos_lines.append(
            f"{symbol}  {pos[i, 0]:.10f} {pos[i, 1]:.10f} {pos[i, 2]:.10f}"
        )

    base_lines = Path(base_cell).read_text().split("\n")

    new_lines = replace_block(
        base_lines, "POSITIONS_ABS", "POSITIONS_(FRAC|ABS)", pos_lines
    )
    new_lines = replace_block(
        new_lines, "LATTICE_CART", "LATTICE_(CART|ABC)", cell_lines
    )

    return new_lines
