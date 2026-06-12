"""
CLI commands for miscellaneous tools.
"""

import click


@click.group("tools")
@click.pass_context
def tools(ctx):
    """Collection of utility tools."""
    _ = ctx


@tools.command("modcell")
@click.argument("base_cell")
@click.argument("other_cell")
def modcell(base_cell, other_cell):
    """Modify the structure of a CELL file using another.

    BASE_CELL is the template cell file. OTHER_CELL is the file
    providing the new structure (any ASE-supported format).
    """
    from ase.io import read

    from airsspy.tools.modcell import modify_cell

    atoms = read(other_cell)
    lines = modify_cell(base_cell, atoms)
    click.echo("\n".join(lines))
