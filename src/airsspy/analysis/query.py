"""
Bridge between the jobflow query layer and analysis DataFrames.

Provides functions that collect results from a ``SearchStore`` into
DataFrames compatible with the analysis utilities in ``collect.py``.
"""

import pandas as pd

from ..jf.store import SearchStore


def collect_results_df(
    store: SearchStore,
    project_name: str,
    norm_mode: str = "per_atom",
) -> pd.DataFrame:
    """
    Collect results from the jobflow store into an analysis-ready DataFrame.

    Args:
        store: A connected ``SearchStore`` instance.
        project_name: Project identifier to query.
        norm_mode: Normalisation mode for energy and volume.
            ``"per_atom"`` (default) or ``"per_formula_unit"``.

    Returns:
        A DataFrame sorted by normalised enthalpy.
    """
    results = store.retrieve_project(project_name)
    if not results:
        return pd.DataFrame()

    records = []
    for r in results:
        entry = r.model_dump()

        # Compute derived fields
        if r.structure is not None:
            entry["chemsys"] = r.structure.composition.chemical_system
            comp = r.structure.composition
            nform = comp.get_reduced_formula_and_factor()[1]
            entry["nform"] = nform
        else:
            entry["chemsys"] = None
            entry["nform"] = 1

        records.append(entry)

    dframe = pd.DataFrame(records)
    if dframe.empty:
        return dframe

    # Normalise energy and volume
    if norm_mode == "per_atom":
        dframe["H"] = dframe["energy"] / dframe["natoms"]
        dframe["V"] = dframe["volume"] / dframe["natoms"]
    else:
        dframe["H"] = dframe["energy"] / dframe["nform"]
        dframe["V"] = dframe["volume"] / dframe["nform"]

    dframe.sort_values("H", inplace=True)
    return dframe
