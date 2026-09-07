"""
Data collection and analysis utilities for AIRSS search results.

Provides functions for collecting RES file data into DataFrames,
reading `ca` command output, combining similar structures via `cryan`,
and computing minsep ranges from ensembles.
"""

import subprocess
from typing import Optional, Union

import numpy as np
import pandas as pd
from pymatgen.entries.computed_entries import ComputedEntry

from ..restools import RESFile, read_res_atoms


def read_stream(stream) -> tuple[list, list]:
    """
    Read from a stream of RES file contents and return lists of
    TitlInfo and Atoms objects.

    Args:
        stream: Iterable of lines from concatenated RES files.

    Returns:
        Tuple of (titl_list, atoms_list).
    """
    lines = []
    atoms_list = []
    titl_list = []
    for line in stream:
        if line.startswith("END"):
            titl, atoms = read_res_atoms(lines)
            titl_list.append(titl)
            atoms_list.append(atoms)
            lines = []
        else:
            lines.append(line)
    return titl_list, atoms_list


def read_ca(lines: list[str]) -> pd.DataFrame:
    """
    Read results from the ``ca`` command into a DataFrame.

    Args:
        lines: String lines as returned by the ``ca`` command.

    Returns:
        A DataFrame with columns for label, pressure, volume,
        enthalpy, spin info, formula, symmetry, etc.
    """
    records = []
    if not lines:
        return pd.DataFrame()
    # Detect format from first non-empty line
    ntokens = None
    for line in lines:
        if line.strip():
            ntokens = len(line.split())
            break
    if ntokens is None:
        return pd.DataFrame()

    for line in lines:
        if not line:
            continue
        tokens = line.split()
        if len(tokens) != ntokens:
            continue
        # If has spin (10 tokens per line)
        if ntokens == 10:
            records.append(
                {
                    "label": tokens[0],
                    "press": float(tokens[1]),
                    "volume": float(tokens[2]),
                    "H": float(tokens[3]),
                    "spin": float(tokens[4]),
                    "aspin": float(tokens[5]),
                    "nform": int(tokens[6]),
                    "formula": tokens[7],
                    "symm": tokens[8],
                    "nseen": int(tokens[9]),
                }
            )
        else:
            records.append(
                {
                    "label": tokens[0],
                    "press": float(tokens[1]),
                    "volume": float(tokens[2]),
                    "H": float(tokens[3]),
                    "nform": int(tokens[4]),
                    "formula": tokens[5],
                    "symm": tokens[6],
                    "nseen": int(tokens[7]),
                }
            )

    dataframe = pd.DataFrame.from_records(records)

    # Fix the H field: cryan outputs relative enthalpy after the first entry
    if len(dataframe) > 1:
        h_col = dataframe.columns.get_loc("H")
        dataframe.iloc[1:, h_col] += dataframe.iloc[0, h_col]
    return dataframe


def collect_res_in_df(
    res_collection: list[RESFile],
    norm_mode: str = "per_atom",
) -> pd.DataFrame:
    """
    Collect a list of RESFile objects into a DataFrame.

    Args:
        res_collection: A collection of RESFile objects.
        norm_mode: Normalisation mode for energy and volume.
            ``"per_atom"`` (default) or ``"per_formula_unit"``.

    Returns:
        A DataFrame with collected data from the RESFile objects.
    """
    records = []
    for res in res_collection:
        entry = {}
        entry.update(res.data)
        entry["formula"] = res.formula
        entry["reduced_formula"] = res.reduced_formula
        entry["nform"] = res.n_formula_units
        entry["res"] = res
        entry["chemsys"] = res.composition.chemical_system if res.composition else None
        records.append(entry)

    dframe = pd.DataFrame(records)
    if dframe.empty:
        return dframe

    # Normalise the energy and volumes
    if norm_mode == "per_atom":
        dframe["H"] = dframe["enthalpy"] / dframe["natoms"]
        dframe["V"] = dframe["volume"] / dframe["natoms"]
    else:
        # Guard against None nform from structures without loaded data
        if dframe["nform"].isna().any():
            dframe["nform"] = dframe["nform"].fillna(1)
        dframe["H"] = dframe["enthalpy"] / dframe["nform"]
        dframe["V"] = dframe["volume"] / dframe["nform"]
    dframe.sort_values("H", inplace=True)

    return dframe


def combine_res_cryan(
    dframe: pd.DataFrame,
    thres: float = 0.1,
    ntop: int = 30,
) -> pd.DataFrame:
    """
    Reduce similar structures using the ``cryan`` command.

    Args:
        dframe: DataFrame with a ``res`` column containing RESFile objects.
        thres: Threshold for combining structures.
        ntop: The number of top structures to be returned.

    Returns:
        A DataFrame of output from the ``cryan`` command.
    """
    lines = []
    for _, row in dframe.iterrows():
        res = row["res"]
        if res.lines:
            lines.extend(res.lines)
        else:
            lines.extend(res.to_res_lines())

    if lines[0].endswith("\n"):
        join_base = ""
    else:
        join_base = "\n"
    inpd = join_base.join(lines)
    cryan_out = subprocess.check_output(
        ["cryan", "-u", str(thres), "-r", "-t", str(ntop), "-l"],
        text=True,
        input=inpd,
    ).split("\n")
    cadf = read_ca(cryan_out)

    return cadf


def get_minsep_range(
    minseps: list[dict[str, float]],
    cap: Optional[tuple[float, float]] = None,
) -> dict[str, list[float]]:
    """
    Create ranged minseps from an ensemble of minsep entries.

    Args:
        minseps: A list of minsep dictionaries (species pair -> distance).
        cap: Optional (min, max) cap for distances.

    Returns:
        A dictionary mapping species pairs to [min, max] ranges.
    """
    base: dict[str, list[float]] = {
        key: [value, value] for key, value in minseps[0].items()
    }
    for minsep in minseps:
        for key, value in minsep.items():
            if key in base:
                existing = base[key]
                if cap and value < cap[0]:
                    existing[0] = cap[0]
                elif cap and value > cap[1]:
                    existing[1] = cap[1]
                elif existing[0] > value:
                    existing[0] = value
                elif existing[1] < value:
                    existing[1] = value
            else:
                base[key] = [value, value]
    return base


def get_entry(
    dataframe: pd.DataFrame,
    pmg_col: str = "pmg_struct",
    label_col: str = "label",
    uuid_col: str = "uuid",
    umap_col: str = "umap",
    xc_col: str = "functional",
    eng_col: str = "energy",
) -> list[ComputedEntry]:
    """
    Create ComputedEntry objects from a DataFrame containing structure data.

    Args:
        dataframe: DataFrame with structure and energy data.
        pmg_col: Column name for pymatgen Structure objects.
        label_col: Column name for structure labels.
        uuid_col: Column name for UUIDs.
        umap_col: Column name for Hubbard U mapping.
        xc_col: Column name for functional labels.
        eng_col: Column name for energy values.

    Returns:
        A list of ComputedEntry objects.
    """
    pd_entries = []
    for idx, row in dataframe.iterrows():
        comp = row[pmg_col].composition
        attrs = {
            "struct_name": row[label_col],
            "entry_type": "MP" if "mp" in str(row.get(label_col, "")) else "AIRSS",
            "structure_uuid": row.get(uuid_col),
            "calc_u": row.get(umap_col),
            "functional": row.get(xc_col),
            "volume": row[pmg_col].volume,
            "dataframe_idx": idx,
        }
        pd_entries.append(ComputedEntry(comp, energy=row[eng_col], parameters=attrs))
    return pd_entries


def export_dataframe_as_res(
    dataframe: pd.DataFrame,
    comment: str = "VASP export",
    extra_comments: Optional[list] = None,
    stress_key: Optional[str] = None,
) -> None:
    """
    Write all structures in a DataFrame into RES format for export.

    Creates an ``exports/`` directory and writes one ``.res`` file per row.

    Args:
        dataframe: DataFrame with structure and energy data.
            Must have ``pmg_struct_relaxed``, ``energy_per_atom``,
            ``volume_per_fu``, ``nform_refine``, and ``label`` columns.
        comment: Comment to include in REM lines.
        extra_comments: Additional REM comment strings.
        stress_key: Column name for stress data (optional).
    """
    from pathlib import Path

    comments = [comment]
    if extra_comments:
        comments.extend(extra_comments)

    export_dir = Path("exports")
    export_dir.mkdir(exist_ok=True)

    for _, row in dataframe.iterrows():
        relaxed = row["pmg_struct_relaxed"]
        res = RESFile(
            relaxed,
            {
                "enthalpy": row["energy_per_atom"] * len(relaxed.sites),
                "volume": row.get("volume_per_fu", relaxed.volume)
                * row.get("nform_refine", 1),
                "pressure": 0.0
                if stress_key is None or stress_key not in row.index
                else row[stress_key],
                "label": row["label"],
                "rem": comments
                + [f"{key} = {row[key]}" for key in row.index if "struct" not in key],
            },
        )
        content = "\n".join(res.to_res_lines())
        (export_dir / f"{row['label']}.res").write_text(content)


def get_pressure_gpa(stress: Union[list, np.ndarray]) -> float:
    """
    Convert a stress tensor to isostatic pressure in GPa.

    Args:
        stress: A 3x3 stress tensor (in kBar units).

    Returns:
        The isostatic pressure in GPa.
    """
    stress = np.asarray(stress)
    # 1 kBar = 0.1 GPa
    return float(np.trace(stress) * 0.1 / 3.0)
