"""
Analysis tools for AIRSS search results.

Provides data collection, DataFrame conversion, convex hull analysis,
and utilities for working with collections of RES files.
"""

from .collect import (
    collect_res_in_df,
    combine_res_cryan,
    export_dataframe_as_res,
    get_entry,
    get_minsep_range,
    get_pressure_gpa,
    read_ca,
    read_stream,
)
from .hull import PlotlyPDPlotter, entry_name, entry_type, make_axis
from .query import collect_results_df

__all__ = [
    "collect_res_in_df",
    "collect_results_df",
    "combine_res_cryan",
    "export_dataframe_as_res",
    "get_entry",
    "get_minsep_range",
    "get_pressure_gpa",
    "PlotlyPDPlotter",
    "read_ca",
    "read_stream",
    "entry_name",
    "entry_type",
    "make_axis",
]
