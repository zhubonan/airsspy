"""
Output document models for AIRSS jobflow jobs.

Defines pydantic models following the atomate2 TaskDoc pattern.
One ``AirssJobDoc`` contains N ``AirssResultDoc`` entries, supporting
multi-structure jobs (e.g., 50 random structures per search job).
"""

from datetime import datetime
from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field
from pymatgen.core import Structure


class RelaxOutcome(str, Enum):
    """Outcome of a single structure relaxation."""

    FINISHED = "finished"
    TIMEDOUT = "timedout"
    ERRORED = "errored"
    FAILED = "failed"
    UNDETERMINED = "undetermined"
    CYCLE_EXCEEDED = "cycle_exceeded"


class AirssResultDoc(BaseModel):
    """Single relaxed structure result from an AIRSS job."""

    struct_name: str
    seed_name: str
    project_name: str

    # Structures
    structure: Optional[Structure] = None
    initial_structure: Optional[Structure] = None

    # Computed properties
    energy: Optional[float] = None
    energy_per_atom: Optional[float] = None
    volume: Optional[float] = None
    pressure: Optional[float] = None
    spin: float = 0.0
    mod_spin: float = 0.0
    symmetry: Optional[str] = None
    formula: Optional[str] = None
    reduced_formula: Optional[str] = None
    natoms: Optional[int] = None

    # Raw RES content
    res_content: Optional[str] = None

    # Computational metadata
    parallel_efficiency: Optional[float] = None
    total_time: Optional[float] = None
    relax_status: RelaxOutcome = RelaxOutcome.FINISHED
    error_message: Optional[str] = None

    # REM block content (raw lines for .res reconstruction)
    rem_lines: Optional[list[str]] = None


class AirssJobDoc(BaseModel):
    """Output document from any AIRSS job (search or relax).

    Contains N results from a single job execution. For
    ``AirssSearchMaker``, this holds N randomly-generated-then-relaxed
    structures. For ``AirssRelaxMaker``, this holds N
    pre-existing-now-relaxed structures.
    """

    project_name: str
    seed_name: str
    job_type: str = "search"  # "search", "relax", "singlepoint"

    # Seed info (present for search jobs)
    seed_content: Optional[str] = None
    seed_hash: Optional[str] = None
    param_content: Optional[str] = None

    results: list[AirssResultDoc] = Field(default_factory=list)

    # Summary
    n_structures: int = 0
    n_finished: int = 0
    n_errored: int = 0
    n_failed: int = 0
    n_timedout: int = 0

    created_on: Optional[datetime] = None
