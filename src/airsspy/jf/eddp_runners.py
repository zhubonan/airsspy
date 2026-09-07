"""Native EDDPotentials.jl runners.

The runners in this module deliberately use a Julia subprocess instead of an
ASE calculator shim.  EDDPotentials therefore owns the optimisation loop and
can use its native calculator workspace, stress handling, and optimisers.  The
Julia bridge writes a temporary RES geometry plus a small JSON result; Python
normalises those results to extxyz so the existing ML result composer can be
reused unchanged.
"""

from __future__ import annotations

import json
import logging
import shlex
import subprocess
from pathlib import Path
from typing import Union, cast

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write as ase_write

logger = logging.getLogger(__name__)

StructureInput = Union[str, Atoms]


def _bridge_path() -> Path:
    """Return the packaged Julia bridge path."""
    return Path(__file__).resolve().parent.parent / "julia" / "eddp_bridge.jl"


def _structure_input_to_atoms(structure_input: StructureInput) -> Atoms:
    """Convert a cell string or copy an ASE structure."""
    if isinstance(structure_input, Atoms):
        return structure_input.copy()

    from .ml_runners import _cell_content_to_atoms

    return _cell_content_to_atoms(structure_input)


def _write_eddp_input(path: Path, atoms: Atoms, label: str) -> None:
    """Write one structure in the RES format consumed by CellBase."""
    from ..convert import _atoms_to_res_lines

    atoms = atoms.copy()
    atoms.calc = None
    atoms.info["label"] = label
    path.write_text("\n".join(_atoms_to_res_lines(atoms)) + "\n")


def _julia_matrix_to_numpy(values, *, transpose: bool = False) -> np.ndarray:
    """Convert a JSON-encoded Julia matrix to a NumPy array."""
    array = np.asarray(values, dtype=float)
    return cast(np.ndarray, array.T if transpose else array)


class _AirssEddpBaseRunner:
    """Shared Julia invocation and result normalisation for EDDP runners."""

    _cleanup_extensions = [
        ".cell",
        ".extxyz",
        ".traj",
        ".res",
        ".err",
        ".eddp-input.res",
        ".eddp-output.res",
        ".eddp-result.json",
        ".eddp.log",
        "-orig.cell",
    ]
    mode = "singlepoint"

    def __init__(
        self,
        model_path: str | Path,
        executable: str = "julia",
        project: str | Path | None = None,
        pressure: float = 0.0,
        timeout: float | None = None,
    ) -> None:
        self.model_path = Path(model_path).expanduser().resolve()
        self.executable = executable
        self.project = (
            Path(project).expanduser().resolve() if project is not None else None
        )
        self.pressure = float(pressure)
        self.timeout = timeout
        self.last_result: dict | None = None
        self.last_command: list[str] | None = None

    def clean_failed(self, struct_name: str) -> None:
        from .runners import clean_files

        clean_files(struct_name, self._cleanup_extensions)

    def _command(
        self,
        struct_name: str,
        input_path: Path,
        output_path: Path,
        result_path: Path,
    ) -> list[str]:
        command = shlex.split(self.executable)
        if self.project is not None:
            command.append(f"--project={self.project}")
        command.extend(
            [
                str(_bridge_path()),
                self.mode,
                str(self.model_path),
                str(input_path),
                str(output_path),
                str(result_path),
                struct_name,
                self._method(),
                str(self._max_steps()),
                str(self._force_tolerance()),
                str(self._stress_tolerance()),
                str(self.pressure),
                str(self._relax_cell()).lower(),
            ]
        )
        return command

    def _method(self) -> str:
        return "tpsd"

    def _max_steps(self) -> int:
        return 0

    def _force_tolerance(self) -> float:
        return 0.0

    def _stress_tolerance(self) -> float:
        return 0.0

    def _relax_cell(self) -> bool:
        return False

    def _normalise_result(
        self,
        struct_name: str,
        output_path: Path,
        result_path: Path,
    ) -> int:
        from ..restools import RESFile

        result = json.loads(result_path.read_text())
        res = RESFile.from_file(str(output_path), include_structure=True)
        atoms = res.atoms
        if atoms is None:
            raise RuntimeError("EDDP output did not contain a structure")

        energy = float(result["energy"])
        forces = _julia_matrix_to_numpy(result["forces"], transpose=True)
        if forces.shape != (len(atoms), 3):
            raise RuntimeError(
                "EDDP force shape mismatch: "
                f"expected {(len(atoms), 3)}, got {forces.shape}"
            )

        calc_results: dict = {"energy": energy, "forces": forces}
        stress_values = result.get("stress")
        if stress_values is not None:
            stress = _julia_matrix_to_numpy(stress_values)
            if stress.shape == (3, 3):
                calc_results["stress"] = stress

        atoms.calc = SinglePointCalculator(atoms, **calc_results)
        atoms.info["label"] = struct_name
        atoms.info["extern_pressure"] = self.pressure
        atoms.info["eddp_model"] = str(self.model_path)
        atoms.info["eddp_pressure_gpa"] = float(result.get("pressure_gpa", 0.0))
        converged = bool(result.get("converged", True))
        if self.mode == "relax":
            atoms.info["relax_converged"] = converged
            atoms.info["relax_status"] = "converged" if converged else "max_steps"
            atoms.info["relax_steps"] = int(result.get("iterations", 0))
            atoms.info["relax_fmax"] = float(result.get("fmax", 0.0))
            atoms.info["relax_smax_gpa"] = float(result.get("smax_gpa", 0.0))

        ase_write(struct_name + ".extxyz", atoms, format="extxyz")
        self.last_result = result
        return 0 if converged else 1

    def run(self, struct_name: str, structure_input: StructureInput) -> int:
        """Run EDDPotentials.jl and write ``<struct_name>.extxyz``."""
        input_path = Path(struct_name + ".eddp-input.res")
        output_path = Path(struct_name + ".eddp-output.res")
        result_path = Path(struct_name + ".eddp-result.json")
        log_path = Path(struct_name + ".eddp.log")
        extxyz_path = Path(struct_name + ".extxyz")

        try:
            atoms = _structure_input_to_atoms(structure_input)
            _write_eddp_input(input_path, atoms, struct_name)
            extxyz_path.unlink(missing_ok=True)
            output_path.unlink(missing_ok=True)
            result_path.unlink(missing_ok=True)

            command = self._command(struct_name, input_path, output_path, result_path)
            self.last_command = command
            completed = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=self.timeout,
                check=False,
            )
            log_path.write_text(
                "STDOUT\n" + completed.stdout + "\nSTDERR\n" + completed.stderr
            )
            if completed.returncode != 0:
                logger.error(
                    "EDDP %s failed for %s with return code %d; see %s",
                    self.mode,
                    struct_name,
                    completed.returncode,
                    log_path,
                )
                return 1
            if not output_path.is_file() or not result_path.is_file():
                logger.error("EDDP did not produce complete output for %s", struct_name)
                return 1
            return self._normalise_result(struct_name, output_path, result_path)
        except Exception:
            logger.error("EDDP %s failed for %s", self.mode, struct_name, exc_info=True)
            return 1


class AirssEddpSinglePointRunner(_AirssEddpBaseRunner):
    """Run native EDDP energy, force, and stress evaluation."""

    mode = "singlepoint"


class AirssEddpRelaxRunner(_AirssEddpBaseRunner):
    """Relax a structure with EDDPotentials.jl native optimisers."""

    mode = "relax"

    def __init__(
        self,
        model_path: str | Path,
        executable: str = "julia",
        project: str | Path | None = None,
        method: str = "tpsd",
        max_steps: int = 500,
        force_tolerance: float = 0.05,
        stress_tolerance_gpa: float = 0.1,
        pressure: float = 0.0,
        relax_cell: bool = True,
        timeout: float | None = None,
    ) -> None:
        super().__init__(
            model_path=model_path,
            executable=executable,
            project=project,
            pressure=pressure,
            timeout=timeout,
        )
        method = method.lower()
        if method not in {"tpsd", "fire"}:
            raise ValueError("EDDP method must be 'tpsd' or 'fire'")
        if max_steps < 1:
            raise ValueError("EDDP max_steps must be >= 1")
        if force_tolerance <= 0.0:
            raise ValueError("EDDP force_tolerance must be > 0")
        if stress_tolerance_gpa <= 0.0:
            raise ValueError("EDDP stress_tolerance_gpa must be > 0")
        self.method = method
        self.max_steps = int(max_steps)
        self.force_tolerance = float(force_tolerance)
        self.stress_tolerance_gpa = float(stress_tolerance_gpa)
        self.relax_cell = bool(relax_cell)

    def _method(self) -> str:
        return self.method

    def _max_steps(self) -> int:
        return self.max_steps

    def _force_tolerance(self) -> float:
        return self.force_tolerance

    def _stress_tolerance(self) -> float:
        return self.stress_tolerance_gpa

    def _relax_cell(self) -> bool:
        return self.relax_cell


def compose_eddp_task_doc(struct_name: str, model_path: str = "") -> dict:
    """Compose an AIRSS result through the common extxyz result pipeline."""
    from .ml_runners import compose_ml_task_doc

    return cast(
        dict,
        compose_ml_task_doc(
            struct_name,
            calculator_spec=model_path,
            calculator_label="EDDP Model",
            metadata_label="EDDP",
        ),
    )
