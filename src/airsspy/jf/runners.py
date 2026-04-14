"""
Execution runners for AIRSS calculations.

Pure computation classes with no jobflow dependency. Each runner handles
one buildcell invocation or one CASTEP relaxation cycle. They are usable
standalone or within jobflow Makers.
"""

import logging
import re
import shutil
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def run_buildcell(
    seed_name: str,
    seed_content: str,
    build_timeout: int = 30,
    write_seed: bool = True,
) -> Optional[dict[str, str]]:
    """
    Run the buildcell executable to generate a random structure.

    Args:
        seed_name: Name of the seed (without extension).
        seed_content: Content of the seed .cell file.
        build_timeout: Timeout in seconds for each buildcell attempt.
        write_seed: Whether to write the seed .cell file to disk.

    Returns:
        Dictionary with ``struct_name``, ``seed_name``, ``seed_hash``,
        ``struct_content`` keys, or None if all attempts timed out.
    """
    from ..casteptools import get_rand_cell_name

    logger.info("Starting random structure generation...")
    attempt = 3
    stdout: Optional[str] = None
    while attempt > 0:
        try:
            proc = subprocess.Popen(
                "buildcell",
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
            out, _ = proc.communicate(seed_content, timeout=build_timeout)
        except subprocess.TimeoutExpired:
            attempt -= 1
            proc.kill()
        else:
            stdout = out
            break

    if attempt <= 0:
        logger.error("Random structure generation timed out")
        return None

    logger.info("Random structure generation completed")
    cell_name = get_rand_cell_name(seed_name)
    struct_name = cell_name.replace(".cell", "")

    if write_seed:
        Path(seed_name + ".cell").write_text(seed_content)

    Path(cell_name).write_text(stdout)
    Path(struct_name + "-orig.cell").write_text(stdout)

    return {
        "struct_name": struct_name,
        "seed_name": seed_name,
        "struct_content": stdout,
    }


class AirssCastepSinglePointRunner:
    """Execute a CASTEP single-point calculation."""

    def __init__(self, executable: str = "castep.mpi") -> None:
        self.executable = executable

    def prepare_inputs(self, struct_name: str, cellinput, paraminput) -> None:
        """Write .cell and .param files to disk.

        Args:
            struct_name: Seed name (without extension).
            cellinput: A CastepInput/CellInput instance or string.
            paraminput: A CastepInput/ParamInput instance or string.
        """
        from castepinput.inputs import CastepInput

        cell_name = struct_name + ".cell"
        param_name = struct_name + ".param"

        if isinstance(cellinput, CastepInput):
            cellinput = cellinput.get_string()
        if isinstance(paraminput, CastepInput):
            paraminput = paraminput.get_string()

        Path(cell_name).write_text(cellinput)
        Path(param_name).write_text(paraminput)

    def run(self, struct_name: str, cellinput, paraminput) -> int:
        """
        Run a single-point calculation.

        Args:
            struct_name: Seed name (without extension).
            cellinput: Cell input (CellInput instance or string).
            paraminput: Param input (ParamInput instance or string).

        Returns:
            0 on success, non-zero on failure.
        """
        paraminput["task"] = "singlepoint"
        self.prepare_inputs(struct_name, cellinput, paraminput)
        output = subprocess.run(self.executable.split() + [struct_name], check=False)
        return output.returncode


class AirssCastepRelaxRunner(AirssCastepSinglePointRunner):
    """
    Execute a cyclic CASTEP geometry optimisation.

    Runs up to ``cycles`` successive relaxations, copying the output
    cell back to the input between cycles. Requires two consecutive
    successful runs to declare convergence.
    """

    def __init__(
        self,
        executable: str = "castep.mpi",
        cycles: int = 4,
        max_fails: int = 2,
        max_iterations: int = 200,
    ) -> None:
        super().__init__(executable=executable)
        self.cycles = cycles
        self.max_fails = max_fails
        self.max_iterations = max_iterations

    def run(self, struct_name: str, cellinput, paraminput) -> int:
        """
        Run cyclic CASTEP relaxation.

        Args:
            struct_name: Seed name (without extension).
            cellinput: Cell input (CellInput instance or string).
            paraminput: Param input (ParamInput instance or string).

        Returns:
            0 if converged, 1 if not converged or failed.
        """
        from castepinput.inputs import CellInput

        paraminput["task"] = "geometryoptimization"
        paraminput["write_cell_structure"] = True

        self.prepare_inputs(struct_name, cellinput, paraminput)
        fail_counter = 0
        cycle = 1
        success_counter = 0
        iter_counter = 0

        while cycle <= self.cycles:
            if fail_counter > self.max_fails:
                return 1
            output = subprocess.run(
                self.executable.split() + [struct_name], check=False
            )
            if output.returncode != 0:
                fail_counter += 1
                continue
            fail_counter = 0

            result = None
            max_iter = 0
            with open(struct_name + ".castep") as fhandle:
                for line in fhandle:
                    match = re.search(r"Geometry optimization ([a-z]+)", line)
                    if match is not None:
                        if match.group(1) == "completed":
                            result = True
                        elif match.group(1) == "failed":
                            result = False
                    match = re.search(r"Finished iteration +(\d+)", line)
                    if match is not None:
                        max_iter = int(match.group(1))

            iter_counter += max_iter
            if result is True:
                success_counter += 1
            if result is False:
                success_counter = 0

            if success_counter >= 2:
                break
            if iter_counter >= self.max_iterations:
                break

            out_cell = CellInput.from_file(struct_name + "-out.cell")
            in_cell = CellInput.from_file(struct_name + ".cell")
            in_cell.set_cell(out_cell.get_cell())
            in_cell.set_positions(*out_cell.get_positions())
            in_cell.save(struct_name + ".cell")
            cycle += 1

        if success_counter >= 2:
            return 0
        return 1


def compose_task_doc(struct_name: str) -> dict:
    """
    Extract results from a completed CASTEP calculation.

    Reads the .castep and .cell files, computes derived properties,
    writes a ``<struct_name>.res`` file to disk, and returns a dictionary suitable for
    constructing an ``AirssResultDoc``.

    Args:
        struct_name: Seed name (without extension).

    Returns:
        Dictionary with energy, structure, volume, formula, etc.
    """
    from ase import Atoms
    from castepinput.inputs import CellInput
    from pymatgen.io.ase import AseAtomsAdaptor

    from ..restools import save_airss_res

    energy = None
    pressure = None
    efficiency = None
    spin = 0.0
    modspin = 0.0
    spin_moms: list[float] = []
    in_spin_group = False
    total_time = None

    castep_file = struct_name + ".castep"
    if Path(castep_file).is_file():
        with open(castep_file) as fhandle:
            for line in fhandle:
                if "NB est. 0K energy" in line:
                    energy = float(line.split()[-2])
                if "Pressure: " in line:
                    pressure = float(line.split()[-2])
                if "Overall parallel efficiency" in line:
                    match = re.search(r"(\d+)%", line)
                    if match:
                        efficiency = float(match.group(1)) / 100.0
                if "Total time" in line and "=" in line:
                    match = re.search(r"=\s*([0-9.]+)", line)
                    if match:
                        total_time = float(match.group(1))
                if "Spin den" in line:
                    spin = float(line.split()[-2])
                if "|Spin den" in line:
                    modspin = float(line.split()[-2])
                if " Total  Charge(e)   Spin(hbar/2)" in line:
                    in_spin_group = True
                    spin_moms = []
                if " Length (A)" in line:
                    in_spin_group = False
                if in_spin_group:
                    tokens = line.split()
                    if re.match(r"^ +[A-Za-z]+ ", line):
                        spin_moms.append(float(tokens[-1]))

    if Path(struct_name + "-out.cell").is_file():
        cell = CellInput.from_file(struct_name + "-out.cell")
    else:
        cell = CellInput.from_file(struct_name + ".cell")
    elements, positions, _tags = cell.get_positions()
    atoms = Atoms(symbols=elements, positions=positions, cell=cell.get_cell(), pbc=True)

    info = {"uid": struct_name, "H": energy}
    save_airss_res(atoms, info, fname=struct_name + ".res", force_write=True)
    structure = AseAtomsAdaptor.get_structure(atoms)
    if spin_moms:
        structure.add_site_property("spin", spin_moms)

    return {
        "structure": structure,
        "volume": structure.volume,
        "reduced_formula": structure.reduced_formula,
        "formula": structure.composition.formula.replace(" ", ""),
        "natoms": len(atoms),
        "label": struct_name,
        "energy": energy,
        "energy_per_atom": energy / len(atoms) if energy else None,
        "spin": spin,
        "mod_spin": modspin,
        "pressure": pressure,
        "parallel_efficiency": efficiency,
        "total_time": total_time,
        "res_content": Path(struct_name + ".res").read_text()
        if Path(struct_name + ".res").is_file()
        else None,
    }


class AirssScriptRelaxRunner:
    """
    Base runner for external AIRSS relaxation scripts (gulp_relax, pp3_relax).

    Calls an external script and checks the output file for success.
    Subclasses must override ``_get_cmd`` and set ``_param_suffix``.
    """

    _param_suffix: str = ".param"

    def __init__(
        self,
        executable: str = "gulp",
        timeout: int = 600,
        max_attempts: int = 3,
    ) -> None:
        self.executable = executable
        self.timeout = timeout
        self.max_attempts = max_attempts

    def _get_cmd(self, struct_name: str) -> list[str]:
        """Construct the shell command. Must be overridden by subclasses."""
        raise NotImplementedError

    def _check_success(self, struct_name: str, stdout: str) -> bool:
        """Check if the relaxation finished successfully."""
        from ..casteptools import gulp_relax_finish_ok

        return gulp_relax_finish_ok(struct_name + ".castep")

    def _prepare_inputs(
        self,
        struct_name: str,
        struct_content: str,
        param_content: str,
        seed_name: Optional[str] = None,
    ) -> None:
        """Write .cell and code-specific param files to disk."""
        Path(struct_name + ".cell").write_text(struct_content)
        Path(struct_name + self._param_suffix).write_text(param_content)

    def run(
        self,
        struct_name: str,
        struct_content: str,
        param_content: str,
        seed_name: Optional[str] = None,
    ) -> int:
        """
        Run relaxation via the external script.

        Args:
            struct_name: Structure name (without extension).
            struct_content: Content of the .cell file.
            param_content: Content of the code-specific param file.
            seed_name: Seed name (needed by GULP for .lib file rename).

        Returns:
            0 on success, 1 on failure.
        """
        self._prepare_inputs(struct_name, struct_content, param_content, seed_name)
        cmd = self._get_cmd(struct_name)

        attempt = 0
        while attempt < self.max_attempts:
            attempt += 1
            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=self.timeout,
                    check=False,
                )
            except subprocess.TimeoutExpired:
                logger.warning(
                    "Relaxation attempt %d/%d timed out for %s",
                    attempt,
                    self.max_attempts,
                    struct_name,
                )
                continue

            if self._check_success(struct_name, result.stdout):
                return 0
            logger.warning(
                "Relaxation attempt %d/%d failed for %s",
                attempt,
                self.max_attempts,
                struct_name,
            )

        return 1


class AirssGulpRelaxRunner(AirssScriptRelaxRunner):
    """
    Runner for GULP relaxation via the external ``gulp_relax`` script.

    Calls ``gulp_relax <exe> <cluster> <pressure> <struct_name>`` and
    checks for success via ``gulp_relax_finish_ok()`` and ``"Volume"``
    in stdout.
    """

    _param_suffix: str = ".lib"

    def __init__(
        self,
        executable: str = "ggulp",
        timeout: int = 600,
        max_attempts: int = 3,
        cluster: bool = False,
        pressure: float = 0.0,
    ) -> None:
        super().__init__(
            executable=executable, timeout=timeout, max_attempts=max_attempts
        )
        self.cluster = cluster
        self.pressure = pressure

    def _get_cmd(self, struct_name: str) -> list[str]:
        return [
            "gulp_relax",
            self.executable,
            str(int(self.cluster)),
            str(self.pressure),
            struct_name,
        ]

    def _check_success(self, struct_name: str, stdout: str) -> bool:
        from ..casteptools import gulp_relax_finish_ok

        return gulp_relax_finish_ok(struct_name + ".castep") and "Volume" in stdout

    def _prepare_inputs(
        self,
        struct_name: str,
        struct_content: str,
        param_content: str,
        seed_name: Optional[str] = None,
    ) -> None:
        super()._prepare_inputs(struct_name, struct_content, param_content, seed_name)
        # gulp_relax looks for <seed_name>.lib, not <struct_name>.lib
        if seed_name is not None and seed_name != struct_name:
            shutil.move(struct_name + ".lib", seed_name + ".lib")


class AirssPp3RelaxRunner(AirssScriptRelaxRunner):
    """
    Runner for pp3 relaxation via the external ``pp3_relax`` script.

    Calls ``pp3_relax <exe> <struct_name>`` and checks for success
    via ``gulp_relax_finish_ok()``.
    """

    _param_suffix: str = ".pp"

    def __init__(
        self,
        executable: str = "pp3",
        timeout: int = 600,
        max_attempts: int = 3,
    ) -> None:
        super().__init__(
            executable=executable, timeout=timeout, max_attempts=max_attempts
        )

    def _get_cmd(self, struct_name: str) -> list[str]:
        return ["pp3_relax", self.executable, struct_name]
