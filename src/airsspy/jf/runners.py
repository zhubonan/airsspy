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

    Runs successive relaxations, copying the output cell back to the
    input between cycles. Requires two consecutive successful runs to
    declare convergence. Stops at max_iterations.
    """

    def __init__(
        self,
        executable: str = "castep.mpi",
        max_fails: int = 2,
        max_iterations: int = 200,
    ) -> None:
        super().__init__(executable=executable)
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
        cycle = 0
        success_counter = 0
        iter_counter = 0

        while iter_counter < self.max_iterations:
            if fail_counter > self.max_fails:
                return 1
            cycle += 1
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

            # Copy -out.cell to .cell for next cycle, stripping ANG
            # (like airss.pl: just "sed -e '/ANG/d' out.cell > cell")
            out_cell = Path(struct_name + "-out.cell")
            if out_cell.is_file():
                out_content = out_cell.read_text()
                out_content = "\n".join(
                    line for line in out_content.splitlines()
                    if line.strip() != "ANG"
                )
                Path(struct_name + ".cell").write_text(out_content)

        return 0 if success_counter >= 2 else 1


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
                    continue
                if " Length (A)" in line:
                    in_spin_group = False
                if in_spin_group:
                    tokens = line.split()
                    if re.match(r"^ +[A-Za-z]+ ", line):
                        spin_moms.append(float(tokens[-1]))

    out_cell_path = Path(struct_name + "-out.cell")
    if out_cell_path.is_file():
        # Strip ANG keyword that castepinput can't parse
        out_content = out_cell_path.read_text()
        out_content_clean = "\n".join(
            line for line in out_content.splitlines() if line.strip() != "ANG"
        )
        tmp_path = Path(struct_name + "-out-tmp.cell")
        tmp_path.write_text(out_content_clean)
        cell = CellInput.from_file(str(tmp_path))
        tmp_path.unlink(missing_ok=True)
    else:
        cell = CellInput.from_file(struct_name + ".cell")
    elements, positions, _tags = cell.get_positions()
    atoms = Atoms(symbols=elements, positions=positions, cell=cell.get_cell(), pbc=True)

    volume = atoms.get_volume()

    # Compute symmetry via spglib
    try:
        import spglib

        sg = spglib.get_spacegroup(
            (atoms.get_cell().array, atoms.get_scaled_positions(), atoms.get_atomic_numbers()),
            symprec=0.1,
        )
        sym = sg.split()[0] if sg else "P1"
    except (ImportError, Exception):
        sym = "P1"

    # Build REM lines from .castep and .cell metadata
    from ..casteptools import build_rem_lines

    rem_lines = build_rem_lines(struct_name)

    info = {
        "uid": struct_name,
        "H": energy if energy else 0.0,
        "P": pressure if pressure is not None else 0.0,
        "V": volume,
        "nat": len(atoms),
        "sym": sym,
        "rem": rem_lines,
    }
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
        "rem_lines": rem_lines,
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


class AirssAbacusRelaxRunner:
    """
    Execute a cyclic ABACUS geometry optimisation.

    Calls the ABACUS binary directly and manages the relaxation loop
    in Python, following the same two-phase pattern as ``abacus_relax``:

    1. Three short rough runs with ``relax_nmax=3``
    2. Full convergence loop until two successive convergences

    Between each ABACUS invocation, the structure is read from
    ``STRU_ION_D``, converted back to .cell format, and fed into the
    next iteration.
    """

    def __init__(
        self,
        executable: str = "abacus",
        max_fails: int = 2,
        max_iterations: int = 200,
        pressure: float = 0.0,
    ) -> None:
        self.executable = executable
        self.max_fails = max_fails
        self.max_iterations = max_iterations
        self.pressure = pressure

    def _set_input_param(self, input_path: str, key: str, value: str) -> None:
        """Set or add a parameter in an ABACUS INPUT file."""
        content = Path(input_path).read_text()
        lines = content.splitlines()
        found = False
        for i, line in enumerate(lines):
            if re.match(rf"^\s*{re.escape(key)}", line):
                lines[i] = f"{key} {value}"
                found = True
                break
        if not found:
            lines.append(f"{key} {value}")
        Path(input_path).write_text("\n".join(lines))

    def _detect_logfile(self, workdir: str, input_path: str) -> str:
        """Detect the ABACUS log file path based on calculation type."""
        from ..abacustools import detect_logfile

        result = detect_logfile(workdir, input_path)
        return result or ""

    def prepare_inputs(
        self,
        struct_name: str,
        cell_content: str,
        input_content: str,
    ) -> None:
        """Write .cell and .INPUT files, convert .cell to STRU.

        Args:
            struct_name: Structure name (without extension).
            cell_content: Content of the .cell file.
            input_content: Content of the ABACUS INPUT file.
        """
        workdir = f"{struct_name}.abacus"
        Path(workdir).mkdir(parents=True, exist_ok=True)

        # Write .cell file
        cell_path = struct_name + ".cell"
        Path(cell_path).write_text(cell_content)

        # Write INPUT file
        input_path = struct_name + ".INPUT"
        Path(input_path).write_text(input_content)

        # Convert .cell to STRU
        from ..abacustools import cell_to_stru

        stru_content = cell_to_stru(cell_content)
        Path(f"{workdir}/STRU").write_text(stru_content)

        # Copy INPUT to workdir
        Path(f"{workdir}/INPUT").write_text(input_content)

    def _run_single(
        self,
        struct_name: str,
        workdir: str,
        input_path: str,
    ) -> Optional[dict]:
        """Run a single ABACUS calculation and parse results.

        Returns:
            Dict with converged, n_steps, energy, pressure, volume
            or None if the calculation crashed.
        """
        out_path = Path(f"{workdir}/abacus_out")
        with open(out_path, "w") as outf:
            proc = subprocess.run(
                self.executable.split(),
                stdout=outf,
                stderr=subprocess.STDOUT,
                cwd=workdir,
                check=False,
            )

        output = out_path.read_text()

        # Check for crash (ABACUS writes "TOTAL  Time" on success)
        if "TOTAL  Time" not in output:
            logger.warning("ABACUS crashed for %s", struct_name)
            for line in output.splitlines()[-5:]:
                logger.info("  | %s", line)
            return None

        logfile = self._detect_logfile(workdir, input_path)
        if not logfile:
            logger.warning("No ABACUS log file found for %s", struct_name)
            return None

        from ..abacustools import parse_abacus_log

        log_data = parse_abacus_log(logfile)
        converged = log_data.get("converged", False)
        n_steps = log_data.get("n_ionic_steps", 0)
        energy = log_data.get("energy")
        logger.info(
            "%s: %s, %d ionic steps, E=%.6f eV",
            struct_name,
            "converged" if converged else "not converged",
            n_steps,
            energy,
        )
        return log_data

    def run(
        self,
        struct_name: str,
        cell_content: str,
        input_content: str,
    ) -> int:
        """
        Run cyclic ABACUS relaxation.

        Phase 1: three rough runs with relax_nmax=3 for initial optimisation.
        Phase 2: full convergence loop until two consecutive convergences
        or max_iterations is reached.

        Args:
            struct_name: Structure name (without extension).
            cell_content: Content of the .cell file.
            input_content: Content of the ABACUS INPUT file.

        Returns:
            0 if converged, 1 if not converged or failed.
        """
        self.prepare_inputs(struct_name, cell_content, input_content)

        workdir = f"{struct_name}.abacus"
        input_path = f"{workdir}/INPUT"

        # Save original relax_nmax if set by user
        user_relax_nmax = None
        for line in Path(input_path).read_text().splitlines():
            m = re.match(r"^\s*relax_nmax\s+(\S+)", line)
            if m:
                user_relax_nmax = m.group(1)
                break

        fail_counter = 0
        success_counter = 0
        iter_counter = 0

        # Phase 1: three rough runs with relax_nmax=3
        logger.info("%s: phase 1 — 3 rough runs (relax_nmax=3)", struct_name)
        self._set_input_param(input_path, "relax_nmax", "3")

        for rough_i in range(3):
            if fail_counter > self.max_fails:
                logger.error(
                    "%s: too many failures (%d), aborting", struct_name, fail_counter
                )
                return 1

            logger.info(
                "%s: rough run %d/3 (iter=%d/%d)",
                struct_name, rough_i + 1, iter_counter, self.max_iterations,
            )
            result = self._run_single(struct_name, workdir, input_path)
            if result is None:
                fail_counter += 1
                continue

            fail_counter = 0
            iter_counter += result.get("n_ionic_steps", 0)

            if result.get("converged"):
                success_counter += 1
            else:
                success_counter = 0

            self._update_cell(struct_name, workdir)

        # Restore original relax_nmax for phase 2
        if user_relax_nmax is None:
            content = Path(input_path).read_text()
            lines = [
                line
                for line in content.splitlines()
                if not re.match(r"^\s*relax_nmax", line)
            ]
            Path(input_path).write_text("\n".join(lines))
        else:
            self._set_input_param(input_path, "relax_nmax", user_relax_nmax)

        # Phase 2: full convergence loop
        if success_counter < 2:
            cycle = 0
            while iter_counter < self.max_iterations:
                if fail_counter > self.max_fails:
                    logger.error(
                        "%s: too many failures (%d), aborting", struct_name, fail_counter
                    )
                    return 1

                cycle += 1
                logger.info(
                    "%s: phase 2 cycle %d (iter=%d/%d, consecutive_ok=%d)",
                    struct_name, cycle, iter_counter, self.max_iterations, success_counter,
                )
                result = self._run_single(struct_name, workdir, input_path)
                if result is None:
                    fail_counter += 1
                    continue

                fail_counter = 0
                iter_counter += result.get("n_ionic_steps", 0)

                if result.get("converged"):
                    success_counter += 1
                else:
                    success_counter = 0

                if success_counter >= 2:
                    break

                self._update_cell(struct_name, workdir)

        converged = success_counter >= 2
        logger.info(
            "%s: finished — %s (%d total iterations)",
            struct_name,
            "converged" if converged else "not converged",
            iter_counter,
        )
        return 0 if converged else 1

    def _update_cell(self, struct_name: str, workdir: str) -> None:
        """Read STRU_ION_D output and update the input .cell file."""
        from ..abacustools import parse_abacus_stru

        stru_path = Path(workdir) / "OUT.ABACUS" / "STRU_ION_D"
        if not stru_path.is_file():
            return

        elements, positions, cell = parse_abacus_stru(str(stru_path))

        # Convert fractional to Cartesian positions
        cart_positions = positions @ cell

        cell_path = struct_name + ".cell"
        if not Path(cell_path).is_file():
            return

        # Rewrite the .cell file with updated lattice and positions,
        # preserving other blocks (SPECIES_POT, KPOINTS, etc.)
        old_content = Path(cell_path).read_text()
        new_lattice = "\n".join(
            f"  {v[0]:.10f}  {v[1]:.10f}  {v[2]:.10f}" for v in cell.tolist()
        )
        new_positions = "\n".join(
            f"{e}  {p[0]:.10f}  {p[1]:.10f}  {p[2]:.10f}"
            for e, p in zip(elements, cart_positions.tolist())
        )

        # Replace LATTICE_CART block
        new_content = re.sub(
            r"%BLOCK\s+LATTICE_CART.*?%ENDBLOCK\s+LATTICE_CART",
            f"%BLOCK LATTICE_CART\n{new_lattice}\n%ENDBLOCK LATTICE_CART",
            old_content,
            count=1,
            flags=re.DOTALL | re.IGNORECASE,
        )

        # Replace POSITIONS_FRAC block
        new_content = re.sub(
            r"%BLOCK\s+POSITIONS_FRAC.*?%ENDBLOCK\s+POSITIONS_FRAC",
            f"%BLOCK POSITIONS_FRAC\n{new_positions}\n%ENDBLOCK POSITIONS_FRAC",
            new_content,
            count=1,
            flags=re.DOTALL | re.IGNORECASE,
        )

        # Remove positions_abs block (not needed, avoids #label pollution)
        new_content = re.sub(
            r"%BLOCK\s+POSITIONS_ABS.*?%ENDBLOCK\s+POSITIONS_ABS\s*",
            "",
            new_content,
            count=1,
            flags=re.DOTALL | re.IGNORECASE,
        )

        Path(cell_path).write_text(new_content)

        # Also regenerate STRU for next ABACUS run
        cell_content = Path(cell_path).read_text()
        from ..abacustools import cell_to_stru

        stru_content = cell_to_stru(cell_content)
        Path(f"{workdir}/STRU").write_text(stru_content)
