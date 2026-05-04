"""
CLI commands for running AIRSS searches locally (non-jobflow, like airss.pl).
"""

import logging
import os
import shutil
import sys
from pathlib import Path

import click

from airsspy.cli.cmd_deploy import SUFFIX_MAP
from airsspy.scheduler import Scheduler

logger = logging.getLogger(__name__)

STOP_FILE_NAME = "stop"

EXE_DEFAULTS = {
    "castep": "castep.mpi",
    "gulp": "ggulp",
    "pp3": "pp3",
    "abacus": "abacus",
}


def _check_stop_file(workdir: Path) -> bool:
    """Return True if the sentinel ``stop`` file exists in *workdir*."""
    return (workdir / STOP_FILE_NAME).exists()


def _clean_failed(struct_name: str, code: str) -> None:
    """Remove intermediate files from a failed relaxation."""
    extensions = {
        "castep": [".castep", ".cell", ".param", "-out.cell", "-orig.cell"],
        "gulp": [".cell", ".lib", ".castep", ".gout", "-orig.cell"],
        "pp3": [".cell", ".pp", ".castep", "-orig.cell"],
        "abacus": [".cell", ".INPUT", "-orig.cell"],
        "ml": [".cell", ".extxyz", ".traj", "-orig.cell"],
    }
    for ext in extensions.get(code, []):
        p = Path(struct_name + ext)
        if p.is_file():
            p.unlink()
    abacus_dir = Path(struct_name + ".abacus")
    if abacus_dir.is_dir():
        shutil.rmtree(abacus_dir, ignore_errors=True)


def _emit_diagnostics(struct_name: str, code: str) -> None:
    """Write diagnostic information for a failed relaxation to stderr."""
    print(f"\n  --- Diagnostics for {struct_name} ---", file=sys.stderr)

    if code == "abacus":
        abacus_out = Path(struct_name + ".abacus") / "abacus_out"
        if abacus_out.is_file():
            lines = abacus_out.read_text().splitlines()
            tail = lines[-20:] if len(lines) > 20 else lines
            for line in tail:
                print(f"  | {line}", file=sys.stderr)

        warning_log = Path(struct_name + ".abacus") / "warning.log"
        if warning_log.is_file():
            content = warning_log.read_text().strip()
            if content:
                print("  Warnings:", file=sys.stderr)
                for line in content.splitlines()[-10:]:
                    print(f"  ! {line}", file=sys.stderr)

    elif code == "castep":
        castep_file = Path(struct_name + ".castep")
        if castep_file.is_file():
            lines = castep_file.read_text().splitlines()
            # Find the last error or warning
            for line in reversed(lines[-30:]):
                low = line.lower()
                if "error" in low or "warning" in low or "failed" in low:
                    print(f"  | {line.strip()}", file=sys.stderr)
                    break

    print(f"  --- End diagnostics ---\n", file=sys.stderr)


def _pack_res_files(workdir: Path, output_name: str = "packed.res") -> Path:
    """Concatenate all ``.res`` files in *workdir* into a single file."""
    res_files = sorted(workdir.glob("*.res"))
    if not res_files:
        raise click.ClickException("No .res files found to pack")
    packed = workdir / output_name
    with open(packed, "w") as out:
        for res_file in res_files:
            out.write(res_file.read_text())
            out.write("\n")
    return packed


def _apply_mpinp(exe: str, code: str, mpinp: int | None) -> str:
    """Prepend ``mpirun`` to *exe* when *mpinp* is set and *code* supports it."""
    if mpinp is None or code not in ("castep", "abacus"):
        return exe
    if mpinp == 0:
        return f"mpirun {exe}"
    return f"mpirun -np {mpinp} {exe}"


def _is_torchsim_model(calculator_spec: str) -> bool:
    """Check if a calculator spec uses torchsim backend format (backend:model).

    TorchSim specs use a simple ``backend:model`` format (e.g. ``mace:medium``)
    while ASE specs use ``module.path:Class`` or ``module.path.Class`` format.
    """
    from airsspy.jf.ml_runners import has_torchsim

    if not has_torchsim():
        return False
    if ":" not in calculator_spec:
        return False
    backend = calculator_spec.split(":", 1)[0]
    # TorchSim backends are single words without dots
    return "." not in backend


def _create_runner(code, exe, max_iterations, cluster, pressure, mpinp=None,
                   calculator_spec=None, calculator_kwargs=None,
                   optimizer="FIRE", fmax=0.05):
    """Create the appropriate relaxation runner for the given *code*."""
    from airsspy.jf.runners import (
        AirssAbacusRelaxRunner,
        AirssCastepRelaxRunner,
        AirssGulpRelaxRunner,
        AirssPp3RelaxRunner,
    )

    exe = _apply_mpinp(exe, code, mpinp)

    if code == "castep":
        return AirssCastepRelaxRunner(
            executable=exe,
            max_iterations=max_iterations,
        )
    elif code == "gulp":
        return AirssGulpRelaxRunner(
            executable=exe,
            cluster=cluster,
            pressure=pressure,
        )
    elif code == "pp3":
        return AirssPp3RelaxRunner(
            executable=exe,
        )
    elif code == "abacus":
        return AirssAbacusRelaxRunner(
            executable=exe,
            max_iterations=max_iterations,
            pressure=pressure,
        )
    elif code == "ml":
        from airsspy.jf.ml_runners import AirssMlRelaxRunner

        return AirssMlRelaxRunner(
            calculator_spec=calculator_spec,
            calculator_kwargs=calculator_kwargs,
            optimizer=optimizer,
            fmax=fmax,
            max_steps=max_iterations,
            pressure=pressure,
        )
    else:
        raise click.ClickException(f"Unknown code: {code}")


def _create_sp_runner(code, exe, calculator_spec=None, calculator_kwargs=None):
    """Create the appropriate single-point runner for the given *code*."""
    if code == "castep":
        from airsspy.jf.runners import AirssCastepSinglePointRunner

        return AirssCastepSinglePointRunner(executable=exe)
    elif code == "abacus":
        from airsspy.jf.runners import AirssAbacusSinglePointRunner

        return AirssAbacusSinglePointRunner(executable=exe)
    elif code == "ml":
        from airsspy.jf.ml_runners import AirssMlSinglePointRunner

        return AirssMlSinglePointRunner(
            calculator_spec=calculator_spec,
            calculator_kwargs=calculator_kwargs,
        )
    else:
        raise click.ClickException(
            f"Single-point not supported for code: {code}"
        )


def _collect_result(struct_name: str, code: str, calculator_spec: str = None) -> None:
    """Write a .res file from completed calculation output."""
    if code == "castep":
        from airsspy.jf.runners import compose_task_doc

        compose_task_doc(struct_name)
    elif code == "abacus":
        from airsspy.abacustools import compose_abacus_task_doc

        compose_abacus_task_doc(struct_name)
    elif code == "ml":
        from airsspy.jf.ml_runners import compose_ml_task_doc

        compose_ml_task_doc(struct_name, calculator_spec=calculator_spec or "")
    else:
        # pp3, gulp — use the external castep2res tool
        import subprocess

        with open(struct_name + ".res", "w") as f:
            subprocess.run(
                ["castep2res", struct_name],
                stdout=f,
                stderr=subprocess.DEVNULL,
                check=False,
            )


@click.group("run")
def run():
    """Run AIRSS searches locally (non-jobflow, like airss.pl)."""


@run.command("search")
@click.option(
    "--seed",
    required=True,
    help="Seed name (<seed>.cell and param file must exist)",
)
@click.option(
    "--nmax",
    default=1000000,
    type=int,
    show_default=True,
    help="Max number of structures to generate",
)
@click.option(
    "--code",
    default="castep",
    show_default=True,
    type=click.Choice(["castep", "gulp", "pp3", "abacus"]),
    help="DFT code to use",
)
@click.option("--exe", default=None, help="Relaxation executable (default: auto)")
@click.option(
    "--workdir",
    default=".",
    show_default=True,
    type=click.Path(),
    help="Working directory for output files",
)
@click.option("--keep", is_flag=True, help="Keep intermediate files from failed runs")
@click.option("--pack", is_flag=True, help="Concatenate .res files into packed.res")
@click.option(
    "--build-only",
    is_flag=True,
    help="Only generate structures, do not relax",
)
@click.option(
    "--pressure",
    default=0.0,
    type=float,
    show_default=True,
    help="External pressure (GPa)",
)
@click.option(
    "--max-iterations",
    default=200,
    type=int,
    show_default=True,
    help="Max total geometry iterations",
)
@click.option(
    "--build-timeout",
    default=60,
    type=int,
    show_default=True,
    help="Buildcell timeout in seconds",
)
@click.option("--cluster", is_flag=True, help="Use cluster boundary conditions (GULP)")
@click.option(
    "--mpinp",
    default=None,
    type=int,
    help="Number of MPI processes. Omit for serial, 0 for mpirun (auto), N for mpirun -np N (castep/abacus only)",
)
@click.option(
    "--walltime-buffer",
    default=300,
    type=int,
    show_default=True,
    help="Seconds before walltime to stop",
)
def run_search(
    seed,
    nmax,
    code,
    exe,
    workdir,
    keep,
    pack,
    build_only,
    pressure,
    max_iterations,
    build_timeout,
    cluster,
    mpinp,
    walltime_buffer,
):
    """Run an AIRSS random structure search locally."""
    from airsspy.jf.runners import run_buildcell

    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    # Read seed cell content
    seed_cell = Path(seed + ".cell")
    if not seed_cell.exists():
        raise click.ClickException(f"Seed cell file not found: {seed_cell}")
    seed_content = seed_cell.read_text()

    # Copy param file to workdir so runners can find it after chdir
    if not build_only:
        param_suffix = SUFFIX_MAP[code]
        param_file = Path(seed + param_suffix)
        if not param_file.exists():
            raise click.ClickException(f"Param file not found: {param_file}")
        if Path(workdir).resolve() != Path().resolve():
            shutil.copy2(param_file, workdir / param_file.name)

    if exe is None:
        exe = EXE_DEFAULTS[code]

    # Detect scheduler for walltime awareness
    sched = Scheduler.get_scheduler()
    if sched is None:
        from airsspy.scheduler import Dummy

        sched = Dummy()

    orig_dir = os.getcwd()
    os.chdir(workdir)

    try:
        n_built = 0
        n_relaxed = 0
        n_failed = 0

        logger.info("Starting AIRSS search: seed=%s, nmax=%d, code=%s", seed, nmax, code)
        logger.info("Working directory: %s", workdir)

        for i in range(1, nmax + 1):
            # Check stop file
            if _check_stop_file(workdir):
                logger.info("Stop file detected. Gracefully terminating.")
                break

            # Check walltime
            remaining = sched.get_remaining_seconds()
            if remaining < walltime_buffer:
                logger.info(
                    "Walltime running low (%ds < %ds buffer). Stopping.",
                    remaining,
                    walltime_buffer,
                )
                break

            # Build
            build_result = run_buildcell(
                seed_name=seed,
                seed_content=seed_content,
                build_timeout=build_timeout,
                write_seed=True,
            )
            if build_result is None:
                logger.info("[%d] Buildcell timed out, skipping", i)
                n_failed += 1
                continue

            struct_name = build_result["struct_name"]
            n_built += 1

            if build_only:
                logger.info("[%d] Built: %s", i, struct_name)
                continue

            # Relax
            runner = _create_runner(
                code, exe, max_iterations, cluster, pressure, mpinp
            )

            try:
                if code == "castep":
                    from castepinput.inputs import ParamInput

                    # Pass raw cell content to avoid castepinput re-serializing
                    # duplicate lattice/position blocks (e.g. both ABC and CART)
                    cellinput = Path(struct_name + ".cell").read_text()
                    paraminput = ParamInput.from_file(seed + param_suffix)
                    rc = runner.run(struct_name, cellinput, paraminput)
                elif code in ("gulp", "pp3"):
                    struct_content = Path(struct_name + ".cell").read_text()
                    param_content = Path(seed + param_suffix).read_text()
                    rc = runner.run(
                        struct_name, struct_content, param_content, seed_name=seed
                    )
                elif code == "abacus":
                    struct_content = Path(struct_name + ".cell").read_text()
                    param_content = Path(seed + param_suffix).read_text()
                    rc = runner.run(struct_name, struct_content, param_content)

                if rc == 0:
                    _collect_result(struct_name, code)
                    n_relaxed += 1
                    logger.info("[%d] Relaxed OK: %s", i, struct_name)
                else:
                    try:
                        _collect_result(struct_name, code)
                        logger.info("[%d] Not converged: %s", i, struct_name)
                    except Exception:
                        logger.info("[%d] Relax FAILED: %s", i, struct_name)
                        _emit_diagnostics(struct_name, code)
                        if not keep:
                            _clean_failed(struct_name, code)
                    n_failed += 1

            except Exception:
                logger.error(
                    "[%d] Relax crashed: %s", i, struct_name, exc_info=True
                )
                _emit_diagnostics(struct_name, code)
                if not keep:
                    _clean_failed(struct_name, code)
                n_failed += 1

        # Summary
        logger.info(
            "Search complete: %d built, %d relaxed, %d failed",
            n_built,
            n_relaxed,
            n_failed,
        )

        if pack and n_relaxed > 0:
            packed = _pack_res_files(workdir)
            logger.info("Packed %d .res files into %s", n_relaxed, packed)

    finally:
        os.chdir(orig_dir)


@run.command("relax")
@click.option("--cell", required=True, help="Glob pattern for cell files to relax")
@click.option(
    "--seed",
    default="USER-RELAX",
    show_default=True,
    help="Seed name for metadata and param file lookup",
)
@click.option(
    "--code",
    default="castep",
    show_default=True,
    type=click.Choice(["castep", "gulp", "pp3", "abacus", "ml"]),
    help="DFT code or ML mode to use",
)
@click.option("--exe", default=None, help="Relaxation executable (default: auto)")
@click.option(
    "--workdir",
    default=".",
    show_default=True,
    type=click.Path(),
    help="Working directory",
)
@click.option("--keep", is_flag=True, help="Keep intermediate files from failed runs")
@click.option("--pack", is_flag=True, help="Concatenate .res files into packed.res")
@click.option(
    "--pressure",
    default=0.0,
    type=float,
    show_default=True,
    help="External pressure (GPa)",
)
@click.option(
    "--max-iterations",
    default=200,
    type=int,
    show_default=True,
    help="Max total geometry iterations",
)
@click.option("--cluster", is_flag=True, help="Use cluster boundary conditions (GULP)")
@click.option(
    "--mpinp",
    default=None,
    type=int,
    help="Number of MPI processes. Omit for serial, 0 for mpirun (auto), N for mpirun -np N (castep/abacus only)",
)
@click.option(
    "--walltime-buffer",
    default=300,
    type=int,
    show_default=True,
    help="Seconds before walltime to stop",
)
@click.option(
    "--calculator",
    "calculator_spec",
    default=None,
    help="Model spec for --code ml. TorchSim: 'mace:medium'; ASE: 'module:Class@model'.",
)
@click.option(
    "--optimizer",
    default="FIRE",
    show_default=True,
    type=click.Choice(["FIRE", "BFGS"]),
    help="ASE optimizer for --code ml",
)
@click.option(
    "--fmax",
    default=0.05,
    type=float,
    show_default=True,
    help="Force convergence threshold (eV/Ang) for --code ml",
)
def run_relax(
    cell,
    seed,
    code,
    exe,
    workdir,
    keep,
    pack,
    pressure,
    max_iterations,
    cluster,
    mpinp,
    walltime_buffer,
    calculator_spec,
    optimizer,
    fmax,
):
    """Relax existing cell files locally."""
    if code == "ml" and not calculator_spec:
        raise click.ClickException("--calculator is required when --code ml")

    workdir = Path(workdir).resolve()
    cell_files = sorted(workdir.glob(cell))
    if not cell_files:
        raise click.ClickException(f"No files matched pattern: {cell}")

    # Read param file (not needed for ML)
    param_file_name = None
    param_suffix = None
    if code != "ml":
        param_suffix = SUFFIX_MAP[code]
        param_file = workdir / (seed + param_suffix)
        if not param_file.exists():
            param_file = workdir / (cell_files[0].stem + param_suffix)
        if not param_file.exists():
            raise click.ClickException(f"Param file not found: {param_file}")
        param_file_name = param_file.stem  # for use after chdir

    if exe is None:
        exe = EXE_DEFAULTS.get(code, "")

    # Detect scheduler
    sched = Scheduler.get_scheduler()
    if sched is None:
        from airsspy.scheduler import Dummy

        sched = Dummy()

    runner = _create_runner(
        code, exe, max_iterations, cluster, pressure, mpinp,
        calculator_spec=calculator_spec,
        optimizer=optimizer,
        fmax=fmax,
    )

    orig_dir = os.getcwd()
    os.chdir(workdir)

    try:
        n_relaxed = 0
        n_failed = 0
        total = len(cell_files)

        # TorchSim batch path: process all structures at once on GPU
        if code == "ml" and _is_torchsim_model(calculator_spec):
            from airsspy.jf.ml_runners import _torchsim_relax_batch

            struct_names = [cp.stem for cp in cell_files]
            cell_contents = [cp.read_text() for cp in cell_files]

            logger.info(
                "TorchSim batch relaxing %d structures with %s",
                total, calculator_spec,
            )
            try:
                batch_results = _torchsim_relax_batch(
                    calculator_spec, struct_names, cell_contents,
                    max_steps=max_iterations,
                    force_tol=fmax,
                    optimizer=optimizer.lower(),
                    scalar_pressure=pressure,
                )
                for sname, rc in batch_results.items():
                    if rc == 0:
                        _collect_result(sname, code, calculator_spec=calculator_spec)
                        n_relaxed += 1
                    else:
                        n_failed += 1
                        _clean_failed(sname, code)
            except Exception:
                logger.error(
                    "TorchSim batch failed", exc_info=True,
                )
                n_failed = total

        else:
            # Standard one-by-one loop (DFT codes or ASE fallback)
            logger.info("Relaxing %d structures with %s", total, code)

            for i, cell_path in enumerate(cell_files, 1):
                struct_name = cell_path.stem

                # Check walltime
                remaining = sched.get_remaining_seconds()
                if remaining < walltime_buffer:
                    logger.info(
                        "Walltime running low (%ds < %ds buffer). Stopping.",
                        remaining,
                        walltime_buffer,
                    )
                    break

                try:
                    if code == "castep":
                        from castepinput.inputs import ParamInput

                        cellinput = cell_path.read_text()
                        paraminput = ParamInput.from_file(
                            param_file_name + param_suffix
                        )
                        rc = runner.run(struct_name, cellinput, paraminput)
                    elif code in ("gulp", "pp3"):
                        struct_content = cell_path.read_text()
                        param_content = Path(
                            param_file_name + param_suffix
                        ).read_text()
                        rc = runner.run(
                            struct_name, struct_content, param_content,
                            seed_name=seed,
                        )
                    elif code == "abacus":
                        struct_content = cell_path.read_text()
                        param_content = Path(
                            param_file_name + param_suffix
                        ).read_text()
                        rc = runner.run(struct_name, struct_content, param_content)
                    elif code == "ml":
                        cell_content = cell_path.read_text()
                        rc = runner.run(struct_name, cell_content)

                    if rc == 0:
                        _collect_result(
                            struct_name, code, calculator_spec=calculator_spec
                        )
                        n_relaxed += 1
                        logger.info("[%d/%d] OK: %s", i, total, struct_name)
                    else:
                        try:
                            _collect_result(
                                struct_name, code,
                                calculator_spec=calculator_spec,
                            )
                            logger.info(
                                "[%d/%d] Not converged: %s", i, total, struct_name
                            )
                        except Exception:
                            logger.info(
                                "[%d/%d] FAILED: %s", i, total, struct_name
                            )
                            _emit_diagnostics(struct_name, code)
                            if not keep:
                                _clean_failed(struct_name, code)
                        n_failed += 1

                except Exception:
                    logger.error(
                        "[%d/%d] Relax crashed: %s",
                        i, total, struct_name, exc_info=True,
                    )
                    _emit_diagnostics(struct_name, code)
                    if not keep:
                        _clean_failed(struct_name, code)
                    n_failed += 1

        logger.info(
            "Relaxation complete: %d/%d succeeded, %d failed",
            n_relaxed,
            total,
            n_failed,
        )

        if pack and n_relaxed > 0:
            packed = _pack_res_files(workdir)
            logger.info("Packed into %s", packed)

    finally:
        os.chdir(orig_dir)


@run.command("sp")
@click.option("--cell", required=True, help="Glob pattern for cell files")
@click.option(
    "--seed",
    default="USER-SP",
    show_default=True,
    help="Seed name for metadata and param file lookup",
)
@click.option(
    "--code",
    default="castep",
    show_default=True,
    type=click.Choice(["castep", "abacus", "ml"]),
    help="Code to use for single-point calculation",
)
@click.option("--exe", default=None, help="Executable (default: auto)")
@click.option(
    "--workdir",
    default=".",
    show_default=True,
    type=click.Path(),
    help="Working directory",
)
@click.option("--keep", is_flag=True, help="Keep intermediate files from failed runs")
@click.option("--pack", is_flag=True, help="Concatenate .res files into packed.res")
@click.option(
    "--walltime-buffer",
    default=300,
    type=int,
    show_default=True,
    help="Seconds before walltime to stop",
)
@click.option(
    "--calculator",
    "calculator_spec",
    default=None,
    help="Model spec for --code ml. TorchSim: 'mace:medium'; ASE: 'module:Class@model'.",
)
def run_sp(cell, seed, code, exe, workdir, keep, pack, walltime_buffer, calculator_spec):
    """Run single-point calculations on existing cell files."""
    if code == "ml" and not calculator_spec:
        raise click.ClickException("--calculator is required when --code ml")

    workdir = Path(workdir).resolve()
    cell_files = sorted(workdir.glob(cell))
    if not cell_files:
        raise click.ClickException(f"No files matched pattern: {cell}")

    # Read param file (not needed for ML)
    param_file_name = None
    param_suffix = None
    if code != "ml":
        param_suffix = SUFFIX_MAP[code]
        param_file = workdir / (seed + param_suffix)
        if not param_file.exists():
            param_file = workdir / (cell_files[0].stem + param_suffix)
        if not param_file.exists():
            raise click.ClickException(f"Param file not found: {param_file}")
        param_file_name = param_file.stem

    if exe is None:
        exe = EXE_DEFAULTS.get(code, "")

    # Detect scheduler
    sched = Scheduler.get_scheduler()
    if sched is None:
        from airsspy.scheduler import Dummy

        sched = Dummy()

    runner = _create_sp_runner(code, exe, calculator_spec=calculator_spec)

    orig_dir = os.getcwd()
    os.chdir(workdir)

    try:
        n_done = 0
        n_failed = 0
        total = len(cell_files)

        # TorchSim batch path
        if code == "ml" and _is_torchsim_model(calculator_spec):
            from airsspy.jf.ml_runners import _torchsim_static_batch

            struct_names = [cp.stem for cp in cell_files]
            cell_contents = [cp.read_text() for cp in cell_files]

            logger.info(
                "TorchSim batch SP on %d structures with %s",
                total, calculator_spec,
            )
            try:
                batch_results = _torchsim_static_batch(
                    calculator_spec, struct_names, cell_contents,
                )
                for sname, rc in batch_results.items():
                    if rc == 0:
                        _collect_result(
                            sname, code, calculator_spec=calculator_spec
                        )
                        n_done += 1
                    else:
                        n_failed += 1
                        _clean_failed(sname, code)
            except Exception:
                logger.error("TorchSim batch SP failed", exc_info=True)
                n_failed = total

        else:
            logger.info(
                "Running single-point on %d structures with %s", total, code
            )

            for i, cell_path in enumerate(cell_files, 1):
                struct_name = cell_path.stem

                # Check walltime
                remaining = sched.get_remaining_seconds()
                if remaining < walltime_buffer:
                    logger.info(
                        "Walltime running low (%ds < %ds buffer). Stopping.",
                        remaining,
                        walltime_buffer,
                    )
                    break

                try:
                    if code == "castep":
                        from castepinput.inputs import ParamInput

                        cellinput = cell_path.read_text()
                        paraminput = ParamInput.from_file(
                            param_file_name + param_suffix
                        )
                        rc = runner.run(struct_name, cellinput, paraminput)
                    elif code == "abacus":
                        struct_content = cell_path.read_text()
                        param_content = Path(
                            param_file_name + param_suffix
                        ).read_text()
                        rc = runner.run(
                            struct_name, struct_content, param_content
                        )
                    elif code == "ml":
                        cell_content = cell_path.read_text()
                        rc = runner.run(struct_name, cell_content)

                    if rc == 0:
                        _collect_result(
                            struct_name, code, calculator_spec=calculator_spec
                        )
                        n_done += 1
                        logger.info("[%d/%d] OK: %s", i, total, struct_name)
                    else:
                        try:
                            _collect_result(
                                struct_name, code,
                                calculator_spec=calculator_spec,
                            )
                            logger.info(
                                "[%d/%d] Not converged: %s",
                                i, total, struct_name,
                            )
                        except Exception:
                            logger.info(
                                "[%d/%d] FAILED: %s", i, total, struct_name
                            )
                            _emit_diagnostics(struct_name, code)
                            if not keep:
                                _clean_failed(struct_name, code)
                        n_failed += 1

                except Exception:
                    logger.error(
                        "[%d/%d] SP crashed: %s",
                        i, total, struct_name, exc_info=True,
                    )
                    _emit_diagnostics(struct_name, code)
                    if not keep:
                        _clean_failed(struct_name, code)
                    n_failed += 1

        logger.info(
            "SP complete: %d/%d succeeded, %d failed",
            n_done,
            total,
            n_failed,
        )

        if pack and n_done > 0:
            packed = _pack_res_files(workdir)
            logger.info("Packed into %s", packed)

    finally:
        os.chdir(orig_dir)
