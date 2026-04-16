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


def _create_runner(code, exe, max_iterations, cluster, pressure):
    """Create the appropriate relaxation runner for the given *code*."""
    from airsspy.jf.runners import (
        AirssAbacusRelaxRunner,
        AirssCastepRelaxRunner,
        AirssGulpRelaxRunner,
        AirssPp3RelaxRunner,
    )

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
    else:
        raise click.ClickException(f"Unknown code: {code}")


def _collect_result(struct_name: str, code: str) -> None:
    """Write a .res file from completed relaxation output."""
    if code == "castep":
        from airsspy.jf.runners import compose_task_doc

        compose_task_doc(struct_name)
    elif code == "abacus":
        from airsspy.abacustools import compose_abacus_task_doc

        compose_abacus_task_doc(struct_name)
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
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        stream=sys.stderr,
    )


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

        click.echo(f"Starting AIRSS search: seed={seed}, nmax={nmax}, code={code}")
        click.echo(f"Working directory: {workdir}")

        for i in range(1, nmax + 1):
            # Check stop file
            if _check_stop_file(workdir):
                click.echo("Stop file detected. Gracefully terminating.")
                break

            # Check walltime
            remaining = sched.get_remaining_seconds()
            if remaining < walltime_buffer:
                click.echo(
                    f"Walltime running low ({remaining}s < {walltime_buffer}s buffer). Stopping."
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
                click.echo(f"  [{i}] Buildcell timed out, skipping")
                n_failed += 1
                continue

            struct_name = build_result["struct_name"]
            n_built += 1

            if build_only:
                click.echo(f"  [{i}] Built: {struct_name}")
                continue

            # Relax
            runner = _create_runner(
                code, exe, max_iterations, cluster, pressure
            )

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
                click.echo(f"  [{i}] Relaxed OK: {struct_name}")
            else:
                try:
                    _collect_result(struct_name, code)
                    click.echo(f"  [{i}] Not converged: {struct_name}")
                except Exception:
                    click.echo(f"  [{i}] Relax FAILED: {struct_name}")
                    _emit_diagnostics(struct_name, code)
                    if not keep:
                        _clean_failed(struct_name, code)
                n_failed += 1

        # Summary
        click.echo(
            f"\nSearch complete: {n_built} built, {n_relaxed} relaxed, {n_failed} failed"
        )

        if pack and n_relaxed > 0:
            packed = _pack_res_files(workdir)
            click.echo(f"Packed {n_relaxed} .res files into {packed}")

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
    type=click.Choice(["castep", "gulp", "pp3", "abacus"]),
    help="DFT code to use",
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
    "--walltime-buffer",
    default=300,
    type=int,
    show_default=True,
    help="Seconds before walltime to stop",
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
    walltime_buffer,
):
    """Relax existing cell files locally."""
    workdir = Path(workdir).resolve()
    cell_files = sorted(workdir.glob(cell))
    if not cell_files:
        raise click.ClickException(f"No files matched pattern: {cell}")

    # Read param file
    param_suffix = SUFFIX_MAP[code]
    param_file = workdir / (seed + param_suffix)
    if not param_file.exists():
        param_file = workdir / (cell_files[0].stem + param_suffix)
    if not param_file.exists():
        raise click.ClickException(f"Param file not found: {param_file}")
    param_file_name = param_file.stem  # for use after chdir

    if exe is None:
        exe = EXE_DEFAULTS[code]

    # Detect scheduler
    sched = Scheduler.get_scheduler()
    if sched is None:
        from airsspy.scheduler import Dummy

        sched = Dummy()

    runner = _create_runner(code, exe, max_iterations, cluster, pressure)

    orig_dir = os.getcwd()
    os.chdir(workdir)

    try:
        n_relaxed = 0
        n_failed = 0
        total = len(cell_files)

        click.echo(f"Relaxing {total} structures with {code}")

        for i, cell_path in enumerate(cell_files, 1):
            struct_name = cell_path.stem

            # Check walltime
            remaining = sched.get_remaining_seconds()
            if remaining < walltime_buffer:
                click.echo(
                    f"Walltime running low ({remaining}s < {walltime_buffer}s buffer). Stopping."
                )
                break

            if code == "castep":
                from castepinput.inputs import ParamInput

                cellinput = cell_path.read_text()
                paraminput = ParamInput.from_file(param_file_name + param_suffix)
                rc = runner.run(struct_name, cellinput, paraminput)
            elif code in ("gulp", "pp3"):
                struct_content = cell_path.read_text()
                param_content = Path(param_file_name + param_suffix).read_text()
                rc = runner.run(
                    struct_name, struct_content, param_content, seed_name=seed
                )
            elif code == "abacus":
                struct_content = cell_path.read_text()
                param_content = Path(param_file_name + param_suffix).read_text()
                rc = runner.run(struct_name, struct_content, param_content)

            if rc == 0:
                _collect_result(struct_name, code)
                n_relaxed += 1
                click.echo(f"  [{i}/{total}] OK: {struct_name}")
            else:
                try:
                    _collect_result(struct_name, code)
                    click.echo(f"  [{i}/{total}] Not converged: {struct_name}")
                except Exception:
                    click.echo(f"  [{i}/{total}] FAILED: {struct_name}")
                    _emit_diagnostics(struct_name, code)
                    if not keep:
                        _clean_failed(struct_name, code)
                n_failed += 1

        click.echo(
            f"\nRelaxation complete: {n_relaxed}/{total} succeeded, {n_failed} failed"
        )

        if pack and n_relaxed > 0:
            packed = _pack_res_files(workdir)
            click.echo(f"Packed into {packed}")

    finally:
        os.chdir(orig_dir)
