"""
CLI commands for running AIRSS searches locally (non-jobflow, like airss.pl).
"""

import logging
import os
import random
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
    "vasp": "vasp_std",
}


def _check_stop_file(workdir: Path) -> bool:
    """Return True if the sentinel ``stop`` file exists in *workdir*."""
    return (workdir / STOP_FILE_NAME).exists()


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

    elif code == "vasp":
        for rel_path in ("vasp.out", "OUTCAR", "OSZICAR"):
            path = Path(struct_name + ".vasp") / rel_path
            if not path.is_file():
                continue
            lines = path.read_text(errors="ignore").splitlines()
            tail = lines[-20:] if len(lines) > 20 else lines
            print(f"  {rel_path}:", file=sys.stderr)
            for line in tail:
                low = line.lower()
                if any(word in low for word in ("error", "warning", "fail", "fatal")):
                    print(f"  | {line.strip()}", file=sys.stderr)
            if tail:
                break

    print("  --- End diagnostics ---\n", file=sys.stderr)


def _pack_res_files(
    workdir: Path,
    output_name: str = "packed.res",
    files: list[Path] | None = None,
) -> Path:
    """Concatenate all ``.res`` files in *workdir* into a single file."""
    res_files = sorted(files) if files is not None else sorted(workdir.glob("*.res"))
    if not res_files:
        raise click.ClickException("No .res files found to pack")
    packed = workdir / output_name
    with open(packed, "w") as out:
        for res_file in res_files:
            out.write(res_file.read_text())
            out.write("\n")
    return packed


def _is_packed_res_input(path: Path) -> bool:
    """Return True if *path* appears to contain more than one RES structure."""
    if path.suffix.lower() != ".res":
        return False
    titl_count = 0
    try:
        with open(path) as handle:
            for line in handle:
                if line.startswith("TITL"):
                    titl_count += 1
                    if titl_count > 1:
                        return True
    except OSError:
        return False
    return False


def _filter_packed_res_inputs(paths: list[Path]) -> list[Path]:
    """Drop packed RES files from per-structure relax/SP input lists."""
    skipped = [path for path in paths if _is_packed_res_input(path)]
    if skipped:
        preview = ", ".join(path.name for path in skipped[:5])
        if len(skipped) > 5:
            preview += ", ..."
        logger.warning("Skipping %d packed RES input(s): %s", len(skipped), preview)
    return [path for path in paths if path not in skipped]


def _cleanup_ml_transients(struct_name: str) -> None:
    """Remove transient ML outputs without deleting input/output RES files."""
    for suffix in (".extxyz", ".traj", ".err"):
        Path(struct_name + suffix).unlink(missing_ok=True)


def _apply_mpinp(exe: str, code: str, mpinp: int | None) -> str:
    """Prepend ``mpirun`` to *exe* when *mpinp* is set and *code* supports it."""
    if mpinp is None or code not in ("castep", "abacus", "vasp"):
        return exe
    if mpinp == 0:
        return f"mpirun {exe}"
    return f"mpirun -np {mpinp} {exe}"


def _parse_formula_key_values(values, parser, option_name: str) -> dict:
    """Parse repeatable CLI ``KEY=value`` options."""
    parsed = {}
    for value in values:
        try:
            key, item = parser(value)
        except ValueError as exc:
            raise click.ClickException(f"Invalid {option_name}: {exc}") from exc
        parsed[key] = item
    return parsed


def _parse_oxidation_state_options(values) -> dict[str, list[int]]:
    """Parse oxidation states from comma-separated ``Element=state`` groups."""
    parsed: dict[str, list[int]] = {}
    for value in values:
        assignments = [item.strip() for item in value.split(",") if item.strip()]
        for assignment in assignments:
            try:
                key, state = assignment.split("=", 1)
                key = key.strip()
                state = state.strip()
                if not key or not state:
                    raise ValueError
                parsed.setdefault(key, []).append(int(state))
            except ValueError as exc:
                raise click.ClickException(
                    f"Invalid --oxidation-state: expected Element=state "
                    f"assignments, got {assignment!r}"
                ) from exc
    return parsed


def _parse_potcar_map_options(values) -> dict[str, str]:
    """Parse repeatable VASP ``Element=symbol`` POTCAR map options."""
    from airsspy.vasptools import parse_potcar_map

    try:
        return parse_potcar_map(values)
    except ValueError as exc:
        raise click.ClickException(f"Invalid --potcar-map: {exc}") from exc


def _resolve_optional_path(path: str | None) -> str | None:
    """Resolve an optional user path before command handlers chdir."""
    return str(Path(path).expanduser().resolve()) if path else None


def _parse_formula_option(value: str) -> list[str]:
    """Parse formulas from a comma-separated CLI value."""
    return [item.strip() for item in value.split(",") if item.strip()]


def _parse_formula_option_callback(ctx, param, value):
    """Click callback for the non-repeatable comma-separated formula option."""
    del ctx, param
    if len(value) > 1:
        raise click.BadParameter("may be specified only once")
    if not value:
        return []
    return _parse_formula_option(value[0])


def _move_pruned_file(candidate, workdir: Path) -> None:
    """Move a rejected candidate's RES file into ``workdir/pruned``."""
    if candidate.res_path is None or not candidate.res_path.exists():
        return
    pruned_dir = workdir / "pruned"
    pruned_dir.mkdir(exist_ok=True)
    target = pruned_dir / candidate.res_path.name
    candidate.res_path.replace(target)


def _strip_structure_blocks(cell_text: str) -> list[str]:
    """Return cell lines with lattice and positions blocks removed."""
    import re

    block_start = re.compile(
        r"^\s*%BLOCK\s+(LATTICE_(?:CART|ABC)|POSITIONS_(?:FRAC|ABS))\b", re.I
    )
    block_end = re.compile(
        r"^\s*%ENDBLOCK\s+(LATTICE_(?:CART|ABC)|POSITIONS_(?:FRAC|ABS))\b", re.I
    )

    lines: list[str] = []
    in_structure_block = False
    for line in cell_text.splitlines():
        if block_start.search(line):
            in_structure_block = True
            continue
        if in_structure_block:
            if block_end.search(line):
                in_structure_block = False
            continue
        lines.append(line)
    return lines


def _read_crud_res_spins(res_lines: list[str], natoms: int) -> list[float]:
    """Read unambiguous per-site spin columns from RES atom lines."""
    spins: list[float] = []
    in_atoms = False
    for line in res_lines:
        tokens = line.split()
        if not tokens:
            continue
        if tokens[0] == "SFAC":
            in_atoms = True
            continue
        if tokens[0] == "END":
            break
        if not in_atoms or not tokens[0][0].isalpha():
            continue

        if len(tokens) == 7 or len(tokens) >= 10:
            try:
                spins.append(float(tokens[6]))
            except ValueError:
                return []
        else:
            return []

    return spins if len(spins) == natoms else []


def _res_to_cell_lines(res_path: Path, root_cell_path: Path) -> list[str]:
    """Build a CASTEP cell file from a RES geometry and root cell settings."""
    from airsspy.restools import read_res_atoms

    res_lines = res_path.read_text().splitlines()
    _, atoms = read_res_atoms(res_lines)
    spins = _read_crud_res_spins(res_lines, len(atoms))

    lines = ["%BLOCK LATTICE_CART"]
    for vec in atoms.cell:
        lines.append(f"{vec[0]:.10f} {vec[1]:.10f} {vec[2]:.10f}")
    lines.append("%ENDBLOCK LATTICE_CART")
    lines.append("%BLOCK POSITIONS_ABS")
    for i, (symbol, pos) in enumerate(
        zip(atoms.get_chemical_symbols(), atoms.positions)
    ):
        line = f"{symbol}  {pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f}"
        if spins:
            line += f" SPIN={spins[i]:.3f}"
        lines.append(line)
    lines.append("%ENDBLOCK POSITIONS_ABS")

    rest = _strip_structure_blocks(root_cell_path.read_text())
    while rest and not rest[0].strip():
        rest.pop(0)
    if rest:
        lines.append("")
        lines.extend(rest)
    return lines


def _claim_crud_job(workdir: Path, num: int) -> Path | None:
    """Atomically move one queued RES file from hopper into *workdir*."""
    hopper = workdir / "hopper"
    files = list(hopper.glob("*-*.res"))
    random.shuffle(files)
    for src in files[:num]:
        dst = workdir / src.name
        if dst.exists():
            continue
        try:
            src.rename(dst)
        except FileNotFoundError:
            continue
        except OSError:
            continue
        if dst.exists():
            return dst
    return None


def _crud_root_from_seed(seed: str) -> str:
    """Return AIRSS root name from a structure label."""
    return seed.split("-", 1)[0]


def _sync_castep_spin_param(cell_path: Path, param_path: Path) -> None:
    """Match CASTEP param spin to per-site SPIN tags in the cell file."""
    spin_total = 0.0
    has_spin = False
    for line in cell_path.read_text().splitlines():
        if "#SPIN=" in line or "SPIN=" not in line:
            continue
        try:
            spin_total += float(line.split("SPIN=", 1)[1].split()[0])
            has_spin = True
        except (ValueError, IndexError):
            continue

    if not has_spin:
        return

    lines = []
    if param_path.exists():
        for line in param_path.read_text().splitlines():
            tokens = line.split()
            if tokens and tokens[0].lower() == "spin":
                continue
            lines.append(line)
    lines.append(f"spin : {spin_total:10.3f}")
    param_path.write_text("\n".join(lines) + "\n")


def _prepare_crud_inputs(seed: str, code: str) -> tuple[str, str | None]:
    """Create per-structure input files for a claimed CRUD job."""
    root = _crud_root_from_seed(seed)
    if code == "ml":
        return root, None

    root_cell = Path(root + ".cell")
    if not root_cell.exists():
        raise click.ClickException(f"Root cell file not found: {root_cell}")

    cell_path = Path(seed + ".cell")
    cell_path.write_text(
        "\n".join(_res_to_cell_lines(Path(seed + ".res"), root_cell)) + "\n"
    )

    param_suffix = SUFFIX_MAP[code]
    root_param = Path(root + param_suffix)
    if not root_param.exists():
        raise click.ClickException(f"Root input file not found: {root_param}")
    seed_param = Path(seed + param_suffix)
    if root_param.resolve() != seed_param.resolve():
        shutil.copy2(root_param, seed_param)
    if code == "castep":
        _sync_castep_spin_param(cell_path, seed_param)
    if code == "vasp":
        root_kpoints = Path(root + ".KPOINTS")
        if root_kpoints.exists():
            shutil.copy2(root_kpoints, seed + ".KPOINTS")
    return root, param_suffix


def _resolve_relax_template_cell(input_path: Path, workdir: Path) -> Path:
    """Find a template .cell file for converting a RES input to CASTEP cell text."""
    root = _crud_root_from_seed(input_path.stem)
    candidates = [workdir / f"{root}.cell"]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    names = ", ".join(str(path.name) for path in candidates)
    raise click.ClickException(
        f"Template .cell file not found for RES input {input_path.name}; "
        f"looked for {names}"
    )


def _prepare_relax_input(
    input_path: Path,
    workdir: Path,
    *,
    write_res_cell: bool = True,
    convert_res_cell: bool = True,
) -> tuple[Path, str, str]:
    """Return cell path/name/content for a relax input, converting RES if needed."""
    suffix = input_path.suffix.lower()
    if suffix == ".cell":
        return input_path, input_path.stem, input_path.read_text()
    if suffix != ".res":
        raise click.ClickException(
            f"Unsupported relax input extension for {input_path.name}; "
            "expected .cell or .res"
        )

    cell_path = input_path.with_suffix(".cell")
    if not convert_res_cell:
        return cell_path, cell_path.stem, ""

    template_cell = _resolve_relax_template_cell(input_path, workdir)
    cell_content = "\n".join(_res_to_cell_lines(input_path, template_cell)) + "\n"
    if write_res_cell:
        cell_path.write_text(cell_content)
    return cell_path, cell_path.stem, cell_content


def _read_res_as_atoms(res_path: Path):
    """Read a RES input directly as ASE Atoms for ML backends."""
    from airsspy.restools import read_res_atoms

    _, atoms = read_res_atoms(res_path.read_text().splitlines())
    return atoms


def _prepare_ml_structure_input(input_path: Path, cell_content: str):
    """Return the preferred ML input representation for a relax/SP structure."""
    if input_path.suffix.lower() == ".res":
        return _read_res_as_atoms(input_path)
    return cell_content


def _resolve_relax_param_file(
    input_path: Path,
    cell_path: Path,
    seed: str,
    workdir: Path,
    param_suffix: str,
) -> Path:
    """Find the parameter file for a relax input."""
    if input_path.suffix.lower() == ".res":
        candidates = [workdir / f"{_crud_root_from_seed(input_path.stem)}{param_suffix}"]
    else:
        candidates = [
            workdir / f"{seed}{param_suffix}",
            workdir / f"{cell_path.stem}{param_suffix}",
        ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    names = ", ".join(str(path.name) for path in candidates)
    raise click.ClickException(f"Param file not found; looked for {names}")


def _run_local_structure_one(
    input_path: Path,
    cell_path: Path,
    struct_name: str,
    cell_content: str,
    code: str,
    runner,
    param_suffix: str | None,
    seed: str,
    workdir: Path,
    *,
    singlepoint: bool = False,
) -> int:
    """Run one local relax/SP structure with an existing runner."""
    if code == "castep":
        from castepinput.inputs import ParamInput

        param_file = _resolve_relax_param_file(
            input_path, cell_path, seed, workdir, param_suffix
        )
        return runner.run(struct_name, cell_content, ParamInput.from_file(param_file))
    if code in ("gulp", "pp3"):
        if singlepoint:
            raise click.ClickException(f"Single-point not supported for code: {code}")
        param_file = _resolve_relax_param_file(
            input_path, cell_path, seed, workdir, param_suffix
        )
        return runner.run(
            struct_name,
            cell_content,
            param_file.read_text(),
            seed_name=_crud_root_from_seed(struct_name),
        )
    if code == "abacus":
        param_file = _resolve_relax_param_file(
            input_path, cell_path, seed, workdir, param_suffix
        )
        return runner.run(struct_name, cell_content, param_file.read_text())
    if code == "vasp":
        param_file = _resolve_relax_param_file(
            input_path, cell_path, seed, workdir, param_suffix
        )
        kpoints_path = param_file.with_suffix(".KPOINTS")
        if not kpoints_path.exists():
            kpoints_path = workdir / f"{seed}.KPOINTS"
        return runner.run(
            struct_name,
            cell_content,
            param_file.read_text(),
            kpoints_path=kpoints_path if kpoints_path.exists() else None,
        )
    if code == "ml":
        ml_input = _prepare_ml_structure_input(input_path, cell_content)
        return runner.run(struct_name, ml_input)
    raise click.ClickException(f"Unknown code: {code}")


def _run_crud_one(
    seed: str,
    code: str,
    runner,
    param_suffix: str | None,
    workdir: Path,
    *,
    singlepoint: bool = False,
) -> int:
    """Run one claimed CRUD structure with an existing local runner."""
    cell_path = Path(seed + ".cell")
    input_path = Path(seed + ".res") if code == "ml" else cell_path
    cell_content = "" if code == "ml" else cell_path.read_text()
    return _run_local_structure_one(
        input_path,
        cell_path,
        seed,
        cell_content,
        code,
        runner,
        param_suffix,
        seed,
        workdir,
        singlepoint=singlepoint,
    )


def _create_task_runner(
    code,
    exe,
    max_iterations,
    cluster,
    pressure,
    mpinp=None,
    calculator_spec=None,
    calculator_kwargs=None,
    optimizer="FIRE",
    fmax=0.05,
    potcar_dir=None,
    potcar_map=None,
    cell_axis_map=None,
    *,
    singlepoint: bool = False,
):
    """Create a relaxation or single-point runner for a local run command."""
    if not singlepoint:
        return _create_runner(
            code,
            exe,
            max_iterations,
            cluster,
            pressure,
            mpinp,
            calculator_spec=calculator_spec,
            calculator_kwargs=calculator_kwargs,
            optimizer=optimizer,
            fmax=fmax,
            potcar_dir=potcar_dir,
            potcar_map=potcar_map,
            cell_axis_map=cell_axis_map,
        )
    exe = _apply_mpinp(exe, code, mpinp)
    return _create_sp_runner(
        code,
        exe,
        calculator_spec=calculator_spec,
        calculator_kwargs=calculator_kwargs,
        pressure=pressure,
        potcar_dir=potcar_dir,
        potcar_map=potcar_map,
        cell_axis_map=cell_axis_map,
    )


def _run_torchsim_batch(
    torchsim_runner,
    struct_names: list[str],
    structures: list,
    code: str,
    calculator_spec: str,
    *,
    singlepoint: bool = False,
    max_iterations: int = 200,
    fmax: float = 0.05,
    optimizer: str = "FIRE",
    pressure: float = 0.0,
) -> tuple[int, int, list[Path]]:
    """Run and collect one torch-sim batch."""
    if singlepoint:
        batch_results = torchsim_runner.static_batch(
            struct_names,
            structures,
            scalar_pressure=pressure,
        )
    else:
        batch_results = torchsim_runner.relax_batch(
            struct_names,
            structures,
            max_steps=max_iterations,
            force_tol=fmax,
            optimizer=optimizer.lower(),
            scalar_pressure=pressure,
        )

    n_done = 0
    n_failed = 0
    collected_res_files: list[Path] = []
    for sname, rc in batch_results.items():
        if rc == 0:
            try:
                _collect_result(sname, code, calculator_spec=calculator_spec)
                collected_res_files.append(Path(f"{sname}.res"))
                n_done += 1
            except Exception:
                logger.error(
                    "TorchSim result collection failed: %s", sname, exc_info=True
                )
                n_failed += 1
                _cleanup_ml_transients(sname)
        else:
            n_failed += 1
            _cleanup_ml_transients(sname)
    return n_done, n_failed, collected_res_files


def _crud_artifacts(seed: str) -> list[Path]:
    """Return files/directories belonging to a CRUD structure."""
    paths = list(Path().glob(seed + ".*"))
    out_cell = Path(seed + "-out.cell")
    if out_cell.exists():
        paths.append(out_cell)
    seen = set()
    unique = []
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        unique.append(path)
    return unique


def _cleanup_crud_artifacts(seed: str) -> None:
    """Remove less useful intermediate files before finalizing outputs."""
    keep_suffixes = {".res", ".cif", ".magres", ".castep", ".odo", ".dos", ".den_fmt"}
    for path in _crud_artifacts(seed):
        if path.is_dir() or path.suffix in keep_suffixes:
            continue
        path.unlink(missing_ok=True)


def _move_crud_artifacts(seed: str, target_dir: Path, keep: bool) -> None:
    """Move structure artifacts into good/bad output directories."""
    if not keep:
        _cleanup_crud_artifacts(seed)
    target_dir.mkdir(parents=True, exist_ok=True)
    for path in _crud_artifacts(seed):
        target = target_dir / path.name
        if target.exists():
            if target.is_dir():
                shutil.rmtree(target)
            else:
                target.unlink()
        shutil.move(str(path), str(target))


def _crud_retryable_failure(seed: str) -> bool:
    """Return whether a failed CASTEP-like job should be cycled."""
    for err_path in Path().glob(seed + "*.err"):
        try:
            if "electronic_minimisation" in err_path.read_text(errors="ignore"):
                return True
        except OSError:
            continue
    return False


def _requeue_crud_job(seed: str, workdir: Path) -> bool:
    """Move a checked-out RES file back into hopper."""
    res_path = Path(seed + ".res")
    if not res_path.exists():
        return False
    hopper = workdir / "hopper"
    hopper.mkdir(exist_ok=True)
    target = hopper / res_path.name
    if target.exists():
        return False
    res_path.rename(target)
    return True


def _get_scheduler_for_walltime():
    """Return a scheduler object for walltime checks, falling back to Dummy."""
    try:
        sched = Scheduler.get_scheduler()
    except Exception as exc:
        logger.warning("Could not detect scheduler for walltime checks: %s", exc)
        sched = None
    if sched is None:
        from airsspy.scheduler import Dummy

        sched = Dummy()
    return sched


def _walltime_remaining_ok(sched, walltime_buffer: int) -> bool:
    """Return whether there is enough walltime, ignoring malformed scheduler data."""
    try:
        remaining = sched.get_remaining_seconds()
    except Exception as exc:
        logger.warning("Could not determine remaining walltime; continuing: %s", exc)
        return True
    if remaining < walltime_buffer:
        logger.info(
            "Walltime running low (%ds < %ds buffer). Stopping.",
            remaining,
            walltime_buffer,
        )
        return False
    return True


def _is_torchsim_model(calculator_spec: str) -> bool:
    """Return whether a model spec should use the torch-sim backend."""
    if calculator_spec.startswith("ase:"):
        return False
    if ":" not in calculator_spec:
        return False
    backend = calculator_spec.split(":", 1)[0]
    return "." not in backend


def _ensure_torchsim_available() -> None:
    """Raise a CLI error if torch-sim is not importable."""
    from airsspy.jf.ml_runners import has_torchsim

    if not has_torchsim():
        raise click.ClickException(
            "torch-sim is required for plain ML model specs. "
            "Install torch-sim or use an explicit ASE fallback spec such as "
            "ase:mace:medium."
        )


def _normalize_ml_ase_spec(model_spec: str) -> str:
    """Convert an explicit ``ase:`` model spec into an ASE calculator spec."""
    if not model_spec.startswith("ase:"):
        return model_spec
    ase_spec = model_spec[4:]
    if ase_spec.startswith("mace:"):
        model_id = ase_spec.split(":", 1)[1]
        return f"mace.calculators:MACECalculator@{model_id}"
    return ase_spec


def _create_runner(
    code,
    exe,
    max_iterations,
    cluster,
    pressure,
    mpinp=None,
    calculator_spec=None,
    calculator_kwargs=None,
    optimizer="FIRE",
    fmax=0.05,
    potcar_dir=None,
    potcar_map=None,
    cell_axis_map=None,
):
    """Create the appropriate relaxation runner for the given *code*."""
    from airsspy.jf.runners import (
        AirssAbacusRelaxRunner,
        AirssCastepRelaxRunner,
        AirssGulpRelaxRunner,
        AirssPp3RelaxRunner,
        AirssVaspRelaxRunner,
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
            cell_axis_map=cell_axis_map,
        )
    elif code == "vasp":
        return AirssVaspRelaxRunner(
            executable=exe,
            max_iterations=max_iterations,
            pressure=pressure,
            potcar_dir=potcar_dir,
            potcar_map=potcar_map,
        )
    elif code == "ml":
        from airsspy.jf.ml_runners import AirssMlRelaxRunner

        return AirssMlRelaxRunner(
            calculator_spec=_normalize_ml_ase_spec(calculator_spec),
            calculator_kwargs=calculator_kwargs,
            optimizer=optimizer,
            fmax=fmax,
            max_steps=max_iterations,
            pressure=pressure,
        )
    else:
        raise click.ClickException(f"Unknown code: {code}")


def _create_sp_runner(
    code,
    exe,
    calculator_spec=None,
    calculator_kwargs=None,
    pressure=0.0,
    potcar_dir=None,
    potcar_map=None,
    cell_axis_map=None,
):
    """Create the appropriate single-point runner for the given *code*."""
    if code == "castep":
        from airsspy.jf.runners import AirssCastepSinglePointRunner

        return AirssCastepSinglePointRunner(executable=exe)
    elif code == "abacus":
        from airsspy.jf.runners import AirssAbacusSinglePointRunner

        return AirssAbacusSinglePointRunner(
            executable=exe,
            pressure=pressure,
            cell_axis_map=cell_axis_map,
        )
    elif code == "vasp":
        from airsspy.jf.runners import AirssVaspSinglePointRunner

        return AirssVaspSinglePointRunner(
            executable=exe,
            pressure=pressure,
            potcar_dir=potcar_dir,
            potcar_map=potcar_map,
        )
    elif code == "ml":
        from airsspy.jf.ml_runners import AirssMlSinglePointRunner

        return AirssMlSinglePointRunner(
            calculator_spec=_normalize_ml_ase_spec(calculator_spec),
            calculator_kwargs=calculator_kwargs,
            pressure=pressure,
        )
    else:
        raise click.ClickException(f"Single-point not supported for code: {code}")


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
    elif code == "vasp":
        from airsspy.vasptools import compose_vasp_task_doc

        metadata = getattr(calculator_spec, "last_metadata", None)
        if isinstance(calculator_spec, dict):
            metadata = calculator_spec
        compose_vasp_task_doc(struct_name, metadata=metadata)
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
    type=click.Choice(["castep", "gulp", "pp3", "abacus", "vasp"]),
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
@click.option(
    "--formula",
    "formulas",
    multiple=True,
    callback=_parse_formula_option_callback,
    help="Comma-separated formulas to sample by injecting #FORMULA.",
)
@click.option(
    "--elements",
    "formula_elements",
    default="",
    help="Comma-separated elements for formula enumeration when --formula is omitted.",
)
@click.option(
    "--max-coeff",
    "formula_max_coeff",
    default=6,
    type=int,
    show_default=True,
    help="Maximum reduced-formula coefficient for formula enumeration.",
)
@click.option(
    "--target-volume",
    "formula_target_volumes",
    multiple=True,
    help="Per-element target volume as Element=value. Repeat as needed.",
)
@click.option(
    "--oxidation-state",
    "formula_oxidation_states",
    multiple=True,
    help="Allowed oxidation states as comma-separated Element=state assignments.",
)
@click.option(
    "--no-charge-neutral",
    "no_formula_charge_neutral",
    is_flag=True,
    help="Disable charge-neutral filtering for formula sampling.",
)
@click.option(
    "--diagnose",
    "formula_diagnose",
    default=0,
    type=int,
    help="Print N sampled buildcell inputs and exit.",
)
@click.option(
    "--volume-minsep-source",
    default="none",
    show_default=True,
    type=click.Choice(["none", "dataset", "baseline", "reference"]),
    help="Auto-configure #VARVOL, #MINSEP, and #NFORM from a volume/minsep estimate.",
)
@click.option(
    "--volume-minsep-dataset",
    default=None,
    type=click.Path(exists=True),
    help="Curated minsep/volume dataset JSON for --volume-minsep-source dataset.",
)
@click.option(
    "--volume-minsep-bundle",
    default=None,
    type=click.Path(exists=True),
    help="Baseline bundle JSON for --volume-minsep-source baseline.",
)
@click.option(
    "--reference-structure",
    "reference_structures",
    multiple=True,
    type=click.Path(exists=True),
    help="Reference structure for --volume-minsep-source reference. Repeat as needed.",
)
@click.option(
    "--volume-scale",
    default=1.1,
    type=float,
    show_default=True,
    help="Scale estimated formula-unit volume before writing #VARVOL.",
)
@click.option(
    "--minsep-scale-low",
    default=0.9,
    type=float,
    show_default=True,
    help="Lower scale for generated #MINSEP pair ranges.",
)
@click.option(
    "--minsep-scale-high",
    default=1.1,
    type=float,
    show_default=True,
    help="Upper scale for generated #MINSEP pair ranges.",
)
@click.option(
    "--max-atoms",
    default=80,
    type=int,
    show_default=True,
    help="Maximum atoms per generated structure for automatic #NFORM.",
)
@click.option(
    "--max-nform",
    default=8,
    type=int,
    show_default=True,
    help="Maximum formula-unit multiplier for automatic #NFORM.",
)
@click.option("--prune", is_flag=True, help="Enable post-relax RSS pruning.")
@click.option(
    "--prune-pool-size",
    default=100,
    type=int,
    show_default=True,
    help="Raw relaxed candidates per pruning pool.",
)
@click.option(
    "--prune-keep-fraction",
    default=0.1,
    type=float,
    show_default=True,
    help="Fraction of each pruning pool to keep.",
)
@click.option(
    "--prune-dedup-tol",
    default=0.1,
    type=float,
    show_default=True,
    help="Fingerprint duplicate tolerance for pruning.",
)
@click.option(
    "--prune-fingerprint-cutoff",
    default=5.0,
    type=float,
    show_default=True,
    help="Fingerprint neighbour cutoff in Angstrom for pruning.",
)
@click.option(
    "--prune-stable-window",
    default=20,
    type=int,
    show_default=True,
    help="Statistics history window for stable pruning flush.",
)
@click.option(
    "--prune-mean-abs-tol",
    default=1e-3,
    type=float,
    show_default=True,
    help="Mean energy-per-atom stability tolerance.",
)
@click.option(
    "--prune-median-abs-tol",
    default=1e-3,
    type=float,
    show_default=True,
    help="Median energy-per-atom stability tolerance.",
)
@click.option(
    "--prune-zweight",
    is_flag=True,
    help="Use element-weighted fingerprint distances when pruning.",
)
@click.option(
    "--prune-keep-rejected",
    is_flag=True,
    help="Leave rejected .res files in place instead of moving to pruned/.",
)
@click.option("--potcar-dir", default=None, help="VASP POTCAR library directory.")
@click.option(
    "--potcar-map",
    "potcar_map_values",
    multiple=True,
    help="VASP POTCAR mapping as Element=symbol. Repeat as needed.",
)
@click.option(
    "--cell-axis-map",
    default=None,
    help="ABACUS-only cell axis permutation, e.g. z:x to make old z become new x.",
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
    formulas,
    formula_elements,
    formula_max_coeff,
    formula_target_volumes,
    formula_oxidation_states,
    no_formula_charge_neutral,
    formula_diagnose,
    volume_minsep_source,
    volume_minsep_dataset,
    volume_minsep_bundle,
    reference_structures,
    volume_scale,
    minsep_scale_low,
    minsep_scale_high,
    max_atoms,
    max_nform,
    prune,
    prune_pool_size,
    prune_keep_fraction,
    prune_dedup_tol,
    prune_fingerprint_cutoff,
    prune_stable_window,
    prune_mean_abs_tol,
    prune_median_abs_tol,
    prune_zweight,
    prune_keep_rejected,
    potcar_dir,
    potcar_map_values,
    cell_axis_map,
):
    """Run an AIRSS random structure search locally."""
    from airsspy.jf.runners import run_buildcell
    from airsspy.search import (
        DEFAULT_ESTIMATE_REMOVE_DIRECTIVES,
        FormulaSamplingOptions,
        RssPruneOptions,
        build_formula_sampling_context,
        candidate_from_res,
        make_seed_text_transform,
        parse_key_float,
        pool_statistics,
        remove_buildcell_directives,
        select_pruned_candidates,
        should_flush_prune_pool,
        validate_prune_options,
    )
    from airsspy.volume_minsep import (
        apply_minsep_headroom,
        build_seed_text_from_estimate,
        lookup_exact_volume_minsep_estimate,
        predict_baseline_volume_minsep_estimate,
        reference_volume_minsep_estimate,
        resolve_nform,
        validate_estimate_options,
    )

    if prune and build_only:
        raise click.ClickException("--prune cannot be used with --build-only")

    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    # Read seed cell content
    seed_cell = Path(seed + ".cell")
    if not seed_cell.exists():
        raise click.ClickException(f"Seed cell file not found: {seed_cell}")
    seed_content = seed_cell.read_text()

    seed_text_transform = None
    use_volume_minsep = volume_minsep_source != "none"
    if use_volume_minsep:
        try:
            validate_estimate_options(
                volume_scale=volume_scale,
                minsep_scale_low=minsep_scale_low,
                minsep_scale_high=minsep_scale_high,
                max_atoms=max_atoms,
                max_nform=max_nform,
            )
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc
        if not formulas and not formula_elements:
            raise click.ClickException(
                "--volume-minsep-source requires --formula or --elements"
            )
        if volume_minsep_source == "dataset" and volume_minsep_dataset is None:
            raise click.ClickException(
                "--volume-minsep-source dataset requires --volume-minsep-dataset"
            )
        if volume_minsep_source == "baseline" and volume_minsep_bundle is None:
            raise click.ClickException(
                "--volume-minsep-source baseline requires --volume-minsep-bundle"
            )
        if volume_minsep_source == "reference" and not reference_structures:
            raise click.ClickException(
                "--volume-minsep-source reference requires --reference-structure"
            )
        volume_minsep_dataset = _resolve_optional_path(volume_minsep_dataset)
        volume_minsep_bundle = _resolve_optional_path(volume_minsep_bundle)
        reference_structures = tuple(
            _resolve_optional_path(path) for path in reference_structures
        )

    if formulas or formula_elements:
        elements = [
            item.strip() for item in formula_elements.split(",") if item.strip()
        ]
        target_volumes = _parse_formula_key_values(
            formula_target_volumes,
            parse_key_float,
            "--target-volume",
        )
        oxidation_states = _parse_oxidation_state_options(formula_oxidation_states)
        seed_text_for_formula_filter = (
            remove_buildcell_directives(seed_content, DEFAULT_ESTIMATE_REMOVE_DIRECTIVES)
            if use_volume_minsep
            else seed_content
        )
        try:
            formula_context = build_formula_sampling_context(
                FormulaSamplingOptions(
                    formulas=formulas,
                    elements=elements,
                    max_coeff=formula_max_coeff,
                    target_atom_volumes=target_volumes,
                    oxidation_states=oxidation_states,
                    require_charge_neutral=not no_formula_charge_neutral,
                ),
                seed_text=seed_text_for_formula_filter,
            )
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc

        reference_atoms_by_formula = None
        if volume_minsep_source == "reference":
            from ase.io import read

            from airsspy.volume_minsep import normalize_formula

            reference_atoms_by_formula = {}
            for path in reference_structures:
                atoms = read(path)
                formula = normalize_formula(atoms.get_chemical_formula())
                reference_atoms_by_formula.setdefault(formula, []).append(atoms)

        def resolve_estimate_for_formula(formula: str):
            if not use_volume_minsep:
                return None
            try:
                if volume_minsep_source == "dataset":
                    return lookup_exact_volume_minsep_estimate(
                        formula=formula,
                        dataset_path=volume_minsep_dataset,
                    )
                if volume_minsep_source == "baseline":
                    return predict_baseline_volume_minsep_estimate(
                        formula=formula,
                        bundle_path=volume_minsep_bundle,
                    )
                if volume_minsep_source == "reference":
                    from airsspy.volume_minsep import normalize_formula

                    reduced_formula = normalize_formula(formula)
                    references = (reference_atoms_by_formula or {}).get(
                        reduced_formula,
                        [],
                    )
                    return reference_volume_minsep_estimate(
                        formula=formula,
                        references=references,
                    )
            except ValueError as exc:
                raise click.ClickException(str(exc)) from exc
            raise click.ClickException(
                f"Unsupported volume/minsep source: {volume_minsep_source}"
            )

        def build_sample_seed(seed_text: str):
            sampled_seed, formula, varvol = formula_context.sample(seed_text)
            estimate = resolve_estimate_for_formula(formula)
            nform = None
            if estimate is not None:
                try:
                    nform = resolve_nform(
                        atoms_per_formula_unit=estimate.atoms_per_formula_unit,
                        max_atoms=max_atoms,
                        max_nform=max_nform,
                    )
                except ValueError as exc:
                    raise click.ClickException(str(exc)) from exc
                sampled_seed = build_seed_text_from_estimate(
                    seed_text,
                    estimate=estimate,
                    volume_scale=volume_scale,
                    minsep_scale_low=minsep_scale_low,
                    minsep_scale_high=minsep_scale_high,
                    nform=nform,
                )
                varvol = estimate.total_volume * volume_scale
            return sampled_seed, formula, varvol, estimate, nform

        if formula_diagnose > 0:
            for i in range(formula_diagnose):
                sampled_seed, formula, varvol, estimate, nform = build_sample_seed(
                    seed_content
                )
                click.echo(f"Sample {i + 1}")
                click.echo("----------------------------------------")
                click.echo("User settings")
                click.echo(f"formula = {formula}")
                if varvol is not None:
                    click.echo(f"varvol = {varvol:g}")
                if estimate is not None:
                    minsep_ranges = apply_minsep_headroom(
                        estimate.canonical_minsep,
                        low_scale=minsep_scale_low,
                        high_scale=minsep_scale_high,
                    )
                    click.echo(f"volume_minsep_source = {volume_minsep_source}")
                    click.echo(f"reduced_formula = {estimate.reduced_formula}")
                    click.echo(f"volume_per_atom = {estimate.volume_per_atom:g}")
                    click.echo(f"total_volume = {estimate.total_volume:g}")
                    click.echo(f"nform = {nform}")
                    click.echo(
                        "minsep = "
                        + ", ".join(
                            f"{pair}={low:g}-{high:g}"
                            for pair, (low, high) in sorted(minsep_ranges.items())
                        )
                    )
                click.echo(f"seedfile = {seed_cell}")
                click.echo("----------------------------------------")
                click.echo("buildcell input")
                click.echo("----------------------------------------")
                click.echo(sampled_seed)
                click.echo("----------------------------------------")
            return

        if use_volume_minsep:

            def seed_text_transform(text: str) -> str:
                sampled_seed, _, _, _, _ = build_sample_seed(text)
                return sampled_seed

        else:
            seed_text_transform = make_seed_text_transform(formula_context)
    elif formula_diagnose > 0:
        raise click.ClickException("--diagnose requires --formula or --elements")

    prune_options = None
    if prune:
        prune_options = RssPruneOptions(
            enabled=True,
            pool_size=prune_pool_size,
            keep_fraction=prune_keep_fraction,
            dedup_tol=prune_dedup_tol,
            fingerprint_cutoff=prune_fingerprint_cutoff,
            stable_window=prune_stable_window,
            mean_abs_tol=prune_mean_abs_tol,
            median_abs_tol=prune_median_abs_tol,
            zweight=prune_zweight,
        )
        try:
            validate_prune_options(prune_options)
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc

    potcar_dir = _resolve_optional_path(potcar_dir)
    potcar_map = _parse_potcar_map_options(potcar_map_values)

    # Copy param file to workdir so runners can find it after chdir
    if not build_only:
        param_suffix = SUFFIX_MAP[code]
        param_file = Path(seed + param_suffix)
        if not param_file.exists():
            raise click.ClickException(f"Param file not found: {param_file}")
        if Path(workdir).resolve() != Path().resolve():
            shutil.copy2(param_file, workdir / param_file.name)
            if code == "vasp":
                kpoints = Path(seed + ".KPOINTS")
                if kpoints.exists():
                    shutil.copy2(kpoints, workdir / kpoints.name)

    if exe is None:
        exe = EXE_DEFAULTS[code]

    # Detect scheduler for walltime awareness
    sched = _get_scheduler_for_walltime()

    orig_dir = os.getcwd()
    os.chdir(workdir)

    try:
        n_built = 0
        n_relaxed = 0
        n_failed = 0
        kept_res_files: list[Path] = []
        prune_pool = []
        prune_stats = []

        def flush_prune_pool(reason: str) -> None:
            nonlocal prune_pool, prune_stats, kept_res_files
            if prune_options is None or not prune_pool:
                return
            kept, rejected = select_pruned_candidates(prune_pool, prune_options)
            kept_res_files.extend(
                candidate.res_path
                for candidate in kept
                if candidate.res_path is not None and candidate.res_path.exists()
            )
            if not prune_keep_rejected:
                for candidate in rejected:
                    _move_pruned_file(candidate, workdir)
            logger.info(
                "Prune flush (%s): kept %d, rejected %d",
                reason,
                len(kept),
                len(rejected),
            )
            prune_pool = []
            prune_stats = []

        logger.info(
            "Starting AIRSS search: seed=%s, nmax=%d, code=%s", seed, nmax, code
        )
        logger.info("Working directory: %s", workdir)

        for i in range(1, nmax + 1):
            # Check stop file
            if _check_stop_file(workdir):
                logger.info("Stop file detected. Gracefully terminating.")
                break

            # Check walltime
            if not _walltime_remaining_ok(sched, walltime_buffer):
                break

            # Build
            build_result = run_buildcell(
                seed_name=seed,
                seed_content=seed_content,
                build_timeout=build_timeout,
                write_seed=True,
                seed_text_transform=seed_text_transform,
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
                code,
                exe,
                max_iterations,
                cluster,
                pressure,
                mpinp,
                potcar_dir=potcar_dir,
                potcar_map=potcar_map,
                cell_axis_map=cell_axis_map,
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
                elif code == "vasp":
                    struct_content = Path(struct_name + ".cell").read_text()
                    incar_content = Path(seed + param_suffix).read_text()
                    kpoints_path = Path(seed + ".KPOINTS")
                    rc = runner.run(
                        struct_name,
                        struct_content,
                        incar_content,
                        kpoints_path=kpoints_path if kpoints_path.exists() else None,
                    )

                if rc == 0:
                    _collect_result(
                        struct_name,
                        code,
                        calculator_spec=runner if code == "vasp" else None,
                    )
                    n_relaxed += 1
                    if prune_options is not None:
                        candidate = candidate_from_res(Path(struct_name + ".res"))
                        prune_pool.append(candidate)
                        prune_stats.append(pool_statistics(prune_pool))
                        reason, should_flush = should_flush_prune_pool(
                            prune_pool,
                            prune_stats,
                            prune_options,
                        )
                        if should_flush:
                            flush_prune_pool(reason)
                    logger.info("[%d] Relaxed OK: %s", i, struct_name)
                else:
                    try:
                        _collect_result(
                            struct_name,
                            code,
                            calculator_spec=runner if code == "vasp" else None,
                        )
                        if prune_options is not None:
                            candidate = candidate_from_res(Path(struct_name + ".res"))
                            prune_pool.append(candidate)
                            prune_stats.append(pool_statistics(prune_pool))
                            reason, should_flush = should_flush_prune_pool(
                                prune_pool,
                                prune_stats,
                                prune_options,
                            )
                            if should_flush:
                                flush_prune_pool(reason)
                        logger.info("[%d] Not converged: %s", i, struct_name)
                    except Exception:
                        logger.info("[%d] Relax FAILED: %s", i, struct_name)
                        _emit_diagnostics(struct_name, code)
                        if not keep:
                            runner.clean_failed(struct_name)
                    n_failed += 1

            except Exception:
                logger.error("[%d] Relax crashed: %s", i, struct_name, exc_info=True)
                _emit_diagnostics(struct_name, code)
                if not keep:
                    runner.clean_failed(struct_name)
                n_failed += 1

        # Summary
        logger.info(
            "Search complete: %d built, %d relaxed, %d failed",
            n_built,
            n_relaxed,
            n_failed,
        )

        if prune_options is not None:
            flush_prune_pool("final")

        if pack and n_relaxed > 0:
            pack_files = kept_res_files if prune_options is not None else None
            packed = _pack_res_files(workdir, files=pack_files)
            logger.info("Packed .res files into %s", packed)

    finally:
        os.chdir(orig_dir)


@run.command("crud")
@click.option(
    "--code",
    default="castep",
    show_default=True,
    type=click.Choice(["castep", "gulp", "pp3", "abacus", "vasp", "ml"]),
    help="DFT code or ML mode to use",
)
@click.option("--exe", default=None, help="Relaxation executable (default: auto)")
@click.option(
    "--workdir",
    default=".",
    show_default=True,
    type=click.Path(),
    help="Queue root containing hopper/",
)
@click.option(
    "--num",
    default=1000,
    type=int,
    show_default=True,
    help="Max number of queued files to inspect when claiming",
)
@click.option("--nostop", is_flag=True, help="Keep polling when no jobs are available")
@click.option("--cycle", is_flag=True, help="Retry electronic minimisation failures")
@click.option("--keep", is_flag=True, help="Keep intermediate files")
@click.option(
    "--singlepoint",
    "--single-point",
    is_flag=True,
    help="Run a single-point calculation instead of a full relaxation",
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
    help=(
        "Model spec for --code ml. Default backend is torch-sim, e.g. "
        "'mace:medium'. Use 'ase:mace:medium' or 'ase:module:Class@model' "
        "for ASE fallback."
    ),
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
@click.option(
    "--device",
    default=None,
    help="Torch device for torch-sim ML runs, e.g. cuda or cpu.",
)
@click.option(
    "--batch-size",
    default=1,
    type=click.IntRange(min=1),
    show_default=True,
    help="Number of structures per torch-sim ML batch.",
)
@click.option("--potcar-dir", default=None, help="VASP POTCAR library directory.")
@click.option(
    "--potcar-map",
    "potcar_map_values",
    multiple=True,
    help="VASP POTCAR mapping as Element=symbol. Repeat as needed.",
)
@click.option(
    "--cell-axis-map",
    default=None,
    help="ABACUS-only cell axis permutation, e.g. z:x to make old z become new x.",
)
def run_crud(
    code,
    exe,
    workdir,
    num,
    nostop,
    cycle,
    keep,
    singlepoint,
    pressure,
    max_iterations,
    cluster,
    mpinp,
    walltime_buffer,
    calculator_spec,
    optimizer,
    fmax,
    device,
    batch_size,
    potcar_dir,
    potcar_map_values,
    cell_axis_map,
):
    """Consume hopper/*.res jobs locally, like crud.pl."""
    if code == "ml" and not calculator_spec:
        raise click.ClickException("--calculator is required when --code ml")
    if singlepoint and code not in ("castep", "abacus", "vasp", "ml"):
        raise click.ClickException(f"Single-point not supported for code: {code}")
    use_torchsim = code == "ml" and _is_torchsim_model(calculator_spec)
    if use_torchsim:
        _ensure_torchsim_available()

    workdir = Path(workdir).resolve()
    hopper = workdir / "hopper"
    if not hopper.exists():
        raise click.ClickException(f"Hopper directory not found: {hopper}")
    workdir.mkdir(parents=True, exist_ok=True)

    if exe is None:
        exe = EXE_DEFAULTS.get(code, "")

    potcar_dir = _resolve_optional_path(potcar_dir)
    potcar_map = _parse_potcar_map_options(potcar_map_values)
    sched = _get_scheduler_for_walltime()
    runner = None
    if not use_torchsim:
        runner = _create_task_runner(
            code,
            exe,
            max_iterations,
            cluster,
            pressure,
            mpinp,
            calculator_spec=calculator_spec,
            optimizer=optimizer,
            fmax=fmax,
            potcar_dir=potcar_dir,
            potcar_map=potcar_map,
            cell_axis_map=cell_axis_map,
            singlepoint=singlepoint,
        )
    torchsim_runner = None

    orig_dir = os.getcwd()
    os.chdir(workdir)

    try:
        n_done = 0
        n_failed = 0
        n_requeued = 0

        logger.info("Starting CRUD worker: workdir=%s, code=%s", workdir, code)

        def handle_failure(seed: str) -> bool:
            logger.error("CRUD failed: %s", seed, exc_info=True)
            _emit_diagnostics(seed, code)
            if not singlepoint and cycle and _crud_retryable_failure(seed):
                if _requeue_crud_job(seed, workdir):
                    _cleanup_crud_artifacts(seed)
                    logger.info("CRUD requeued: %s", seed)
                    return True
            _move_crud_artifacts(seed, workdir / "bad_castep", keep=True)
            return False

        def finalize_success(seed: str, rc: int) -> None:
            _move_crud_artifacts(seed, workdir / "good_castep", keep)
            if rc == 0:
                logger.info("CRUD OK: %s", seed)
            else:
                logger.info("CRUD not converged but collected: %s", seed)

        while True:
            if Path("STOP_CRUD").exists():
                logger.info("STOP_CRUD detected. Gracefully terminating.")
                break

            if not _walltime_remaining_ok(sched, walltime_buffer):
                break

            claimed = _claim_crud_job(workdir, num)
            if claimed is None:
                if not nostop:
                    break
                import time

                time.sleep(random.random())
                continue

            seed = claimed.stem

            try:
                _, param_suffix = _prepare_crud_inputs(seed, code)
                if use_torchsim:
                    claimed_seeds = [seed]
                    while len(claimed_seeds) < batch_size:
                        extra = _claim_crud_job(workdir, num)
                        if extra is None:
                            break
                        extra_seed = extra.stem
                        try:
                            _prepare_crud_inputs(extra_seed, code)
                        except Exception:
                            if handle_failure(extra_seed):
                                n_requeued += 1
                            else:
                                n_failed += 1
                            continue
                        claimed_seeds.append(extra_seed)

                    if torchsim_runner is None:
                        from airsspy.jf.ml_runners import TorchSimRunner

                        torchsim_runner = TorchSimRunner(
                            calculator_spec, device=device
                        )
                    structures = [
                        _read_res_as_atoms(Path(sname + ".res"))
                        for sname in claimed_seeds
                    ]
                    try:
                        done, failed, collected = _run_torchsim_batch(
                            torchsim_runner,
                            claimed_seeds,
                            structures,
                            code,
                            calculator_spec,
                            singlepoint=singlepoint,
                            max_iterations=max_iterations,
                            fmax=fmax,
                            optimizer=optimizer,
                            pressure=pressure,
                        )
                    except Exception:
                        logger.error(
                            "CRUD TorchSim batch failed for structures %s",
                            ", ".join(claimed_seeds),
                            exc_info=True,
                        )
                        for sname in claimed_seeds:
                            _move_crud_artifacts(
                                sname, workdir / "bad_castep", keep=True
                            )
                        n_failed += len(claimed_seeds)
                        continue

                    collected_names = {path.stem for path in collected}
                    for sname in claimed_seeds:
                        if sname in collected_names:
                            finalize_success(sname, 0)
                        else:
                            _move_crud_artifacts(
                                sname, workdir / "bad_castep", keep=True
                            )
                    n_done += done
                    n_failed += failed
                    continue

                rc = _run_crud_one(
                    seed,
                    code,
                    runner,
                    param_suffix,
                    workdir,
                    singlepoint=singlepoint,
                )

                try:
                    _collect_result(
                        seed,
                        code,
                        calculator_spec=runner if code == "vasp" else calculator_spec,
                    )
                except Exception as exc:
                    if rc == 0:
                        raise
                    raise RuntimeError(
                        "calculation did not produce collectable output"
                    ) from exc

                finalize_success(seed, rc)
                n_done += 1

            except Exception:
                if handle_failure(seed):
                    n_requeued += 1
                else:
                    n_failed += 1

        logger.info(
            "CRUD complete: %d collected, %d failed, %d requeued",
            n_done,
            n_failed,
            n_requeued,
        )

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
    type=click.Choice(["castep", "gulp", "pp3", "abacus", "vasp", "ml"]),
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
    "--singlepoint",
    "--single-point",
    is_flag=True,
    help="Run a single-point calculation instead of a full relaxation",
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
    help=(
        "Model spec for --code ml. Default backend is torch-sim, e.g. "
        "'mace:medium'. Use 'ase:mace:medium' or 'ase:module:Class@model' "
        "for ASE fallback."
    ),
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
@click.option(
    "--device",
    default=None,
    help="Torch device for torch-sim ML runs, e.g. cuda or cpu.",
)
@click.option(
    "--batch-size",
    default=1,
    type=click.IntRange(min=1),
    show_default=True,
    help="Number of structures per torch-sim ML batch.",
)
@click.option(
    "--debug",
    is_flag=True,
    help="Print detailed run-relax debug logging.",
)
@click.option("--potcar-dir", default=None, help="VASP POTCAR library directory.")
@click.option(
    "--potcar-map",
    "potcar_map_values",
    multiple=True,
    help="VASP POTCAR mapping as Element=symbol. Repeat as needed.",
)
@click.option(
    "--cell-axis-map",
    default=None,
    help="ABACUS-only cell axis permutation, e.g. z:x to make old z become new x.",
)
def run_relax(
    cell,
    seed,
    code,
    exe,
    workdir,
    keep,
    pack,
    singlepoint,
    pressure,
    max_iterations,
    cluster,
    mpinp,
    walltime_buffer,
    calculator_spec,
    optimizer,
    fmax,
    device,
    batch_size,
    debug,
    potcar_dir,
    potcar_map_values,
    cell_axis_map,
):
    """Relax existing cell files locally."""
    if debug:
        logging.getLogger().setLevel(logging.DEBUG)
        for handler in logging.getLogger().handlers:
            handler.setLevel(logging.DEBUG)

    if code == "ml" and not calculator_spec:
        raise click.ClickException("--calculator is required when --code ml")
    if singlepoint and code not in ("castep", "abacus", "vasp", "ml"):
        raise click.ClickException(f"Single-point not supported for code: {code}")
    potcar_dir = _resolve_optional_path(potcar_dir)
    potcar_map = _parse_potcar_map_options(potcar_map_values)

    workdir = Path(workdir).resolve()
    cell_files = sorted(workdir.glob(cell))
    if not cell_files:
        raise click.ClickException(f"No files matched pattern: {cell}")
    cell_files = _filter_packed_res_inputs(cell_files)
    if not cell_files:
        raise click.ClickException(f"No single-structure inputs matched pattern: {cell}")
    use_torchsim = code == "ml" and _is_torchsim_model(calculator_spec)
    if use_torchsim:
        _ensure_torchsim_available()
    logger.debug(
        "run relax resolved %d single-structure input(s): %s",
        len(cell_files),
        ", ".join(path.name for path in cell_files[:10])
        + (" ..." if len(cell_files) > 10 else ""),
    )
    relax_inputs = [
        (
            input_path,
            *_prepare_relax_input(
                input_path,
                workdir,
                write_res_cell=code == "castep",
                convert_res_cell=code != "ml",
            ),
        )
        for input_path in cell_files
    ]

    # Read param file (not needed for ML)
    param_suffix = None
    if code != "ml":
        param_suffix = SUFFIX_MAP[code]
        for input_path, cell_path, _, _ in relax_inputs:
            _resolve_relax_param_file(
                input_path, cell_path, seed, workdir, param_suffix
            )

    if exe is None:
        exe = EXE_DEFAULTS.get(code, "")

    # Detect scheduler
    sched = _get_scheduler_for_walltime()

    runner = None
    if code != "ml" or not use_torchsim:
        logger.debug(
            "Creating %s runner for code=%s exe=%s backend=%s",
            "single-point" if singlepoint else "relax",
            code,
            exe,
            "ase" if code == "ml" else "external",
        )
        runner = _create_task_runner(
            code,
            exe,
            max_iterations,
            cluster,
            pressure,
            mpinp,
            calculator_spec=calculator_spec,
            optimizer=optimizer,
            fmax=fmax,
            potcar_dir=potcar_dir,
            potcar_map=potcar_map,
            cell_axis_map=cell_axis_map,
            singlepoint=singlepoint,
        )

    orig_dir = os.getcwd()
    os.chdir(workdir)

    try:
        n_relaxed = 0
        n_failed = 0
        collected_res_files: list[Path] = []
        total = len(relax_inputs)

        # TorchSim batch path: process all structures at once on GPU
        if use_torchsim:
            from airsspy.jf.ml_runners import TorchSimRunner

            struct_names = [struct_name for _, _, struct_name, _ in relax_inputs]
            structures = [
                _prepare_ml_structure_input(input_path, cell_content)
                for input_path, _, _, cell_content in relax_inputs
            ]
            logger.debug(
                "Prepared ML structures: total_atoms=%d, first_batch=%s",
                sum(len(structure) for structure in structures),
                ", ".join(struct_names[:batch_size]),
            )

            mode_label = "SP" if singlepoint else "relaxing"
            logger.info(
                "TorchSim batch %s %d structures with %s (batch size %d)",
                mode_label,
                total,
                calculator_spec,
                batch_size,
            )
            torchsim_runner = TorchSimRunner(calculator_spec, device=device)
            for start in range(0, total, batch_size):
                end = min(start + batch_size, total)
                batch_names = struct_names[start:end]
                batch_structures = structures[start:end]
                logger.info(
                    "TorchSim batch [%d-%d/%d]",
                    start + 1,
                    end,
                    total,
                )
                try:
                    done, failed, collected = _run_torchsim_batch(
                        torchsim_runner,
                        batch_names,
                        batch_structures,
                        code,
                        calculator_spec,
                        singlepoint=singlepoint,
                        max_iterations=max_iterations,
                        fmax=fmax,
                        optimizer=optimizer,
                        pressure=pressure,
                    )
                    collected_res_files.extend(workdir / path.name for path in collected)
                    n_relaxed += done
                    n_failed += failed
                except Exception:
                    logger.error(
                        "TorchSim batch failed for structures %s",
                        ", ".join(batch_names),
                        exc_info=True,
                    )
                    n_failed += len(batch_names)
                    continue

        else:
            # Standard one-by-one loop (DFT codes or ASE fallback)
            backend_label = "ASE ML" if code == "ml" else code
            logger.info(
                "%s %d structures with %s",
                "Running single-point on" if singlepoint else "Relaxing",
                total,
                backend_label,
            )

            for i, (input_path, cell_path, struct_name, cell_content) in enumerate(
                relax_inputs, 1
            ):

                # Check walltime
                if not _walltime_remaining_ok(sched, walltime_buffer):
                    break

                try:
                    rc = _run_local_structure_one(
                        input_path,
                        cell_path,
                        struct_name,
                        cell_content,
                        code,
                        runner,
                        param_suffix,
                        seed,
                        workdir,
                        singlepoint=singlepoint,
                    )

                    if rc == 0:
                        _collect_result(
                            struct_name,
                            code,
                            calculator_spec=runner if code == "vasp" else calculator_spec,
                        )
                        collected_res_files.append(workdir / f"{struct_name}.res")
                        n_relaxed += 1
                        logger.info("[%d/%d] OK: %s", i, total, struct_name)
                    else:
                        try:
                            _collect_result(
                                struct_name,
                                code,
                                calculator_spec=runner
                                if code == "vasp"
                                else calculator_spec,
                            )
                            collected_res_files.append(workdir / f"{struct_name}.res")
                            logger.info(
                                "[%d/%d] Not converged: %s", i, total, struct_name
                            )
                        except Exception:
                            logger.info("[%d/%d] FAILED: %s", i, total, struct_name)
                            _emit_diagnostics(struct_name, code)
                            if not keep:
                                runner.clean_failed(struct_name)
                        n_failed += 1

                except Exception:
                    logger.error(
                        "[%d/%d] %s crashed: %s",
                        i,
                        total,
                        "SP" if singlepoint else "Relax",
                        struct_name,
                        exc_info=True,
                    )
                    _emit_diagnostics(struct_name, code)
                    if not keep:
                        runner.clean_failed(struct_name)
                    n_failed += 1

        logger.info(
            "%s complete: %d/%d succeeded, %d failed",
            "SP" if singlepoint else "Relaxation",
            n_relaxed,
            total,
            n_failed,
        )

        if pack and n_relaxed > 0:
            packed = _pack_res_files(workdir, files=collected_res_files)
            logger.info("Packed into %s", packed)

        if use_torchsim and n_relaxed == 0 and n_failed > 0:
            raise click.ClickException(
                f"TorchSim {'single-point' if singlepoint else 'relaxation'} "
                "failed for all matched structures; "
                "see verbose log above for the failing batch."
            )

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
    type=click.Choice(["castep", "abacus", "vasp", "ml"]),
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
    help=(
        "Model spec for --code ml. Default backend is torch-sim, e.g. "
        "'mace:medium'. Use 'ase:mace:medium' or 'ase:module:Class@model' "
        "for ASE fallback."
    ),
)
@click.option(
    "--device",
    default=None,
    help="Torch device for torch-sim ML runs, e.g. cuda or cpu.",
)
@click.option(
    "--batch-size",
    default=1,
    type=click.IntRange(min=1),
    show_default=True,
    help="Number of structures per torch-sim ML batch.",
)
@click.option(
    "--pressure",
    default=0.0,
    type=float,
    show_default=True,
    help="External pressure (GPa)",
)
@click.option("--potcar-dir", default=None, help="VASP POTCAR library directory.")
@click.option(
    "--potcar-map",
    "potcar_map_values",
    multiple=True,
    help="VASP POTCAR mapping as Element=symbol. Repeat as needed.",
)
@click.option(
    "--cell-axis-map",
    default=None,
    help="ABACUS-only cell axis permutation, e.g. z:x to make old z become new x.",
)
def run_sp(
    cell,
    seed,
    code,
    exe,
    workdir,
    keep,
    pack,
    walltime_buffer,
    calculator_spec,
    device,
    batch_size,
    pressure,
    potcar_dir,
    potcar_map_values,
    cell_axis_map,
):
    """Run single-point calculations on existing cell files."""
    if code == "ml" and not calculator_spec:
        raise click.ClickException("--calculator is required when --code ml")

    workdir = Path(workdir).resolve()
    cell_files = sorted(workdir.glob(cell))
    if not cell_files:
        raise click.ClickException(f"No files matched pattern: {cell}")
    cell_files = _filter_packed_res_inputs(cell_files)
    if not cell_files:
        raise click.ClickException(f"No single-structure inputs matched pattern: {cell}")
    if code not in ("ml", "vasp") and any(
        path.suffix.lower() == ".res" for path in cell_files
    ):
        raise click.ClickException(
            "RES input for run sp is currently supported only with --code ml or vasp"
        )
    sp_inputs = None
    if code in ("ml", "vasp"):
        sp_inputs = [
            (
                input_path,
                *_prepare_relax_input(
                    input_path,
                    workdir,
                    write_res_cell=code == "vasp",
                    convert_res_cell=input_path.suffix.lower() != ".res"
                    or code == "vasp",
                ),
            )
            for input_path in cell_files
        ]
    use_torchsim = code == "ml" and _is_torchsim_model(calculator_spec)
    if use_torchsim:
        _ensure_torchsim_available()

    # Read param file (not needed for ML)
    param_file_name = None
    param_suffix = None
    if code != "ml":
        param_suffix = SUFFIX_MAP[code]
        if code == "vasp" and sp_inputs is not None:
            for input_path, cell_path, _, _ in sp_inputs:
                _resolve_relax_param_file(
                    input_path, cell_path, seed, workdir, param_suffix
                )
        else:
            param_file = workdir / (seed + param_suffix)
            if not param_file.exists():
                param_file = workdir / (cell_files[0].stem + param_suffix)
            if not param_file.exists():
                raise click.ClickException(f"Param file not found: {param_file}")
            param_file_name = param_file.stem

    if exe is None:
        exe = EXE_DEFAULTS.get(code, "")

    potcar_dir = _resolve_optional_path(potcar_dir)
    potcar_map = _parse_potcar_map_options(potcar_map_values)

    # Detect scheduler
    sched = _get_scheduler_for_walltime()

    runner = None
    if code != "ml" or not use_torchsim:
        runner = _create_sp_runner(
            code,
            exe,
            calculator_spec=calculator_spec,
            pressure=pressure,
            potcar_dir=potcar_dir,
            potcar_map=potcar_map,
            cell_axis_map=cell_axis_map,
        )

    orig_dir = os.getcwd()
    os.chdir(workdir)

    try:
        n_done = 0
        n_failed = 0
        collected_res_files: list[Path] = []
        total = len(sp_inputs) if sp_inputs is not None else len(cell_files)

        # TorchSim batch path
        if use_torchsim:
            from airsspy.jf.ml_runners import TorchSimRunner

            struct_names = [struct_name for _, _, struct_name, _ in sp_inputs]
            structures = [
                _prepare_ml_structure_input(input_path, cell_content)
                for input_path, _, _, cell_content in sp_inputs
            ]

            logger.info(
                "TorchSim batch SP on %d structures with %s (batch size %d)",
                total,
                calculator_spec,
                batch_size,
            )
            torchsim_runner = TorchSimRunner(calculator_spec, device=device)
            for start in range(0, total, batch_size):
                end = min(start + batch_size, total)
                batch_names = struct_names[start:end]
                batch_structures = structures[start:end]
                logger.info(
                    "TorchSim batch SP [%d-%d/%d]",
                    start + 1,
                    end,
                    total,
                )
                try:
                    batch_results = torchsim_runner.static_batch(
                        batch_names,
                        batch_structures,
                        scalar_pressure=pressure,
                    )
                except Exception:
                    logger.error(
                        "TorchSim batch SP failed for structures %s",
                        ", ".join(batch_names),
                        exc_info=True,
                    )
                    n_failed += len(batch_names)
                    continue

                for sname, rc in batch_results.items():
                    if rc == 0:
                        try:
                            _collect_result(
                                sname,
                                code,
                                calculator_spec=calculator_spec,
                            )
                            collected_res_files.append(workdir / f"{sname}.res")
                            n_done += 1
                        except Exception:
                            logger.error(
                                "TorchSim SP result collection failed: %s",
                                sname,
                                exc_info=True,
                            )
                            n_failed += 1
                            _cleanup_ml_transients(sname)
                    else:
                        n_failed += 1
                        _cleanup_ml_transients(sname)

        else:
            logger.info("Running single-point on %d structures with %s", total, code)

            loop_inputs = (
                sp_inputs
                if sp_inputs is not None
                else [(cell_path, cell_path, cell_path.stem, None) for cell_path in cell_files]
            )
            for i, (input_path, cell_path, struct_name, cell_content) in enumerate(
                loop_inputs, 1
            ):

                # Check walltime
                if not _walltime_remaining_ok(sched, walltime_buffer):
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
                        param_content = Path(param_file_name + param_suffix).read_text()
                        rc = runner.run(struct_name, struct_content, param_content)
                    elif code == "vasp":
                        struct_content = cell_content or cell_path.read_text()
                        param_file = _resolve_relax_param_file(
                            input_path, cell_path, seed, workdir, param_suffix
                        )
                        kpoints_path = param_file.with_suffix(".KPOINTS")
                        rc = runner.run(
                            struct_name,
                            struct_content,
                            param_file.read_text(),
                            kpoints_path=kpoints_path if kpoints_path.exists() else None,
                        )
                    elif code == "ml":
                        ml_input = _prepare_ml_structure_input(
                            input_path, cell_content
                        )
                        rc = runner.run(struct_name, ml_input)

                    if rc == 0:
                        _collect_result(
                            struct_name,
                            code,
                            calculator_spec=runner if code == "vasp" else calculator_spec,
                        )
                        collected_res_files.append(workdir / f"{struct_name}.res")
                        n_done += 1
                        logger.info("[%d/%d] OK: %s", i, total, struct_name)
                    else:
                        try:
                            _collect_result(
                                struct_name,
                                code,
                                calculator_spec=runner
                                if code == "vasp"
                                else calculator_spec,
                            )
                            collected_res_files.append(workdir / f"{struct_name}.res")
                            logger.info(
                                "[%d/%d] Not converged: %s",
                                i,
                                total,
                                struct_name,
                            )
                        except Exception:
                            logger.info("[%d/%d] FAILED: %s", i, total, struct_name)
                            _emit_diagnostics(struct_name, code)
                            if not keep:
                                runner.clean_failed(struct_name)
                        n_failed += 1

                except Exception:
                    logger.error(
                        "[%d/%d] SP crashed: %s",
                        i,
                        total,
                        struct_name,
                        exc_info=True,
                    )
                    _emit_diagnostics(struct_name, code)
                    if not keep:
                        runner.clean_failed(struct_name)
                    n_failed += 1

        logger.info(
            "SP complete: %d/%d succeeded, %d failed",
            n_done,
            total,
            n_failed,
        )

        if pack and n_done > 0:
            packed = _pack_res_files(workdir, files=collected_res_files)
            logger.info("Packed into %s", packed)

        if (
            code == "ml" and use_torchsim and n_done == 0 and n_failed > 0
        ):
            raise click.ClickException(
                "TorchSim single-point failed for all matched structures; "
                "see verbose log above for the failing batch."
            )

    finally:
        os.chdir(orig_dir)
