"""Utilities for local VASP AIRSS runs."""

from __future__ import annotations

import gzip
import hashlib
import importlib
import os
import re
import shutil
import tempfile
from pathlib import Path
from typing import Any

CONTROL_INPUT_SET = "AIRSSPY_VASP_INPUT_SET"
CONTROL_PREFIX = "AIRSSPY_"
EV_PER_ANG3_TO_GPA = 160.21766208


def parse_potcar_map(values: tuple[str, ...] | list[str]) -> dict[str, str]:
    """Parse ``Element=symbol`` POTCAR map options."""
    mapping: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"expected Element=symbol, got {value!r}")
        element, symbol = (part.strip() for part in value.split("=", 1))
        if not element or not symbol:
            raise ValueError(f"expected Element=symbol, got {value!r}")
        mapping[element] = symbol
    return mapping


def parse_seed_incar(incar_content: str):
    """Parse seed INCAR content with pymatgen."""
    from pymatgen.io.vasp.inputs import Incar

    return Incar.from_str(incar_content)


def _normalize_control_key(key: str) -> str:
    return key.upper().replace("-", "_")


def split_incar_settings(incar) -> tuple[dict[str, Any], dict[str, Any]]:
    """Split airsspy control keys from VASP INCAR settings."""
    control: dict[str, Any] = {}
    user: dict[str, Any] = {}
    for key, value in dict(incar).items():
        normalized = _normalize_control_key(str(key))
        if normalized.startswith(CONTROL_PREFIX):
            control[normalized] = value
        else:
            user[str(key)] = value
    return control, user


def _coerce_input_set_name(value: Any, default: str) -> str:
    if value is None:
        return default
    if isinstance(value, (list, tuple)) and value:
        value = value[0]
    name = str(value).strip()
    return name or default


def resolve_input_set_class(name: str):
    """Resolve a pymatgen VASP input set class."""
    if ":" in name:
        module_name, class_name = name.split(":", 1)
    elif "." in name:
        module_name, class_name = name.rsplit(".", 1)
    else:
        module_name = "pymatgen.io.vasp.sets"
        class_name = name
    module = importlib.import_module(module_name)
    try:
        cls = getattr(module, class_name)
    except AttributeError as exc:
        for attr_name in dir(module):
            if attr_name.casefold() == class_name.casefold():
                return getattr(module, attr_name)
        raise ValueError(f"Unknown VASP input set: {name}") from exc
    return cls


def structure_from_cell_text(cell_text: str):
    """Convert CASTEP/AIRSS cell text to a pymatgen Structure."""
    from ase import Atoms
    from castepinput.inputs import CellInput
    from pymatgen.io.ase import AseAtomsAdaptor

    with tempfile.NamedTemporaryFile("w", suffix=".cell", delete=False) as handle:
        handle.write(cell_text)
        tmp_path = Path(handle.name)
    try:
        cell = CellInput.from_file(str(tmp_path))
    finally:
        tmp_path.unlink(missing_ok=True)
    elements, positions, _tags = cell.get_positions()
    atoms = Atoms(symbols=elements, positions=positions, cell=cell.get_cell(), pbc=True)
    return AseAtomsAdaptor.get_structure(atoms)


def structure_from_res(path: str | Path):
    """Read a single RES file as a pymatgen Structure."""
    from pymatgen.io.ase import AseAtomsAdaptor

    from .restools import read_res_atoms

    _, atoms = read_res_atoms(Path(path).read_text().splitlines())
    return AseAtomsAdaptor.get_structure(atoms)


def _potcar_root(potcar_dir: str | None = None) -> Path | None:
    if potcar_dir:
        return Path(potcar_dir)
    for env_name in ("AIRSSPY_POTCAR_DIR",):
        env = os.environ.get(env_name)
        if env:
            return Path(env)
    try:
        import pymatgen.core as pmg_core
    except Exception:
        pmg_core = None
    for setting_name in ("PMG_VASP_PSP_DIR", "PMG_VAP_PSP_DIR"):
        value = _pymatgen_setting(pmg_core, setting_name)
        if value:
            return Path(str(value))
    for env_name in ("PMG_VASP_PSP_DIR", "VASP_PSP_DIR"):
        env = os.environ.get(env_name)
        if env:
            return Path(env)
    return None


def _pymatgen_setting(pmg_core, setting_name: str) -> Any:
    """Read pymatgen settings, including .pmgrc.yml compatibility."""
    if pmg_core is None:
        return None
    settings = getattr(pmg_core, "SETTINGS", {})
    value = settings.get(setting_name)
    if value:
        return value

    load_settings = getattr(pmg_core, "_load_pmg_settings", None)
    if not callable(load_settings):
        return None

    candidate_paths = []
    for attr in ("SETTINGS_FILE", "OLD_SETTINGS_FILE"):
        path = getattr(pmg_core, attr, None)
        if path:
            candidate_paths.append(Path(path).with_suffix(".yml"))

    original_config = os.environ.get("PMG_CONFIG_FILE")
    try:
        for path in candidate_paths:
            if not path.is_file():
                continue
            os.environ["PMG_CONFIG_FILE"] = str(path)
            try:
                value = load_settings().get(setting_name)
            except Exception:
                value = None
            if value:
                return value
    finally:
        if original_config is None:
            os.environ.pop("PMG_CONFIG_FILE", None)
        else:
            os.environ["PMG_CONFIG_FILE"] = original_config
    return None


def _unique_species(structure) -> list[str]:
    symbols: list[str] = []
    for site in structure:
        symbol = str(site.specie.symbol)
        if symbol not in symbols:
            symbols.append(symbol)
    return symbols


def _input_set_potcar_symbol_map(input_set, fallback_structure) -> dict[str, str]:
    """Return element-to-POTCAR symbols resolved by a pymatgen input set."""
    poscar = getattr(input_set, "poscar", None)
    structure = _input_set_potcar_structure(input_set, fallback_structure)
    elements = list(getattr(poscar, "site_symbols", None) or _unique_species(structure))
    potcar_symbols = getattr(input_set, "potcar_symbols", None)
    if callable(potcar_symbols):
        potcar_symbols = potcar_symbols()
    if potcar_symbols is None or len(potcar_symbols) != len(elements):
        return {}
    return {
        str(element): str(symbol)
        for element, symbol in zip(elements, potcar_symbols)
        if str(element) and str(symbol)
    }


def _find_potcar_file(root: Path, symbol: str) -> Path:
    candidates = [
        root / symbol / "POTCAR",
        root / symbol / "POTCAR.gz",
        root / f"{symbol}.POTCAR",
        root / f"{symbol}.POTCAR.gz",
        root / f"POTCAR.{symbol}",
        root / f"POTCAR.{symbol}.gz",
        root / "potpaw_PBE" / symbol / "POTCAR",
        root / "potpaw_PBE" / symbol / "POTCAR.gz",
        root / "potpaw_PBE" / f"{symbol}.POTCAR",
        root / "potpaw_PBE" / f"{symbol}.POTCAR.gz",
        root / "potpaw_PBE" / f"POTCAR.{symbol}",
        root / "potpaw_PBE" / f"POTCAR.{symbol}.gz",
        root / "POT_GGA_PAW_PBE" / symbol / "POTCAR",
        root / "POT_GGA_PAW_PBE" / symbol / "POTCAR.gz",
        root / "POT_GGA_PAW_PBE" / f"{symbol}.POTCAR",
        root / "POT_GGA_PAW_PBE" / f"{symbol}.POTCAR.gz",
        root / "POT_GGA_PAW_PBE" / f"POTCAR.{symbol}",
        root / "POT_GGA_PAW_PBE" / f"POTCAR.{symbol}.gz",
        root / "POT_PAW_PBE_64" / f"POTCAR.{symbol}",
        root / "POT_PAW_PBE_64" / f"POTCAR.{symbol}.gz",
        root / "POT_GGA_PAW_PBE_54" / f"POTCAR.{symbol}",
        root / "POT_GGA_PAW_PBE_54" / f"POTCAR.{symbol}.gz",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"POTCAR for {symbol!r} not found in {root}; looked for "
        + ", ".join(str(path) for path in candidates)
    )


def _sha256(path: Path) -> str:
    hsh = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hsh.update(chunk)
    return hsh.hexdigest()


def _read_potcar_bytes(path: Path) -> bytes:
    """Read a POTCAR file, transparently handling pymatgen-style gzip files."""
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            return handle.read()
    return path.read_bytes()


def assemble_potcar(
    structure,
    output_path: str | Path,
    *,
    potcar_dir: str | None = None,
    potcar_map: dict[str, str] | None = None,
    resolved_symbols: dict[str, str] | None = None,
) -> list[dict[str, str]]:
    """Assemble a concatenated POTCAR from a simple local potential library."""
    root = _potcar_root(potcar_dir)
    if root is None:
        raise FileNotFoundError("No POTCAR directory supplied")
    potcar_map = potcar_map or {}
    resolved_symbols = resolved_symbols or {}
    metadata: list[dict[str, str]] = []
    with open(output_path, "wb") as out:
        for element in _unique_species(structure):
            symbol = resolved_symbols.get(element, potcar_map.get(element, element))
            source = _find_potcar_file(root, symbol)
            content = _read_potcar_bytes(source)
            out.write(content)
            if not content.endswith(b"\n"):
                out.write(b"\n")
            metadata.append(
                {
                    "element": element,
                    "symbol": symbol,
                    "source": str(source),
                    "sha256": _sha256(source),
                }
            )
    return metadata


def _input_set_potcar_structure(input_set, fallback_structure):
    """Return the structure ordering used by the generated POSCAR."""
    poscar = getattr(input_set, "poscar", None)
    return getattr(poscar, "structure", fallback_structure)


def _write_vasp_input_objects(input_set, workdir: Path) -> None:
    """Write pymatgen-generated VASP input objects except POTCAR."""
    input_set.poscar.write_file(str(workdir / "POSCAR"))
    input_set.incar.write_file(str(workdir / "INCAR"))
    if getattr(input_set, "kpoints", None) is not None:
        input_set.kpoints.write_file(str(workdir / "KPOINTS"))


def _patch_incar(path: Path, settings: dict[str, Any]) -> None:
    from pymatgen.io.vasp.inputs import Incar

    incar = Incar.from_file(str(path))
    for key, value in settings.items():
        incar[key] = value
    incar.write_file(str(path))


def prepare_vasp_inputs(
    struct_name: str,
    structure,
    incar_content: str,
    *,
    mode: str,
    pressure: float = 0.0,
    potcar_dir: str | None = None,
    potcar_map: dict[str, str] | None = None,
    kpoints_path: str | Path | None = None,
    workdir: str | Path | None = None,
) -> dict[str, Any]:
    """Prepare VASP input files for a structure and return preparation metadata."""
    from pymatgen.io.vasp.inputs import Kpoints

    incar = parse_seed_incar(incar_content)
    control, user_incar = split_incar_settings(incar)

    user_incar["PSTRESS"] = pressure * 10.0
    if mode == "sp":
        user_incar["NSW"] = 0
        user_incar["IBRION"] = -1

    default_set = "MPStaticSet" if mode == "sp" else "MPRelaxSet"
    input_set_name = _coerce_input_set_name(control.get(CONTROL_INPUT_SET), default_set)
    input_set_cls = resolve_input_set_class(input_set_name)
    input_set_name = f"{input_set_cls.__module__}:{input_set_cls.__name__}" if (
        ":" in input_set_name or "." in input_set_name
    ) else input_set_cls.__name__

    kwargs: dict[str, Any] = {"user_incar_settings": user_incar}
    if potcar_map:
        kwargs["user_potcar_settings"] = potcar_map

    explicit_kpoints = None
    if kpoints_path and Path(kpoints_path).is_file():
        explicit_kpoints = Kpoints.from_file(str(kpoints_path))

    input_set = input_set_cls(structure, **kwargs)

    workdir_path = Path(workdir) if workdir is not None else Path(f"{struct_name}.vasp")
    workdir_path.mkdir(parents=True, exist_ok=True)

    potcar_metadata: list[dict[str, str]] = []
    root = _potcar_root(potcar_dir)
    if root is not None:
        _write_vasp_input_objects(input_set, workdir_path)
        potcar_metadata = assemble_potcar(
            _input_set_potcar_structure(input_set, structure),
            workdir_path / "POTCAR",
            potcar_dir=str(root),
            potcar_map=potcar_map,
            resolved_symbols=_input_set_potcar_symbol_map(input_set, structure),
        )
    else:
        try:
            input_set.write_input(str(workdir_path), make_dir_if_not_present=True)
        except TypeError:
            input_set.write_input(str(workdir_path))

    if explicit_kpoints is not None:
        shutil.copy2(str(kpoints_path), workdir_path / "KPOINTS")
    _patch_incar(workdir_path / "INCAR", user_incar)

    return {
        "workdir": workdir_path,
        "input_set": input_set_name,
        "incar_overrides": sorted(user_incar),
        "potcars": potcar_metadata,
        "kpoints_source": str(kpoints_path) if explicit_kpoints is not None else "input-set",
    }


def _read_structure_from_vasp(workdir: Path):
    from pymatgen.io.vasp.inputs import Poscar

    for name in ("CONTCAR", "POSCAR"):
        path = workdir / name
        if path.is_file():
            return Poscar.from_file(str(path)).structure
    raise FileNotFoundError(f"No CONTCAR or POSCAR found in {workdir}")


def _parse_vasprun(workdir: Path) -> dict[str, Any]:
    data: dict[str, Any] = {}
    path = workdir / "vasprun.xml"
    if not path.is_file():
        return data
    try:
        from pymatgen.io.vasp.outputs import Vasprun

        vasprun = Vasprun(str(path), parse_potcar_file=False)
        data["energy"] = getattr(vasprun, "final_energy", None)
        data["converged"] = bool(getattr(vasprun, "converged", False))
        ionic_steps = getattr(vasprun, "ionic_steps", None)
        if ionic_steps is not None:
            data["ionic_steps"] = len(ionic_steps)
        if getattr(vasprun, "final_structure", None) is not None:
            data["structure"] = vasprun.final_structure
    except Exception as exc:
        data["parse_error"] = str(exc)
    return data


def _parse_text_outputs(workdir: Path) -> dict[str, Any]:
    data: dict[str, Any] = {}
    outcar = workdir / "OUTCAR"
    if outcar.is_file():
        text = outcar.read_text(errors="ignore")
        matches = re.findall(r"enthalpy\s+is\s+(-?\d+(?:\.\d+)?)", text)
        if matches:
            data["enthalpy"] = float(matches[-1])
        pressure = re.findall(r"pressure\s+=?\s*(-?\d+(?:\.\d+)?)", text, re.I)
        if pressure:
            data["pressure"] = float(pressure[-1]) * 0.1
    oszicar = workdir / "OSZICAR"
    if oszicar.is_file():
        text = oszicar.read_text(errors="ignore")
        matches = re.findall(r"\bE0=\s*(-?\d+(?:\.\d+)?)", text)
        if matches:
            data["energy"] = float(matches[-1])
    return data


def _parse_external_pressure_gpa(workdir: Path) -> float | None:
    """Parse VASP PSTRESS from INCAR as external pressure in GPa."""
    incar = workdir / "INCAR"
    if not incar.is_file():
        return None
    text = incar.read_text(errors="ignore")
    for line in text.splitlines():
        stripped = re.split(r"[#!]", line, maxsplit=1)[0].strip()
        if not stripped:
            continue
        match = re.match(r"^PSTRESS\s*=?\s*([-+0-9.eEdD]+)", stripped, re.I)
        if match:
            return float(match.group(1).replace("D", "E").replace("d", "e")) / 10.0
    return None


def _enthalpy_from_energy_pressure_volume(
    energy: float,
    pressure: float,
    volume: float,
) -> float:
    """Return E + P*V with pressure in GPa and volume in Angstrom^3."""
    return energy + pressure * volume / EV_PER_ANG3_TO_GPA


def build_vasp_rem_lines(struct_name: str, metadata: dict[str, Any] | None = None) -> list[str]:
    """Build REM metadata lines for VASP results."""
    metadata = metadata or {}
    lines = [f"VASP input set {metadata.get('input_set', 'unknown')}"]
    if metadata.get("kpoints_source"):
        lines.append(f"KPOINTS {metadata['kpoints_source']}")
    for item in metadata.get("potcars", []):
        lines.append(
            "POTCAR {element} {symbol} sha256={sha256}".format(
                element=item["element"],
                symbol=item["symbol"],
                sha256=item["sha256"],
            )
        )
    outcar = Path(f"{struct_name}.vasp") / "OUTCAR"
    if outcar.is_file():
        for line in outcar.read_text(errors="ignore").splitlines():
            if "vasp." in line.lower() or "executed on" in line.lower():
                lines.append(line.strip())
                break
    return lines


def compose_vasp_task_doc(
    struct_name: str,
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Collect VASP outputs and write an AIRSS RES file."""
    from pymatgen.io.ase import AseAtomsAdaptor

    from .restools import save_airss_res

    workdir = Path(f"{struct_name}.vasp")
    parsed = _parse_vasprun(workdir)
    text_parsed = _parse_text_outputs(workdir)

    structure = parsed.get("structure")
    if structure is None:
        structure = _read_structure_from_vasp(workdir)

    raw_energy = parsed.get("energy")
    text_energy = text_parsed.get("energy")
    outcar_enthalpy = text_parsed.get("enthalpy")
    energy = raw_energy if raw_energy is not None else text_energy
    if energy is None:
        energy = outcar_enthalpy
    if energy is None:
        raise ValueError(f"No VASP energy found for {struct_name}")
    pressure = _parse_external_pressure_gpa(workdir)
    if pressure is None:
        pressure = text_parsed.get("pressure", 0.0)
    enthalpy = (
        outcar_enthalpy
        if outcar_enthalpy is not None
        else _enthalpy_from_energy_pressure_volume(energy, pressure, structure.volume)
    )

    atoms = AseAtomsAdaptor.get_atoms(structure)
    try:
        import spglib

        sg = spglib.get_spacegroup(
            (
                atoms.get_cell().array,
                atoms.get_scaled_positions(),
                atoms.get_atomic_numbers(),
            ),
            symprec=0.1,
        )
        sym = sg.split()[0] if sg else "P1"
    except Exception:
        sym = "P1"

    rem_lines = build_vasp_rem_lines(struct_name, metadata)
    info = {
        "uid": struct_name,
        "H": enthalpy,
        "P": pressure,
        "V": structure.volume,
        "nat": len(structure),
        "sym": sym,
        "rem": rem_lines,
    }
    save_airss_res(atoms, info, fname=struct_name + ".res", force_write=True)
    return {
        "structure": structure,
        "volume": structure.volume,
        "reduced_formula": structure.composition.reduced_formula,
        "formula": structure.composition.formula.replace(" ", ""),
        "natoms": len(structure),
        "label": struct_name,
        "energy": energy,
        "energy_per_atom": energy / len(structure) if energy is not None else None,
        "spin": 0.0,
        "mod_spin": 0.0,
        "pressure": pressure,
        "parallel_efficiency": None,
        "total_time": None,
        "res_content": Path(struct_name + ".res").read_text()
        if Path(struct_name + ".res").is_file()
        else None,
        "rem_lines": rem_lines,
        "converged": parsed.get("converged"),
    }
