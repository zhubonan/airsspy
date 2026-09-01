"""
ML interatomic potential runners.

Two orchestration drivers:

- **torch-sim** (default): GPU-accelerated batch processing via ``torch_sim``.
  Handles all structures in a single batch for maximum throughput.
- **ASE** (explicit): Uses any ASE-compatible calculator with ASE optimizers.
  Processes structures one at a time; works without torchsim.

The ``torch_sim`` package is an optional dependency (``pip install airsspy[ml]``).
Select the ASE driver explicitly with an ``ase:`` calculator specification.

Model specification format (torchsim)::

    backend:model_id

For example ``mace:medium``, ``torch-sim:mace:medium``,
``mace:/path/to/model.pt``, or ``sevennet:sevennet-mf-ompa``. The optional
``torch-sim:`` prefix explicitly selects this framework.

ASE calculator specification format::

    module.path:ClassName@model

For example ``mace.calculators:MACECalculator@medium``.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import logging
import os
import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read as ase_read
from ase.io import write as ase_write

logger = logging.getLogger(__name__)

# 1 eV/Ang^3 = 160.21766208 GPa
EV_PER_ANG3_TO_GPA = 160.21766208
StructureInput = str | Atoms
_SYMMETRIX_INTERNAL_PREFIX = "symmetrix:Symmetrix@"
_SYMMETRIX_GENERIC_ASE_PREFIX = "ase:symmetrix:Symmetrix@"
_SYMMETRIX_ASE_PREFIXES = ("ase:symmetrix:", "ase:symmetrics:")


def _enthalpy_from_energy_pressure_volume(
    energy: float | None,
    pressure: float | None,
    volume: float | None,
) -> float | None:
    """Return E + P*V with pressure in GPa and volume in Angstrom^3."""
    if energy is None:
        return None
    if pressure is None or volume is None:
        return energy
    return energy + pressure * volume / EV_PER_ANG3_TO_GPA


def _resolve_calculator(calculator_spec: str, **kwargs):
    """Import and instantiate an ASE calculator from a spec string.

    Args:
        calculator_spec: Calculator specification in the format
            ``module.path:ClassName@model``. The ``@model`` part is optional
            and is passed as the first positional argument to the constructor.
        **kwargs: Additional keyword arguments passed to the constructor.

    Returns:
        An instance of the calculator class.

    Raises:
        ValueError: If the spec cannot be parsed.
    """
    calculator_spec = _normalize_symmetrix_calculator_spec(calculator_spec)
    if calculator_spec.startswith(_SYMMETRIX_INTERNAL_PREFIX):
        return _resolve_symmetrix_calculator(calculator_spec, **kwargs)

    # Split off the @model suffix
    model = None
    if "@" in calculator_spec:
        calculator_spec, model = calculator_spec.rsplit("@", 1)

    # Split module:class
    if ":" in calculator_spec:
        module_path, class_name = calculator_spec.rsplit(":", 1)
    else:
        # Try treating the last dot-separated segment as the class
        parts = calculator_spec.rsplit(".", 1)
        if len(parts) == 2:
            module_path, class_name = parts
        else:
            raise ValueError(
                f"Cannot parse calculator spec: {calculator_spec!r}. "
                "Use 'module.path:ClassName' or 'module.path.ClassName'."
            )

    module = importlib.import_module(module_path)
    cls = getattr(module, class_name)

    if model is not None:
        return cls(model, **kwargs)
    return cls(**kwargs)


def _normalize_symmetrix_calculator_spec(calculator_spec: str) -> str:
    """Normalize public Symmetrix ASE specs to the internal calculator form."""
    if calculator_spec.startswith(_SYMMETRIX_GENERIC_ASE_PREFIX):
        return _normalize_symmetrix_calculator_spec(calculator_spec[len("ase:") :])

    if calculator_spec.startswith(_SYMMETRIX_INTERNAL_PREFIX):
        if len(calculator_spec) == len(_SYMMETRIX_INTERNAL_PREFIX):
            raise ValueError("A Symmetrix MACE model name is required.")
        return calculator_spec

    if calculator_spec.startswith("symmetrix:"):
        parts = calculator_spec.split(":", 2)
        if len(parts) != 3 or parts[1] != "mace" or not parts[2]:
            raise ValueError(
                "Symmetrix ML specs must use "
                "'ase:symmetrix:<mace-model>' or the legacy "
                "'symmetrix:mace:<model>' form."
            )
        return f"{_SYMMETRIX_INTERNAL_PREFIX}{parts[2]}"

    for prefix in _SYMMETRIX_ASE_PREFIXES:
        if calculator_spec.startswith(prefix):
            model_id = calculator_spec[len(prefix) :]
            if not model_id:
                raise ValueError(
                    "Symmetrix ASE specs must use "
                    "'ase:symmetrix:<mace-model>'."
                )
            return f"{_SYMMETRIX_INTERNAL_PREFIX}{model_id}"
    return calculator_spec


def _resolve_symmetrix_calculator(calculator_spec: str, **kwargs):
    """Instantiate the explicit Symmetrix MACE ASE calculator backend."""
    spec = calculator_spec[len("symmetrix:") :]
    if "@" not in spec:
        raise ValueError(
            "Cannot parse Symmetrix calculator spec. "
            "Use 'symmetrix:Symmetrix@<mace-model>'."
        )
    class_name, model_id = spec.rsplit("@", 1)
    if class_name != "Symmetrix":
        raise ValueError(
            f"Unsupported Symmetrix calculator: {class_name!r}. "
            "Expected 'symmetrix:Symmetrix@<mace-model>'."
        )

    try:
        module = importlib.import_module("symmetrix")
    except ImportError as exc:
        raise ImportError(
            "Symmetrix is required for 'ase:symmetrix:<mace-model>' specs. "
            "Install the symmetrix Python package or use another ML backend."
        ) from exc

    calc_kwargs = {"dtype": "float64", "use_kokkos": True}
    calc_kwargs.update(kwargs)
    model_id = _normalize_symmetrix_mace_model_name(model_id, calc_kwargs)
    model_file = _symmetrix_model_file(model_id, calc_kwargs)
    return module.Symmetrix(model_file, **calc_kwargs)


def _normalize_symmetrix_mace_model_name(model_id: str, calc_kwargs: dict) -> str:
    """Decode branded MACE-MH model names into checkpoint and head arguments."""
    model_name, separator, head = model_id.partition(":")
    prefix = "mace-mh-"
    generation = model_name.lower().removeprefix(prefix)
    if not model_name.lower().startswith(prefix) or not generation.isdigit():
        return model_id

    if separator:
        if not head:
            raise ValueError("A MACE-MH model head is required after ':'.")
        calc_kwargs.setdefault("head", head)
    return f"mh-{generation}"


def _resolve_mace_model_file(model_id: str) -> Path | str:
    """Resolve a MACE model id or local path for torch-sim-compatible names."""
    model_path = Path(model_id).expanduser()
    if model_path.is_file():
        return model_path

    from mace.calculators.foundations_models import download_mace_mp_checkpoint

    return Path(download_mace_mp_checkpoint(model=model_id))


def _symmetrix_model_file(model_id: str, calc_kwargs: dict) -> Path | str:
    """Return a cached Symmetrix JSON file for MACE models when possible."""
    model_file = _resolve_mace_model_file(model_id)
    extract_kwargs = {
        key: calc_kwargs[key]
        for key in ("species", "head", "num_spline_points")
        if key in calc_kwargs
    }
    if not extract_kwargs:
        return model_file

    extract_mace_data = importlib.import_module(
        "symmetrix.extract_mace_data"
    ).extract_mace_data

    cache_dir = _symmetrix_cache_dir()
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / (
        _symmetrix_cache_key(model_id, model_file, extract_kwargs) + ".json"
    )
    if not cache_path.is_file():
        data = extract_mace_data(model_file, **extract_kwargs)
        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile(
                "w",
                suffix=".json",
                prefix=cache_path.stem + "-",
                dir=cache_dir,
                delete=False,
                encoding="utf-8",
            ) as handle:
                json.dump(data, handle)
                tmp_path = Path(handle.name)
            tmp_path.replace(cache_path)
        finally:
            if tmp_path is not None and tmp_path.exists():
                tmp_path.unlink()
    return cache_path


def _symmetrix_cache_dir() -> Path:
    """Return the cache directory for converted Symmetrix JSON models."""
    env_path = os.environ.get("AIRSSPY_SYMMETRIX_CACHE")
    if env_path:
        return Path(env_path).expanduser()
    return (
        Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
        / "airsspy"
        / "symmetrix"
    )


def _symmetrix_cache_key(model_id: str, model_file: Path | str, extract_kwargs: dict) -> str:
    """Build a stable cache key for a Symmetrix model conversion."""
    model_path = Path(model_file).expanduser()
    identity: dict[str, object] = {"model_id": model_id}
    if model_path.is_file():
        stat = model_path.stat()
        identity.update(
            {
                "path": str(model_path.resolve()),
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            }
        )
    else:
        identity["path"] = str(model_file)
    identity["extract_kwargs"] = extract_kwargs
    payload = json.dumps(identity, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _cell_to_atoms(cell_path: str):
    """Read a .cell file into an ASE Atoms object."""
    atoms = ase_read(cell_path)
    atoms.pbc = True
    return atoms


def _cell_content_to_atoms(cell_content: str):
    """Read CASTEP cell content into ASE Atoms without leaving a .cell file."""
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            "w", suffix=".cell", prefix="airsspy-ml-", delete=False
        ) as handle:
            handle.write(cell_content)
            tmp_path = handle.name
        return _cell_to_atoms(tmp_path)
    finally:
        if tmp_path is not None:
            Path(tmp_path).unlink(missing_ok=True)


def _structure_input_to_atoms(structure_input: StructureInput):
    """Return ASE Atoms from either cell text or an existing Atoms object."""
    if isinstance(structure_input, Atoms):
        atoms = structure_input.copy()
        atoms.pbc = True
        return atoms
    return _cell_content_to_atoms(structure_input)


def _calculator_kwargs_for_atoms(
    calculator_spec: str,
    calculator_kwargs: dict | None,
    atoms: Atoms,
) -> dict:
    """Return calculator kwargs augmented with structure-specific metadata."""
    kwargs = dict(calculator_kwargs or {})
    normalized_spec = _normalize_symmetrix_calculator_spec(calculator_spec)
    if (
        normalized_spec.startswith(_SYMMETRIX_INTERNAL_PREFIX)
        and "species" not in kwargs
    ):
        kwargs["species"] = sorted({int(z) for z in atoms.get_atomic_numbers()})
    return kwargs


def _get_pressure_gpa(atoms) -> float:
    """Extract scalar pressure in GPa from an ASE Atoms with stress.

    ASE stress is in eV/Ang^3 (Voigt notation).
    Hydrostatic pressure = -trace(stress) / 3, converted to GPa.
    """
    try:
        stress = atoms.get_stress()  # eV/Ang^3
        pressure_ev_ang3 = -(stress[0] + stress[1] + stress[2]) / 3.0
        return pressure_ev_ang3 * EV_PER_ANG3_TO_GPA
    except Exception:
        return 0.0


def _normalize_static_forces(forces):
    """Return forces for one static-batch structure."""
    if forces is None:
        return None
    if getattr(forces, "ndim", 0) == 3 and forces.shape[0] == 1:
        return forces[0]
    return forces


def _normalize_static_stress(stress, index: int):
    """Return stress for one static-batch structure."""
    if stress is None:
        return None
    if getattr(stress, "ndim", 0) == 3:
        if stress.shape[0] == 1:
            return stress[0]
        return stress[index]
    return stress


def _as_numpy_array(value):
    """Convert torch/array-like values from torch-sim to numpy."""
    if hasattr(value, "detach"):
        value = value.detach().cpu().numpy()
    return np.asarray(value)


def _static_prop(props: dict, *names: str):
    """Return the first available torch-sim static property."""
    for name in names:
        if name in props:
            return props[name]
    names_text = ", ".join(names)
    raise KeyError(f"torch-sim static output missing expected key(s): {names_text}")


class AirssMlSinglePointRunner:
    """Run a single-point energy/forces/stress calculation using an ASE calculator."""

    _cleanup_extensions = [
        ".cell",
        ".extxyz",
        ".traj",
        "-orig.cell",
        ".res",
        ".err",
    ]

    def __init__(
        self,
        calculator_spec: str,
        calculator_kwargs: dict | None = None,
        pressure: float = 0.0,
    ) -> None:
        self.calculator_spec = calculator_spec
        self.calculator_kwargs = calculator_kwargs or {}
        self.pressure = pressure

    def clean_failed(self, struct_name: str) -> None:
        from .runners import clean_files

        clean_files(struct_name, self._cleanup_extensions)

    def run(self, struct_name: str, structure_input: StructureInput) -> int:
        """Attach calculator, compute energy/forces/stress, save results.

        Args:
            struct_name: Structure name (without extension).
            structure_input: Content of a .cell file or an ASE Atoms object.

        Returns:
            0 on success, 1 on failure.
        """
        try:
            atoms = _structure_input_to_atoms(structure_input)
            calc = _resolve_calculator(
                self.calculator_spec,
                **_calculator_kwargs_for_atoms(
                    self.calculator_spec,
                    self.calculator_kwargs,
                    atoms,
                ),
            )
            atoms.calc = calc

            energy = atoms.get_potential_energy()
            forces = atoms.get_forces()

            # Stress may not be available for all calculators
            try:
                stress = atoms.get_stress()
            except Exception:
                stress = None

            # Store results as SinglePointCalculator for downstream use
            sp_kwargs = {"energy": energy, "forces": forces}
            if stress is not None:
                sp_kwargs["stress"] = stress
            atoms.calc = SinglePointCalculator(atoms, **sp_kwargs)
            atoms.info["extern_pressure"] = self.pressure

            # Write output with attached calculator results
            ase_write(struct_name + ".extxyz", atoms, format="extxyz")

            logger.info("ML single-point OK: %s, E=%.6f eV", struct_name, energy)
            return 0
        except Exception:
            logger.error("ML single-point failed for %s", struct_name, exc_info=True)
            return 1


class AirssMlRelaxRunner:
    """Run geometry optimisation using an ASE calculator and ASE optimizer."""

    _cleanup_extensions = [
        ".cell",
        ".extxyz",
        ".traj",
        "-orig.cell",
        ".res",
        ".err",
    ]

    def __init__(
        self,
        calculator_spec: str,
        calculator_kwargs: dict | None = None,
        optimizer: str = "FIRE",
        fmax: float = 0.05,
        max_steps: int = 500,
        pressure: float = 0.0,
    ) -> None:
        self.calculator_spec = calculator_spec
        self.calculator_kwargs = calculator_kwargs or {}
        self.optimizer = optimizer
        self.fmax = fmax
        self.max_steps = max_steps
        self.pressure = pressure

    def clean_failed(self, struct_name: str) -> None:
        from .runners import clean_files

        clean_files(struct_name, self._cleanup_extensions)

    def run(self, struct_name: str, structure_input: StructureInput) -> int:
        """Attach calculator, run ASE optimizer, save results.

        Args:
            struct_name: Structure name (without extension).
            structure_input: Content of a .cell file or an ASE Atoms object.

        Returns:
            0 if converged, 1 if not converged or failed.
        """
        from ase.optimize import BFGS, FIRE

        try:
            atoms = _structure_input_to_atoms(structure_input)
            calc = _resolve_calculator(
                self.calculator_spec,
                **_calculator_kwargs_for_atoms(
                    self.calculator_spec,
                    self.calculator_kwargs,
                    atoms,
                ),
            )
            atoms.calc = calc

            # Choose optimizer
            opt_map = {"FIRE": FIRE, "BFGS": BFGS}
            opt_cls = opt_map.get(self.optimizer.upper())
            if opt_cls is None:
                raise ValueError(f"Unknown optimizer: {self.optimizer}")

            if self.pressure > 0.0:
                from ase.filters import ExpCellFilter

                # scalar_pressure sign: positive = compress (in eV/Ang^3)
                # GPa -> eV/Ang^3: divide by 160.21766208
                sp = self.pressure / EV_PER_ANG3_TO_GPA
                ecf = ExpCellFilter(atoms, scalar_pressure=sp)
                dyn = opt_cls(ecf, trajectory=struct_name + ".traj")
            else:
                dyn = opt_cls(atoms, trajectory=struct_name + ".traj")

            converged = False
            try:
                converged = bool(dyn.run(fmax=self.fmax, steps=self.max_steps))
            except Exception:
                logger.warning("Optimizer failed for %s", struct_name)
                return 1

            relax_status = "converged" if converged else "max_steps"
            atoms.info["relax_converged"] = converged
            atoms.info["relax_status"] = relax_status
            atoms.info["relax_steps"] = getattr(dyn, "nsteps", None)
            atoms.info["relax_fmax"] = self.fmax
            atoms.info["relax_max_steps"] = self.max_steps
            atoms.info["extern_pressure"] = self.pressure

            # Store final results as SinglePointCalculator
            try:
                final_energy = atoms.get_potential_energy()
                final_forces = atoms.get_forces()
                sp_kwargs = {"energy": final_energy, "forces": final_forces}
                try:
                    sp_kwargs["stress"] = atoms.get_stress()
                except Exception:
                    pass
                atoms.calc = SinglePointCalculator(atoms, **sp_kwargs)
            except Exception:
                pass

            # Write output
            ase_write(struct_name + ".extxyz", atoms, format="extxyz")

            logger.info(
                "ML relax %s: %s",
                struct_name,
                "converged" if converged else "not converged",
            )
            return 0
        except Exception:
            logger.error("ML relaxation failed for %s", struct_name, exc_info=True)
            return 1


def compose_ml_task_doc(struct_name: str, calculator_spec: str = "") -> dict:
    """Extract results from a completed ML calculation.

    Reads the ``.extxyz`` output file (with SinglePointCalculator attached),
    creates an ASE Atoms object, writes a ``.res`` file, and returns a
    dictionary suitable for constructing an ``AirssResultDoc``.

    Args:
        struct_name: Structure name (without extension).
        calculator_spec: The calculator spec string (for REM metadata).

    Returns:
        Dictionary with energy, structure, volume, formula, etc.
    """
    from pymatgen.io.ase import AseAtomsAdaptor

    from ..restools import save_airss_res

    extxyz_path = struct_name + ".extxyz"
    atoms = ase_read(extxyz_path)

    energy = None
    pressure = None
    forces = None

    if atoms.calc is not None:
        try:
            energy = atoms.get_potential_energy()
        except Exception:
            pass
        try:
            forces = atoms.get_forces()
        except Exception:
            pass
        pressure = atoms.info.get("extern_pressure", atoms.info.get("pressure"))
        if pressure is None:
            pressure = _get_pressure_gpa(atoms)
        else:
            pressure = float(pressure)

    volume = atoms.get_volume()

    # Compute symmetry via spglib
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
    except (ImportError, Exception):
        sym = "P1"

    # REM lines for ML calculation
    rem_lines = ["", f"ML Calculator {calculator_spec}"]
    relax_status = atoms.info.get("relax_status")
    if relax_status is not None:
        rem_lines.append(f"ML Relax status {relax_status}")
    if "relax_converged" in atoms.info:
        rem_lines.append(f"ML Relax converged {bool(atoms.info['relax_converged'])}")
    if atoms.info.get("relax_steps") is not None:
        rem_lines.append(f"ML Relax steps {atoms.info['relax_steps']}")
    rem_lines.append("")

    enthalpy = _enthalpy_from_energy_pressure_volume(energy, pressure, volume)

    info = {
        "uid": struct_name,
        "P": pressure if pressure is not None else 0.0,
        "V": volume,
        "H": enthalpy if enthalpy is not None else 0.0,
        "nat": len(atoms),
        "sym": sym,
        "rem": rem_lines,
    }

    # Build atom annotations from forces if available
    atom_annotations = None
    if forces is not None:
        atom_annotations = []
        for f in forces:
            atom_annotations.append(f"{f[0]:.6f} {f[1]:.6f} {f[2]:.6f}")

    save_airss_res(
        atoms,
        info,
        fname=struct_name + ".res",
        force_write=True,
        atom_annotations=atom_annotations,
    )

    structure = AseAtomsAdaptor.get_structure(atoms)

    return {
        "structure": structure,
        "volume": structure.volume,
        "reduced_formula": structure.reduced_formula,
        "formula": structure.composition.formula.replace(" ", ""),
        "natoms": len(atoms),
        "label": struct_name,
        "energy": energy,
        "energy_per_atom": energy / len(atoms) if energy else None,
        "pressure": pressure,
        "total_time": None,
        "relax_converged": atoms.info.get("relax_converged"),
        "relax_status": atoms.info.get("relax_status"),
        "relax_steps": atoms.info.get("relax_steps"),
        "res_content": Path(struct_name + ".res").read_text()
        if Path(struct_name + ".res").is_file()
        else None,
        "rem_lines": rem_lines,
    }


# ---------------------------------------------------------------------------
# TorchSim batch runners (optional dependency)
# ---------------------------------------------------------------------------


def has_torchsim() -> bool:
    """Check whether ``torch_sim`` is importable."""
    try:
        import torch_sim  # noqa: F401

        return True
    except Exception as exc:
        logger.debug("torch_sim is unavailable: %s", exc)
        return False


class TorchSimRunner:
    """Reusable torch-sim model context for chunked ML runs."""

    def __init__(self, model_spec: str, *, device: str | None = None) -> None:
        import torch

        self.model_spec = _normalize_torchsim_model_spec(model_spec)
        self.device = (
            torch.device(device)
            if device
            else (
                torch.device("cuda")
                if torch.cuda.is_available()
                else torch.device("cpu")
            )
        )
        self.dtype = torch.float32 if self.device.type == "cuda" else torch.float64
        self.model = _load_torchsim_model(
            self.model_spec,
            device=self.device,
            dtype=self.dtype,
        )

    def relax_batch(
        self,
        struct_names: list[str],
        structures: list[StructureInput],
        *,
        max_steps: int = 300,
        force_tol: float = 0.05,
        optimizer: str = "fire",
        cell_filter: str = "frechet",
        convergence_mode: str = "force_stress",
        scalar_pressure: float = 0.0,
    ) -> dict[str, int]:
        """Relax one batch of structures using the loaded model."""
        import torch_sim as ts

        results: dict[str, int] = {}

        opt_map = {
            "fire": ts.Optimizer.fire,
            "lbfgs": ts.Optimizer.lbfgs,
            "bfgs": ts.Optimizer.bfgs,
            "gradient_descent": ts.Optimizer.gradient_descent,
        }
        opt = opt_map[optimizer.lower()]

        cf_map = {
            "frechet": ts.CellFilter.frechet,
            "unit": ts.CellFilter.unit,
        }
        cf = cf_map[cell_filter.lower()]

        conv_fn = _resolve_torchsim_convergence(
            convergence_mode,
            force_tol=force_tol,
        )

        atoms_list = [_structure_input_to_atoms(structure) for structure in structures]

        state = ts.io.atoms_to_state(
            atoms_list,
            device=self.device,
            dtype=self.dtype,
        )

        converged_state = ts.optimize(
            system=state,
            model=self.model,
            optimizer=opt,
            convergence_fn=conv_fn,
            max_steps=max_steps,
            init_kwargs={
                "cell_filter": cf,
                "scalar_pressure": scalar_pressure / EV_PER_ANG3_TO_GPA,
            },
        )

        final_atoms_list = ts.io.state_to_atoms(converged_state)
        outputs = self.model.forward(converged_state)

        energy_values = outputs["energy"].detach().cpu().reshape(-1).tolist()
        force_values = outputs["forces"].detach().cpu().numpy()
        stress_tensor = outputs.get("stress")
        stress_values = (
            stress_tensor.detach().cpu().numpy() if stress_tensor is not None else None
        )

        force_offset = 0
        for index, (name, atoms) in enumerate(zip(struct_names, final_atoms_list)):
            natoms = len(atoms)
            calc_kwargs: dict = {
                "energy": float(energy_values[index]),
                "forces": force_values[force_offset : force_offset + natoms],
            }
            force_offset += natoms
            if stress_values is not None and index < len(stress_values):
                calc_kwargs["stress"] = stress_values[index]
            atoms.calc = SinglePointCalculator(atoms, **calc_kwargs)
            atoms.info["extern_pressure"] = scalar_pressure

            ase_write(name + ".extxyz", atoms, format="extxyz")
            results[name] = 0

        return results

    def static_batch(
        self,
        struct_names: list[str],
        structures: list[StructureInput],
        *,
        scalar_pressure: float = 0.0,
    ) -> dict[str, int]:
        """Run one static batch using the loaded model."""
        import torch_sim as ts

        results: dict[str, int] = {}
        atoms_list = [_structure_input_to_atoms(structure) for structure in structures]
        state = ts.io.atoms_to_state(
            atoms_list,
            device=self.device,
            dtype=self.dtype,
        )
        props_list = ts.static(system=state, model=self.model)

        final_atoms_list = ts.io.state_to_atoms(state)
        for index, (name, atoms) in enumerate(zip(struct_names, final_atoms_list)):
            props = props_list[index]
            energy_prop = _static_prop(props, "energy", "potential_energy")
            energy = float(_as_numpy_array(energy_prop).reshape(-1)[0])
            forces = _as_numpy_array(props["forces"]) if "forces" in props else None
            stress = _as_numpy_array(props["stress"]) if "stress" in props else None

            calc_kwargs: dict = {"energy": energy}
            if forces is not None:
                calc_kwargs["forces"] = _normalize_static_forces(forces)
            if stress is not None:
                calc_kwargs["stress"] = _normalize_static_stress(stress, index)

            atoms.calc = SinglePointCalculator(atoms, **calc_kwargs)
            atoms.info["extern_pressure"] = scalar_pressure
            ase_write(name + ".extxyz", atoms, format="extxyz")
            results[name] = 0

        return results


def _torchsim_relax_batch(
    model_spec: str,
    struct_names: list[str],
    structures: list[StructureInput],
    *,
    device: str | None = None,
    max_steps: int = 300,
    force_tol: float = 0.05,
    optimizer: str = "fire",
    cell_filter: str = "frechet",
    convergence_mode: str = "force_stress",
    scalar_pressure: float = 0.0,
) -> dict[str, int]:
    """Relax a batch of structures using torchsim."""
    return TorchSimRunner(model_spec, device=device).relax_batch(
        struct_names,
        structures,
        max_steps=max_steps,
        force_tol=force_tol,
        optimizer=optimizer,
        cell_filter=cell_filter,
        convergence_mode=convergence_mode,
        scalar_pressure=scalar_pressure,
    )


def _torchsim_static_batch(
    model_spec: str,
    struct_names: list[str],
    structures: list[StructureInput],
    *,
    device: str | None = None,
    scalar_pressure: float = 0.0,
) -> dict[str, int]:
    """Run single-point calculations on a batch of structures using torchsim."""
    return TorchSimRunner(model_spec, device=device).static_batch(
        struct_names,
        structures,
        scalar_pressure=scalar_pressure,
    )


def _normalize_torchsim_model_spec(model_spec: str) -> str:
    """Remove an explicit framework prefix and validate a torch-sim model spec."""
    if model_spec.startswith("torch-sim:"):
        model_spec = model_spec[len("torch-sim:") :]

    if ":" not in model_spec:
        raise ValueError(
            "Torch-sim model specs must use 'backend:<model>', for example "
            "'torch-sim:mace:medium'."
        )
    backend, model_id = model_spec.split(":", 1)
    if not backend or not model_id:
        raise ValueError(
            "Torch-sim model specs must use 'backend:<model>', for example "
            "'torch-sim:mace:medium'."
        )
    return model_spec


def _load_torchsim_model(model_spec: str, *, device, dtype):
    """Load a torchsim model from a ``backend:model_id`` spec string.

    Supported backends: mace, fairchem, mattersim, sevennet, orb,
    graphpes, metatomic, nequix.
    """
    backend, model_id = model_spec.split(":", 1)

    if backend == "mace":
        from pathlib import Path as _P

        from torch_sim.models.mace import MaceModel

        model_path = _P(model_id)
        if model_path.is_file():
            return MaceModel(
                model=model_path,
                device=device,
                dtype=dtype,
                compute_forces=True,
                compute_stress=True,
            )
        from mace.calculators.foundations_models import mace_mp

        raw_model = mace_mp(
            model=model_id,
            return_raw_model=True,
            default_dtype=str(dtype).removeprefix("torch."),
            device=str(device),
        )
        return MaceModel(
            model=raw_model,
            device=device,
            dtype=dtype,
            compute_forces=True,
            compute_stress=True,
        )

    if backend == "fairchem":
        from torch_sim.models.fairchem import FairChemModel

        return FairChemModel(
            model=model_id,
            device=device,
            dtype=dtype,
            compute_stress=True,
        )

    if backend == "mattersim":
        from mattersim.forcefield import Potential
        from torch_sim.models.mattersim import MatterSimModel

        raw_model = Potential.from_checkpoint(
            load_path=model_id,
            model_name="m3gnet",
            device=str(device),
            load_training_state=False,
        )
        return MatterSimModel(model=raw_model, device=device, dtype=dtype)

    if backend == "sevennet":
        from torch_sim.models.sevennet import SevenNetModel

        return SevenNetModel(model=model_id, device=device, dtype=dtype)

    if backend == "orb":
        from torch_sim.models.orb import OrbModel

        return OrbModel(model=model_id, device=device)

    if backend == "graphpes":
        from torch_sim.models.graphpes_framework import GraphPESWrapper

        return GraphPESWrapper(
            model=model_id,
            device=device,
            dtype=dtype,
            compute_forces=True,
            compute_stress=True,
        )

    if backend == "metatomic":
        from torch_sim.models.metatomic import MetatomicModel

        return MetatomicModel(
            model=model_id,
            device=device,
            compute_forces=True,
            compute_stress=True,
        )

    if backend == "nequix":
        from torch_sim.models.nequix import NequixModel

        return NequixModel(model=model_id, device=device, dtype=dtype)

    raise ValueError(
        f"Unsupported torchsim backend: {backend!r}. "
        "Expected: mace, fairchem, mattersim, sevennet, orb, "
        "graphpes, metatomic, nequix."
    )


def _resolve_torchsim_convergence(mode: str, *, force_tol: float = 0.05):
    """Build a convergence function for torchsim.optimize()."""
    import torch_sim as ts

    if mode == "force":
        return ts.generate_force_convergence_fn(force_tol=force_tol)
    if mode == "force_stress":
        stress_tol_ev = 1.0 / EV_PER_ANG3_TO_GPA  # 1 GPa in eV/Ang^3

        def convergence_fn(state, _last_energy=None):
            force_conv = ts.system_wise_max_force(state) < force_tol
            if state.stress is None:
                raise ValueError("Stress required for force_stress convergence.")
            stress_conv = state.stress.abs().amax(dim=(1, 2)) < stress_tol_ev
            return force_conv & stress_conv

        return convergence_fn
    if mode == "energy":
        return ts.generate_energy_convergence_fn()

    raise ValueError(f"Unsupported convergence mode: {mode!r}")
