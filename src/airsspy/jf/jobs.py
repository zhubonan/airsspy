"""
Jobflow Makers for AIRSS searches and relaxations.

Both ``AirssSearchMaker`` and ``AirssRelaxMaker`` produce the same
``AirssJobDoc`` output type, supporting multi-structure jobs for
high-throughput scenarios.
"""

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from jobflow import Maker, Response, job
from pymatgen.core import Structure

from .documents import AirssJobDoc, AirssResultDoc, RelaxOutcome
from .runners import (
    AirssAbacusRelaxRunner,
    AirssCastepRelaxRunner,
    AirssGulpRelaxRunner,
    AirssPp3RelaxRunner,
    compose_task_doc,
    run_buildcell,
)

logger = logging.getLogger(__name__)


def _get_hash(content: str) -> str:
    """Return a short hash of a string."""
    import hashlib

    return hashlib.md5(content.encode()).hexdigest()[:8]


@dataclass
class AirssSearchMaker(Maker):
    """
    Run N build+relax cycles as a single jobflow job.

    Generates N random structures from a seed using buildcell,
    then relaxes each one with CASTEP. All results are collected
    into a single ``AirssJobDoc``.
    """

    name: str = "airss search"
    n_structures: int = 50
    executable: str = "castep.mpi"
    build_timeout: int = 60
    cycles: int = 4
    max_fails: int = 2
    max_iterations: int = 200
    write_seed: bool = True
    stop_if_not_converged: bool = False
    code: str = "castep"

    @job
    def make(
        self,
        seed_name: str,
        seed_content: str,
        paraminput,
        project_name: str,
    ) -> Response:
        """
        Generate N random structures and relax them.

        Args:
            seed_name: Name of the seed (without extension).
            seed_content: Content of the seed .cell file.
            paraminput: ParamInput instance for CASTEP.
            project_name: Project identifier for grouping results.
        """
        results: list[AirssResultDoc] = []
        n_failed = 0

        for i in range(self.n_structures):
            logger.info("Building structure %d/%d", i + 1, self.n_structures)

            struct_name = None
            runner = None

            try:
                build_output = run_buildcell(
                    seed_name,
                    seed_content,
                    build_timeout=self.build_timeout,
                    write_seed=self.write_seed,
                )
                if build_output is None:
                    logger.warning("Buildcell timed out for structure %d", i + 1)
                    continue

                struct_name = build_output["struct_name"]

                if self.code == "castep":
                    from castepinput.inputs import CellInput

                    cellinput = CellInput.from_file(struct_name + ".cell")

                    runner = AirssCastepRelaxRunner(
                        executable=self.executable,
                        max_fails=self.max_fails,
                        max_iterations=self.max_iterations,
                    )
                    return_code = runner.run(struct_name, cellinput, paraminput)
                elif self.code == "gulp":
                    struct_content = Path(struct_name + ".cell").read_text()
                    param_content = (
                        paraminput.get_string()
                        if hasattr(paraminput, "get_string")
                        else str(paraminput)
                    )
                    runner = AirssGulpRelaxRunner(executable=self.executable)
                    return_code = runner.run(
                        struct_name, struct_content, param_content, seed_name=seed_name
                    )
                elif self.code == "pp3":
                    struct_content = Path(struct_name + ".cell").read_text()
                    param_content = (
                        paraminput.get_string()
                        if hasattr(paraminput, "get_string")
                        else str(paraminput)
                    )
                    runner = AirssPp3RelaxRunner(executable=self.executable)
                    return_code = runner.run(
                        struct_name, struct_content, param_content, seed_name=seed_name
                    )
                elif self.code == "abacus":
                    from ..abacustools import compose_abacus_task_doc

                    struct_content = Path(struct_name + ".cell").read_text()
                    param_content = (
                        paraminput.get_string()
                        if hasattr(paraminput, "get_string")
                        else str(paraminput)
                    )
                    runner = AirssAbacusRelaxRunner(
                        executable=self.executable,
                        max_iterations=self.max_iterations,
                    )
                    return_code = runner.run(struct_name, struct_content, param_content)
                else:
                    raise ValueError(f"Unknown code: {self.code}")

                if return_code == 0:
                    relax_status = RelaxOutcome.FINISHED
                else:
                    relax_status = RelaxOutcome.ERRORED

                if self.code == "abacus":
                    from ..abacustools import compose_abacus_task_doc

                    task_doc = compose_abacus_task_doc(struct_name)
                else:
                    task_doc = compose_task_doc(struct_name)
                result_doc = AirssResultDoc(
                    struct_name=struct_name,
                    seed_name=seed_name,
                    project_name=project_name,
                    structure=task_doc.get("structure"),
                    energy=task_doc.get("energy"),
                    energy_per_atom=task_doc.get("energy_per_atom"),
                    volume=task_doc.get("volume"),
                    pressure=task_doc.get("pressure"),
                    spin=task_doc.get("spin", 0.0),
                    mod_spin=task_doc.get("mod_spin", 0.0),
                    symmetry=task_doc.get("symmetry"),
                    formula=task_doc.get("formula"),
                    reduced_formula=task_doc.get("reduced_formula"),
                    natoms=task_doc.get("natoms"),
                    res_content=task_doc.get("res_content"),
                    parallel_efficiency=task_doc.get("parallel_efficiency"),
                    total_time=task_doc.get("total_time"),
                    relax_status=relax_status,
                    rem_lines=task_doc.get("rem_lines"),
                )
                results.append(result_doc)

            except Exception as e:
                logger.error(
                    "Structure %d/%d failed with exception: %s",
                    i + 1,
                    self.n_structures,
                    e,
                    exc_info=True,
                )
                if struct_name and runner:
                    try:
                        runner.clean_failed(struct_name)
                    except Exception:
                        pass
                n_failed += 1
                results.append(
                    AirssResultDoc(
                        struct_name=struct_name or f"unknown-{i}",
                        seed_name=seed_name,
                        project_name=project_name,
                        relax_status=RelaxOutcome.FAILED,
                        error_message=str(e),
                    )
                )
                continue

        n_finished = sum(1 for r in results if r.relax_status == RelaxOutcome.FINISHED)
        n_errored = sum(1 for r in results if r.relax_status == RelaxOutcome.ERRORED)

        search_doc = AirssJobDoc(
            project_name=project_name,
            seed_name=seed_name,
            job_type="search",
            seed_content=seed_content,
            seed_hash=_get_hash(seed_content),
            results=results,
            n_structures=len(results),
            n_finished=n_finished,
            n_errored=n_errored,
            n_failed=n_failed,
        )

        stop = False
        if (
            self.stop_if_not_converged
            and len(results) > 0
            and n_errored == len(results)
        ):
            stop = True
        return Response(stop_children=stop, output=search_doc)


@dataclass
class AirssRelaxMaker(Maker):
    """
    Relax one or more provided structures as a single jobflow job.

    Accepts lists of structures, names, and cell inputs. All structures
    share the same param input and project/seed metadata.
    """

    name: str = "airss relax"
    executable: str = "castep.mpi"
    cycles: int = 4
    max_fails: int = 2
    max_iterations: int = 200
    stop_if_not_converged: bool = False
    code: str = "castep"

    @job
    def make(
        self,
        structures: list[Structure],
        struct_names: list[str],
        cellinputs: list,
        paraminput,
        project_name: str,
        seed_name: str,
    ) -> Response:
        """
        Relax N structures in a single job.

        Args:
            structures: List of pymatgen Structure objects.
            struct_names: Corresponding structure names.
            cellinputs: Corresponding CellInput instances.
            paraminput: Shared ParamInput instance.
            project_name: Project identifier.
            seed_name: Seed name for metadata.
        """
        results: list[AirssResultDoc] = []
        n_failed = 0

        for structure, struct_name, cellinput in zip(
            structures, struct_names, cellinputs
        ):
            runner = None
            try:
                if self.code == "castep":
                    runner = AirssCastepRelaxRunner(
                        executable=self.executable,
                        max_fails=self.max_fails,
                        max_iterations=self.max_iterations,
                    )

                    cellinput.set_positions(
                        [str(elem) for elem in structure.species],
                        structure.cart_coords,
                    )
                    cellinput.set_cell(structure.lattice.matrix)
                    return_code = runner.run(struct_name, cellinput, paraminput)
                elif self.code == "gulp":
                    struct_content = (
                        cellinput.get_string()
                        if hasattr(cellinput, "get_string")
                        else str(cellinput)
                    )
                    param_content = (
                        paraminput.get_string()
                        if hasattr(paraminput, "get_string")
                        else str(paraminput)
                    )
                    runner = AirssGulpRelaxRunner(executable=self.executable)
                    return_code = runner.run(
                        struct_name, struct_content, param_content, seed_name=seed_name
                    )
                elif self.code == "pp3":
                    struct_content = (
                        cellinput.get_string()
                        if hasattr(cellinput, "get_string")
                        else str(cellinput)
                    )
                    param_content = (
                        paraminput.get_string()
                        if hasattr(paraminput, "get_string")
                        else str(paraminput)
                    )
                    runner = AirssPp3RelaxRunner(executable=self.executable)
                    return_code = runner.run(
                        struct_name, struct_content, param_content, seed_name=seed_name
                    )
                elif self.code == "abacus":
                    from ..abacustools import compose_abacus_task_doc

                    struct_content = (
                        cellinput.get_string()
                        if hasattr(cellinput, "get_string")
                        else str(cellinput)
                    )
                    param_content = (
                        paraminput.get_string()
                        if hasattr(paraminput, "get_string")
                        else str(paraminput)
                    )
                    runner = AirssAbacusRelaxRunner(
                        executable=self.executable,
                        max_iterations=self.max_iterations,
                    )
                    return_code = runner.run(struct_name, struct_content, param_content)
                else:
                    raise ValueError(f"Unknown code: {self.code}")

                relax_status = (
                    RelaxOutcome.FINISHED if return_code == 0 else RelaxOutcome.ERRORED
                )

                if self.code == "abacus":
                    from ..abacustools import compose_abacus_task_doc

                    task_doc = compose_abacus_task_doc(struct_name)
                else:
                    task_doc = compose_task_doc(struct_name)
                result_doc = AirssResultDoc(
                    struct_name=struct_name,
                    seed_name=seed_name,
                    project_name=project_name,
                    structure=task_doc.get("structure"),
                    energy=task_doc.get("energy"),
                    energy_per_atom=task_doc.get("energy_per_atom"),
                    volume=task_doc.get("volume"),
                    pressure=task_doc.get("pressure"),
                    spin=task_doc.get("spin", 0.0),
                    mod_spin=task_doc.get("mod_spin", 0.0),
                    symmetry=task_doc.get("symmetry"),
                    formula=task_doc.get("formula"),
                    reduced_formula=task_doc.get("reduced_formula"),
                    natoms=task_doc.get("natoms"),
                    res_content=task_doc.get("res_content"),
                    parallel_efficiency=task_doc.get("parallel_efficiency"),
                    total_time=task_doc.get("total_time"),
                    relax_status=relax_status,
                    rem_lines=task_doc.get("rem_lines"),
                )
                results.append(result_doc)

            except Exception as e:
                logger.error(
                    "Structure %s failed with exception: %s",
                    struct_name,
                    e,
                    exc_info=True,
                )
                if runner:
                    try:
                        runner.clean_failed(struct_name)
                    except Exception:
                        pass
                n_failed += 1
                results.append(
                    AirssResultDoc(
                        struct_name=struct_name,
                        seed_name=seed_name,
                        project_name=project_name,
                        relax_status=RelaxOutcome.FAILED,
                        error_message=str(e),
                    )
                )
                continue

        n_finished = sum(1 for r in results if r.relax_status == RelaxOutcome.FINISHED)
        n_errored = sum(1 for r in results if r.relax_status == RelaxOutcome.ERRORED)

        relax_doc = AirssJobDoc(
            project_name=project_name,
            seed_name=seed_name,
            job_type="relax",
            results=results,
            n_structures=len(results),
            n_finished=n_finished,
            n_errored=n_errored,
            n_failed=n_failed,
        )

        stop = False
        if (
            self.stop_if_not_converged
            and len(results) > 0
            and n_errored == len(results)
        ):
            stop = True
        return Response(stop_children=stop, output=relax_doc)


@dataclass
class AirssValidateMaker(Maker):
    """Validate that required AIRSS executables are installed."""

    name: str = "airss validate"
    additional_exes: tuple = ()
    required_exes: tuple = ("buildcell", "castep_relax", "castep2res")

    @job
    def make(self) -> Optional[Response]:
        """Check that all required executables are on PATH."""
        exes_to_check = self.additional_exes + self.required_exes
        not_found = []
        for exe_name in exes_to_check:
            try:
                subprocess.run(["which", exe_name], check=True)
            except subprocess.CalledProcessError:
                not_found.append(exe_name)

        if not_found:
            logger.error(
                "AIRSS installation incomplete, executables not found: %s",
                not_found,
            )
            return Response(stop_jobflow=True)
        return None
