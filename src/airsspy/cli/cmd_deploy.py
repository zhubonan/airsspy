"""
CLI commands for deploying AIRSS searches and relaxations.
"""

from pathlib import Path

import click

SUFFIX_MAP = {
    "castep": ".param",
    "gulp": ".lib",
    "pp3": ".pp",
    "abacus": ".INPUT",
}


@click.group("deploy")
def deploy():
    """Deploy search/relaxation jobs."""


@deploy.command("search")
@click.option(
    "--seed",
    required=True,
    help="Name of the seed (<seed>.cell and <seed>.param must exist)",
)
@click.option("--project", required=True, help="Project name")
@click.option("--num", required=True, type=int, help="Number of structures to search")
@click.option("--exe", default="castep.mpi", show_default=True, help="Executable name")
@click.option(
    "--cycles",
    default=4,
    type=int,
    show_default=True,
    help="Max relaxation cycles per structure",
)
@click.option(
    "--max-iterations",
    default=200,
    type=int,
    show_default=True,
    help="Max total geometry iterations",
)
@click.option(
    "--n-structures", default=50, type=int, show_default=True, help="Structures per job"
)
@click.option(
    "--build-timeout",
    default=60,
    type=int,
    show_default=True,
    help="Buildcell timeout in seconds",
)
@click.option(
    "--code", default="castep", show_default=True, help="Code to use: castep, gulp, pp3"
)
@click.option("--dryrun", is_flag=True, help="Print configuration without submitting")
@click.pass_context
def deploy_search(
    ctx,
    seed,
    project,
    num,
    exe,
    cycles,
    max_iterations,
    n_structures,
    build_timeout,
    code,
    dryrun,
):
    """Deploy an AIRSS search using jobflow."""
    from castepinput.inputs import ParamInput

    seed_content = Path(seed + ".cell").read_text()
    param_content = Path(seed + SUFFIX_MAP[code]).read_text()

    # Override executable for GULP/pp3/abacus if still set to a CASTEP default
    if code == "gulp" and "castep" in exe:
        exe = "ggulp"
    elif code == "pp3" and "castep" in exe:
        exe = "pp3"
    elif code == "abacus" and "castep" in exe:
        exe = "abacus"

    if dryrun:
        click.echo(f"Project: {project}, Seed: {seed}")
        click.echo(f"Structures requested: {num}")
        click.echo(f"Structures per job: {n_structures}")
        click.echo(f"Executable: {exe}, Code: {code}")
        click.echo(f"Cycles: {cycles}, Max iterations: {max_iterations}")
        click.echo(f"Build timeout: {build_timeout}s")
        click.echo(f"\nSeed content:\n{seed_content}")
        click.echo(f"\nParam content:\n{param_content}")
        return

    from jobflow import Flow, JobStore, run_locally
    from maggma.stores import MongoStore

    from airsspy.jf.jobs import AirssSearchMaker

    n_jobs = (num + n_structures - 1) // n_structures
    click.echo(f"Creating {n_jobs} job(s) for {num} structures")

    maker = AirssSearchMaker(
        n_structures=n_structures,
        executable=exe,
        cycles=cycles,
        max_iterations=max_iterations,
        build_timeout=build_timeout,
        code=code,
    )

    jobs = []
    for i in range(n_jobs):
        remaining = num - i * n_structures
        if remaining < n_structures:
            maker_i = AirssSearchMaker(
                n_structures=remaining,
                executable=exe,
                cycles=cycles,
                max_iterations=max_iterations,
                build_timeout=build_timeout,
                code=code,
            )
        else:
            maker_i = maker

        job = maker_i.make(
            seed_name=seed,
            seed_content=seed_content,
            paraminput=ParamInput.from_file(seed + SUFFIX_MAP[code]),
            project_name=project,
        )
        jobs.append(job)

    flow = Flow(jobs)

    docs_store = MongoStore(
        database=ctx.obj["db_name"],
        collection_name="jobs",
        host=ctx.obj["db_host"],
        port=ctx.obj["db_port"],
    )
    store = JobStore(docs_store=docs_store)

    click.echo("Running flow...")
    run_locally(flow, store=store)
    click.echo("Done")


@deploy.command("relax")
@click.option(
    "--seed", default="USER-RELAX", show_default=True, help="Seed name for metadata"
)
@click.option("--project", required=True, help="Project name")
@click.option(
    "--cell",
    required=True,
    help="Glob pattern for cell/structure files",
)
@click.option("--param", required=True, type=click.Path(exists=True), help="Param file")
@click.option("--base-cell", help="Base cell file (structures are overlaid onto this)")
@click.option("--exe", default="castep.mpi", show_default=True, help="Executable name")
@click.option("--cycles", default=4, type=int, show_default=True)
@click.option("--max-iterations", default=200, type=int, show_default=True)
@click.option("--code", default="castep", show_default=True)
@click.option("--dryrun", is_flag=True)
@click.pass_context
def deploy_relax(
    ctx,
    seed,
    project,
    cell,
    param,
    base_cell,
    exe,
    cycles,
    max_iterations,
    code,
    dryrun,
):
    """Deploy relaxation of existing structures using jobflow."""
    from ase.io import read as ase_read
    from castepinput.inputs import CellInput, ParamInput
    from pymatgen.io.ase import AseAtomsAdaptor

    from airsspy.jf.jobs import AirssRelaxMaker

    # Override executable for GULP/pp3/abacus if still set to a CASTEP default
    if code == "gulp" and "castep" in exe:
        exe = "ggulp"
    elif code == "pp3" and "castep" in exe:
        exe = "pp3"
    elif code == "abacus" and "castep" in exe:
        exe = "abacus"

    paraminput = ParamInput.from_file(param)

    cell_files = list(Path(".").glob(cell))
    if not cell_files:
        raise click.ClickException(f"No files matched pattern: {cell}")

    structures = []
    struct_names = []
    cellinputs = []

    for cell_path in cell_files:
        struct_name = cell_path.stem
        struct_names.append(struct_name)

        if base_cell is not None:
            atoms = ase_read(str(cell_path))
            from airsspy.tools.modcell import modify_cell

            cell_lines = modify_cell(base_cell, atoms)
            cell_content = "\n".join(cell_lines)
            ci = CellInput.from_string(cell_content)
        elif cell_path.suffix == ".cell":
            cell_content = cell_path.read_text()
            ci = CellInput.from_string(cell_content)
        else:
            atoms = ase_read(str(cell_path))
            pmg_struct = AseAtomsAdaptor.get_structure(atoms)
            structures.append(pmg_struct)
            ci = CellInput()
            ci.set_cell(atoms.cell.array)
            ci.set_positions(atoms.get_chemical_symbols(), atoms.positions)

        cellinputs.append(ci)

        if not structures:
            atoms = ase_read(str(cell_path))
            pmg_struct = AseAtomsAdaptor.get_structure(atoms)
            structures.append(pmg_struct)

    if len(structures) != len(struct_names):
        # Re-read all structures
        structures = []
        for cell_path in cell_files:
            atoms = ase_read(str(cell_path))
            structures.append(AseAtomsAdaptor.get_structure(atoms))

    if dryrun:
        click.echo(f"Project: {project}, Seed: {seed}")
        click.echo(f"Structures: {len(struct_names)}")
        for name in struct_names:
            click.echo(f"  {name}")
        click.echo(f"Executable: {exe}, Code: {code}")
        click.echo(f"Cycles: {cycles}, Max iterations: {max_iterations}")
        return

    from jobflow import Flow, JobStore, run_locally
    from maggma.stores import MongoStore

    maker = AirssRelaxMaker(
        executable=exe,
        cycles=cycles,
        max_iterations=max_iterations,
        code=code,
    )
    job = maker.make(
        structures=structures,
        struct_names=struct_names,
        cellinputs=cellinputs,
        paraminput=paraminput,
        project_name=project,
        seed_name=seed,
    )

    flow = Flow([job])

    docs_store = MongoStore(
        database=ctx.obj["db_name"],
        collection_name="jobs",
        host=ctx.obj["db_host"],
        port=ctx.obj["db_port"],
    )
    store = JobStore(docs_store=docs_store)

    click.echo("Running flow...")
    run_locally(flow, store=store)
    click.echo("Done")
