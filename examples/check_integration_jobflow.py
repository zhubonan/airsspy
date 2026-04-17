"""
Integration test: jobflow search with real buildcell + mocked CASTEP + real MongoDB.

Tests the full pipeline:
  1. AirssSearchMaker creates a jobflow job
  2. run_buildcell generates real random structures
  3. CASTEP relaxation is mocked (would be too expensive)
  4. Results are stored in MongoDB via JobStore
  5. SearchStore can retrieve and query the results
"""

import sys
import tempfile
import os
from pathlib import Path
from unittest.mock import MagicMock, patch

# ── Seed file for a simple Si2 search ──
SEED_CELL = """%BLOCK LATTICE_CART
ANGSTROM
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Si 0.0 0.0 0.0
Si 0.5 0.5 0.5
%ENDBLOCK POSITIONS_FRAC

%BLOCK SPECIES
Si
%ENDBLOCK SPECIES

%BLOCK IONS
NUMAT 2
%ENDBLOCK IONS

%BLOCK CELL
SYMOPS 48
VOL 60-100
%ENDBLOCK CELL
"""

PARAM_CONTENT = """task: geometryoptimization
cut_off_energy: 200 eV
"""


def main():
    from castepinput.inputs import ParamInput
    from jobflow import Flow, JobStore, run_locally
    from maggma.stores import MongoStore

    from airsspy.jf.jobs import AirssSearchMaker
    from airsspy.jf.store import SearchStore

    tmpdir = tempfile.mkdtemp(prefix="airsspy_test_")
    old_cwd = os.getcwd()
    os.chdir(tmpdir)
    print(f"Working in: {tmpdir}")

    try:
        # 1. Write seed and param files
        seed_name = "Si2-test"
        Path(f"{seed_name}.cell").write_text(SEED_CELL)
        Path(f"{seed_name}.param").write_text(PARAM_CONTENT)

        seed_content = SEED_CELL
        Path(f"{seed_name}.param").write_text(PARAM_CONTENT)
        paraminput = ParamInput.from_file(f"{seed_name}.param")
        project_name = "test-jobflow-integration"

        # 2. Create the search maker
        n_structures = 3
        maker = AirssSearchMaker(n_structures=n_structures, code="castep")

        job = maker.make(
            seed_name=seed_name,
            seed_content=seed_content,
            paraminput=paraminput,
            project_name=project_name,
        )

        # 3. Set up MongoDB JobStore
        docs_store = MongoStore(
            database="airss",
            collection_name="jobs",
            host="localhost",
            port=27017,
        )
        store = JobStore(docs_store=docs_store)

        # 4. Run with real buildcell but mocked CASTEP relaxation
        with (
            patch("airsspy.jf.jobs.AirssCastepRelaxRunner") as mock_runner_cls,
            patch("airsspy.jf.jobs.compose_task_doc") as mock_compose,
        ):
            # Mock runner: pretend relaxation succeeded
            mock_runner = MagicMock()
            mock_runner.run.return_value = 0
            mock_runner_cls.return_value = mock_runner

            # Mock compose_task_doc: return plausible data
            mock_compose.return_value = {
                "structure": None,
                "volume": 80.0,
                "reduced_formula": "Si",
                "formula": "Si2",
                "natoms": 2,
                "energy": -10.0,
                "energy_per_atom": -5.0,
                "pressure": 0.1,
                "spin": 0.0,
                "mod_spin": 0.0,
                "symmetry": "(Fd-3m)",
                "res_content": "TITL Si2-test-001 Fd-3m\nCELL ...",
                "parallel_efficiency": 0.9,
                "total_time": 50.0,
            }

            print(f"\n=== Running {n_structures}-structure search ===")
            flow = Flow([job])
            responses = run_locally(flow, store=store, create_folders=False)

        # 5. Check the local response
        output = responses[job.uuid][1].output
        print(f"\n=== Local response ===")
        print(f"  project:  {output.project_name}")
        print(f"  seed:     {output.seed_name}")
        print(f"  job_type: {output.job_type}")
        print(f"  n_structures: {output.n_structures}")
        print(f"  n_finished:   {output.n_finished}")
        print(f"  n_errored:   {output.n_errored}")
        print(f"  n_failed:    {output.n_failed}")

        assert output.n_structures == n_structures, (
            f"Expected {n_structures} structures, got {output.n_structures}"
        )
        assert output.n_finished == n_structures, (
            f"Expected {n_structures} finished, got {output.n_finished}"
        )

        for i, r in enumerate(output.results):
            print(f"  [{i}] {r.struct_name}  status={r.relax_status}  "
                  f"E/A={r.energy_per_atom}  V={r.volume}")

        # 6. Query from MongoDB via SearchStore
        print(f"\n=== Querying MongoDB ===")
        with SearchStore(database="airss", host="localhost", port=27017) as ss:
            projects = ss.list_projects()
            print(f"  Projects: {projects}")
            assert project_name in projects, (
                f"Project '{project_name}' not found in {projects}"
            )

            results = ss.retrieve_project(project_name)
            print(f"  Retrieved {len(results)} result docs from DB")

            df = ss.retrieve_project_df(project_name)
            print(f"  DataFrame shape: {df.shape}")
            if not df.empty:
                print(f"  Columns: {list(df.columns)}")
                print(df[["struct_name", "relax_status", "energy_per_atom", "volume"]].to_string())

            seeds = ss.list_seeds(project_name)
            print(f"  Seeds: {seeds}")

            counts = ss.show_struct_counts(project_name)
            print(f"\n  Struct counts:")
            print(counts.to_string())

        print("\n=== ALL TESTS PASSED ===")
        return 0

    except Exception:
        import traceback
        traceback.print_exc()
        return 1
    finally:
        os.chdir(old_cwd)
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
