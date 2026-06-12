"""Tests for jf/documents module."""


from airsspy.jf.documents import AirssJobDoc, AirssResultDoc, RelaxOutcome


def test_relax_outcome_enum():
    """Test RelaxOutcome enum values."""
    assert RelaxOutcome.FINISHED == "finished"
    assert RelaxOutcome.ERRORED == "errored"
    assert RelaxOutcome.TIMEDOUT == "timedout"
    assert RelaxOutcome.UNDETERMINED == "undetermined"
    assert RelaxOutcome.CYCLE_EXCEEDED == "cycle_exceeded"


def test_airss_result_doc_defaults():
    """Test AirssResultDoc default values."""
    doc = AirssResultDoc(
        struct_name="test-123",
        seed_name="Si",
        project_name="test_project",
    )
    assert doc.struct_name == "test-123"
    assert doc.seed_name == "Si"
    assert doc.project_name == "test_project"
    assert doc.structure is None
    assert doc.energy is None
    assert doc.spin == 0.0
    assert doc.mod_spin == 0.0
    assert doc.relax_status == RelaxOutcome.FINISHED


def test_airss_result_doc_with_data():
    """Test AirssResultDoc with computed properties."""
    doc = AirssResultDoc(
        struct_name="Si-abc123",
        seed_name="Si",
        project_name="silicon",
        energy=-42.5,
        energy_per_atom=-10.625,
        volume=40.0,
        pressure=0.1,
        natoms=4,
        formula="Si4",
        reduced_formula="Si",
        symmetry="(Fd-3m)",
        relax_status=RelaxOutcome.FINISHED,
    )
    assert doc.energy == -42.5
    assert doc.energy_per_atom == -10.625
    assert doc.volume == 40.0
    assert doc.natoms == 4
    assert doc.symmetry == "(Fd-3m)"


def test_airss_job_doc_empty():
    """Test AirssJobDoc with no results."""
    doc = AirssJobDoc(
        project_name="test",
        seed_name="Si",
    )
    assert doc.project_name == "test"
    assert doc.seed_name == "Si"
    assert doc.job_type == "search"
    assert doc.results == []
    assert doc.n_structures == 0
    assert doc.n_finished == 0
    assert doc.n_errored == 0


def test_airss_job_doc_with_results():
    """Test AirssJobDoc with multiple results."""
    r1 = AirssResultDoc(
        struct_name="Si-001",
        seed_name="Si",
        project_name="test",
        energy=-40.0,
        natoms=4,
        relax_status=RelaxOutcome.FINISHED,
    )
    r2 = AirssResultDoc(
        struct_name="Si-002",
        seed_name="Si",
        project_name="test",
        energy=-38.0,
        natoms=4,
        relax_status=RelaxOutcome.FINISHED,
    )
    r3 = AirssResultDoc(
        struct_name="Si-003",
        seed_name="Si",
        project_name="test",
        relax_status=RelaxOutcome.ERRORED,
    )

    doc = AirssJobDoc(
        project_name="test",
        seed_name="Si",
        results=[r1, r2, r3],
        n_structures=3,
        n_finished=2,
        n_errored=1,
    )
    assert len(doc.results) == 3
    assert doc.n_structures == 3
    assert doc.n_finished == 2
    assert doc.n_errored == 1


def test_airss_job_doc_seed_info():
    """Test AirssJobDoc seed metadata."""
    doc = AirssJobDoc(
        project_name="test",
        seed_name="Si",
        seed_content="%BLOCK LATTICE_CART\n%ENDBLOCK LATTICE_CART",
        seed_hash="abc12345",
        param_content="task: geometryoptimization",
    )
    assert doc.seed_content is not None
    assert doc.seed_hash == "abc12345"
    assert doc.param_content is not None


def test_airss_job_doc_serialization():
    """Test that documents can be serialised and deserialised."""
    r = AirssResultDoc(
        struct_name="Si-001",
        seed_name="Si",
        project_name="test",
        energy=-40.0,
        natoms=4,
    )
    doc = AirssJobDoc(
        project_name="test",
        seed_name="Si",
        results=[r],
        n_structures=1,
        n_finished=1,
    )

    # Round-trip through dict
    data = doc.model_dump()
    doc2 = AirssJobDoc(**data)
    assert doc2.project_name == "test"
    assert len(doc2.results) == 1
    assert doc2.results[0].energy == -40.0
