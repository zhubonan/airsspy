"""Tests for SearchStore query layer."""

from datetime import datetime, timedelta
from unittest.mock import patch

import pytest

mongomock = pytest.importorskip("mongomock")

from airsspy.jf.documents import AirssJobDoc, AirssResultDoc, RelaxOutcome  # noqa: E402


def _make_job_doc(**overrides):
    defaults = {
        "project_name": "test_proj",
        "seed_name": "Si",
        "job_type": "search",
        "n_structures": 10,
        "n_finished": 8,
        "n_errored": 2,
        "results": [
            AirssResultDoc(
                struct_name="Si-001",
                seed_name="Si",
                project_name="test_proj",
                energy=-10.0,
                natoms=4,
                relax_status=RelaxOutcome.FINISHED,
            ),
            AirssResultDoc(
                struct_name="Si-002",
                seed_name="Si",
                project_name="test_proj",
                energy=-8.0,
                natoms=4,
                relax_status=RelaxOutcome.ERRORED,
            ),
        ],
    }
    defaults.update(overrides)
    return AirssJobDoc(**defaults)


def _make_store_doc(job_doc):
    """Wrap an AirssJobDoc in the dict format that MongoStore returns."""
    return {"output": job_doc.model_dump(), "uuid": "test-uuid"}


@patch("airsspy.jf.store.MongoStore")
class TestSearchStoreList:
    """Test project and seed listing methods."""

    def test_list_projects(self, MockMongoStore):
        mock_store = MockMongoStore.return_value
        mock_store.distinct.return_value = ["proj_b", "proj_a"]

        from airsspy.jf.store import SearchStore

        store = SearchStore()
        projects = store.list_projects()

        assert projects == ["proj_a", "proj_b"]
        mock_store.distinct.assert_called_once_with("output.project_name")

    def test_list_seeds_with_project(self, MockMongoStore):
        mock_store = MockMongoStore.return_value
        mock_store.distinct.return_value = ["Si", "C"]

        from airsspy.jf.store import SearchStore

        store = SearchStore()
        seeds = store.list_seeds(project_name="test_proj")

        assert seeds == ["C", "Si"]
        mock_store.distinct.assert_called_once_with(
            "output.seed_name", {"output.project_name": "test_proj"}
        )

    def test_list_seeds_no_filter(self, MockMongoStore):
        mock_store = MockMongoStore.return_value
        mock_store.distinct.return_value = ["Si", "C"]

        from airsspy.jf.store import SearchStore

        store = SearchStore()
        seeds = store.list_seeds()

        assert seeds == ["C", "Si"]
        mock_store.distinct.assert_called_once_with("output.seed_name", {})


@patch("airsspy.jf.store.MongoStore")
class TestSearchStoreRetrieve:
    """Test result retrieval methods."""

    def test_retrieve_project(self, MockMongoStore):
        job_doc = _make_job_doc()
        mock_store = MockMongoStore.return_value
        mock_store.query.return_value = [_make_store_doc(job_doc)]

        from airsspy.jf.store import SearchStore

        store = SearchStore()
        results = store.retrieve_project("test_proj")

        assert len(results) == 2
        assert results[0].struct_name == "Si-001"
        assert results[1].struct_name == "Si-002"
        mock_store.query.assert_called_once_with(criteria={"output.project_name": "test_proj"})

    def test_retrieve_project_bad_doc(self, MockMongoStore):
        mock_store = MockMongoStore.return_value
        mock_store.query.return_value = [{"output": {"bad_key": "value"}, "uuid": "bad"}]

        from airsspy.jf.store import SearchStore

        store = SearchStore()
        results = store.retrieve_project("bad_proj")

        assert results == []

    def test_retrieve_project_empty(self, MockMongoStore):
        mock_store = MockMongoStore.return_value
        mock_store.query.return_value = []

        from airsspy.jf.store import SearchStore

        store = SearchStore()
        results = store.retrieve_project("empty_proj")

        assert results == []

    def test_retrieve_project_df(self, MockMongoStore):
        job_doc = _make_job_doc()
        mock_store = MockMongoStore.return_value
        mock_store.query.return_value = [_make_store_doc(job_doc)]

        from airsspy.jf.store import SearchStore

        store = SearchStore()
        df = store.retrieve_project_df("test_proj")

        assert len(df) == 2
        assert "structure" not in df.columns
        assert "initial_structure" not in df.columns
        assert "energy" in df.columns

    def test_retrieve_project_df_empty(self, MockMongoStore):
        mock_store = MockMongoStore.return_value
        mock_store.query.return_value = []

        from airsspy.jf.store import SearchStore

        store = SearchStore()
        df = store.retrieve_project_df("empty_proj")

        assert len(df) == 0


@patch("airsspy.jf.store.MongoStore")
class TestSearchStoreSummary:
    """Test summary methods."""

    def test_show_struct_counts(self, MockMongoStore):
        job_doc = _make_job_doc()
        mock_store = MockMongoStore.return_value
        mock_store.query.return_value = [_make_store_doc(job_doc)]

        from airsspy.jf.store import SearchStore

        store = SearchStore()
        df = store.show_struct_counts("test_proj")

        assert len(df) == 1
        assert list(df.columns) == [
            "project_name",
            "seed_name",
            "job_type",
            "n_structures",
            "n_finished",
            "n_errored",
        ]
        assert df.iloc[0]["n_structures"] == 10
        assert df.iloc[0]["n_finished"] == 8

    def test_show_struct_counts_empty(self, MockMongoStore):
        mock_store = MockMongoStore.return_value
        mock_store.query.return_value = []

        from airsspy.jf.store import SearchStore

        store = SearchStore()
        df = store.show_struct_counts("empty")

        assert len(df) == 0

    def test_throughput_summary(self, MockMongoStore):
        now = datetime.now()
        yesterday = now - timedelta(days=1)

        job_doc_1 = _make_job_doc(
            n_structures=5,
            n_finished=4,
            created_on=yesterday.replace(hour=10),
        )
        job_doc_2 = _make_job_doc(
            project_name="other_proj",
            n_structures=3,
            n_finished=2,
            created_on=yesterday.replace(hour=15),
        )
        mock_store = MockMongoStore.return_value
        mock_store.query.return_value = [
            _make_store_doc(job_doc_1),
            _make_store_doc(job_doc_2),
        ]

        from airsspy.jf.store import SearchStore

        store = SearchStore()
        df = store.throughput_summary(past_days=2)

        assert len(df) == 1
        day_key = yesterday.strftime("%Y-%m-%d")
        assert df.iloc[0]["date"] == day_key
        assert df.iloc[0]["n_structures"] == 8
        assert df.iloc[0]["n_finished"] == 6

    def test_throughput_summary_empty(self, MockMongoStore):
        mock_store = MockMongoStore.return_value
        mock_store.query.return_value = []

        from airsspy.jf.store import SearchStore

        store = SearchStore()
        df = store.throughput_summary()

        assert len(df) == 0
