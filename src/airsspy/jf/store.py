"""
Query layer for AIRSS search results stored in jobflow JobStore.

Provides ``SearchStore`` — a high-level interface for querying,
filtering, and summarising AIRSS results without direct MongoDB
knowledge.
"""

import logging
from datetime import datetime, timedelta
from typing import Optional

import pandas as pd
from maggma.stores import MongoStore

from .documents import AirssJobDoc, AirssResultDoc

logger = logging.getLogger(__name__)


class SearchStore:
    """Query interface for AIRSS results in a jobflow JobStore.

    Wraps a ``MongoStore`` (maggma) and provides project-oriented
    query methods that flatten multi-structure job documents into
    individual result records.
    """

    def __init__(
        self,
        database: str = "airss",
        host: str = "localhost",
        port: int = 27017,
        collection: str = "jobs",
        username: Optional[str] = None,
        password: Optional[str] = None,
        **kwargs,
    ) -> None:
        """
        Connect to the jobflow store.

        Args:
            database: MongoDB database name.
            host: MongoDB host.
            port: MongoDB port.
            collection: Collection name for jobflow jobs.
            username: Optional MongoDB username.
            password: Optional MongoDB password.
            **kwargs: Additional arguments passed to MongoStore.
        """
        self.docs_store = MongoStore(
            database=database,
            collection_name=collection,
            host=host,
            port=port,
            username=username,
            password=password,
            **kwargs,
        )

    def connect(self) -> None:
        """Open the connection to the MongoDB store."""
        self.docs_store.connect()

    def close(self) -> None:
        """Close the connection."""
        self.docs_store.close()

    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, *args):
        self.close()

    def retrieve_project(self, project_name: str) -> list[AirssResultDoc]:
        """
        Get all results for a project, flattened from multi-structure jobs.

        Args:
            project_name: The project identifier.

        Returns:
            A list of ``AirssResultDoc`` instances.
        """
        results = []
        criteria = {"output.project_name": project_name}
        for doc in self.docs_store.query(criteria=criteria):
            output = doc.get("output", {})
            try:
                job_doc = AirssJobDoc(**output)
                results.extend(job_doc.results)
            except (TypeError, KeyError, ValueError):
                logger.warning("Failed to parse job document %s", doc.get("uuid"), exc_info=True)
        return results

    def retrieve_project_df(
        self, project_name: str, **filters
    ) -> pd.DataFrame:
        """
        Get all results for a project as a DataFrame.

        Args:
            project_name: The project identifier.
            **filters: Additional MongoDB query filters.

        Returns:
            A DataFrame with one row per result.
        """
        results = self.retrieve_project(project_name)
        if not results:
            return pd.DataFrame()
        records = [r.model_dump() for r in results]

        df = pd.DataFrame(records)
        # Drop columns with unserialisable structure objects for display
        struct_cols = ["structure", "initial_structure"]
        for col in struct_cols:
            if col in df.columns:
                df = df.drop(columns=[col])
        return df

    def list_projects(self) -> list[str]:
        """List all distinct project names in the store."""
        return sorted(self.docs_store.distinct("output.project_name"))

    def list_seeds(self, project_name: Optional[str] = None) -> list[str]:
        """
        List seed names, optionally filtered by project.

        Args:
            project_name: If given, filter seeds to this project.
        """
        criteria: dict = {}
        if project_name:
            criteria["output.project_name"] = project_name
        return sorted(self.docs_store.distinct("output.seed_name", criteria))

    def show_struct_counts(self, project_name: Optional[str] = None) -> pd.DataFrame:
        """
        Summary of structures per project/seed.

        Args:
            project_name: If given, filter to this project.

        Returns:
            DataFrame with columns project_name, seed_name,
            n_structures, n_finished, n_errored.
        """
        criteria: dict = {}
        if project_name:
            criteria["output.project_name"] = project_name

        records = []
        for doc in self.docs_store.query(criteria=criteria):
            output = doc.get("output", {})
            try:
                job_doc = AirssJobDoc(**output)
            except (TypeError, KeyError, ValueError):
                continue
            records.append({
                "project_name": job_doc.project_name,
                "seed_name": job_doc.seed_name,
                "job_type": job_doc.job_type,
                "n_structures": job_doc.n_structures,
                "n_finished": job_doc.n_finished,
                "n_errored": job_doc.n_errored,
            })

        if not records:
            return pd.DataFrame()
        return pd.DataFrame(records)

    def throughput_summary(
        self, past_days: int = 2, project_name: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Throughput over recent days.

        Args:
            past_days: Number of days to look back.
            project_name: If given, filter to this project.

        Returns:
            DataFrame with date, n_structures, n_finished columns.
        """
        cutoff = datetime.now() - timedelta(days=past_days)
        criteria: dict = {"output.created_on": {"$gte": cutoff.isoformat()}}
        if project_name:
            criteria["output.project_name"] = project_name

        daily: dict[str, dict[str, int]] = {}
        for doc in self.docs_store.query(criteria=criteria):
            output = doc.get("output", {})
            try:
                job_doc = AirssJobDoc(**output)
            except (TypeError, KeyError, ValueError):
                continue
            created = job_doc.created_on
            if created is None:
                continue
            day_key = created.strftime("%Y-%m-%d")
            if day_key not in daily:
                daily[day_key] = {"n_structures": 0, "n_finished": 0}
            daily[day_key]["n_structures"] += job_doc.n_structures
            daily[day_key]["n_finished"] += job_doc.n_finished

        if not daily:
            return pd.DataFrame()
        records = [
            {"date": k, **v} for k, v in sorted(daily.items())
        ]
        return pd.DataFrame(records)
