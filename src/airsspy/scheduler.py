"""
Job scheduler interface for Slurm, SGE, and local execution.

Provides a unified interface for detecting the current scheduler
environment, querying job metadata (CPUs, remaining walltime),
and managing job arrays.
"""

import logging
import os
import re
import subprocess
import tempfile
from datetime import datetime, timedelta, timezone
from typing import ClassVar, Optional

logger = logging.getLogger(__name__)


class Scheduler:
    """Base class for job scheduler interfaces."""

    def __init__(self) -> None:
        self._job_id: Optional[str] = None
        self._ncpus: Optional[int] = None

    def get_n_cpus(self) -> Optional[int]:
        """Return the number of CPUs allocated to this job."""
        raise NotImplementedError

    @property
    def user_name(self) -> str:
        """Return the name of the current user."""
        return os.environ["USER"]

    def get_remaining_seconds(self) -> int:
        """Get the remaining walltime in seconds."""
        raise NotImplementedError

    @property
    def is_in_job(self) -> bool:
        """Return whether we are inside a scheduler job."""
        return self.job_id is not None

    @property
    def job_id(self) -> Optional[str]:
        """Return the job ID."""
        raise NotImplementedError

    @classmethod
    def get_scheduler(cls) -> Optional["Scheduler"]:
        """
        Detect and return a scheduler instance for the current environment.

        Tries Slurm, SGE, then Dummy. Returns None if not in any job
        and Dummy is not appropriate.
        """
        for trial in [Slurm, SGE, Dummy]:
            obj = trial()
            if obj.is_in_job:
                return obj
        return None


class Dummy(Scheduler):
    """Dummy scheduler for local execution."""

    DEFAULT_REMAINING_TIME = 3600 * 24 * 30  # 30 days

    def __init__(self) -> None:
        super().__init__()
        self._job_id = "0"

    def get_n_cpus(self) -> int:
        return 4

    @property
    def job_id(self) -> str:
        return self._job_id

    def get_remaining_seconds(self) -> int:
        """Get the remaining time. Defaults to 30 days."""
        return self.DEFAULT_REMAINING_TIME

    @property
    def is_in_job(self) -> bool:
        return True


class SGE(Scheduler):
    """Scheduler object for Sun Grid Engine (SGE)."""

    def __init__(self) -> None:
        super().__init__()
        self._start_time: Optional[datetime] = None
        if self.is_in_job:
            self._read_task_info()

    @property
    def job_id(self) -> Optional[str]:
        if self._job_id is None:
            job_id = os.environ.get("JOB_ID")
            task_id = os.environ.get("SGE_TASK_ID")
            if job_id:
                if task_id and task_id != "undefined":
                    job_id = job_id + "." + task_id
                    logger.warning("WARNING: REMAINING TIME IS NOT CORRECT FOR TASK ARRAY")
                self._job_id = job_id
        return self._job_id

    @property
    def is_in_job(self) -> bool:
        return "JOB_ID" in os.environ

    def _read_task_info(self) -> None:
        """Read detailed task information from qstat."""
        raw_data = subprocess.check_output(
            ["qstat", "-j", str(self.job_id)], universal_newlines=True
        )
        raw_data = raw_data.split("\n")
        task_info: dict[str, str] = {}
        for line in raw_data[1:]:
            try:
                key, value = line.split(":", maxsplit=1)
            except ValueError:
                continue
            task_info[key.strip()] = value.strip()
        self._task_info = task_info

    def get_n_cpus(self) -> Optional[int]:
        """Get the number of CPUs from the NSLOTS environment variable."""
        nslots = os.environ.get("NSLOTS")
        if nslots:
            return int(nslots)
        return None

    def get_max_run_seconds(self) -> Optional[int]:
        """Return the maximum run time in seconds."""
        rlist = self._task_info["hard resource_list"]
        match = re.search(r"h_rt=(\d+)", rlist)
        if match:
            return int(match.group(1))
        return None

    def get_end_time(self) -> Optional[datetime]:
        """Return the expected finish time of this job."""
        start = self.get_start_time()
        max_run = self.get_max_run_seconds()
        if start and max_run:
            return start + timedelta(seconds=max_run)
        return None

    def get_start_time(self) -> Optional[datetime]:
        """Return the start time of this job."""
        output = subprocess.check_output(
            ["qstat", "-j", str(self.job_id), "-xml"], universal_newlines=True
        )
        match = re.search(r"<JAT_start_time>(.+)</JAT_start_time>", output)
        if match:
            raw = match.group(1)
            time_int = int(raw)
            start_time = datetime.fromtimestamp(time_int, tz=timezone.utc)
            self._start_time = start_time
        return self._start_time

    def get_remaining_seconds(self) -> int:
        """Return the remaining walltime in seconds."""
        end_time = self.get_end_time()
        if end_time is None:
            return 0
        tdelta = end_time - datetime.now().astimezone()
        return int(tdelta.total_seconds())


class Slurm(Scheduler):
    """Slurm scheduler interface."""

    _task_info: ClassVar[Optional[dict[str, str]]] = None
    _warning: ClassVar[int] = 0

    def __init__(self) -> None:
        super().__init__()
        self.task_info: dict[str, str] = {}
        if Slurm._task_info is None:
            self._read_task_info()
            Slurm._task_info = self.task_info
        else:
            self.task_info = Slurm._task_info

    @property
    def is_in_job(self) -> bool:
        """Whether we are inside a Slurm job."""
        return "SLURM_JOB_ID" in os.environ

    @property
    def job_id(self) -> Optional[str]:
        if self._job_id is None:
            self._job_id = os.environ.get("SLURM_JOB_ID")
        return self._job_id

    def _read_task_info(self) -> None:
        """Extract job information from scontrol."""
        if not self.is_in_job:
            if Slurm._warning == 0:
                logger.debug("NOT STARTED FROM SLURM")
                Slurm._warning += 1
            self.task_info = {}
            return

        sinfo_dict: dict[str, str] = {}
        with tempfile.TemporaryFile(mode="w+") as tmp_file:
            subprocess.run(
                ["scontrol", "show", f"jobid={self.job_id}"],
                check=True,
                stdout=tmp_file,
            )
            tmp_file.seek(0)
            for line in tmp_file:
                for pair in line.split():
                    pair_s = pair.split("=", maxsplit=2)
                    if len(pair_s) == 2:
                        sinfo_dict[pair_s[0]] = pair_s[1]
                    elif len(pair_s) == 1:
                        sinfo_dict[pair_s[0]] = ""

        Slurm._task_info = sinfo_dict
        self.task_info = sinfo_dict

    def get_end_time(self) -> Optional[datetime]:
        """Return the end time of this job."""
        if self.task_info:
            end_time = datetime.strptime(
                self.task_info["EndTime"], "%Y-%m-%dT%H:%M:%S"
            )
            return end_time.astimezone()
        return None

    def get_remaining_seconds(self) -> int:
        """Return the remaining walltime in seconds."""
        end_time = self.get_end_time()
        if end_time is None:
            return 0
        return int((end_time - datetime.now().astimezone()).total_seconds())

    def get_user_name(self) -> Optional[str]:
        """Parse the user name from task info."""
        if self.task_info:
            match = re.match(r"([A-Za-z]+[0-9]+)\([0-9]+\)", self.task_info["UserId"])
            if match:
                return match.group(1)
        return None

    def get_n_cpus(self) -> Optional[str]:
        """Return number of CPUs allocated."""
        return self.task_info.get("NumCPUs")

    def get_array_id(self) -> Optional[str]:
        """Return the array job ID, or None."""
        return self.task_info.get("ArrayJobId")

    def get_array_task_id(self) -> Optional[str]:
        """Return the array task ID, or None."""
        return self.task_info.get("ArrayTaskId")

    def get_array_job_id(self) -> Optional[str]:
        """Return the array job ID as ``{array_id}_{task_id}``."""
        task = self.task_info
        if not task:
            return None
        try:
            return f"{task['ArrayJobId']}_{task['ArrayTaskId']}"
        except KeyError:
            return None

    def hold_array(self, array_num: Optional[str] = None) -> None:
        """Hold the array this job belongs to."""
        array_id = self.get_array_job_id()
        if array_num:
            array = array_num
        elif array_id:
            array = array_id.split("_")[0]
        else:
            array = None
        if array is not None:
            proc = subprocess.run(["scontrol", "hold", str(array)], check=False)
            if proc.returncode == 0:
                logger.info("Successfully held array %s", array)
            else:
                logger.error("Cannot hold array %s", array)
        else:
            logger.error("Cannot find the array to hold")

    def hold_all_pd_arrays(self, user_name: Optional[str] = None) -> None:
        """Hold all pending arrays."""
        arrays = self.get_pd_arrays(user_name)
        if not arrays:
            logger.warning("No array found to hold")
            return
        logger.info("Holding all pending arrays")
        array_list = ",".join(arrays)
        subprocess.run(["scontrol", "hold", array_list], check=True)

    def release_all_pd_arrays(self, user_name: Optional[str] = None) -> None:
        """Release all pending arrays."""
        arrays = self.get_pd_arrays(user_name)
        if not arrays:
            logger.warning("No array found to release")
            return
        array_list = ",".join(arrays)
        subprocess.run(["scontrol", "release", array_list], check=True)

    def get_running_jobs(self, user_name: Optional[str] = None) -> Optional[list[str]]:
        """Return a list of running job IDs for the current user."""
        return self._get_id_of_state("R", "%A %t", user_name)

    def get_pd_arrays(self, user_name: Optional[str] = None) -> Optional[list[str]]:
        """Return a list of pending array job IDs."""
        return self._get_id_of_state("PD", "%F %t", user_name)

    def _get_id_of_state(
        self, criteria: str, fmt_str: str, user_name: Optional[str] = None
    ) -> Optional[list[str]]:
        """Get IDs for jobs matching a given state."""
        if not self.task_info and not user_name:
            return None

        if not user_name:
            user_name = self.get_user_name()
        if not user_name:
            return None

        res = tempfile.TemporaryFile(mode="w+")
        subprocess.run(
            ["squeue", "-u", user_name, f"-o{fmt_str}"],
            check=True,
            stdout=res,
        )
        res.seek(0)
        task_ids: list[str] = []
        for line in res:
            sline = line.split()
            if len(sline) >= 2 and sline[1] == criteria:
                task_ids.append(sline[0])
        res.close()
        return task_ids

    def __bool__(self) -> bool:
        return bool(self.task_info)

    def __getitem__(self, key: str) -> Optional[str]:
        if self.task_info:
            return self.task_info[key]
        return None

    def __contains__(self, key: str) -> bool:
        return key in self.task_info
