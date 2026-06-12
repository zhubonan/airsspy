"""Tests for scheduler module."""

import os

import pytest

from airsspy.scheduler import SGE, Dummy, Scheduler, Slurm


def test_scheduler_get_scheduler_returns_dummy():
    """When not in Slurm or SGE, get_scheduler returns None or Dummy."""
    # In a non-job environment, get_scheduler tries Slurm, SGE, Dummy.
    # Dummy always returns is_in_job=True, so it will be returned.
    sched = Scheduler.get_scheduler()
    # Since we're not in a real job, but Dummy claims True
    assert sched is not None
    assert isinstance(sched, Dummy)


def test_dummy_scheduler():
    """Test Dummy scheduler properties."""
    dummy = Dummy()
    assert dummy.is_in_job is True
    assert dummy.job_id == "0"
    assert dummy.get_n_cpus() == 4
    assert dummy.get_remaining_seconds() == Dummy.DEFAULT_REMAINING_TIME


def test_slurm_not_in_job():
    """Test Slurm when not in a job."""
    slurm = Slurm()
    assert slurm.is_in_job is False
    assert slurm.job_id is None
    assert bool(slurm) is False


def test_sge_not_in_job():
    """Test SGE when not in a job."""
    sge = SGE()
    assert sge.is_in_job is False
    assert sge.job_id is None


def test_slurm_get_end_time_not_in_job():
    """Test that get_end_time returns None when not in a job."""
    slurm = Slurm()
    assert slurm.get_end_time() is None


def test_slurm_get_remaining_seconds_not_in_job():
    """Test remaining seconds is 0 when not in a job."""
    slurm = Slurm()
    assert slurm.get_remaining_seconds() == 0


def test_scheduler_base_not_implemented():
    """Test that base Scheduler raises NotImplementedError."""
    base = Scheduler()
    with pytest.raises(NotImplementedError):
        base.get_n_cpus()
    with pytest.raises(NotImplementedError):
        base.get_remaining_seconds()
    with pytest.raises(NotImplementedError):
        _ = base.job_id


def test_user_name():
    """Test user_name property reads from environment."""
    dummy = Dummy()
    assert dummy.user_name == os.environ.get("USER")
