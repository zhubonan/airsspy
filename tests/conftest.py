# -*- coding: utf-8 -*-
###########################################################################
# airss-ase                                                               #
# Copyright (C) 2019  Bonan Zhu                                           #
#                                                                         #
# This program is free software; you can redistribute it and/or modify    #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation; either version 2 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# This program is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License along #
# with this program; if not, write to the Free Software Foundation, Inc., #
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             #
###########################################################################
"""
Test configuration
"""
import shutil
import sys
import os
import tempfile
from pathlib import Path
from unittest import mock

# Add src directory to Python path
src_path = Path(__file__).parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from ase import Atoms
from tempfile import mkstemp
import pytest


@pytest.fixture
def memory_jobstore():
    """In-memory JobStore for jobflow tests."""
    from jobflow import JobStore
    from maggma.stores import MemoryStore

    store = JobStore(MemoryStore(), additional_stores={"data": MemoryStore()})
    store.connect()
    return store


@pytest.fixture(autouse=True)
def mock_jobflow_settings(memory_jobstore):
    """Mock jobflow settings to use an in-memory JobStore."""
    from jobflow.settings import JobflowSettings

    settings = JobflowSettings(JOB_STORE=memory_jobstore)
    with mock.patch("jobflow.SETTINGS", settings):
        yield


@pytest.fixture(autouse=True)
def clean_dir():
    """Run each test in a fresh temporary working directory."""
    old_cwd = os.getcwd()
    new_path = tempfile.mkdtemp()
    os.chdir(new_path)
    yield
    os.chdir(old_cwd)
    shutil.rmtree(new_path)


@pytest.fixture
def al_atoms():
    return Atoms(
        "Al2", cell=[2, 2, 2], positions=[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], pbc=True
    )


@pytest.fixture
def tmpfile():
    fname = mkstemp()[1]
    yield fname
    os.remove(fname)


def pytest_addoption(parser):
    parser.addoption(
        "--run-e2e",
        action="store_true",
        default=False,
        help="Run end-to-end tests that call CASTEP/ABACUS executables",
    )


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--run-e2e"):
        skip_e2e = pytest.mark.skip(reason="Needs --run-e2e option to run")
        for item in items:
            if "e2e" in item.keywords:
                item.add_marker(skip_e2e)
