"""Tests for structure shaking helpers."""

import pytest

from airsspy.shake import collect_shake_input_paths


def test_collect_shake_input_paths_combines_positional_and_filelist(tmp_path):
    """Positional files and --filelist entries are combined in order."""
    first = tmp_path / "first.res"
    second = tmp_path / "second.res"
    first.write_text("TITL first 0 1 0 0 0 1 (P1) n - 1\nEND\n")
    second.write_text("TITL second 0 1 0 0 0 1 (P1) n - 1\nEND\n")
    filelist = tmp_path / "inputs.txt"
    filelist.write_text("# comments are allowed\nsecond.res\n\n")

    paths = collect_shake_input_paths(
        [str(first)],
        filelist=str(filelist),
        base_dir=tmp_path,
    )

    assert paths == [first.resolve(), second.resolve()]


def test_collect_shake_input_paths_rejects_missing_file(tmp_path):
    """Missing filelist entries should fail before starting a relaxation."""
    filelist = tmp_path / "inputs.txt"
    filelist.write_text("missing.res\n")

    with pytest.raises(FileNotFoundError, match="missing.res"):
        collect_shake_input_paths([], filelist=str(filelist), base_dir=tmp_path)


def test_collect_shake_input_paths_resolves_filelist_entries_from_filelist_dir(tmp_path):
    """Relative --filelist entries are resolved beside the filelist itself."""
    inputs = tmp_path / "nested"
    inputs.mkdir()
    source = inputs / "source.res"
    source.write_text("TITL source 0 1 0 0 0 1 (P1) n - 1\nEND\n")
    filelist = inputs / "inputs.txt"
    filelist.write_text("source.res\n")

    paths = collect_shake_input_paths([], filelist=filelist, base_dir=tmp_path)

    assert paths == [source.resolve()]
