"""Tests for CLI commands."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

from ase import Atoms
from ase.io import write
from click.testing import CliRunner

from airsspy.cli import cmd_run
from airsspy.cli.main import cli


def test_cli_help():
    """Test the main CLI help output."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "AIRSS" in result.output
    assert "deploy" in result.output
    assert "db" in result.output
    assert "check" in result.output
    assert "tools" in result.output


def test_cli_version():
    """Test --version flag."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert "0.1.4" in result.output


def test_check_airss():
    """Test 'check airss' subcommand."""
    runner = CliRunner()
    result = runner.invoke(cli, ["check", "airss", "--help"])
    assert result.exit_code == 0


def test_check_scheduler():
    """Test 'check scheduler --help'."""
    runner = CliRunner()
    result = runner.invoke(cli, ["check", "scheduler", "--help"])
    assert result.exit_code == 0


def test_deploy_search_dryrun(tmp_path, monkeypatch):
    """Test 'deploy search' dryrun prints config."""
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path):
        # Create dummy seed files
        with open("Si.cell", "w") as f:
            f.write(
                "%BLOCK LATTICE_CART\n5.43 0 0\n0 5.43 0\n0 0 5.43\n%ENDBLOCK LATTICE_CART\n"
            )
        with open("Si.param", "w") as f:
            f.write("task: geometryoptimization\ncut_off_energy: 300\n")

        result = runner.invoke(
            cli,
            [
                "deploy",
                "search",
                "--seed",
                "Si",
                "--project",
                "test",
                "--num",
                "10",
                "--dryrun",
            ],
        )
        assert result.exit_code == 0
        assert "Project: test, Seed: Si" in result.output
        assert "Structures requested: 10" in result.output
        assert "task: geometryoptimization" in result.output


def test_deploy_search_vasp_dryrun_uses_incar(tmp_path):
    """VASP deploy search dryrun reads INCAR text and selects VASP executable."""
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path):
        Path("Si.cell").write_text(
            "%BLOCK LATTICE_CART\n5.43 0 0\n0 5.43 0\n0 0 5.43\n"
            "%ENDBLOCK LATTICE_CART\n"
        )
        Path("Si.INCAR").write_text("ENCUT = 400\n")

        result = runner.invoke(
            cli,
            [
                "deploy",
                "search",
                "--seed",
                "Si",
                "--project",
                "test",
                "--num",
                "1",
                "--code",
                "vasp",
                "--dryrun",
            ],
        )

    assert result.exit_code == 0
    assert "Executable: vasp_std, Code: vasp" in result.output
    assert "ENCUT = 400" in result.output


def test_deploy_search_rejects_unknown_code():
    """Deploy validates code choices before indexing suffix maps."""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "deploy",
            "search",
            "--seed",
            "Si",
            "--project",
            "test",
            "--num",
            "1",
            "--code",
            "unknown",
        ],
    )

    assert result.exit_code != 0
    assert "Invalid value for '--code'" in result.output


def test_db_commands_help():
    """Test all db subcommands have valid help."""
    runner = CliRunner()
    for cmd in [
        "list-projects",
        "list-seeds",
        "summary",
        "throughput",
        "retrieve-project",
    ]:
        result = runner.invoke(cli, ["db", cmd, "--help"])
        assert result.exit_code == 0, f"db {cmd} --help failed"


def test_tools_modcell_help():
    """Test tools modcell help."""
    runner = CliRunner()
    result = runner.invoke(cli, ["tools", "modcell", "--help"])
    assert result.exit_code == 0
    assert "BASE_CELL" in result.output


def test_rank_help():
    """Test 'rank --help'."""
    runner = CliRunner()
    result = runner.invoke(cli, ["rank", "--help"])
    assert result.exit_code == 0
    assert "enthalpy" in result.output.lower()
    assert "-h" in result.output
    assert "-?" in result.output
    assert "--help" in result.output
    assert "-r, --rank" in result.output
    assert "--denergy" in result.output
    assert "--not_relative" in result.output
    assert "--long" in result.output
    assert "-fu, --formula-unit" in result.output
    assert "-sn, --speciesnumber" in result.output
    assert "-in, --ionsnumber" in result.output
    assert "-dr, --distance" in result.output


def test_rank_cryan_help_aliases():
    """Test cryan-style help aliases."""
    runner = CliRunner()
    for alias in ("-h", "-?"):
        result = runner.invoke(cli, ["rank", alias])
        assert result.exit_code == 0
        assert "Usage: " in result.output


def test_rank_from_stdin():
    """Test ranking structures piped via stdin."""
    runner = CliRunner()
    packed_res = (
        "TITL Si-001 -0.05 40.0 -42.500000 0 0 4 (Fd-3m) n - 1\n"
        "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "Si     1  0.25 0.25 0.25 1.0\n"
        "Si     1  0.5  0.5  0.5  1.0\n"
        "Si     1  0.75 0.75 0.75 1.0\n"
        "END\n"
        "TITL Si-002 -0.05 42.0 -43.200000 0 0 4 (Pm-3m) n - 1\n"
        "CELL 1.0  5.50 5.50 5.50 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "Si     1  0.5  0.5  0.5  1.0\n"
        "Si     1  0.25 0.25 0.25 1.0\n"
        "Si     1  0.75 0.75 0.75 1.0\n"
        "END\n"
    )
    result = runner.invoke(cli, ["rank"], input=packed_res)
    assert result.exit_code == 0
    # Most stable (Si-002, lower enthalpy) should appear first
    assert "Si-002" in result.output


def test_rank_from_file(tmp_path):
    """Test ranking structures from a file argument."""
    runner = CliRunner()
    res_file = tmp_path / "test.res"
    res_file.write_text(
        "TITL Si-001 -0.05 40.0 -42.500000 0 0 4 (Fd-3m) n - 1\n"
        "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "Si     1  0.25 0.25 0.25 1.0\n"
        "Si     1  0.5  0.5  0.5  1.0\n"
        "Si     1  0.75 0.75 0.75 1.0\n"
        "END\n"
    )
    result = runner.invoke(cli, ["rank", str(res_file)])
    assert result.exit_code == 0
    assert "Si-001" in result.output


def test_rank_res_file_drops_raw_lines_when_not_merging(tmp_path, monkeypatch):
    """Plain ranking avoids retaining full RES blocks in memory."""
    from airsspy.ranking import read_res_file as real_read_res_file

    seen_keep_raw = []

    def fake_read_res_file(path, keep_raw=True):
        seen_keep_raw.append(keep_raw)
        return real_read_res_file(path, keep_raw=keep_raw)

    monkeypatch.setattr("airsspy.ranking.read_res_file", fake_read_res_file)
    res_file = tmp_path / "test.res"
    res_file.write_text(_single_atom_res("Si-001", -1.0))

    runner = CliRunner()
    result = runner.invoke(cli, ["rank", str(res_file)])

    assert result.exit_code == 0
    assert seen_keep_raw == [False]


def test_rank_keeps_raw_lines_for_unite(tmp_path, monkeypatch):
    """Fingerprint merging still keeps raw RES blocks."""
    from airsspy.ranking import read_res_file as real_read_res_file

    seen_keep_raw = []

    def fake_read_res_file(path, keep_raw=True):
        seen_keep_raw.append(keep_raw)
        return real_read_res_file(path, keep_raw=keep_raw)

    monkeypatch.setattr("airsspy.ranking.read_res_file", fake_read_res_file)
    res_file = tmp_path / "test.res"
    res_file.write_text(_single_atom_res("Si-001", -1.0))

    runner = CliRunner()
    result = runner.invoke(cli, ["rank", "-u", "0.1", str(res_file)])

    assert result.exit_code == 0
    assert seen_keep_raw == [True]


def test_rank_summary_mode():
    """Test -s summary flag."""
    runner = CliRunner()
    packed_res = (
        "TITL Si-001 -0.05 40.0 -42.500000 0 0 4 (Fd-3m) n - 1\n"
        "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "Si     1  0.25 0.25 0.25 1.0\n"
        "Si     1  0.5  0.5  0.5  1.0\n"
        "Si     1  0.75 0.75 0.75 1.0\n"
        "END\n"
        "TITL Si-002 -0.05 42.0 -43.200000 0 0 4 (Pm-3m) n - 1\n"
        "CELL 1.0  5.50 5.50 5.50 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "Si     1  0.5  0.5  0.5  1.0\n"
        "Si     1  0.25 0.25 0.25 1.0\n"
        "Si     1  0.75 0.75 0.75 1.0\n"
        "END\n"
    )
    result = runner.invoke(cli, ["rank", "-s"], input=packed_res)
    assert result.exit_code == 0
    # Summary should show only 1 Si structure (most stable)
    assert "Number of structures" in result.output
    assert "Number of compositions" in result.output


def test_rank_no_structures():
    """Test rank with empty stdin and no files."""
    runner = CliRunner()
    result = runner.invoke(cli, ["rank"], input="")
    assert result.exit_code == 0
    assert "No structures found" in result.output


def test_rank_formula_filter():
    """Test -f formula filter."""
    runner = CliRunner()
    packed_res = (
        "TITL Si-001 -0.05 40.0 -42.500000 0 0 4 (Fd-3m) n - 1\n"
        "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "Si     1  0.25 0.25 0.25 1.0\n"
        "Si     1  0.5  0.5  0.5  1.0\n"
        "Si     1  0.75 0.75 0.75 1.0\n"
        "END\n"
        "TITL SiO2-001 0.0 45.0 -80.500000 0 0 6 (P1) n - 1\n"
        "CELL 1.0  4.0 5.0 6.0 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si O\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "O      2  0.5  0.5  0.5  1.0\n"
        "O      2  0.3  0.3  0.3  1.0\n"
        "Si     1  0.6  0.6  0.6  1.0\n"
        "O      2  0.1  0.1  0.1  1.0\n"
        "O      2  0.8  0.8  0.8  1.0\n"
        "END\n"
    )
    result = runner.invoke(cli, ["rank", "-f", "SiO2"], input=packed_res)
    assert result.exit_code == 0
    assert "SiO2-001" in result.output
    assert "Si-001" not in result.output


def test_rank_formula_filter_accepts_reordered_formula():
    """Test -f matches formulas independent of element order."""
    runner = CliRunner()
    packed_res = _res_block("SiO2-001", -80.5, ["Si", "O", "O"])
    result = runner.invoke(cli, ["rank", "-f", "O2Si"], input=packed_res)
    assert result.exit_code == 0
    assert "SiO2-001" in result.output


def _single_atom_res(label, energy, element="Si"):
    return (
        f"TITL {label} 0.0 10.0 {energy:.6f} 0 0 1 (P1) n - 1\n"
        "CELL 1.0  3.0 3.0 3.0 90.0 90.0 90.0\n"
        "LATT -1\n"
        f"SFAC {element}\n"
        f"{element}     1  0.0  0.0  0.0  1.0\n"
        "END\n"
    )


def _res_block(label, energy, symbols):
    species = []
    for sym in symbols:
        if sym not in species:
            species.append(sym)

    lines = [
        f"TITL {label} 0.0 10.0 {energy:.6f} 0 0 {len(symbols)} (P1) n - 1",
        "CELL 1.0  3.0 3.0 3.0 90.0 90.0 90.0",
        "LATT -1",
        "SFAC " + " ".join(species),
    ]
    for i, sym in enumerate(symbols):
        sfac_index = species.index(sym) + 1
        coord = i / 10.0
        lines.append(f"{sym}     {sfac_index}  {coord:.1f}  0.0  0.0  1.0")
    lines.append("END")
    return "\n".join(lines) + "\n"


def test_rank_cryan_rank_alias_behaves_like_default():
    """Test -r is accepted as a no-op rank compatibility flag."""
    runner = CliRunner()
    packed_res = _single_atom_res("Si-001", -1.0)

    default = runner.invoke(cli, ["rank"], input=packed_res)
    compat = runner.invoke(cli, ["rank", "-r"], input=packed_res)

    assert default.exit_code == 0
    assert compat.exit_code == 0
    assert default.stdout == compat.stdout


def test_rank_cryan_option_aliases_are_accepted():
    """Test common cryan-style aliases parse successfully."""
    runner = CliRunner()
    packed_res = _single_atom_res("Si-001", -1.0)

    result = runner.invoke(
        cli,
        [
            "rank",
            "--not_relative",
            "--long",
            "--denergy",
            "1.0",
            "--elementlist",
            "Si,O",
            "-dr",
            "4.0",
        ],
        input=packed_res,
    )

    assert result.exit_code == 0
    assert "Si-001" in result.stdout


def test_rank_cryan_formula_unit_filter():
    """Test -fu filters by exact number of formula units."""
    runner = CliRunner()
    packed_res = (
        _res_block("Si4", -4.0, ["Si", "Si", "Si", "Si"])
        + _res_block("Si2", -2.0, ["Si", "Si"])
        + _res_block("SiO2", -6.0, ["Si", "O", "O"])
    )

    result = runner.invoke(cli, ["rank", "-fu", "4"], input=packed_res)

    assert result.exit_code == 0
    assert "Si4" in result.stdout
    assert "Si2" not in result.stdout
    assert "SiO2" not in result.stdout


def test_rank_cryan_species_and_ions_filters():
    """Test -sn and -in filters from cryan-style options."""
    runner = CliRunner()
    packed_res = (
        _res_block("Si4", -4.0, ["Si", "Si", "Si", "Si"])
        + _res_block("Si2", -2.0, ["Si", "Si"])
        + _res_block("SiO2", -6.0, ["Si", "O", "O"])
    )

    species = runner.invoke(cli, ["rank", "-sn", "2"], input=packed_res)
    ions = runner.invoke(cli, ["rank", "-in", "-3"], input=packed_res)

    assert species.exit_code == 0
    assert "SiO2" in species.stdout
    assert "Si4" not in species.stdout
    assert "Si2" not in species.stdout
    assert ions.exit_code == 0
    assert "Si2" in ions.stdout
    assert "SiO2" in ions.stdout
    assert "Si4" not in ions.stdout


def test_rank_pathological_prune_hides_rejected_from_stdout():
    """Test trimmed-MAD pathological pruning filters rank output."""
    try:
        runner = CliRunner(mix_stderr=False)
    except TypeError:
        runner = CliRunner()
    packed_res = "".join(
        [
            _single_atom_res("Si-pathological", -10.0),
            _single_atom_res("Si-low", -5.2),
            _single_atom_res("Si-a", -5.1),
            _single_atom_res("Si-b", -5.0),
            _single_atom_res("Si-c", -4.9),
            _single_atom_res("Si-d", -4.8),
        ]
    )

    result = runner.invoke(
        cli,
        [
            "rank",
            "--prune-pathological",
            "--pathology-tail-fraction",
            "1.0",
            "--pathology-min-tail-size",
            "3",
        ],
        input=packed_res,
    )

    assert result.exit_code == 0
    assert "Si-pathological" not in result.stdout
    assert "Si-low" in result.stdout
    assert "Pathological prune: kept 5, rejected 1" in result.stderr
    assert "Rejected pathological structures: Si-pathological" in result.stderr


def test_rank_pathological_prune_absent_keeps_existing_output():
    """Test rank output is unchanged unless pathological pruning is enabled."""
    runner = CliRunner()
    packed_res = "".join(
        [
            _single_atom_res("Si-pathological", -10.0),
            _single_atom_res("Si-low", -5.2),
            _single_atom_res("Si-a", -5.1),
        ]
    )

    result = runner.invoke(cli, ["rank"], input=packed_res)

    assert result.exit_code == 0
    assert "Si-pathological" in result.stdout
    stderr = result.stderr if result.stderr_bytes is not None else ""
    assert "Pathological prune" not in stderr


def test_rank_maxwell_collapses_duplicate_compositions():
    """Test -m emits one row per reduced composition."""
    runner = CliRunner()
    packed_res = (
        "TITL Si-001 0.0 10.0 -5.000000 0 0 1 (P1) n - 1\n"
        "CELL 1.0  3.0 3.0 3.0 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "END\n"
        "TITL O2-001 0.0 20.0 -4.000000 0 0 2 (P1) n - 1\n"
        "CELL 1.0  3.0 3.0 3.0 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC O\n"
        "O      1  0.0  0.0  0.0  1.0\n"
        "O      1  0.5  0.5  0.5  1.0\n"
        "END\n"
        "TITL SiO-001 0.0 20.0 -20.000000 0 0 2 (P1) n - 1\n"
        "CELL 1.0  3.0 3.0 3.0 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si O\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "O      2  0.5  0.5  0.5  1.0\n"
        "END\n"
        "TITL SiO-002 0.0 21.0 -18.000000 0 0 2 (P1) n - 3\n"
        "CELL 1.0  3.1 3.1 3.1 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si O\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "O      2  0.5  0.5  0.5  1.0\n"
        "END\n"
    )
    result = runner.invoke(cli, ["rank", "-m", "-el", "Si,O"], input=packed_res)
    assert result.exit_code == 0
    lines = [line for line in result.output.splitlines() if "SiO-001" in line]
    assert len(lines) == 1
    assert lines[0].split()[-1] == "4"
    assert "SiO-002" not in result.output


def test_rank_delta_e_filter():
    """Test -de delta-e filter."""
    runner = CliRunner()
    packed_res = (
        "TITL Si-001 -0.05 40.0 -42.500000 0 0 4 (Fd-3m) n - 1\n"
        "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "Si     1  0.25 0.25 0.25 1.0\n"
        "Si     1  0.5  0.5  0.5  1.0\n"
        "Si     1  0.75 0.75 0.75 1.0\n"
        "END\n"
        "TITL Si-002 -0.05 42.0 -30.000000 0 0 4 (Pm-3m) n - 1\n"
        "CELL 1.0  5.50 5.50 5.50 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "Si     1  0.5  0.5  0.5  1.0\n"
        "Si     1  0.25 0.25 0.25 1.0\n"
        "Si     1  0.75 0.75 0.75 1.0\n"
        "END\n"
    )
    # delta_e=0.5 should filter out Si-002 which is much higher energy
    result = runner.invoke(cli, ["rank", "-de", "0.5"], input=packed_res)
    assert result.exit_code == 0
    assert "Si-001" in result.output
    # Si-002 should be filtered out
    lines = [line for line in result.output.split("\n") if "Si-002" in line]
    assert len(lines) == 0


def test_rank_extxyz_file(tmp_path):
    """Test reading an extxyz file."""
    runner = CliRunner()
    xyz_file = tmp_path / "test.xyz"
    xyz_content = (
        "2\n"
        'Lattice="5.43 0.0 0.0 0.0 5.43 0.0 0.0 0.0 5.43" '
        "Properties=species:S:1:pos:R:3 energy=-10.625 label=Si-001\n"
        "Si 0.0 0.0 0.0\n"
        "Si 0.25 0.25 0.25\n"
    )
    xyz_file.write_text(xyz_content)
    result = runner.invoke(cli, ["rank", str(xyz_file)])
    assert result.exit_code == 0
    assert "Si-001" in result.output


def test_rank_formula_filter_no_match():
    """Test -f with a formula that matches nothing."""
    runner = CliRunner()
    packed_res = (
        "TITL Si-001 -0.05 40.0 -42.500000 0 0 4 (Fd-3m) n - 1\n"
        "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "Si     1  0.25 0.25 0.25 1.0\n"
        "Si     1  0.5  0.5  0.5  1.0\n"
        "Si     1  0.75 0.75 0.75 1.0\n"
        "END\n"
    )
    result = runner.invoke(cli, ["rank", "-f", "NaCl"], input=packed_res)
    assert result.exit_code == 0
    assert "No structures remaining after filtering" in result.output


def test_rank_unite_mode():
    """Test -u unite mode for merging similar structures."""
    runner = CliRunner()
    # Two structures with identical coordinates but different energies
    packed_res = (
        "TITL Si-001 -0.05 40.0 -42.500000 0 0 4 (Fd-3m) n - 1\n"
        "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0000000000  0.0000000000  0.0000000000 1.0\n"
        "Si     1  0.2500000000  0.2500000000  0.2500000000 1.0\n"
        "Si     1  0.5000000000  0.5000000000  0.5000000000 1.0\n"
        "Si     1  0.7500000000  0.7500000000  0.7500000000 1.0\n"
        "END\n"
        "TITL Si-002 -0.05 40.0 -43.000000 0 0 4 (Fd-3m) n - 1\n"
        "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0000000000  0.0000000000  0.0000000000 1.0\n"
        "Si     1  0.2500000000  0.2500000000  0.2500000000 1.0\n"
        "Si     1  0.5000000000  0.5000000000  0.5000000000 1.0\n"
        "Si     1  0.7500000000  0.7500000000  0.7500000000 1.0\n"
        "END\n"
    )
    result = runner.invoke(cli, ["rank", "-u", "0.1"], input=packed_res)
    assert result.exit_code == 0
    assert "Merging similar structures" in result.output
    assert "After merging" in result.output


def test_convert_help():
    """Test 'convert --help'."""
    runner = CliRunner()
    result = runner.invoke(cli, ["convert", "--help"])
    assert result.exit_code == 0
    assert "INPUT_PATH" in result.output
    assert "OUTPUT_PATH" in result.output


def test_convert_res_to_xyz(tmp_path):
    """Test converting packed .res to extxyz."""
    runner = CliRunner()
    packed = tmp_path / "packed.res"
    xyz_out = tmp_path / "out.xyz"
    packed.write_text(
        "TITL Si-001 -0.05 40.0 -42.500000 0 0 4 (Fd-3m) n - 1\n"
        "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "Si     1  0.25 0.25 0.25 1.0\n"
        "END\n"
        "TITL Si-002 -0.05 42.0 -43.200000 0 0 4 (Pm-3m) n - 1\n"
        "CELL 1.0  5.50 5.50 5.50 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "Si     1  0.5  0.5  0.5  1.0\n"
        "END\n"
    )
    result = runner.invoke(cli, ["convert", str(packed), str(xyz_out)])
    assert result.exit_code == 0
    assert "Converted 2 structures" in result.output
    assert xyz_out.exists()


def test_convert_extract_by_label(tmp_path):
    """Test extracting a single structure by label."""
    runner = CliRunner()
    packed = tmp_path / "packed.res"
    out = tmp_path / "Si-002.res"
    packed.write_text(
        "TITL Si-001 -0.05 40.0 -42.5 0 0 4 (Fd-3m) n - 1\n"
        "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "END\n"
        "TITL Si-002 -0.05 42.0 -43.2 0 0 4 (Pm-3m) n - 1\n"
        "CELL 1.0  5.50 5.50 5.50 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.5  0.5  0.5  1.0\n"
        "END\n"
    )
    result = runner.invoke(cli, ["convert", "-l", "Si-002", str(packed), str(out)])
    assert result.exit_code == 0
    assert "Extracted" in result.output
    content = out.read_text()
    assert "Si-002" in content
    assert "Si-001" not in content


def test_convert_res_to_res_directory_handles_missing_final_end(tmp_path):
    """RES directory conversion splits raw blocks without full structure parsing."""
    runner = CliRunner()
    packed = tmp_path / "packed.res"
    out_dir = tmp_path / "res_out"
    packed.write_text(
        _single_atom_res("Si-001", -1.0)
        + _single_atom_res("Si-002", -2.0).replace("END\n", "")
    )

    result = runner.invoke(cli, ["convert", str(packed), str(out_dir)])

    assert result.exit_code == 0
    assert "Unpacked 2 structures" in result.output
    assert (out_dir / "Si-001.res").exists()
    assert (out_dir / "Si-002.res").read_text().rstrip().endswith("END")


def test_convert_extract_not_found(tmp_path):
    """Test extracting a nonexistent label."""
    runner = CliRunner()
    packed = tmp_path / "packed.res"
    out = tmp_path / "out.res"
    packed.write_text(
        "TITL Si-001 -0.05 40.0 -42.5 0 0 4 (Fd-3m) n - 1\n"
        "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si     1  0.0  0.0  0.0  1.0\n"
        "END\n"
    )
    result = runner.invoke(cli, ["convert", "-l", "missing", str(packed), str(out)])
    assert result.exit_code == 1
    assert "not found" in result.output


def test_convert_xyz_to_res(tmp_path):
    """Test converting extxyz to unpacked .res files."""
    runner = CliRunner()
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointCalculator

    xyz_in = tmp_path / "in.xyz"
    out_dir = tmp_path / "res_output"
    atoms = Atoms(
        "Si2", positions=[[0, 0, 0], [1.3, 1.3, 1.3]], cell=[5.43] * 3, pbc=True
    )
    atoms.info["label"] = "Si-test"
    atoms.info["pressure"] = 0.0
    atoms.info["spin"] = 0.0
    atoms.info["spin_abs"] = 0.0
    atoms.info["symm"] = "P1"
    atoms.info["copies"] = 1
    calc = SinglePointCalculator(atoms, energy=-42.5)
    atoms.calc = calc
    from ase.io import write

    write(str(xyz_in), atoms, format="extxyz")

    result = runner.invoke(cli, ["convert", str(xyz_in), str(out_dir)])
    assert result.exit_code == 0
    assert "Converted 1 structures" in result.output
    assert (out_dir / "Si-test.res").exists()


def test_pack_and_unpack_res_roundtrip(tmp_path):
    runner = CliRunner()
    src_dir = tmp_path / "src"
    out_dir = tmp_path / "out"
    src_dir.mkdir()
    (src_dir / "a.res").write_text(_single_atom_res("Si-001", -1.0))
    (src_dir / "b.res").write_text(_single_atom_res("Si-002", -2.0).rstrip("\n"))
    packed = tmp_path / "packed.res"

    pack_result = runner.invoke(cli, ["pack", "--from-dir", str(src_dir), str(packed)])
    unpack_result = runner.invoke(cli, ["unpack", str(packed), str(out_dir)])

    assert pack_result.exit_code == 0
    assert unpack_result.exit_code == 0
    assert (out_dir / "Si-001.res").exists()
    assert (out_dir / "Si-002.res").exists()


def test_unpack_extxyz_writes_each_structure(tmp_path):
    from ase.io import write

    runner = CliRunner()
    xyz = tmp_path / "packed.xyz"
    out_dir = tmp_path / "xyz_out"
    atoms_a = Atoms("Si", positions=[[0, 0, 0]], cell=[3, 3, 3], pbc=True)
    atoms_a.info["label"] = "Si-a"
    atoms_b = Atoms("Si", positions=[[0, 0, 0]], cell=[3, 3, 3], pbc=True)
    atoms_b.info["label"] = "Si-b"
    write(str(xyz), [atoms_a, atoms_b], format="extxyz")

    result = runner.invoke(cli, ["unpack", str(xyz), str(out_dir)])

    assert result.exit_code == 0
    assert (out_dir / "Si-a.xyz").exists()
    assert (out_dir / "Si-b.xyz").exists()


def test_run_help():
    """Test 'run --help'."""
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--help"])
    assert result.exit_code == 0
    assert "search" in result.output
    assert "relax" in result.output
    assert "crud" in result.output


def test_all_run_workflows_expose_the_complete_backend_matrix():
    """Search, relax, CRUD, and SP share the same complete backend choices."""
    runner = CliRunner()
    expected = "castep|gulp|pp3|abacus|vasp|ml"

    for command in ("search", "relax", "crud", "sp"):
        result = runner.invoke(cli, ["run", command, "--help"])
        assert result.exit_code == 0
        assert expected in result.output.replace("\n", "")


def test_run_search_help():
    """Test 'run search --help'."""
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "search", "--help"])
    assert result.exit_code == 0
    assert "--seed" in result.output
    assert "--nmax" in result.output
    assert "--build-only" in result.output
    assert "--code" in result.output
    assert "--formula" in result.output
    assert "--elements" in result.output
    assert "--max-coeff" in result.output
    assert "--max-num-atoms" in result.output
    assert "--composition-ratio" in result.output
    assert "--oxidation-state" in result.output
    assert "--volume-minsep-source" in result.output
    assert "--formula-elements" not in result.output
    assert "--prune" in result.output
    assert "--cell-axis-map" in result.output
    assert "--calculator" in result.output
    assert "--device" in result.output


def test_run_search_ml_uses_default_torchsim_driver():
    """Random search can relax generated cells through the default ML driver."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("Si.cell").write_text("#SPECIES=Si\n#NATOM=1\n")
        Path("Si-001.cell").write_text(
            "%BLOCK LATTICE_CART\n"
            "3 0 0\n0 3 0\n0 0 3\n"
            "%ENDBLOCK LATTICE_CART\n"
            "%BLOCK POSITIONS_ABS\n"
            "Si 0 0 0\n"
            "%ENDBLOCK POSITIONS_ABS\n"
        )
        fake_torchsim = MagicMock()
        fake_torchsim.relax_batch.return_value = {"Si-001": 0}
        build_result = {
            "struct_name": "Si-001",
            "seed_name": "Si",
            "struct_content": Path("Si-001.cell").read_text(),
        }
        with patch("airsspy.jf.runners.run_buildcell", return_value=build_result):
            with patch("airsspy.cli.cmd_run._ensure_torchsim_available"):
                with patch(
                    "airsspy.jf.ml_runners.TorchSimRunner",
                    return_value=fake_torchsim,
                ) as torchsim_cls:
                    with patch("airsspy.cli.cmd_run._collect_result") as collect:
                        result = runner.invoke(
                            cli,
                            [
                                "run",
                                "search",
                                "--seed",
                                "Si",
                                "--nmax",
                                "1",
                                "--code",
                                "ml",
                                "--calculator",
                                "mace:medium",
                                "--device",
                                "cpu",
                            ],
                        )

    assert result.exit_code == 0
    torchsim_cls.assert_called_once_with("mace:medium", device="cpu")
    fake_torchsim.relax_batch.assert_called_once()
    assert fake_torchsim.relax_batch.call_args.args[0] == ["Si-001"]
    collect.assert_called_once_with(
        "Si-001", "ml", calculator_spec="mace:medium"
    )


def test_run_search_ml_supports_explicit_ase_driver():
    """Random search routes explicit ASE calculator specs through the ASE runner."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("Si.cell").write_text("#SPECIES=Si\n#NATOM=1\n")
        Path("Si-001.cell").write_text("cell content\n")
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        build_result = {
            "struct_name": "Si-001",
            "seed_name": "Si",
            "struct_content": "cell content\n",
        }
        with patch("airsspy.jf.runners.run_buildcell", return_value=build_result):
            with patch(
                "airsspy.cli.cmd_run._create_runner", return_value=fake_runner
            ) as create:
                with patch("airsspy.cli.cmd_run._collect_result"):
                    result = runner.invoke(
                        cli,
                        [
                            "run",
                            "search",
                            "--seed",
                            "Si",
                            "--nmax",
                            "1",
                            "--code",
                            "ml",
                            "--calculator",
                            "ase:my.module:Calculator",
                        ],
                    )

    assert result.exit_code == 0
    assert create.call_args.kwargs["calculator_spec"] == "ase:my.module:Calculator"
    fake_runner.run.assert_called_once_with("Si-001", "cell content\n")


def test_run_crud_help():
    """Test 'run crud --help'."""
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "crud", "--help"])
    assert result.exit_code == 0
    assert "--workdir" in result.output
    assert "--nostop" in result.output
    assert "--cycle" in result.output
    assert "--singlepoint" in result.output
    assert "--device" in result.output
    assert "--batch-size" in result.output
    assert "--pack" not in result.output


def test_crud_claim_is_single_owner(tmp_path):
    """A hopper job can be checked out only once."""
    hopper = tmp_path / "hopper"
    hopper.mkdir()
    (hopper / "Si-001.res").write_text("res")

    first = cmd_run._claim_crud_job(tmp_path, 1000)
    second = cmd_run._claim_crud_job(tmp_path, 1000)

    assert first == tmp_path / "Si-001.res"
    assert second is None
    assert (tmp_path / "Si-001.res").exists()
    assert not (hopper / "Si-001.res").exists()


def test_crud_res_to_cell_preserves_root_settings(tmp_path):
    """CRUD cell reconstruction replaces geometry and keeps root settings."""
    root_cell = tmp_path / "Si.cell"
    root_cell.write_text(
        "%BLOCK LATTICE_CART\n"
        "1 0 0\n0 1 0\n0 0 1\n"
        "%ENDBLOCK LATTICE_CART\n"
        "%BLOCK POSITIONS_FRAC\n"
        "Si 0 0 0\n"
        "%ENDBLOCK POSITIONS_FRAC\n"
        "kpoints_mp_grid : 2 2 2\n"
    )
    res = tmp_path / "Si-001.res"
    res.write_text(
        "TITL Si-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
        "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si 1 0.5000000000000 0.5000000000000 0.5000000000000 1.0\n"
        "END\n"
    )

    lines = cmd_run._res_to_cell_lines(res, root_cell)
    text = "\n".join(lines)

    assert "%BLOCK LATTICE_CART" in text
    assert "5.0000000000 0.0000000000 0.0000000000" in text
    assert "Si  2.5000000000 2.5000000000 2.5000000000" in text
    assert "kpoints_mp_grid : 2 2 2" in text
    assert text.count("%BLOCK LATTICE_CART") == 1
    assert "POSITIONS_FRAC" not in text


def test_crud_prepare_inputs_preserves_res_spins(tmp_path, monkeypatch):
    """Magnetic RES inputs keep SPIN tags and update CASTEP total spin."""
    monkeypatch.chdir(tmp_path)
    Path("Fe.cell").write_text(
        "%BLOCK LATTICE_CART\n"
        "1 0 0\n0 1 0\n0 0 1\n"
        "%ENDBLOCK LATTICE_CART\n"
        "%BLOCK POSITIONS_ABS\n"
        "Fe 0 0 0\n"
        "Fe 0.5 0.5 0.5\n"
        "%ENDBLOCK POSITIONS_ABS\n"
    )
    Path("Fe.param").write_text("task : geometryoptimization\nspin : 0\n")
    Path("Fe-001.res").write_text(
        "TITL Fe-001 0.000 125.000 -1.0000 2.00 2.00 2 (P1) n - 1\n"
        "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
        "LATT -1\n"
        "SFAC Fe\n"
        "Fe 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0 1.5\n"
        "Fe 1 0.5000000000000 0.5000000000000 0.5000000000000 1.0 0.5\n"
        "END\n"
    )

    cmd_run._prepare_crud_inputs("Fe-001", "castep")

    cell_text = Path("Fe-001.cell").read_text()
    param_text = Path("Fe-001.param").read_text()
    assert "SPIN=1.500" in cell_text
    assert "SPIN=0.500" in cell_text
    assert "spin :      2.000" in param_text
    assert "spin : 0" not in param_text


def test_crud_prepare_inputs_does_not_treat_forces_as_spins(tmp_path, monkeypatch):
    """Force-only RES columns must not become CASTEP SPIN tags."""
    monkeypatch.chdir(tmp_path)
    Path("Si.cell").write_text(
        "%BLOCK LATTICE_CART\n"
        "1 0 0\n0 1 0\n0 0 1\n"
        "%ENDBLOCK LATTICE_CART\n"
    )
    Path("Si.param").write_text("task : geometryoptimization\n")
    Path("Si-001.res").write_text(
        "TITL Si-001 0.000 125.000 -1.0000 0.00 0.00 2 (P1) n - 1\n"
        "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0 0.10 -0.20 0.30\n"
        "Si 1 0.5000000000000 0.5000000000000 0.5000000000000 1.0 -0.10 0.20 -0.30\n"
        "END\n"
    )

    cmd_run._prepare_crud_inputs("Si-001", "castep")

    cell_text = Path("Si-001.cell").read_text()
    param_text = Path("Si-001.param").read_text()
    assert "SPIN=" not in cell_text
    assert "spin :" not in param_text


def test_run_crud_processes_claimed_ml_job():
    """CRUD ML checks out and runs a RES job without candidate cell files."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("hopper").mkdir()
        Path("Si.cell").write_text("kpoints_mp_grid : 1 1 1\n")
        Path("hopper/Si-001.res").write_text(
            "TITL Si-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch("airsspy.cli.cmd_run._create_runner", return_value=fake_runner):
            with patch("airsspy.cli.cmd_run._collect_result") as collect:
                result = runner.invoke(
                    cli,
                    [
                        "run",
                        "crud",
                        "--code",
                        "ml",
                        "--calculator",
                        "ase:dummy:model",
                        "--keep",
                    ],
                )
                assert result.exit_code == 0
                assert collect.call_count == 1
                assert Path("good_castep/Si-001.res").exists()
                assert not Path("good_castep/Si-001.cell").exists()
                assert not Path("hopper/Si-001.res").exists()
                fake_runner.run.assert_called_once()
                assert isinstance(fake_runner.run.call_args.args[1], Atoms)

def test_run_crud_processes_claimed_vasp_job():
    """CRUD VASP converts RES input and passes INCAR/KPOINTS to the runner."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("hopper").mkdir()
        Path("Si.cell").write_text("kpoints_mp_grid : 1 1 1\n")
        Path("Si.INCAR").write_text("ENCUT = 400\n")
        Path("Si.KPOINTS").write_text("explicit kpoints\n")
        Path("hopper/Si-001.res").write_text(
            "TITL Si-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch("airsspy.cli.cmd_run._create_task_runner", return_value=fake_runner):
            with patch("airsspy.cli.cmd_run._collect_result") as collect:
                result = runner.invoke(
                    cli,
                    ["run", "crud", "--code", "vasp", "--keep"],
                )

        assert result.exit_code == 0
        fake_runner.run.assert_called_once()
        assert fake_runner.run.call_args.args[:3] == (
            "Si-001",
            Path("good_castep/Si-001.cell").read_text(),
            "ENCUT = 400\n",
        )
        kpoints_path = fake_runner.run.call_args.kwargs["kpoints_path"]
        assert kpoints_path.name == "Si-001.KPOINTS"
        assert Path("good_castep/Si-001.KPOINTS").read_text() == "explicit kpoints\n"
        collect.assert_called_once()


def test_create_runner_passes_max_iterations_to_vasp():
    """VASP relaxations use the same restart iteration budget as CASTEP."""
    with patch("airsspy.jf.runners.AirssVaspRelaxRunner") as runner_cls:
        cmd_run._create_runner(
            "vasp",
            "vasp_std",
            7,
            False,
            100.0,
            potcar_dir="/potcars",
            potcar_map={"Si": "Si"},
        )

    runner_cls.assert_called_once_with(
        executable="vasp_std",
        max_iterations=7,
        pressure=100.0,
        potcar_dir="/potcars",
        potcar_map={"Si": "Si"},
    )


def test_create_sp_runner_supports_gulp_and_pp3():
    """The SP factory constructs native no-relax runners for both engines."""
    from airsspy.jf.runners import (
        AirssGulpSinglePointRunner,
        AirssPp3SinglePointRunner,
    )

    gulp = cmd_run._create_sp_runner("gulp", "ggulp", pressure=4.0)
    pp3 = cmd_run._create_sp_runner("pp3", "pp3")

    assert isinstance(gulp, AirssGulpSinglePointRunner)
    assert gulp.pressure == 4.0
    assert isinstance(pp3, AirssPp3SinglePointRunner)


def test_run_crud_ml_torchsim_batches_claimed_jobs():
    """CRUD ML supports the same torch-sim model path as run relax."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("hopper").mkdir()
        res = (
            "TITL {label} 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        for label in ("Si-001", "Si-002"):
            Path(f"hopper/{label}.res").write_text(res.format(label=label))
        fake_torchsim = MagicMock()
        captured = {}

        def fake_relax_batch(names, structures, **kwargs):
            captured["names"] = names
            captured["structures"] = structures
            captured["kwargs"] = kwargs
            for name in names:
                Path(f"{name}.extxyz").write_text("")
            return dict.fromkeys(names, 0)

        fake_torchsim.relax_batch.side_effect = fake_relax_batch
        with patch("airsspy.cli.cmd_run._is_torchsim_model", return_value=True):
            with patch("airsspy.cli.cmd_run._ensure_torchsim_available"):
                with patch(
                    "airsspy.jf.ml_runners.TorchSimRunner",
                    return_value=fake_torchsim,
                ) as torchsim_cls:
                    with patch("airsspy.cli.cmd_run._collect_result") as collect:
                        result = runner.invoke(
                            cli,
                            [
                                "run",
                                "crud",
                                "--code",
                                "ml",
                                "--calculator",
                                "mace:medium",
                                "--device",
                                "cuda",
                                "--batch-size",
                                "2",
                            ],
                        )

        assert result.exit_code == 0
        torchsim_cls.assert_called_once_with("mace:medium", device="cuda")
        assert set(captured["names"]) == {"Si-001", "Si-002"}
        assert all(isinstance(struct, Atoms) for struct in captured["structures"])
        assert captured["kwargs"]["optimizer"] == "fire"
        assert collect.call_count == 2
        assert Path("good_castep/Si-001.res").exists()
        assert Path("good_castep/Si-002.res").exists()


def test_run_crud_ml_torchsim_singlepoint_uses_static_batch():
    """CRUD ML single-point uses torch-sim static batching."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("hopper").mkdir()
        res = (
            "TITL {label} 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        for label in ("Si-001", "Si-002"):
            Path(f"hopper/{label}.res").write_text(res.format(label=label))
        fake_torchsim = MagicMock()
        fake_torchsim.static_batch.return_value = {"Si-001": 0, "Si-002": 0}
        with patch("airsspy.cli.cmd_run._is_torchsim_model", return_value=True):
            with patch("airsspy.cli.cmd_run._ensure_torchsim_available"):
                with patch(
                    "airsspy.jf.ml_runners.TorchSimRunner",
                    return_value=fake_torchsim,
                ):
                    with patch("airsspy.cli.cmd_run._collect_result") as collect:
                        result = runner.invoke(
                            cli,
                            [
                                "run",
                                "crud",
                                "--code",
                                "ml",
                                "--calculator",
                                "mace:medium",
                                "--singlepoint",
                                "--batch-size",
                                "2",
                            ],
                        )

        assert result.exit_code == 0
        fake_torchsim.static_batch.assert_called_once()
        assert fake_torchsim.relax_batch.call_count == 0
        assert collect.call_count == 2


def test_run_crud_singlepoint_uses_sp_runner():
    """CRUD can run a single-point calculation instead of full relaxation."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("hopper").mkdir()
        Path("Si.cell").write_text("kpoints_mp_grid : 1 1 1\n")
        Path("Si.param").write_text("task : geometryoptimization\n")
        Path("hopper/Si-001.res").write_text(
            "TITL Si-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch("airsspy.cli.cmd_run._create_sp_runner", return_value=fake_runner):
            with patch("airsspy.cli.cmd_run._collect_result") as collect:
                result = runner.invoke(
                    cli,
                    [
                        "run",
                        "crud",
                        "--code",
                        "castep",
                        "--singlepoint",
                        "--keep",
                    ],
                )

        assert result.exit_code == 0
        fake_runner.run.assert_called_once()
        paraminput = fake_runner.run.call_args.args[2]
        assert paraminput["task"] == "geometryoptimization"
        collect.assert_called_once()


def test_ml_model_backend_resolution():
    """ML model specs default to torch-sim and use ase: for the ASE driver."""
    assert cmd_run._is_torchsim_model("mace:medium")
    assert cmd_run._is_torchsim_model("torch-sim:mace:medium")
    assert not cmd_run._is_torchsim_model("ase:mace:medium")
    assert cmd_run._is_symmetrix_model("symmetrix:mace:medium")
    assert cmd_run._is_symmetrix_model("ase:symmetrics:medium")
    assert cmd_run._is_symmetrix_model("ase:symmetrix:medium")
    assert not cmd_run._is_torchsim_model("symmetrix:mace:medium")
    assert not cmd_run._is_torchsim_model("ase:symmetrics:medium")
    assert (
        cmd_run._normalize_ml_ase_spec("ase:mace:medium")
        == "mace.calculators:MACECalculator@medium"
    )
    assert (
        cmd_run._normalize_ml_ase_spec("symmetrix:mace:medium")
        == "symmetrix:Symmetrix@medium"
    )
    assert (
        cmd_run._normalize_ml_ase_spec("ase:symmetrics:MACE-MH-1:matpes_r2scan")
        == "symmetrix:Symmetrix@MACE-MH-1:matpes_r2scan"
    )
    assert (
        cmd_run._normalize_ml_ase_spec("ase:symmetrix:medium-mpa-0")
        == "symmetrix:Symmetrix@medium-mpa-0"
    )
    assert (
        cmd_run._normalize_ml_ase_spec("ase:symmetrix:Symmetrix@mh-1")
        == "symmetrix:Symmetrix@mh-1"
    )
    assert (
        cmd_run._normalize_ml_ase_spec("ase:my.module:Calc@model")
        == "my.module:Calc@model"
    )


def test_symmetrix_rejects_non_mace_models():
    """Symmetrix backend is intentionally limited to MACE models."""
    try:
        cmd_run._normalize_ml_ase_spec("symmetrix:sevennet:sevennet-mf-ompa")
    except ValueError as exc:
        assert "ase:symmetrix:<mace-model>" in str(exc)
    else:
        raise AssertionError("Expected Symmetrix non-MACE model to fail")


def test_plain_ml_model_fails_clearly_without_torchsim():
    """Plain ML model specs do not fall back to ASE if torch-sim is missing."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LiTaOCl.cell").write_text("kpoints_mp_grid : 1 1 1\n")
        Path("LiTaOCl-001.res").write_text(
            "TITL LiTaOCl-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        with patch("airsspy.jf.ml_runners.has_torchsim", return_value=False):
            result = runner.invoke(
                cli,
                [
                    "run",
                    "relax",
                    "--cell",
                    "*.res",
                    "--code",
                    "ml",
                    "--calculator",
                    "mace:medium",
                ],
            )

        assert result.exit_code != 0
        assert "torch-sim is required" in result.output


def test_run_relax_ase_symmetrics_uses_ase_runner_without_torchsim():
    """Symmetrix MACE specs use ASE runner routing, not torch-sim batching."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LiTaOCl.cell").write_text("kpoints_mp_grid : 1 1 1\n")
        Path("LiTaOCl-001.res").write_text(
            "TITL LiTaOCl-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch("airsspy.jf.ml_runners.has_torchsim", return_value=False):
            with patch("airsspy.cli.cmd_run._create_runner", return_value=fake_runner):
                with patch("airsspy.cli.cmd_run._collect_result") as collect:
                    result = runner.invoke(
                        cli,
                        [
                            "run",
                            "relax",
                            "--cell",
                            "*.res",
                            "--code",
                            "ml",
                            "--calculator",
                            "ase:symmetrics:medium-mpa-0",
                            "--keep",
                        ],
                    )

        assert result.exit_code == 0
        fake_runner.run.assert_called_once()
        collect.assert_called_once()


def test_run_sp_symmetrix_uses_ase_runner_without_torchsim():
    """Symmetrix MACE single-points use the normal ML SP runner path."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LiTaOCl-001.res").write_text(
            "TITL LiTaOCl-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch("airsspy.jf.ml_runners.has_torchsim", return_value=False):
            with patch("airsspy.cli.cmd_run._create_sp_runner", return_value=fake_runner):
                with patch("airsspy.cli.cmd_run._collect_result") as collect:
                    result = runner.invoke(
                        cli,
                        [
                            "run",
                            "sp",
                            "--cell",
                            "*.res",
                            "--code",
                            "ml",
                            "--calculator",
                            "symmetrix:mace:medium",
                        ],
                    )

    assert result.exit_code == 0
    fake_runner.run.assert_called_once()
    collect.assert_called_once()


def test_run_crud_symmetrix_uses_ase_runner_without_torchsim():
    """Symmetrix MACE CRUD work is handled by the non-batch ML runner."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("hopper").mkdir()
        Path("hopper/LiTaOCl-001.res").write_text(
            "TITL LiTaOCl-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch("airsspy.jf.ml_runners.has_torchsim", return_value=False):
            with patch("airsspy.cli.cmd_run._create_task_runner", return_value=fake_runner):
                with patch("airsspy.cli.cmd_run._collect_result") as collect:
                    result = runner.invoke(
                        cli,
                        [
                            "run",
                            "crud",
                            "--code",
                            "ml",
                            "--calculator",
                            "symmetrix:mace:medium",
                        ],
                    )

    assert result.exit_code == 0
    fake_runner.run.assert_called_once()
    collect.assert_called_once()


def test_run_relax_accepts_res_input_for_ase_ml():
    """Explicit ase: model uses the ASE driver with in-memory RES parsing."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LiTaOCl.cell").write_text("kpoints_mp_grid : 1 1 1\n")
        Path("LiTaOCl-001.res").write_text(
            "TITL LiTaOCl-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch("airsspy.cli.cmd_run._is_torchsim_model", return_value=False):
            with patch("airsspy.cli.cmd_run._create_runner", return_value=fake_runner):
                with patch("airsspy.cli.cmd_run._collect_result") as collect:
                    result = runner.invoke(
                        cli,
                        [
                            "run",
                            "relax",
                            "--cell",
                            "*.res",
                            "--code",
                            "ml",
                            "--calculator",
                            "ase:mace:medium",
                            "--keep",
                        ],
                    )

        assert result.exit_code == 0
        assert not Path("LiTaOCl-001.cell").exists()
        fake_runner.run.assert_called_once()
        assert fake_runner.run.call_args.args[0] == "LiTaOCl-001"
        assert isinstance(fake_runner.run.call_args.args[1], Atoms)
        assert fake_runner.run.call_args.args[1].get_chemical_formula() == "Si"
        collect.assert_called_once()


def test_run_relax_torchsim_res_input_passes_device_without_cell_side_effect():
    """TorchSim RES relaxation passes device and avoids writing converted cells."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LiTaOCl.cell").write_text("kpoints_mp_grid : 1 1 1\n")
        Path("LiTaOCl-001.res").write_text(
            "TITL LiTaOCl-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        fake_runner = MagicMock()
        captured = {}

        def fake_batch(*args, **kwargs):
            captured["args"] = args
            captured["kwargs"] = kwargs
            return {"LiTaOCl-001": 0}

        fake_torchsim = MagicMock()
        fake_torchsim.relax_batch.side_effect = fake_batch
        with patch("airsspy.cli.cmd_run._ensure_torchsim_available"):
            with patch("airsspy.cli.cmd_run._create_runner", return_value=fake_runner):
                with patch(
                    "airsspy.jf.ml_runners.TorchSimRunner",
                    return_value=fake_torchsim,
                ) as torchsim_cls:
                    with patch("airsspy.cli.cmd_run._collect_result") as collect:
                        result = runner.invoke(
                            cli,
                            [
                                "run",
                                "relax",
                                "--cell",
                                "*.res",
                                "--code",
                                "ml",
                                "--calculator",
                                "torch-sim:mace:medium-mpa-0",
                                "--device",
                                "cuda",
                                "--keep",
                            ],
                        )

        assert result.exit_code == 0
        torchsim_cls.assert_called_once_with(
            "torch-sim:mace:medium-mpa-0", device="cuda"
        )
        assert not Path("LiTaOCl-001.cell").exists()
        assert fake_torchsim.relax_batch.call_count == 1
        assert captured["args"][0] == ["LiTaOCl-001"]
        assert isinstance(captured["args"][1][0], Atoms)
        assert captured["args"][1][0].get_chemical_formula() == "Si"
        collect.assert_called_once()


def test_run_relax_torchsim_chunks_and_skips_packed_res():
    """TorchSim globbed RES relaxation skips packed files and batches candidates."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LiTaOCl.cell").write_text("kpoints_mp_grid : 1 1 1\n")
        res = (
            "TITL {label} 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        for label in ("LiTaOCl-001", "LiTaOCl-002", "LiTaOCl-003"):
            Path(f"{label}.res").write_text(res.format(label=label))
        Path("packed.res").write_text(
            res.format(label="packed-001") + res.format(label="packed-002")
        )
        fake_runner = MagicMock()
        calls = []

        def fake_batch(names, structures, **kwargs):
            calls.append((names, structures, kwargs))
            return dict.fromkeys(names, 0)

        fake_torchsim = MagicMock()
        fake_torchsim.relax_batch.side_effect = fake_batch
        with patch("airsspy.cli.cmd_run._ensure_torchsim_available"):
            with patch("airsspy.cli.cmd_run._create_runner", return_value=fake_runner):
                with patch(
                    "airsspy.jf.ml_runners.TorchSimRunner",
                    return_value=fake_torchsim,
                ):
                    with patch("airsspy.cli.cmd_run._collect_result") as collect:
                        result = runner.invoke(
                            cli,
                            [
                                "run",
                                "relax",
                                "--cell",
                                "*.res",
                                "--code",
                                "ml",
                                "--calculator",
                                "mace:medium-mpa-0",
                                "--device",
                                "cuda",
                                "--batch-size",
                                "2",
                            ],
                        )

        assert result.exit_code == 0
        assert [call[0] for call in calls] == [
            ["LiTaOCl-001", "LiTaOCl-002"],
            ["LiTaOCl-003"],
        ]
        assert all(isinstance(struct, Atoms) for call in calls for struct in call[1])
        assert fake_torchsim.relax_batch.call_count == 2
        assert collect.call_count == 3
        assert not Path("packed.cell").exists()


def test_run_relax_accepts_res_input_for_castep():
    """Relax command converts RES inputs and loads param files for CASTEP."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LiTaOCl.cell").write_text("kpoints_mp_grid : 1 1 1\n")
        Path("LiTaOCl.param").write_text("task : geometryoptimization\n")
        Path("LiTaOCl-001.res").write_text(
            "TITL LiTaOCl-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.2500000000000 0.2500000000000 0.2500000000000 1.0\n"
            "END\n"
        )
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch("airsspy.cli.cmd_run._create_runner", return_value=fake_runner):
            with patch("airsspy.cli.cmd_run._collect_result") as collect:
                result = runner.invoke(
                    cli,
                    [
                        "run",
                        "relax",
                        "--cell",
                        "*.res",
                        "--seed",
                        "IGNORED",
                        "--code",
                        "castep",
                        "--keep",
                    ],
                )

        assert result.exit_code == 0
        assert Path("LiTaOCl-001.cell").exists()
        cell_text = Path("LiTaOCl-001.cell").read_text()
        assert "Si  1.2500000000 1.2500000000 1.2500000000" in cell_text
        fake_runner.run.assert_called_once()
        assert fake_runner.run.call_args.args[0] == "LiTaOCl-001"
        assert "%BLOCK LATTICE_CART" in fake_runner.run.call_args.args[1]
        collect.assert_called_once()


def test_run_relax_singlepoint_uses_sp_runner():
    """Relax command can dispatch existing structures through SP runners."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("Si-001.cell").write_text(
            "%BLOCK LATTICE_CART\n"
            "1 0 0\n0 1 0\n0 0 1\n"
            "%ENDBLOCK LATTICE_CART\n"
            "%BLOCK POSITIONS_ABS\n"
            "Si 0 0 0\n"
            "%ENDBLOCK POSITIONS_ABS\n"
        )
        Path("Si-001.param").write_text("task : geometryoptimization\n")
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch("airsspy.cli.cmd_run._create_sp_runner", return_value=fake_runner):
            with patch("airsspy.cli.cmd_run._collect_result") as collect:
                result = runner.invoke(
                    cli,
                    [
                        "run",
                        "relax",
                        "--cell",
                        "*.cell",
                        "--code",
                        "castep",
                        "--singlepoint",
                    ],
                )

        assert result.exit_code == 0
        fake_runner.run.assert_called_once()
        assert fake_runner.run.call_args.args[0] == "Si-001"
        collect.assert_called_once()


def test_run_relax_accepts_res_input_for_vasp():
    """Relax command converts RES inputs and dispatches VASP with root INCAR."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LiTaOCl.cell").write_text("kpoints_mp_grid : 1 1 1\n")
        Path("LiTaOCl.INCAR").write_text("ENCUT = 400\n")
        Path("LiTaOCl-001.res").write_text(
            "TITL LiTaOCl-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.2500000000000 0.2500000000000 0.2500000000000 1.0\n"
            "END\n"
        )
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch("airsspy.cli.cmd_run._create_task_runner", return_value=fake_runner):
            with patch("airsspy.cli.cmd_run._collect_result") as collect:
                result = runner.invoke(
                    cli,
                    [
                        "run",
                        "relax",
                        "--cell",
                        "*.res",
                        "--seed",
                        "IGNORED",
                        "--code",
                        "vasp",
                        "--potcar-map",
                        "Si=Si_GW",
                        "--keep",
                    ],
                )

        assert result.exit_code == 0
        fake_runner.run.assert_called_once()
        assert fake_runner.run.call_args.args[0] == "LiTaOCl-001"
        assert fake_runner.run.call_args.args[2] == "ENCUT = 400\n"
        collect.assert_called_once()


def test_run_sp_accepts_res_input_for_vasp():
    """Single-point VASP supports RES input using the root INCAR."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("Si.cell").write_text("kpoints_mp_grid : 1 1 1\n")
        Path("Si.INCAR").write_text("ENCUT = 400\n")
        Path("Si-001.res").write_text(
            "TITL Si-001 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.2500000000000 0.2500000000000 0.2500000000000 1.0\n"
            "END\n"
        )
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch("airsspy.cli.cmd_run._create_sp_runner", return_value=fake_runner):
            with patch("airsspy.cli.cmd_run._collect_result") as collect:
                result = runner.invoke(
                    cli,
                    [
                        "run",
                        "sp",
                        "--cell",
                        "*.res",
                        "--seed",
                        "IGNORED",
                        "--code",
                        "vasp",
                    ],
                )

        assert result.exit_code == 0
        fake_runner.run.assert_called_once()
        assert fake_runner.run.call_args.args[0] == "Si-001"
        assert fake_runner.run.call_args.args[2] == "ENCUT = 400\n"
        collect.assert_called_once()


def test_run_sp_accepts_res_input_for_all_other_external_backends():
    """Standalone SP converts RES input for CASTEP, GULP, PP3, and ABACUS."""
    runner = CliRunner()
    input_files = {
        "castep": (".param", "task : singlepoint\n"),
        "gulp": (".lib", "species\n"),
        "pp3": (".pp", "1.0 1.0\n"),
        "abacus": (".INPUT", "calculation scf\n"),
    }
    for code, (suffix, param_content) in input_files.items():
        with runner.isolated_filesystem():
            Path("Si.cell").write_text("kpoints_mp_grid : 1 1 1\n")
            Path("Si" + suffix).write_text(param_content)
            Path("Si-001.res").write_text(
                "TITL Si-001 0.000 27.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
                "CELL 1.0 3.000000 3.000000 3.000000 90.000000 90.000000 90.000000\n"
                "LATT -1\n"
                "SFAC Si\n"
                "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
                "END\n"
            )
            fake_runner = MagicMock()
            fake_runner.run.return_value = 0
            with patch(
                "airsspy.cli.cmd_run._create_sp_runner", return_value=fake_runner
            ):
                with patch("airsspy.cli.cmd_run._collect_result") as collect:
                    result = runner.invoke(
                        cli,
                        [
                            "run",
                            "sp",
                            "--cell",
                            "*.res",
                            "--seed",
                            "IGNORED",
                            "--code",
                            code,
                        ],
                    )

            assert result.exit_code == 0, (code, result.output, result.exception)
            fake_runner.run.assert_called_once()
            assert fake_runner.run.call_args.args[0] == "Si-001"
            assert "%BLOCK LATTICE_CART" in fake_runner.run.call_args.args[1]
            collect.assert_called_once()


def test_run_sp_res_input_uses_explicit_seed_as_lookup_fallback():
    """An explicit seed supplies templates and parameters for arbitrary RES labels."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("Si.cell").write_text("kpoints_mp_grid : 1 1 1\n")
        Path("Si.param").write_text("task : singlepoint\n")
        Path("candidate.res").write_text(
            "TITL candidate 0.000 27.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
            "CELL 1.0 3.000000 3.000000 3.000000 90.000000 90.000000 90.000000\n"
            "LATT -1\n"
            "SFAC Si\n"
            "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
            "END\n"
        )
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch(
            "airsspy.cli.cmd_run._create_sp_runner", return_value=fake_runner
        ):
            with patch("airsspy.cli.cmd_run._collect_result"):
                result = runner.invoke(
                    cli,
                    [
                        "run",
                        "sp",
                        "--cell",
                        "candidate.res",
                        "--seed",
                        "Si",
                        "--code",
                        "castep",
                    ],
                )

    assert result.exit_code == 0
    fake_runner.run.assert_called_once()
    assert "%BLOCK LATTICE_CART" in fake_runner.run.call_args.args[1]


def test_collect_result_rejects_invalid_castep2res_output_without_clobbering():
    """Failed external conversion preserves an existing RES input and raises."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("Si-001.res").write_text("original input\n")
        Path("Si-001.castep").write_text("current calculation output\n")
        failed = MagicMock(returncode=0)
        with patch("subprocess.run", return_value=failed):
            try:
                cmd_run._collect_result("Si-001", "pp3")
            except RuntimeError as exc:
                assert "valid RES" in str(exc)
            else:
                raise AssertionError("Expected invalid castep2res output to fail")

        assert Path("Si-001.res").read_text() == "original input\n"
        assert not Path("Si-001.res.tmp").exists()


def test_collect_result_rejects_castep2res_cell_only_fallback():
    """A failed external calculation cannot become a zero-energy RES via .cell."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("Si-001.cell").write_text("fresh calculation input\n")
        Path("Si-001.res").write_text("original input\n")
        with patch("subprocess.run") as convert:
            try:
                cmd_run._collect_result("Si-001", "gulp")
            except RuntimeError as exc:
                assert "No CASTEP-like result" in str(exc)
            else:
                raise AssertionError("Expected missing calculation output to fail")

        convert.assert_not_called()
        assert Path("Si-001.res").read_text() == "original input\n"


def test_failed_relax_and_sp_cleanup_restore_original_res_input():
    """Default cleanup must not delete source RES files accepted by the CLI."""
    runner = CliRunner()
    original_res = (
        "TITL candidate 0.000 27.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
        "CELL 1.0 3.000000 3.000000 3.000000 90.000000 90.000000 90.000000\n"
        "LATT -1\n"
        "SFAC Si\n"
        "Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
        "END\n"
    )

    for command, factory in (
        ("relax", "_create_runner"),
        ("sp", "_create_sp_runner"),
    ):
        with runner.isolated_filesystem():
            Path("Si.cell").write_text("kpoints_mp_grid : 1 1 1\n")
            Path("Si.lib").write_text("species\n")
            Path("candidate.res").write_text(original_res)
            fake_runner = MagicMock()
            fake_runner.run.return_value = 1
            fake_runner.clean_failed.side_effect = (
                lambda name: Path(name + ".res").unlink(missing_ok=True)
            )

            with patch(
                f"airsspy.cli.cmd_run.{factory}", return_value=fake_runner
            ):
                with patch(
                    "airsspy.cli.cmd_run._collect_result",
                    side_effect=RuntimeError("no current result"),
                ):
                    result = runner.invoke(
                        cli,
                        [
                            "run",
                            command,
                            "--cell",
                            "candidate.res",
                            "--seed",
                            "Si",
                            "--code",
                            "gulp",
                        ],
                    )

            assert result.exit_code == 0, (command, result.output, result.exception)
            assert Path("candidate.res").read_text() == original_res
            fake_runner.clean_failed.assert_called_once_with("candidate")


def test_run_search_formula_diagnose():
    """Test formula diagnosis prints rewritten seed text and exits."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("Si.cell").write_text("#SPECIES=Si\n#NATOM=2\n#SLACK=0.25\n")
        result = runner.invoke(
            cli,
            [
                "run",
                "search",
                "--seed",
                "Si",
                "--nmax",
                "1",
                "--build-only",
                "--formula",
                "Si2",
                "--diagnose",
                "1",
            ],
        )

    assert result.exit_code == 0
    assert "Sample 1" in result.output
    assert "User settings" in result.output
    assert "buildcell input" in result.output
    assert "seedfile = Si.cell" in result.output
    assert "#FORMULA=Si" in result.output
    assert "#SPECIES=Si" not in result.output
    assert "#SLACK=0.25" in result.output


def test_run_search_volume_minsep_dataset_diagnose():
    """Test dataset-backed volume/minsep diagnosis prints generated directives."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("SiO.cell").write_text(
            "#SPECIES=Si,O\n#FORMULA=Si\n#VARVOL=999\n#MINSEP=9\n#NFORM=1\n#SLACK=0.25\n"
        )
        dataset = Path("curated.json")
        dataset.write_text(
            json.dumps(
                [
                    {
                        "material_id": "mp-sio2",
                        "reduced_formula": "SiO2",
                        "chemical_formula": "SiO2",
                        "composition": {"Si": 1.0, "O": 2.0},
                        "volume": 45.0,
                        "minsep": {
                            "Si-O": 1.6,
                            "O-Si": 1.6,
                            "O-O": 2.5,
                            "Si-Si": 3.0,
                        },
                        "energy_above_hull": 0.0,
                    }
                ]
            ),
            encoding="utf-8",
        )

        result = runner.invoke(
            cli,
            [
                "run",
                "search",
                "--seed",
                "SiO",
                "--build-only",
                "--formula",
                "O2Si",
                "--volume-minsep-source",
                "dataset",
                "--volume-minsep-dataset",
                str(dataset),
                "--volume-scale",
                "1.1",
                "--max-atoms",
                "12",
                "--max-nform",
                "3",
                "--diagnose",
                "1",
            ],
        )

    assert result.exit_code == 0
    assert "volume_minsep_source = dataset" in result.output
    assert "volume_per_atom = 15" in result.output
    assert "#FORMULA=SiO2" in result.output
    assert "#VARVOL=49.5" in result.output
    assert "#MINSEP=0.5-1 O-O=2.25-2.75 O-Si=1.44-1.76 Si-Si=2.7-3.3" in result.output
    assert "#NFORM={2,3}" in result.output
    assert "#SPECIES=" not in result.output
    assert "#VARVOL=999" not in result.output
    assert "#MINSEP=9" not in result.output
    assert "#NFORM=1" not in result.output
    assert "#SLACK=0.25" in result.output


def test_run_search_volume_minsep_diagnose_ignores_superseded_seed_constraints():
    """Estimator mode ignores old NATOM/NFORM constraints it replaces."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("SiO.cell").write_text("#SPECIES=Si,O\n#NATOM=1-2\n#NFORM=1\n")
        dataset = Path("curated.json")
        dataset.write_text(
            json.dumps(
                [
                    {
                        "material_id": "mp-sio2",
                        "reduced_formula": "SiO2",
                        "chemical_formula": "SiO2",
                        "composition": {"Si": 1.0, "O": 2.0},
                        "volume": 45.0,
                        "minsep": {
                            "Si-O": 1.6,
                            "O-Si": 1.6,
                            "O-O": 2.5,
                            "Si-Si": 3.0,
                        },
                    }
                ]
            ),
            encoding="utf-8",
        )

        result = runner.invoke(
            cli,
            [
                "run",
                "search",
                "--seed",
                "SiO",
                "--build-only",
                "--formula",
                "SiO2",
                "--volume-minsep-source",
                "dataset",
                "--volume-minsep-dataset",
                str(dataset),
                "--max-atoms",
                "12",
                "--diagnose",
                "1",
            ],
        )

    assert result.exit_code == 0
    assert "#FORMULA=SiO2" in result.output
    assert "#NATOM=1-2" not in result.output
    assert "#NFORM={2,3,4}" in result.output


def test_run_search_reference_source_filters_references_by_sampled_formula():
    """Reference source can accept references for multiple sampled formulas."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("mix.cell").write_text("#SPECIES=Si,O,Na,Cl\n")
        write(
            "sio2.xyz",
            Atoms(
                symbols=["Si", "O", "O"],
                positions=[[0, 0, 0], [1.6, 0, 0], [0, 2.5, 0]],
                cell=[5, 5, 5],
                pbc=True,
            ),
        )
        write(
            "nacl.xyz",
            Atoms(
                symbols=["Na", "Cl"],
                positions=[[0, 0, 0], [2.8, 0, 0]],
                cell=[5, 5, 5],
                pbc=True,
            ),
        )

        result = runner.invoke(
            cli,
            [
                "run",
                "search",
                "--seed",
                "mix",
                "--build-only",
                "--formula",
                "SiO2",
                "--volume-minsep-source",
                "reference",
                "--reference-structure",
                "sio2.xyz",
                "--reference-structure",
                "nacl.xyz",
                "--diagnose",
                "1",
            ],
        )

    assert result.exit_code == 0
    assert "volume_minsep_source = reference" in result.output
    assert "#FORMULA=SiO2" in result.output
    assert "#MINSEP=" in result.output


def test_tools_volume_minsep_curate_and_train_baseline(tmp_path):
    """Test volume/minsep offline CLI curation and baseline training."""
    runner = CliRunner()
    raw_path = tmp_path / "raw.json"
    mp_docs_path = tmp_path / "mp_docs.json"
    curated_path = tmp_path / "curated.json"
    output_dir = tmp_path / "artifacts"
    raw_path.write_text(
        json.dumps(
            [
                {
                    "material_id": "toy-1",
                    "reduced_formula": "SiO2",
                    "chemical_formula": "SiO2",
                    "composition": {"Si": 1.0, "O": 2.0},
                    "volume": 45.0,
                    "minsep": {
                        "Si-O": 1.60,
                        "O-Si": 1.60,
                        "O-O": 2.55,
                        "Si-Si": 3.05,
                    },
                },
                {
                    "material_id": "toy-2",
                    "reduced_formula": "SiO2",
                    "chemical_formula": "SiO2",
                    "composition": {"Si": 2.0, "O": 4.0},
                    "volume": 91.2,
                    "minsep": {
                        "Si-O": 1.58,
                        "O-Si": 1.58,
                        "O-O": 2.50,
                        "Si-Si": 3.10,
                    },
                },
                {
                    "material_id": "toy-3",
                    "reduced_formula": "MgO",
                    "chemical_formula": "MgO",
                    "composition": {"Mg": 1.0, "O": 1.0},
                    "volume": 22.4,
                    "minsep": {
                        "Mg-O": 1.95,
                        "O-Mg": 1.95,
                        "Mg-Mg": 3.02,
                        "O-O": 3.02,
                    },
                },
                {
                    "material_id": "toy-4",
                    "reduced_formula": "NaCl",
                    "chemical_formula": "NaCl",
                    "composition": {"Na": 1.0, "Cl": 1.0},
                    "volume": 34.0,
                    "minsep": {
                        "Na-Cl": 2.81,
                        "Cl-Na": 2.81,
                        "Na-Na": 3.95,
                        "Cl-Cl": 3.95,
                    },
                },
            ]
        ),
        encoding="utf-8",
    )
    mp_docs_path.write_text(
        json.dumps(
            [
                {"material_id": "toy-1", "energy_above_hull": 0.03},
                {"material_id": "toy-2", "energy_above_hull": 0.0},
                {"material_id": "toy-3", "energy_above_hull": 0.0},
                {"material_id": "toy-4", "energy_above_hull": 0.0},
            ]
        ),
        encoding="utf-8",
    )

    curate = runner.invoke(
        cli,
        [
            "tools",
            "volume-minsep-curate-dataset",
            "--dataset",
            str(raw_path),
            "--mp-docs",
            str(mp_docs_path),
            "--output",
            str(curated_path),
        ],
    )

    assert curate.exit_code == 0
    assert curated_path.exists()
    curated_rows = json.loads(curated_path.read_text(encoding="utf-8"))
    sio2_row = next(row for row in curated_rows if row["reduced_formula"] == "SiO2")
    assert sio2_row["material_id"] == "toy-2"

    train = runner.invoke(
        cli,
        [
            "tools",
            "volume-minsep-train-baseline",
            "--dataset",
            str(raw_path),
            "--output-dir",
            str(output_dir),
            "--volume-alpha-grid",
            "1e-4,1e-2",
            "--minsep-alpha-grid",
            "1e-4,1e-2",
        ],
    )

    assert train.exit_code == 0
    assert (output_dir / "baseline" / "baseline_bundle.json").exists()


def test_run_search_formula_accepts_comma_separated_values():
    """Test --formula accepts one comma-separated formula list."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("SiO.cell").write_text("#SPECIES=Si,O\n#NFORM=1\n")
        result = runner.invoke(
            cli,
            [
                "run",
                "search",
                "--seed",
                "SiO",
                "--build-only",
                "--formula",
                "SiO2,SiO",
                "--diagnose",
                "1",
            ],
        )

    assert result.exit_code == 0
    assert "#FORMULA=" in result.output


def test_run_search_formula_repeated_use_rejected():
    """Test --formula is not repeatable."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("SiO.cell").write_text("#SPECIES=Si,O\n#NFORM=1\n")
        result = runner.invoke(
            cli,
            [
                "run",
                "search",
                "--seed",
                "SiO",
                "--build-only",
                "--formula",
                "SiO2",
                "--formula",
                "SiO",
                "--diagnose",
                "1",
            ],
        )

    assert result.exit_code != 0
    assert "may be specified only once" in result.output


def test_run_search_formula_diagnose_combined_oxidation_states():
    """Test oxidation states can be supplied as one comma-separated option."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LTO.cell").write_text("#SPECIES=Li,Ti,O\n#NATOM=4-24\n#NFORM=1-2\n")
        result = runner.invoke(
            cli,
            [
                "run",
                "search",
                "--seed",
                "LTO",
                "--nmax",
                "1",
                "--build-only",
                "--elements",
                "Li,Ti,O",
                "--max-coeff",
                "3",
                "--oxidation-state",
                "Li=1,Ti=4,O=-2",
                "--diagnose",
                "1",
            ],
        )

    assert result.exit_code == 0
    assert "#FORMULA=" in result.output


def test_run_search_formula_diagnose_atom_budget_and_composition_ratio():
    """Test atom-budget formula sampling is exposed through run search."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LiNa.cell").write_text("#SPECIES=Li,Na\n")
        result = runner.invoke(
            cli,
            [
                "run",
                "search",
                "--seed",
                "LiNa",
                "--nmax",
                "1",
                "--build-only",
                "--elements",
                "Li,Na",
                "--max-num-atoms",
                "6",
                "--composition-ratio",
                "1=0,2=1",
                "--diagnose",
                "1",
            ],
        )

    assert result.exit_code == 0
    assert "formula_arity = 2" in result.output
    assert "formula_counts_by_arity = 1:2, 2:11" in result.output
    assert "#FORMULA=Li\n" not in result.output
    assert "#FORMULA=Na\n" not in result.output


def test_run_search_formula_composition_ratio_rejects_unavailable_arity():
    """Test composition ratio fails clearly when it cannot sample anything."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LiNa.cell").write_text("#SPECIES=Li,Na\n")
        result = runner.invoke(
            cli,
            [
                "run",
                "search",
                "--seed",
                "LiNa",
                "--nmax",
                "1",
                "--build-only",
                "--elements",
                "Li,Na",
                "--max-num-atoms",
                "6",
                "--composition-ratio",
                "3=1",
                "--diagnose",
                "1",
            ],
        )

    assert result.exit_code != 0
    assert "composition ratio" in result.output


def test_run_search_formula_oxidation_state_invalid_assignment():
    """Test invalid combined oxidation-state assignments fail clearly."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LTO.cell").write_text("#SPECIES=Li,Ti,O\n")
        result = runner.invoke(
            cli,
            [
                "run",
                "search",
                "--seed",
                "LTO",
                "--build-only",
                "--elements",
                "Li,Ti,O",
                "--oxidation-state",
                "Li=1,Ti",
                "--diagnose",
                "1",
            ],
        )

    assert result.exit_code != 0
    assert "Invalid --oxidation-state" in result.output


def test_run_search_old_formula_prefixed_alias_rejected():
    """Test old formula-prefixed aliases are not accepted."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("LTO.cell").write_text("#SPECIES=Li,Ti,O\n")
        result = runner.invoke(
            cli,
            [
                "run",
                "search",
                "--seed",
                "LTO",
                "--build-only",
                "--formula-elements",
                "Li,Ti,O",
            ],
        )

    assert result.exit_code != 0
    assert "No such option" in result.output


def test_run_search_formula_transform_passed_to_buildcell():
    """Test formula sampling passes a seed transform into run_buildcell."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("Si.cell").write_text("#SPECIES=Si\n#NATOM=2\n")
        with patch("airsspy.jf.runners.run_buildcell") as mock_buildcell:
            mock_buildcell.return_value = {
                "struct_name": "Si-001",
                "seed_name": "Si",
                "struct_content": "",
            }
            result = runner.invoke(
                cli,
                [
                    "run",
                    "search",
                    "--seed",
                    "Si",
                    "--nmax",
                    "1",
                    "--build-only",
                    "--formula",
                    "Si2",
                ],
            )

    assert result.exit_code == 0
    transform = mock_buildcell.call_args.kwargs["seed_text_transform"]
    assert transform is not None
    transformed = transform("#SPECIES=Si\n#NATOM=2\n")
    assert "#FORMULA=Si" in transformed
    assert "#SPECIES=Si" not in transformed


def test_run_search_volume_minsep_transform_resolves_paths_before_chdir():
    """Test estimator dataset paths still work after run_search enters workdir."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("SiO.cell").write_text("#SPECIES=Si,O\n#NFORM=1\n")
        Path("work").mkdir()
        dataset = Path("curated.json")
        dataset.write_text(
            json.dumps(
                [
                    {
                        "material_id": "mp-sio2",
                        "reduced_formula": "SiO2",
                        "chemical_formula": "SiO2",
                        "composition": {"Si": 1.0, "O": 2.0},
                        "volume": 45.0,
                        "minsep": {
                            "Si-O": 1.6,
                            "O-Si": 1.6,
                            "O-O": 2.5,
                            "Si-Si": 3.0,
                        },
                    }
                ]
            ),
            encoding="utf-8",
        )
        with patch("airsspy.jf.runners.run_buildcell") as mock_buildcell:
            mock_buildcell.return_value = {
                "struct_name": "SiO-001",
                "seed_name": "SiO",
                "struct_content": "",
            }
            result = runner.invoke(
                cli,
                [
                    "run",
                    "search",
                    "--seed",
                    "SiO",
                    "--nmax",
                    "1",
                    "--build-only",
                    "--workdir",
                    "work",
                    "--formula",
                    "SiO2",
                    "--volume-minsep-source",
                    "dataset",
                    "--volume-minsep-dataset",
                    str(dataset),
                ],
            )

            assert result.exit_code == 0
            transform = mock_buildcell.call_args.kwargs["seed_text_transform"]
            transformed = transform("#SPECIES=Si,O\n#NFORM=1\n")
            assert "#FORMULA=SiO2" in transformed
            assert "#MINSEP=" in transformed


def test_run_search_prune_rejects_build_only():
    """Test pruning is rejected for build-only searches."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("Si.cell").write_text("#SPECIES=Si\n")
        result = runner.invoke(
            cli,
            ["run", "search", "--seed", "Si", "--build-only", "--prune"],
        )

    assert result.exit_code != 0
    assert "--prune cannot be used with --build-only" in result.output


def test_walltime_check_ignores_scheduler_parse_errors(caplog):
    """Malformed scheduler walltime data should not abort a run loop."""

    class BrokenScheduler:
        def get_remaining_seconds(self):
            raise ValueError("bad scheduler time")

    with caplog.at_level("WARNING"):
        assert cmd_run._walltime_remaining_ok(BrokenScheduler(), 300)
    assert "Could not determine remaining walltime" in caplog.text


def test_walltime_check_stops_when_buffer_exceeded():
    """Low remaining walltime should still stop the run loop."""

    class LowTimeScheduler:
        def get_remaining_seconds(self):
            return 10

    assert not cmd_run._walltime_remaining_ok(LowTimeScheduler(), 300)


def test_run_relax_help():
    """Test 'run relax --help'."""
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "relax", "--help"])
    assert result.exit_code == 0
    assert "--cell" in result.output
    assert "--cell-axis-map" in result.output


def test_run_sp_help():
    """Test 'run sp --help'."""
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "sp", "--help"])
    assert result.exit_code == 0
    assert "--cell" in result.output
    assert "--cell-axis-map" in result.output


def test_run_search_missing_seed():
    """Test 'run search' fails without seed files."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(
            cli, ["run", "search", "--seed", "NonExistent", "--nmax", "1"]
        )
        assert result.exit_code != 0
        assert "not found" in result.output.lower()


def test_run_search_missing_param():
    """Test 'run search' fails when param file is missing."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        from pathlib import Path

        Path("Si.cell").write_text(
            "%BLOCK LATTICE_CART\n5.43 0 0\n0 5.43 0\n0 0 5.43\n"
            "%ENDBLOCK LATTICE_CART\n"
            "%BLOCK POSITIONS_FRAC\nSi 0.0 0.0 0.0\n%ENDBLOCK POSITIONS_FRAC\n"
        )
        result = runner.invoke(cli, ["run", "search", "--seed", "Si", "--nmax", "1"])
        assert result.exit_code != 0
        assert "not found" in result.output.lower()


def test_run_search_build_only_missing_seed():
    """Test 'run search --build-only' fails without seed cell file."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(
            cli,
            ["run", "search", "--seed", "Missing", "--nmax", "1", "--build-only"],
        )
        assert result.exit_code != 0
        assert "not found" in result.output.lower()


def test_run_relax_no_files():
    """Test 'run relax' fails when no files match pattern."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli, ["run", "relax", "--cell", "*.cell"])
        assert result.exit_code != 0
        assert "No files matched" in result.output
