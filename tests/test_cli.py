"""Tests for CLI commands."""

from click.testing import CliRunner

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
            f.write("%BLOCK LATTICE_CART\n5.43 0 0\n0 5.43 0\n0 0 5.43\n%ENDBLOCK LATTICE_CART\n")
        with open("Si.param", "w") as f:
            f.write("task: geometryoptimization\ncut_off_energy: 300\n")

        result = runner.invoke(
            cli,
            [
                "deploy",
                "search",
                "--seed", "Si",
                "--project", "test",
                "--num", "10",
                "--dryrun",
            ],
        )
        assert result.exit_code == 0
        assert "Project: test, Seed: Si" in result.output
        assert "Structures requested: 10" in result.output
        assert "task: geometryoptimization" in result.output


def test_db_commands_help():
    """Test all db subcommands have valid help."""
    runner = CliRunner()
    for cmd in ["list-projects", "list-seeds", "summary", "throughput", "retrieve-project"]:
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
    lines = [l for l in result.output.split("\n") if "Si-002" in l]
    assert len(lines) == 0


def test_rank_extxyz_file(tmp_path):
    """Test reading an extxyz file."""
    runner = CliRunner()
    xyz_file = tmp_path / "test.xyz"
    xyz_content = (
        '2\n'
        'Lattice="5.43 0.0 0.0 0.0 5.43 0.0 0.0 0.0 5.43" '
        'Properties=species:S:1:pos:R:3 energy=-10.625 label=Si-001\n'
        'Si 0.0 0.0 0.0\n'
        'Si 0.25 0.25 0.25\n'
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
    atoms = Atoms("Si2", positions=[[0, 0, 0], [1.3, 1.3, 1.3]], cell=[5.43] * 3, pbc=True)
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


def test_run_help():
    """Test 'run --help'."""
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--help"])
    assert result.exit_code == 0
    assert "search" in result.output
    assert "relax" in result.output


def test_run_search_help():
    """Test 'run search --help'."""
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "search", "--help"])
    assert result.exit_code == 0
    assert "--seed" in result.output
    assert "--nmax" in result.output
    assert "--build-only" in result.output
    assert "--code" in result.output


def test_run_relax_help():
    """Test 'run relax --help'."""
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "relax", "--help"])
    assert result.exit_code == 0
    assert "--cell" in result.output


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
        result = runner.invoke(
            cli, ["run", "search", "--seed", "Si", "--nmax", "1"]
        )
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
