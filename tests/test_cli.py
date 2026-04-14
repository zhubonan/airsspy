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
