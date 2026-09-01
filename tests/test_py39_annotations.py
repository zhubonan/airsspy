from pathlib import Path

FUTURE_IMPORT = "from __future__ import annotations"
FILES_WITH_PEP604_ANNOTATIONS = [
    "src/airsspy/casteptools.py",
    "src/airsspy/jf/runners.py",
    "src/airsspy/jf/ml_runners.py",
    "src/airsspy/tools/modcell.py",
    "src/airsspy/cli/cmd_run.py",
    "src/airsspy/cli/cmd_deploy.py",
    "src/airsspy/cli/cmd_rank.py",
    "tests/test_volume_minsep.py",
    "tests/test_cli.py",
    "tests/test_rem_extraction.py",
    "tests/test_jf_runners_castep.py",
    "tests/test_ranking.py",
    "tests/test_run_crud_e2e.py",
    "tests/test_modcell.py",
]


def test_pep604_annotations_are_safe_on_python39():
    root = Path(__file__).resolve().parents[1]
    missing = []
    for relative_path in FILES_WITH_PEP604_ANNOTATIONS:
        text = (root / relative_path).read_text()
        if "|" in text and FUTURE_IMPORT not in text:
            missing.append(relative_path)
    assert not missing
