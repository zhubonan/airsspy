"""Tests for VASP input helpers."""

from pathlib import Path
from types import SimpleNamespace

import pytest

from airsspy import vasptools


def test_parse_potcar_map():
    assert vasptools.parse_potcar_map(("Si=Si_GW", "O=O_s")) == {
        "Si": "Si_GW",
        "O": "O_s",
    }


def test_parse_potcar_map_rejects_bad_value():
    with pytest.raises(ValueError):
        vasptools.parse_potcar_map(("Si",))


def test_split_incar_settings_strips_control_keys():
    control, user = vasptools.split_incar_settings(
        {
            "ENCUT": 520,
            "AIRSSPY_VASP_INPUT_SET": "MPRelaxSet",
            "airsspy_custom": "ignored",
        }
    )
    assert control["AIRSSPY_VASP_INPUT_SET"] == "MPRelaxSet"
    assert control["AIRSSPY_CUSTOM"] == "ignored"
    assert user == {"ENCUT": 520}


def test_resolve_input_set_short_name(monkeypatch):
    cls = object()
    monkeypatch.setattr(
        vasptools.importlib,
        "import_module",
        lambda name: SimpleNamespace(MPRelaxSet=cls),
    )
    assert vasptools.resolve_input_set_class("MPRelaxSet") is cls


def test_resolve_input_set_short_name_case_insensitive(monkeypatch):
    cls = object()
    monkeypatch.setattr(
        vasptools.importlib,
        "import_module",
        lambda name: SimpleNamespace(MPStaticSet=cls),
    )
    assert vasptools.resolve_input_set_class("Mpstaticset") is cls


def test_assemble_potcar_uses_map_and_hash(tmp_path):
    potcars = tmp_path / "potcars"
    (potcars / "Si_GW").mkdir(parents=True)
    (potcars / "O_s").mkdir()
    (potcars / "Si_GW" / "POTCAR").write_text("si-potcar\n")
    (potcars / "O_s" / "POTCAR").write_text("o-potcar\n")

    structure = [
        SimpleNamespace(specie=SimpleNamespace(symbol="Si")),
        SimpleNamespace(specie=SimpleNamespace(symbol="O")),
        SimpleNamespace(specie=SimpleNamespace(symbol="Si")),
    ]
    metadata = vasptools.assemble_potcar(
        structure,
        tmp_path / "POTCAR",
        potcar_dir=str(potcars),
        potcar_map={"Si": "Si_GW", "O": "O_s"},
    )

    assert (tmp_path / "POTCAR").read_text() == "si-potcar\no-potcar\n"
    assert [item["symbol"] for item in metadata] == ["Si_GW", "O_s"]
    assert all("sha256" in item for item in metadata)


def test_assemble_potcar_prefers_resolved_symbols(tmp_path):
    potcars = tmp_path / "potcars"
    (potcars / "Si_sv").mkdir(parents=True)
    (potcars / "O_override").mkdir()
    (potcars / "Si_sv" / "POTCAR").write_text("si-sv\n")
    (potcars / "O_override" / "POTCAR").write_text("o-override\n")

    structure = [
        SimpleNamespace(specie=SimpleNamespace(symbol="Si")),
        SimpleNamespace(specie=SimpleNamespace(symbol="O")),
    ]
    metadata = vasptools.assemble_potcar(
        structure,
        tmp_path / "POTCAR",
        potcar_dir=str(potcars),
        potcar_map={"O": "O_override"},
        resolved_symbols={"Si": "Si_sv"},
    )

    assert (tmp_path / "POTCAR").read_text() == "si-sv\no-override\n"
    assert [item["symbol"] for item in metadata] == ["Si_sv", "O_override"]


def test_assemble_potcar_finds_common_pymatgen_layout(tmp_path):
    potcars = tmp_path / "potcars"
    (potcars / "potpaw_PBE" / "Si").mkdir(parents=True)
    (potcars / "potpaw_PBE" / "Si" / "POTCAR").write_text("si-potcar\n")

    structure = [SimpleNamespace(specie=SimpleNamespace(symbol="Si"))]
    vasptools.assemble_potcar(structure, tmp_path / "POTCAR", potcar_dir=str(potcars))

    assert (tmp_path / "POTCAR").read_text() == "si-potcar\n"


def test_prepare_vasp_inputs_uses_configurable_input_set(monkeypatch, tmp_path):
    calls = {}

    class FakeIncar(dict):
        @classmethod
        def from_str(cls, text):
            del text
            return cls({"AIRSSPY_VASP_INPUT_SET": "FakeSet", "ENCUT": 400})

        @classmethod
        def from_file(cls, path):
            content = Path(path).read_text()
            assert "AIRSSPY" not in content
            return cls({"ENCUT": 400})

        def write_file(self, path):
            Path(path).write_text("\n".join(f"{k} = {v}" for k, v in self.items()))

    class FakePoscar:
        structure = [
            SimpleNamespace(specie=SimpleNamespace(symbol="O")),
            SimpleNamespace(specie=SimpleNamespace(symbol="Si")),
        ]

        def write_file(self, path):
            Path(path).write_text("POSCAR\n")

    class FakeKpoints:
        @classmethod
        def from_file(cls, path):
            calls["kpoints_path"] = str(path)
            return cls()

        def write_file(self, path):
            Path(path).write_text("KPOINTS\n")

    class FakeSet:
        def __init__(self, structure, **kwargs):
            calls["structure"] = structure
            calls["kwargs"] = kwargs
            self.incar = FakeIncar(kwargs["user_incar_settings"])
            self.poscar = FakePoscar()
            self.kpoints = FakeKpoints()

    monkeypatch.setattr(vasptools, "parse_seed_incar", lambda content: FakeIncar.from_str(content))
    monkeypatch.setattr(vasptools, "resolve_input_set_class", lambda name: FakeSet)
    def fake_assemble(structure, *args, **kwargs):
        calls["potcar_structure"] = structure
        return []

    monkeypatch.setattr(vasptools, "assemble_potcar", fake_assemble)
    monkeypatch.setattr(vasptools, "_patch_incar", lambda path, settings: None)

    kpoints = tmp_path / "seed.KPOINTS"
    kpoints.write_text("Automatic mesh\n0\nGamma\n1 1 1\n0 0 0\n")
    potcars = tmp_path / "potcars"
    potcars.mkdir()

    metadata = vasptools.prepare_vasp_inputs(
        "Si-001",
        object(),
        "fake",
        mode="sp",
        pressure=5.0,
        potcar_dir=str(potcars),
        potcar_map={"Si": "Si_GW"},
        kpoints_path=kpoints,
        workdir=tmp_path / "run",
    )

    settings = calls["kwargs"]["user_incar_settings"]
    assert settings["ENCUT"] == 400
    assert settings["PSTRESS"] == pytest.approx(50.0)
    assert settings["NSW"] == 0
    assert settings["IBRION"] == -1
    assert calls["kwargs"]["user_potcar_settings"] == {"Si": "Si_GW"}
    assert [site.specie.symbol for site in calls["potcar_structure"]] == ["O", "Si"]
    assert metadata["input_set"] == "FakeSet"
    assert metadata["kpoints_source"] == str(kpoints)


def test_compose_vasp_task_doc_requires_energy(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "Si-001.vasp").mkdir()
    monkeypatch.setattr(vasptools, "_parse_vasprun", lambda workdir: {})
    monkeypatch.setattr(vasptools, "_parse_text_outputs", lambda workdir: {})
    monkeypatch.setattr(
        vasptools,
        "_read_structure_from_vasp",
        lambda workdir: SimpleNamespace(volume=1.0, composition=SimpleNamespace()),
    )

    with pytest.raises(ValueError, match="No VASP energy"):
        vasptools.compose_vasp_task_doc("Si-001")
