from pathlib import Path

import pytest

from plasmid_detection import config

REPO_CONFIG = Path(__file__).resolve().parents[1] / "config" / "config.yaml"


def write(tmp_path, text):
    path = tmp_path / "config.yaml"
    path.write_text(text)
    return path


def test_resolves_nested_references(tmp_path):
    cfg = config.load(write(tmp_path, "a: /x\nb: ${a}/y\nc:\n  d: ${b}/z\n"))
    assert cfg["c"]["d"] == "/x/y/z"


def test_circular_reference_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="circular"):
        config.load(write(tmp_path, "a: ${b}\nb: ${a}\n"))


def test_unknown_reference_is_an_error(tmp_path):
    with pytest.raises(KeyError, match="nope"):
        config.load(write(tmp_path, "a: ${nope}\n"))


def test_lists_resolve_and_flatten_space_separated(tmp_path):
    cfg = config.load(write(tmp_path, 'a: /x\nb: [5, "${a}/y"]\n'))
    assert cfg["b"] == [5, "/x/y"]
    assert dict(config.flatten(cfg))["B"] == "5 /x/y"


def test_flatten_names():
    assert dict(config.flatten({"db": {"gtdb": {"dir": "/g"}}, "x": 1})) == {
        "DB_GTDB_DIR": "/g",
        "X": 1,
    }


def test_repo_config_resolves_fully():
    flat = dict(config.flatten(config.load(REPO_CONFIG)))
    assert not [k for k, v in flat.items() if "${" in str(v)]
    # Slurm reads SLURM_* from the environment; the config must never set one
    assert not [k for k in flat if k.startswith("SLURM_")]
    assert flat["DB_GTDB_RELEASE"] == "r226"
    assert flat["REFERENCE_HUMAN_ACCESSION"] == "GCF_000001405.40"
