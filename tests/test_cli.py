import json

import pytest

from sangerflow.cli import main


def test_cli_version(capsys):
    with pytest.raises(SystemExit, match="0"):
        main(["--version"])
    assert "sangerflow 0.1.0" in capsys.readouterr().out


def test_cli_inspect_reports_json(monkeypatch, capsys, tmp_path):
    from sangerflow.models import ReadData

    path = tmp_path / "read.ab1"
    path.touch()
    read = ReadData("read", path, "A" * 100, [30] * 100)
    monkeypatch.setattr("sangerflow.cli.read_abi", lambda ignored: read)
    assert main(["inspect", str(path)]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["status"] == "PASS"
    assert result["trimmed_length"] == 100
