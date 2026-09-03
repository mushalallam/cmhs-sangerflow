from pathlib import Path

import pytest

from sangerflow.pairing import auto_pair, infer_sample_and_direction, read_sample_sheet


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("GN20-4-VRK1-F_A09_001.ab1", ("GN20-4-VRK1", "forward")),
        ("sample_17_R_G08_014.ab1", ("sample_17", "reverse")),
        ("15.10F_003_C04.ab1", ("15.10", "forward")),
        ("case.reverse.ab1", ("case", "reverse")),
    ],
)
def test_infer_sample_and_direction(name, expected):
    assert infer_sample_and_direction(Path(name)) == expected


def test_auto_pair_requires_complete_pairs(tmp_path):
    (tmp_path / "sample_F_01.ab1").touch()
    with pytest.raises(ValueError, match="missing a reverse"):
        auto_pair(tmp_path)


def test_auto_pair_rejects_duplicates(tmp_path):
    (tmp_path / "sample_F_01.ab1").touch()
    (tmp_path / "sample_F_02.ab1").touch()
    (tmp_path / "sample_R_01.ab1").touch()
    with pytest.raises(ValueError, match="Multiple forward"):
        auto_pair(tmp_path)


def test_sample_sheet_is_explicit(tmp_path):
    (tmp_path / "f.ab1").touch()
    (tmp_path / "r.ab1").touch()
    sheet = tmp_path / "samples.csv"
    sheet.write_text("sample,forward,reverse\npatient-1,f.ab1,r.ab1\n", encoding="utf-8")
    samples = read_sample_sheet(sheet)
    assert samples[0].sample == "patient-1"
    assert samples[0].forward == (tmp_path / "f.ab1").resolve()
