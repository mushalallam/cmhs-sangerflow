from pathlib import Path

import pytest

from sangerflow.models import ReadData
from sangerflow.quality import mott_trim_bounds, quality_summary, trim_and_mask


def test_mott_trims_low_quality_ends():
    qualities = [2] * 10 + [35] * 100 + [2] * 10
    start, end = mott_trim_bounds(qualities)
    assert (start, end) == (10, 110)


def test_mott_rejects_invalid_cutoff():
    with pytest.raises(ValueError, match="between 0 and 1"):
        mott_trim_bounds([20], 1.0)


def test_trim_masks_internal_low_quality_base():
    read = ReadData(
        "read", Path("read.ab1"), "A" * 50 + "C" + "G" * 50, [30] * 50 + [10] + [30] * 50
    )
    trim_and_mask(read, min_length=50, mask_below=20)
    assert read.trimmed_sequence[50] == "N"
    assert quality_summary(read)["status"] == "PASS"


def test_short_read_fails():
    read = ReadData("short", Path("short.ab1"), "A" * 20, [30] * 20)
    trim_and_mask(read, min_length=80)
    assert not read.trimmed_sequence
    assert quality_summary(read)["status"] == "FAIL_SHORT"


def test_missing_peak_data_is_json_compatible():
    read = ReadData("read", Path("read.ab1"), "A" * 100, [30] * 100)
    trim_and_mask(read, min_length=50)
    assert quality_summary(read)["median_secondary_ratio"] is None
