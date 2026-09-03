from pathlib import Path

from sangerflow.alignment import build_consensus, reverse_complement_read
from sangerflow.models import ReadData


def read(name: str, sequence: str, qualities: list[int]) -> ReadData:
    return ReadData(
        name=name,
        path=Path(f"{name}.ab1"),
        sequence=sequence,
        qualities=qualities,
        trimmed_sequence=sequence,
        trimmed_qualities=qualities,
    )


def test_reverse_complement_preserves_quality_orientation():
    value = read("reverse", "ACGTAA", [10, 11, 12, 13, 14, 15])
    reverse_complement_read(value)
    assert value.trimmed_sequence == "TTACGT"
    assert value.trimmed_qualities == [15, 14, 13, 12, 11, 10]


def test_consensus_uses_high_quality_disagreement():
    forward = read("forward", "ACGT", [30, 30, 35, 30])
    reverse = read("reverse", "ACAT", [30, 30, 10, 30])
    result = build_consensus(forward, reverse, quality_delta=8)
    assert result.sequence == "ACGT"
    assert result.conflicts == 1
    assert result.overlap_bases == 4
    assert result.overlap_identity == 0.75


def test_consensus_uses_iupac_when_qualities_are_similar():
    forward = read("forward", "ACGT", [30] * 4)
    reverse = read("reverse", "ACAT", [30] * 4)
    result = build_consensus(forward, reverse)
    assert result.sequence == "ACRT"
    assert result.ambiguous_bases == 1


def test_consensus_reports_unrelated_pair_overlap_identity():
    forward = read("forward", "A" * 100, [30] * 100)
    reverse = read("reverse", "T" * 100, [30] * 100)
    result = build_consensus(forward, reverse)
    assert result.overlap_bases == 0
    assert result.overlap_identity is None
