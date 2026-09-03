from types import SimpleNamespace

from sangerflow.abi import extract_peak_evidence


def test_extract_peak_evidence_accepts_biopython_string_tags():
    raw = {
        "FWO_1": b"GATC",
        "PLOC2": (1, 2),
        "DATA9": (0, 100, 5),
        "DATA10": (0, 20, 120),
        "DATA11": (0, 10, 10),
        "DATA12": (0, 5, 5),
    }
    record = SimpleNamespace(
        seq="GA",
        annotations={"abif_raw": raw},
        letter_annotations={"phred_quality": [30, 31]},
    )
    evidence = extract_peak_evidence(record)
    assert len(evidence) == 2
    assert evidence[0].primary_base == "G"
    assert evidence[0].secondary_base == "A"
    assert evidence[0].secondary_ratio == 0.2
