from pathlib import Path

from sangerflow.evidence import enrich_variant, raw_base_index, variant_consensus_index
from sangerflow.models import (
    ConsensusResult,
    PeakEvidence,
    ReadData,
    ReferenceAlignment,
    Variant,
)


def alignment() -> ReferenceAlignment:
    return ReferenceAlignment("AACCGG", "AATCGG", 1, 6, 6, 5, 5 / 6, 6, 1.0)


def read(name: str, reverse: bool = False) -> ReadData:
    peak = PeakEvidence(2, 20, "T", 35, "T", 500, "C", 50, 0.1, "Y")
    return ReadData(
        name,
        Path(f"{name}.ab1"),
        "AATCGG",
        [35] * 6,
        peak_evidence=[peak],
        trim_start=0,
        trim_end=6,
        trimmed_sequence="AATCGG",
        trimmed_qualities=[35] * 6,
        is_reverse_complemented=reverse,
    )


def test_variant_maps_to_consensus_index():
    variant = Variant("sample", 3, "C", "T", "SNV", 35)
    assert variant_consensus_index(variant, alignment()) == 2


def test_enrich_variant_reports_bidirectional_support():
    variant = Variant("sample", 3, "C", "T", "SNV", 35)
    consensus = ConsensusResult("AATCGG", [35] * 6, "AATCGG", "AATCGG", "AATCGG", 0, 0)
    enriched, centers = enrich_variant(
        variant, alignment(), consensus, read("forward"), read("reverse", reverse=True)
    )
    assert enriched.strand_support == "F+R"
    assert enriched.forward_quality == 35
    assert enriched.forward_peak_ratio == 0.1
    assert centers == {"forward": 2, "reverse": 3}


def test_reverse_read_index_maps_back_to_raw_trace():
    assert raw_base_index(read("reverse", reverse=True), 2) == 3
