from sangerflow.variants import call_variants


def test_snv_call():
    variants, aligned_ref, aligned_query = call_variants("sample", "AACCGG", "AATCGG", [35] * 6)
    assert len(variants) == 1
    assert (variants[0].position, variants[0].reference, variants[0].alternate) == (3, "C", "T")
    assert len(aligned_ref) == len(aligned_query)


def test_low_quality_filter():
    variants, _, _ = call_variants("sample", "AACCGG", "AATCGG", [35, 35, 8, 35, 35, 35])
    assert variants[0].filter == "LowQual"


def test_ambiguous_call_is_symbolic():
    variants, _, _ = call_variants("sample", "AACCGG", "AARCGG", [35] * 6)
    assert variants[0].alternate == "<AMBIG>"
    assert variants[0].filter == "Ambiguous"


def test_insertion_and_deletion_calls():
    inserted, _, _ = call_variants("sample", "AACCGG", "AACTCGG", [35] * 7)
    deleted, _, _ = call_variants("sample", "AACCGG", "AACGG", [35] * 5)
    assert any(variant.kind == "INS" for variant in inserted)
    assert any(variant.kind == "DEL" for variant in deleted)


def test_unsequenced_reference_flanks_are_not_deletions():
    variants, aligned_ref, aligned_query = call_variants(
        "sample", "TTTTAACCGGTTTT", "AACCGG", [35] * 6
    )
    assert variants == []
    assert aligned_ref == "AACCGG"
    assert aligned_query == "AACCGG"


def test_position_in_long_reference_uses_reference_coordinates():
    variants, _, _ = call_variants("sample", "TTTTAACCGGTTTT", "AATCGG", [35] * 6)
    assert len(variants) == 1
    assert variants[0].position == 7
