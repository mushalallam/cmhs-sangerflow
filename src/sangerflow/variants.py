"""Reference alignment and conservative small-variant extraction."""

from __future__ import annotations

from statistics import fmean

from .alignment import IUPAC_BASES, gapped_alignment
from .models import ReferenceAlignment, Variant


def _trim_uncovered_ends(
    aligned_reference: str, aligned_consensus: str
) -> tuple[str, str, int, int]:
    """Remove terminal alignment gaps, which represent noncoverage rather than variants."""
    start = 0
    reference_offset = 0
    query_offset = 0
    while start < len(aligned_reference) and (
        aligned_reference[start] == "-" or aligned_consensus[start] == "-"
    ):
        reference_offset += aligned_reference[start] != "-"
        query_offset += aligned_consensus[start] != "-"
        start += 1
    end = len(aligned_reference)
    while end > start and (aligned_reference[end - 1] == "-" or aligned_consensus[end - 1] == "-"):
        end -= 1
    return (
        aligned_reference[start:end],
        aligned_consensus[start:end],
        reference_offset,
        query_offset,
    )


def analyze_reference(
    sample: str,
    reference: str,
    consensus: str,
    qualities: list[int],
    *,
    min_quality: int = 20,
) -> tuple[list[Variant], ReferenceAlignment]:
    reference = reference.upper()
    consensus = consensus.upper()
    aligned_ref, aligned_consensus = gapped_alignment(reference, consensus)
    aligned_ref, aligned_consensus, reference_offset, query_offset = _trim_uncovered_ends(
        aligned_ref, aligned_consensus
    )
    if not aligned_ref:
        raise ValueError("Consensus does not align to the reference; verify sample and reference")
    variants: list[Variant] = []
    ref_pos = reference_offset
    query_pos = query_offset
    index = 0

    while index < len(aligned_ref):
        ref_base = aligned_ref[index]
        alt_base = aligned_consensus[index]
        if ref_base != "-" and alt_base != "-":
            ref_pos += 1
            quality = qualities[query_pos] if query_pos < len(qualities) else 0
            query_pos += 1
            if ref_base != alt_base:
                if alt_base in "ACGT":
                    filt = "PASS" if quality >= min_quality else "LowQual"
                    variants.append(
                        Variant(sample, ref_pos, ref_base, alt_base, "SNV", quality, filt)
                    )
                else:
                    bases = IUPAC_BASES.get(alt_base, frozenset("ACGT"))
                    note = f"consensus={alt_base};possible_bases={''.join(sorted(bases))}"
                    variants.append(
                        Variant(
                            sample,
                            ref_pos,
                            ref_base,
                            "<AMBIG>",
                            "AMBIGUOUS",
                            quality,
                            "Ambiguous",
                            note,
                        )
                    )
            index += 1
            continue

        if ref_base == "-":
            inserted: list[str] = []
            inserted_q: list[int] = []
            while index < len(aligned_ref) and aligned_ref[index] == "-":
                inserted.append(aligned_consensus[index])
                inserted_q.append(qualities[query_pos] if query_pos < len(qualities) else 0)
                query_pos += 1
                index += 1
            anchor_pos = max(1, ref_pos)
            anchor = reference[anchor_pos - 1]
            alt = anchor + "".join(inserted) if ref_pos > 0 else "".join(inserted) + anchor
            quality = fmean(inserted_q) if inserted_q else 0.0
            variants.append(
                Variant(
                    sample,
                    anchor_pos,
                    anchor,
                    alt,
                    "INS",
                    quality,
                    "PASS" if quality >= min_quality else "LowQual",
                )
            )
            continue

        start_ref = ref_pos
        deleted: list[str] = []
        while index < len(aligned_consensus) and aligned_consensus[index] == "-":
            deleted.append(aligned_ref[index])
            ref_pos += 1
            index += 1
        if start_ref > 0:
            anchor_pos = start_ref
            anchor = reference[anchor_pos - 1]
            ref_allele = anchor + "".join(deleted)
            alt_allele = anchor
        else:
            anchor_pos = 1
            next_base = reference[ref_pos] if ref_pos < len(reference) else "N"
            ref_allele = "".join(deleted) + next_base
            alt_allele = next_base
        flank_q = qualities[query_pos] if query_pos < len(qualities) else 0
        variants.append(
            Variant(
                sample,
                anchor_pos,
                ref_allele,
                alt_allele,
                "DEL",
                flank_q,
                "PASS" if flank_q >= min_quality else "LowQual",
            )
        )
    aligned_columns = len(aligned_ref)
    matches = sum(
        ref == query and ref in "ACGT"
        for ref, query in zip(aligned_ref, aligned_consensus, strict=True)
    )
    reference_start = reference_offset + 1
    reference_end = reference_offset + sum(base != "-" for base in aligned_ref)
    consensus_bases = sum(base != "-" for base in aligned_consensus)
    result = ReferenceAlignment(
        aligned_reference=aligned_ref,
        aligned_consensus=aligned_consensus,
        reference_start=reference_start,
        reference_end=reference_end,
        aligned_columns=aligned_columns,
        matches=matches,
        identity=matches / aligned_columns,
        consensus_bases=consensus_bases,
        consensus_coverage=consensus_bases / len(consensus),
    )
    return variants, result


def call_variants(
    sample: str,
    reference: str,
    consensus: str,
    qualities: list[int],
    *,
    min_quality: int = 20,
) -> tuple[list[Variant], str, str]:
    """Return calls and gapped sequences; retained as the compact public API."""
    variants, alignment = analyze_reference(
        sample, reference, consensus, qualities, min_quality=min_quality
    )
    return variants, alignment.aligned_reference, alignment.aligned_consensus
