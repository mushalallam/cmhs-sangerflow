"""Map reference variants back to contributing reads and trace peaks."""

from __future__ import annotations

from dataclasses import replace

from .alignment import IUPAC_BASES
from .models import ConsensusResult, ReadData, ReferenceAlignment, Variant


def variant_consensus_index(variant: Variant, alignment: ReferenceAlignment) -> int | None:
    """Return the ungapped consensus index nearest a variant's reference anchor."""
    reference_position = alignment.reference_start - 1
    consensus_index = 0
    fallback: int | None = None
    for reference_base, consensus_base in zip(
        alignment.aligned_reference, alignment.aligned_consensus, strict=True
    ):
        if reference_base != "-":
            reference_position += 1
        if consensus_base == "-":
            continue
        if reference_position == variant.position:
            if variant.kind == "INS" and reference_base == "-":
                return consensus_index
            fallback = consensus_index
            if variant.kind != "INS":
                return consensus_index
        consensus_index += 1
    return fallback


def read_index_at_consensus(
    aligned_read: str, consensus_index: int
) -> tuple[int, str] | None:
    read_index = 0
    for index, base in enumerate(aligned_read):
        if index == consensus_index:
            return (read_index, base) if base != "-" else None
        if base != "-":
            read_index += 1
    return None


def raw_base_index(read: ReadData, oriented_index: int) -> int:
    if read.is_reverse_complemented:
        return read.trim_end - 1 - oriented_index
    return read.trim_start + oriented_index


def _read_values(
    read: ReadData | None, aligned_read: str, consensus_index: int
) -> tuple[str, int | None, float | None, int | None]:
    if read is None:
        return "", None, None, None
    mapped = read_index_at_consensus(aligned_read, consensus_index)
    if mapped is None:
        return "", None, None, None
    oriented_index, base = mapped
    quality = (
        read.trimmed_qualities[oriented_index]
        if oriented_index < len(read.trimmed_qualities)
        else None
    )
    raw_index = raw_base_index(read, oriented_index)
    peak_by_index = {peak.index: peak for peak in read.peak_evidence}
    peak = peak_by_index.get(raw_index)
    ratio = round(peak.secondary_ratio, 4) if peak else None
    return base, quality, ratio, raw_index


def enrich_variant(
    variant: Variant,
    alignment: ReferenceAlignment,
    consensus: ConsensusResult,
    forward: ReadData | None,
    reverse: ReadData | None,
) -> tuple[Variant, dict[str, int]]:
    consensus_index = variant_consensus_index(variant, alignment)
    if consensus_index is None:
        return variant, {}
    forward_base, forward_quality, forward_ratio, forward_raw = _read_values(
        forward, consensus.forward_aligned, consensus_index
    )
    reverse_base, reverse_quality, reverse_ratio, reverse_raw = _read_values(
        reverse, consensus.reverse_aligned, consensus_index
    )
    centers = {
        direction: index
        for direction, index in (("forward", forward_raw), ("reverse", reverse_raw))
        if index is not None
    }
    support = "not_assessed"
    if variant.kind in {"SNV", "AMBIGUOUS"}:
        consensus_base = consensus.sequence[consensus_index]
        allowed = IUPAC_BASES.get(consensus_base, frozenset())
        directions = []
        if forward_base in allowed:
            directions.append("F")
        if reverse_base in allowed:
            directions.append("R")
        support = "+".join(directions) or "none"
    return (
        replace(
            variant,
            strand_support=support,
            forward_quality=forward_quality,
            reverse_quality=reverse_quality,
            forward_peak_ratio=forward_ratio,
            reverse_peak_ratio=reverse_ratio,
        ),
        centers,
    )
