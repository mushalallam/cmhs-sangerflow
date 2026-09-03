"""Pairwise alignment helpers and quality-aware consensus generation."""

from __future__ import annotations

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

from .models import ConsensusResult, ReadData

IUPAC_BASES = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
BASES_IUPAC = {bases: code for code, bases in IUPAC_BASES.items()}


def reverse_complement_read(read: ReadData) -> ReadData:
    read.trimmed_sequence = str(Seq(read.trimmed_sequence).reverse_complement())
    read.trimmed_qualities = list(reversed(read.trimmed_qualities))
    return read


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -3.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    aligner.open_end_gap_score = -2.0
    aligner.extend_end_gap_score = -0.5
    return aligner


def gapped_alignment(target: str, query: str) -> tuple[str, str]:
    alignment = make_aligner().align(target, query)[0]
    coordinates = alignment.coordinates
    target_parts: list[str] = []
    query_parts: list[str] = []
    for index in range(coordinates.shape[1] - 1):
        t0, t1 = int(coordinates[0, index]), int(coordinates[0, index + 1])
        q0, q1 = int(coordinates[1, index]), int(coordinates[1, index + 1])
        t_len, q_len = t1 - t0, q1 - q0
        if t_len and q_len:
            target_parts.append(target[t0:t1])
            query_parts.append(query[q0:q1])
        elif t_len:
            target_parts.append(target[t0:t1])
            query_parts.append("-" * t_len)
        elif q_len:
            target_parts.append("-" * q_len)
            query_parts.append(query[q0:q1])
    return "".join(target_parts), "".join(query_parts)


def _choose_base(
    base_a: str, qa: int, base_b: str, qb: int, quality_delta: int
) -> tuple[str, int, bool]:
    if base_a == base_b:
        return base_a, max(qa, qb), False
    if base_a == "N":
        return base_b, qb, False
    if base_b == "N":
        return base_a, qa, False
    if abs(qa - qb) >= quality_delta:
        return (base_a, qa, True) if qa > qb else (base_b, qb, True)
    union = IUPAC_BASES.get(base_a, frozenset("ACGT")) | IUPAC_BASES.get(base_b, frozenset("ACGT"))
    return BASES_IUPAC.get(frozenset(union), "N"), min(qa, qb), True


def build_consensus(
    forward: ReadData | None,
    reverse: ReadData | None,
    *,
    quality_delta: int = 8,
) -> ConsensusResult:
    if not forward and not reverse:
        raise ValueError("At least one read is required")
    if forward and not forward.trimmed_sequence:
        raise ValueError(f"Forward read {forward.name} failed trimming")
    if reverse and not reverse.trimmed_sequence:
        raise ValueError(f"Reverse read {reverse.name} failed trimming")
    if forward is None:
        assert reverse is not None
        sequence = reverse.trimmed_sequence
        return ConsensusResult(
            sequence,
            reverse.trimmed_qualities,
            "-" * len(sequence),
            sequence,
            sequence,
            0,
            sequence.count("N"),
        )
    if reverse is None:
        sequence = forward.trimmed_sequence
        return ConsensusResult(
            sequence,
            forward.trimmed_qualities,
            sequence,
            "-" * len(sequence),
            sequence,
            0,
            sequence.count("N"),
        )

    aligned_f, aligned_r = gapped_alignment(forward.trimmed_sequence, reverse.trimmed_sequence)
    result: list[str] = []
    qualities: list[int] = []
    displayed: list[str] = []
    f_index = r_index = conflicts = overlap_bases = overlap_matches = 0
    for base_f, base_r in zip(aligned_f, aligned_r, strict=True):
        qf = forward.trimmed_qualities[f_index] if base_f != "-" else 0
        qr = reverse.trimmed_qualities[r_index] if base_r != "-" else 0
        if base_f != "-":
            f_index += 1
        if base_r != "-":
            r_index += 1
        if base_f == "-":
            base, quality, conflict = base_r, qr, False
        elif base_r == "-":
            base, quality, conflict = base_f, qf, False
        else:
            overlap_bases += 1
            overlap_matches += int(base_f == base_r and base_f in "ACGT")
            base, quality, conflict = _choose_base(base_f, qf, base_r, qr, quality_delta)
        result.append(base)
        displayed.append(base)
        qualities.append(quality)
        conflicts += int(conflict)
    sequence = "".join(result)
    return ConsensusResult(
        sequence=sequence,
        qualities=qualities,
        forward_aligned=aligned_f,
        reverse_aligned=aligned_r,
        consensus_aligned="".join(displayed),
        conflicts=conflicts,
        ambiguous_bases=sum(base not in "ACGT" for base in sequence),
        overlap_bases=overlap_bases,
        overlap_identity=overlap_matches / overlap_bases if overlap_bases else None,
    )
