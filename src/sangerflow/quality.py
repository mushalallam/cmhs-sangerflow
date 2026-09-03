"""Quality trimming, masking, and summary statistics."""

from __future__ import annotations

import statistics
from collections.abc import Sequence

from .models import ReadData


def mott_trim_bounds(qualities: Sequence[int], error_cutoff: float = 0.05) -> tuple[int, int]:
    """Return the highest-scoring Mott segment as a half-open interval.

    Each base contributes ``error_cutoff - 10**(-Q/10)``. The maximum scoring
    contiguous segment retains the region whose estimated error probabilities
    are consistently below the configured cutoff.
    """
    if not qualities:
        return 0, 0
    if not 0 < error_cutoff < 1:
        raise ValueError("error_cutoff must be between 0 and 1")

    best_score = 0.0
    running_score = 0.0
    best_start = 0
    best_end = 0
    candidate_start = 0

    for index, quality in enumerate(qualities):
        score = error_cutoff - 10 ** (-max(0, quality) / 10)
        running_score += score
        if running_score < 0:
            running_score = 0.0
            candidate_start = index + 1
        elif running_score > best_score:
            best_score = running_score
            best_start = candidate_start
            best_end = index + 1
    return best_start, best_end


def trim_and_mask(
    read: ReadData,
    *,
    error_cutoff: float = 0.05,
    mask_below: int = 20,
    min_length: int = 80,
) -> ReadData:
    start, end = mott_trim_bounds(read.qualities, error_cutoff)
    sequence = read.sequence[start:end].upper()
    qualities = read.qualities[start:end]
    masked = "".join(
        base if quality >= mask_below else "N"
        for base, quality in zip(sequence, qualities, strict=True)
    )
    read.trim_start = start
    read.trim_end = end
    read.trimmed_sequence = masked if len(masked) >= min_length else ""
    read.trimmed_qualities = qualities if len(masked) >= min_length else []
    return read


def quality_summary(read: ReadData) -> dict[str, int | float | str]:
    raw = read.qualities
    trimmed = read.trimmed_qualities
    sequence = read.trimmed_sequence
    ratios = [peak.secondary_ratio for peak in read.peak_evidence]
    return {
        "read": read.name,
        "file": str(read.path),
        "raw_length": read.raw_length,
        "trim_start": read.trim_start,
        "trim_end": read.trim_end,
        "trimmed_length": read.trimmed_length,
        "mean_q_raw": round(statistics.fmean(raw), 2) if raw else 0.0,
        "mean_q_trimmed": round(statistics.fmean(trimmed), 2) if trimmed else 0.0,
        "q20_fraction": round(sum(q >= 20 for q in trimmed) / len(trimmed), 4) if trimmed else 0.0,
        "q30_fraction": round(sum(q >= 30 for q in trimmed) / len(trimmed), 4) if trimmed else 0.0,
        "n_fraction": round(sequence.count("N") / len(sequence), 4) if sequence else 1.0,
        "median_secondary_ratio": round(statistics.median(ratios), 4) if ratios else None,
        "status": "PASS" if sequence else "FAIL_SHORT",
    }
