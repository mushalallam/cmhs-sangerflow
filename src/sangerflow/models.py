"""Typed data models shared by the SangerFlow pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True)
class PeakEvidence:
    index: int
    trace_position: int
    called_base: str
    quality: int
    primary_base: str
    primary_signal: int
    secondary_base: str
    secondary_signal: int
    secondary_ratio: float
    suggested_iupac: str


@dataclass(slots=True)
class ReadData:
    name: str
    path: Path
    sequence: str
    qualities: list[int]
    peak_evidence: list[PeakEvidence] = field(default_factory=list)
    trace_channels: dict[str, list[int]] = field(default_factory=dict)
    peak_locations: list[int] = field(default_factory=list)
    trim_start: int = 0
    trim_end: int = 0
    trimmed_sequence: str = ""
    trimmed_qualities: list[int] = field(default_factory=list)

    @property
    def raw_length(self) -> int:
        return len(self.sequence)

    @property
    def trimmed_length(self) -> int:
        return len(self.trimmed_sequence)


@dataclass(frozen=True, slots=True)
class SampleInput:
    sample: str
    forward: Path | None = None
    reverse: Path | None = None


@dataclass(slots=True)
class ConsensusResult:
    sequence: str
    qualities: list[int]
    forward_aligned: str
    reverse_aligned: str
    consensus_aligned: str
    conflicts: int
    ambiguous_bases: int
    overlap_bases: int = 0
    overlap_identity: float | None = None


@dataclass(frozen=True, slots=True)
class ReferenceAlignment:
    aligned_reference: str
    aligned_consensus: str
    reference_start: int
    reference_end: int
    aligned_columns: int
    matches: int
    identity: float
    consensus_bases: int
    consensus_coverage: float


@dataclass(frozen=True, slots=True)
class Variant:
    sample: str
    position: int
    reference: str
    alternate: str
    kind: str
    quality: float
    filter: str = "PASS"
    note: str = ""
