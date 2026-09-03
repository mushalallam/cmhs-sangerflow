"""ABI/AB1 chromatogram parsing and peak-level evidence extraction."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from .models import PeakEvidence, ReadData

IUPAC_PAIR = {
    frozenset(("A", "G")): "R",
    frozenset(("C", "T")): "Y",
    frozenset(("G", "C")): "S",
    frozenset(("A", "T")): "W",
    frozenset(("G", "T")): "K",
    frozenset(("A", "C")): "M",
}


def _raw_value(raw: dict, *keys: bytes):
    for key in keys:
        if key in raw:
            return raw[key]
        text_key = key.decode("ascii")
        if text_key in raw:
            return raw[text_key]
    return None


def extract_trace_data(record) -> tuple[dict[str, list[int]], list[int]]:
    raw = record.annotations.get("abif_raw", {})
    order_value = _raw_value(raw, b"FWO_1") or b"GATC"
    if isinstance(order_value, bytes):
        order = order_value.decode("ascii", errors="replace")
    else:
        order = str(order_value)
    order = order[:4].upper()
    if set(order) != set("ACGT"):
        return {}, []

    locations = _raw_value(raw, b"PLOC2", b"PLOC1")
    if locations is None:
        return {}, []
    channels: dict[str, list[int]] = {}
    for offset, base in enumerate(order, start=9):
        values = _raw_value(raw, f"DATA{offset}".encode())
        if values is None:
            return {}, []
        channels[base] = list(values)
    return channels, list(locations)


def extract_peak_evidence(record) -> list[PeakEvidence]:
    channels, locations = extract_trace_data(record)
    if not channels or not locations:
        return []

    sequence = str(record.seq).upper()
    qualities = list(record.letter_annotations.get("phred_quality", []))
    result: list[PeakEvidence] = []
    for index, (called, position) in enumerate(zip(sequence, locations, strict=False)):
        signals = sorted(
            (
                (base, channel[position])
                for base, channel in channels.items()
                if position < len(channel)
            ),
            key=lambda item: item[1],
            reverse=True,
        )
        if len(signals) < 2:
            continue
        (primary_base, primary), (secondary_base, secondary) = signals[:2]
        ratio = secondary / primary if primary else 0.0
        suggested = IUPAC_PAIR.get(frozenset((primary_base, secondary_base)), "N")
        result.append(
            PeakEvidence(
                index=index,
                trace_position=int(position),
                called_base=called,
                quality=qualities[index] if index < len(qualities) else 0,
                primary_base=primary_base,
                primary_signal=int(primary),
                secondary_base=secondary_base,
                secondary_signal=int(secondary),
                secondary_ratio=ratio,
                suggested_iupac=suggested,
            )
        )
    return result


def read_abi(path: Path) -> ReadData:
    path = Path(path)
    try:
        record = SeqIO.read(path, "abi")
    except Exception as exc:
        raise ValueError(f"Unable to parse ABI chromatogram {path}: {exc}") from exc
    qualities = list(record.letter_annotations.get("phred_quality", []))
    sequence = str(record.seq).upper()
    if not qualities or len(sequence) != len(qualities):
        raise ValueError(f"ABI file lacks complete base-call qualities: {path}")
    channels, locations = extract_trace_data(record)
    return ReadData(
        name=path.stem,
        path=path.resolve(),
        sequence=sequence,
        qualities=qualities,
        peak_evidence=extract_peak_evidence(record),
        trace_channels=channels,
        peak_locations=locations,
    )
