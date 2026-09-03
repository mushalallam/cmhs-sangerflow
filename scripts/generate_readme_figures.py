"""Generate deterministic synthetic figures embedded in the project README."""

from __future__ import annotations

import math
from pathlib import Path

from sangerflow.models import ReadData
from sangerflow.output import (
    write_qc_dashboard_svg,
    write_quality_svg,
    write_trace_window_svg,
)

OUTPUT = Path(__file__).resolve().parents[1] / "docs" / "images"


def synthetic_read() -> ReadData:
    sequence_bases = list("ACGTTGCAAGTCGATCGTACATGCCGATTCGAGTACGTA")
    center = len(sequence_bases) // 2
    sequence_bases[center] = "T"
    sequence = "".join(sequence_bases)
    peak_locations = [18 + index * 18 for index in range(len(sequence))]
    channel_length = peak_locations[-1] + 20
    channels = {
        base: [6 + int(2 * math.sin(index / 13)) for index in range(channel_length)]
        for base in "ACGT"
    }
    qualities = []
    for index, (base, location) in enumerate(zip(sequence, peak_locations, strict=True)):
        quality = max(12, min(58, 55 - abs(center - index) // 2))
        qualities.append(quality)
        amplitude = 930 if index != center else 820
        for trace_index in range(max(0, location - 16), min(channel_length, location + 17)):
            channels[base][trace_index] += int(
                amplitude * math.exp(-((trace_index - location) ** 2) / 28)
            )
        if index == center:
            for trace_index in range(location - 16, location + 17):
                channels["C"][trace_index] += int(
                    360 * math.exp(-((trace_index - location) ** 2) / 32)
                )
    return ReadData(
        name="SYNTHETIC_SAMPLE_forward",
        path=Path("SYNTHETIC_SAMPLE_forward.ab1"),
        sequence=sequence,
        qualities=qualities,
        trace_channels=channels,
        peak_locations=peak_locations,
        trim_start=3,
        trim_end=len(sequence) - 3,
        trimmed_sequence=sequence[3:-3],
        trimmed_qualities=qualities[3:-3],
    )


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    read = synthetic_read()
    center = len(read.sequence) // 2
    write_trace_window_svg(
        OUTPUT / "example-variant-evidence.svg",
        read,
        center,
        "Synthetic example · SNV C>T at DEMO_REF:123",
    )
    write_quality_svg(OUTPUT / "example-quality-profile.svg", read)
    rows = [
        {
            "sample": f"SYNTHETIC-{index + 1:02d}",
            "direction": direction,
            "q20_fraction": fraction,
            "trimmed_length": length,
            "status": status,
        }
        for index, (direction, fraction, length, status) in enumerate(
            [
                ("forward", 0.98, 412, "PASS"),
                ("reverse", 0.96, 398, "PASS"),
                ("forward", 0.91, 355, "PASS"),
                ("reverse", 0.87, 341, "PASS"),
                ("forward", 0.74, 188, "PASS"),
                ("reverse", 0.00, 0, "FAIL_SHORT"),
            ]
        )
    ]
    write_qc_dashboard_svg(OUTPUT / "example-batch-qc.svg", rows)


if __name__ == "__main__":
    main()
