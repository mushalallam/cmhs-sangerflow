"""Generate a patient-free synthetic result gallery for first-time users."""

from __future__ import annotations

import html
import math
from pathlib import Path

from .models import ReadData
from .output import write_qc_dashboard_svg, write_quality_svg, write_trace_window_svg


def _synthetic_read() -> ReadData:
    sequence = "ACGTTGCAAGTCGATCGTACTTGCCGATTCGAGTACGTA"
    center = len(sequence) // 2
    locations = [18 + index * 18 for index in range(len(sequence))]
    channel_length = locations[-1] + 20
    channels = {
        base: [6 + int(2 * math.sin(index / 13)) for index in range(channel_length)]
        for base in "ACGT"
    }
    qualities: list[int] = []
    for index, (base, location) in enumerate(zip(sequence, locations, strict=True)):
        qualities.append(max(12, min(58, 55 - abs(center - index) // 2)))
        amplitude = 820 if index == center else 930
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
        peak_locations=locations,
        trim_start=3,
        trim_end=len(sequence) - 3,
        trimmed_sequence=sequence[3:-3],
        trimmed_qualities=qualities[3:-3],
    )


def create_demo_results(output_dir: Path) -> Path:
    """Create a compact synthetic report and return its HTML path."""
    output_dir = Path(output_dir)
    traces = output_dir / "traces"
    reports = output_dir / "reports"
    traces.mkdir(parents=True)
    reports.mkdir()
    read = _synthetic_read()
    write_quality_svg(traces / "quality.svg", read)
    write_trace_window_svg(
        traces / "variant.svg",
        read,
        len(read.sequence) // 2,
        "Synthetic example · SNV C>T at DEMO_REF:123",
    )
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
    write_qc_dashboard_svg(traces / "batch_qc.svg", rows)
    title = "CMHS SangerFlow — synthetic demonstration"
    report = reports / "report.html"
    report.write_text(
        "<!doctype html><meta charset='utf-8'>"
        f"<title>{html.escape(title)}</title>"
        "<style>body{font:16px system-ui;max-width:1200px;margin:35px auto;padding:0 20px;"
        "color:#17324d}img{width:100%;border:1px solid #d9e2ea;border-radius:10px;"
        "margin:8px 0 28px}.note{background:#eaf5ff;padding:14px;border-radius:8px}</style>"
        f"<h1>{html.escape(title)}</h1>"
        "<p class='note'><strong>Demonstration only.</strong> All sequences and measurements "
        "on this page are deterministic synthetic examples; no patient or laboratory data "
        "are included.</p><h2>Variant-centred chromatogram evidence</h2>"
        "<img src='../traces/variant.svg' alt='Synthetic chromatogram'>"
        "<h2>Per-read Phred quality</h2>"
        "<img src='../traces/quality.svg' alt='Synthetic quality profile'>"
        "<h2>Batch Q20 dashboard</h2>"
        "<img src='../traces/batch_qc.svg' alt='Synthetic batch dashboard'>",
        encoding="utf-8",
    )
    return report
