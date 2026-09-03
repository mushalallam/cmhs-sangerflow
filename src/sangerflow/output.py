"""Output writers for FASTA/FASTQ, tables, VCF, alignments, and HTML reports."""

from __future__ import annotations

import csv
import html
import json
from collections.abc import Iterable, Sequence
from pathlib import Path

from .models import ConsensusResult, ReadData, Variant

TRACE_COLORS = {"A": "#16883e", "C": "#2468d8", "G": "#20252a", "T": "#d9363e"}


def safe_name(value: str) -> str:
    cleaned = "".join(char if char.isalnum() or char in "._-" else "_" for char in value)
    return cleaned.strip("._") or "sample"


def write_fasta(path: Path, name: str, sequence: str, width: int = 80) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), width):
            handle.write(sequence[start : start + width] + "\n")


def write_multi_fasta(path: Path, records: Iterable[tuple[str, str]], width: int = 80) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")


def write_fastq(path: Path, name: str, sequence: str, qualities: Sequence[int]) -> None:
    encoded = "".join(chr(min(93, max(0, int(q))) + 33) for q in qualities)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"@{name}\n{sequence}\n+\n{encoded}\n")


def write_alignment(path: Path, sample: str, target: str, query: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">reference\n{target}\n>{sample}\n{query}\n")


def write_pair_alignment(path: Path, sample: str, consensus: ConsensusResult) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            f">{sample}_forward\n{consensus.forward_aligned}\n"
            f">{sample}_reverse_rc\n{consensus.reverse_aligned}\n"
            f">{sample}_consensus\n{consensus.consensus_aligned}\n"
        )


def write_tsv(path: Path, rows: Iterable[dict], fieldnames: Sequence[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_peak_table(path: Path, reads: Iterable[ReadData]) -> None:
    fields = [
        "read",
        "file",
        "base_index",
        "trace_position",
        "called_base",
        "quality",
        "primary_base",
        "primary_signal",
        "secondary_base",
        "secondary_signal",
        "secondary_ratio",
        "suggested_iupac",
        "within_trim",
    ]
    rows = []
    for read in reads:
        for peak in read.peak_evidence:
            rows.append(
                {
                    "read": read.name,
                    "file": str(read.path),
                    "base_index": peak.index + 1,
                    "trace_position": peak.trace_position,
                    "called_base": peak.called_base,
                    "quality": peak.quality,
                    "primary_base": peak.primary_base,
                    "primary_signal": peak.primary_signal,
                    "secondary_base": peak.secondary_base,
                    "secondary_signal": peak.secondary_signal,
                    "secondary_ratio": f"{peak.secondary_ratio:.4f}",
                    "suggested_iupac": peak.suggested_iupac,
                    "within_trim": read.trim_start <= peak.index < read.trim_end,
                }
            )
    write_tsv(path, rows, fields)


def write_trace_svg(path: Path, read: ReadData) -> bool:
    """Write a dependency-free four-channel chromatogram view."""
    if set(read.trace_channels) != set("ACGT") or not read.peak_locations:
        return False
    width, height = 1200, 360
    left, right, top, bottom = 55, 20, 35, 45
    plot_width = width - left - right
    plot_height = height - top - bottom
    channel_length = min(len(values) for values in read.trace_channels.values())
    if channel_length < 2:
        return False
    values = sorted(
        value for channel in read.trace_channels.values() for value in channel if value >= 0
    )
    scale = values[min(len(values) - 1, int(len(values) * 0.995))] if values else 1
    scale = max(1, scale)
    step = max(1, channel_length // plot_width)

    def x_position(index: int) -> float:
        return left + (index / (channel_length - 1)) * plot_width

    def polyline(channel: list[int]) -> str:
        points = []
        for index in range(0, channel_length, step):
            signal = min(max(channel[index], 0), scale)
            y = top + plot_height - signal / scale * plot_height
            points.append(f"{x_position(index):.1f},{y:.1f}")
        return " ".join(points)

    trim_start_index = min(read.trim_start, len(read.peak_locations) - 1)
    trim_end_index = max(read.trim_start, min(read.trim_end - 1, len(read.peak_locations) - 1))
    trim_start_position = read.peak_locations[trim_start_index]
    trim_end_position = read.peak_locations[trim_end_index]
    traces = "".join(
        f'<polyline points="{polyline(read.trace_channels[base])}" stroke="{TRACE_COLORS[base]}"/>'
        for base in "ACGT"
    )
    legend = " ".join(f'<tspan fill="{TRACE_COLORS[base]}">{base}</tspan>' for base in "ACGT")
    retained_width = max(0, x_position(trim_end_position) - x_position(trim_start_position))
    trim_x = x_position(trim_start_position)
    baseline_y = top + plot_height
    title = html.escape(read.name)
    document = f"""<svg xmlns="http://www.w3.org/2000/svg"
viewBox="0 0 {width} {height}" role="img">
<title>{title} chromatogram</title>
<rect width="100%" height="100%" fill="white"/>
<rect x="{trim_x:.1f}" y="{top}" width="{retained_width:.1f}"
height="{plot_height}" fill="#e7f5ee"/>
<g fill="none" stroke-width="1.15" stroke-linejoin="round">{traces}</g>
<line x1="{left}" y1="{baseline_y}" x2="{width - right}" y2="{baseline_y}"
stroke="#87949d"/>
<text x="{left}" y="22" font-family="system-ui" font-size="14">
{title} — green area retained after trimming</text>
<text x="{width - 105}" y="22" font-family="system-ui" font-size="14">{legend}</text>
<text x="{left}" y="{height - 12}" font-family="system-ui" font-size="12">trace position 0</text>
<text x="{width - 150}" y="{height - 12}" font-family="system-ui" font-size="12">
{channel_length - 1}</text>
</svg>"""
    path.write_text(document, encoding="utf-8")
    return True


def write_vcf(
    path: Path, reference_name: str, reference_length: int, variants: Iterable[Variant]
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.3\n")
        handle.write("##source=SangerFlow\n")
        handle.write(f"##contig=<ID={reference_name},length={reference_length}>\n")
        handle.write('##FILTER=<ID=LowQual,Description="Consensus quality below threshold">\n')
        handle.write('##FILTER=<ID=Ambiguous,Description="Ambiguous IUPAC consensus call">\n')
        handle.write('##ALT=<ID=AMBIG,Description="Ambiguous consensus base">\n')
        handle.write('##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample identifier">\n')
        handle.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">\n')
        handle.write('##INFO=<ID=NOTE,Number=1,Type=String,Description="Additional evidence">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for variant in sorted(variants, key=lambda value: (value.position, value.sample)):
            quality = f"{variant.quality:.1f}"
            note = variant.note.replace(";", ",").replace(" ", "_") or "."
            info = f"SAMPLE={safe_name(variant.sample)};TYPE={variant.kind};NOTE={note}"
            handle.write(
                f"{reference_name}\t{variant.position}\t.\t{variant.reference}\t"
                f"{variant.alternate}\t{quality}\t{variant.filter}\t{info}\n"
            )


def write_html_report(
    path: Path,
    run: dict,
    samples: list[dict],
    variants: list[dict],
    trace_files: list[dict] | None = None,
) -> None:
    def table(rows: list[dict]) -> str:
        if not rows:
            return "<p>None</p>"
        headers = list(rows[0])
        head = "".join(f"<th>{html.escape(str(value))}</th>" for value in headers)
        body = "".join(
            "<tr>"
            + "".join(f"<td>{html.escape(str(row.get(key, '')))}</td>" for key in headers)
            + "</tr>"
            for row in rows
        )
        return (
            "<div class='table-wrap'><table><thead><tr>"
            f"{head}</tr></thead><tbody>{body}</tbody></table></div>"
        )

    payload = html.escape(json.dumps(run, indent=2, sort_keys=True))
    trace_cards = (
        "".join(
            "<article><h3>"
            f"{html.escape(str(trace['sample']))} — {html.escape(str(trace['direction']))}"
            "</h3>"
            f"<img loading='lazy' src='../traces/{html.escape(str(trace['file']))}' "
            f"alt='{html.escape(str(trace['sample']))} chromatogram'></article>"
            for trace in (trace_files or [])
        )
        or "<p>Raw trace channels were unavailable.</p>"
    )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>SangerFlow report</title>
<style>
body{{font:15px system-ui,sans-serif;color:#17212b;background:#f5f7fa;margin:0}}
main{{max-width:1200px;margin:32px auto;padding:0 20px}} h1,h2{{color:#123c55}}
.card{{background:white;border:1px solid #d9e2e8;border-radius:10px;padding:20px;margin:18px 0}}
.table-wrap{{overflow:auto}} table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{border-bottom:1px solid #e6ecef;padding:8px;text-align:left;white-space:nowrap}}
th{{background:#eef5f7}}
pre{{white-space:pre-wrap;background:#102a37;color:#e7f5f7;padding:14px;border-radius:7px}}
.notice{{border-left:5px solid #d98b1d;padding-left:12px}}
img{{display:block;width:100%;border:1px solid #d9e2e8;background:white}}
</style></head><body><main><h1>SangerFlow analysis report</h1>
<div class="card notice"><strong>Research use:</strong> calls require chromatogram review and
independent validation before clinical use.</div>
<section class="card"><h2>Samples</h2>{table(samples)}</section>
<section class="card"><h2>Variants</h2>{table(variants)}</section>
<section class="card"><h2>Chromatograms</h2>{trace_cards}</section>
<section class="card"><h2>Run metadata</h2><pre>{payload}</pre></section>
</main></body></html>"""
    path.write_text(document, encoding="utf-8")
