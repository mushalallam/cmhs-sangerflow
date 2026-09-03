"""Output writers for FASTA/FASTQ, tables, VCF, alignments, and HTML reports."""

from __future__ import annotations

import csv
import html
import json
from collections.abc import Iterable, Sequence
from pathlib import Path

from .models import ConsensusResult, ReadData, Variant

TRACE_COLORS = {"A": "#16883e", "C": "#2468d8", "G": "#20252a", "T": "#d9363e"}
COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "N": "N",
}


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
<rect width="{width}" height="{height}" fill="white"/>
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


def write_quality_svg(path: Path, read: ReadData) -> bool:
    if not read.qualities:
        return False
    width, height = 1200, 300
    left, right, top, bottom = 55, 20, 35, 45
    plot_width = width - left - right
    plot_height = height - top - bottom
    maximum = max(60, max(read.qualities))

    def x_position(index: int) -> float:
        denominator = max(1, len(read.qualities) - 1)
        return left + index / denominator * plot_width

    points = " ".join(
        f"{x_position(index):.1f},{top + plot_height - quality / maximum * plot_height:.1f}"
        for index, quality in enumerate(read.qualities)
    )
    trim_x = x_position(read.trim_start)
    trim_width = max(0, x_position(max(read.trim_start, read.trim_end - 1)) - trim_x)
    q20_y = top + plot_height - 20 / maximum * plot_height
    title = html.escape(read.name)
    document = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}"
role="img"><title>{title} base-quality profile</title>
<rect width="{width}" height="{height}" fill="white"/>
<rect x="{trim_x:.1f}" y="{top}" width="{trim_width:.1f}" height="{plot_height}"
fill="#e7f5ee"/>
<line x1="{left}" y1="{q20_y:.1f}" x2="{width - right}" y2="{q20_y:.1f}"
stroke="#d97706" stroke-dasharray="6 4"/>
<polyline points="{points}" fill="none" stroke="#235789" stroke-width="1.5"/>
<line x1="{left}" y1="{top + plot_height}" x2="{width - right}"
y2="{top + plot_height}" stroke="#87949d"/>
<text x="{left}" y="22" font-family="system-ui" font-size="14">{title} — Phred quality</text>
<text x="{left + 5}" y="{q20_y - 5:.1f}" font-family="system-ui" font-size="12"
fill="#9a5800">Q20</text>
<text x="{left}" y="{height - 12}" font-family="system-ui" font-size="12">base 1</text>
<text x="{width - 95}" y="{height - 12}" font-family="system-ui" font-size="12">
base {len(read.qualities)}</text></svg>"""
    path.write_text(document, encoding="utf-8")
    return True


def write_trace_window_svg(
    path: Path,
    read: ReadData,
    center_base_index: int,
    label: str,
    *,
    reverse: bool = False,
    half_window: int = 12,
) -> bool:
    """Write a reference-oriented chromatogram window around one evidence base."""
    if set(read.trace_channels) != set("ACGT") or not read.peak_locations:
        return False
    if not 0 <= center_base_index < len(read.peak_locations):
        return False
    first = max(0, center_base_index - half_window)
    last = min(len(read.peak_locations), center_base_index + half_window + 1)
    peak_positions = read.peak_locations[first:last]
    if len(peak_positions) < 2:
        return False
    spacing = max(1, int((peak_positions[-1] - peak_positions[0]) / (len(peak_positions) - 1)))
    sample_start = max(0, peak_positions[0] - spacing)
    channel_length = min(len(values) for values in read.trace_channels.values())
    sample_end = min(channel_length - 1, peak_positions[-1] + spacing)
    if sample_end <= sample_start:
        return False
    width, height = 1200, 420
    left, right, top, plot_height = 55, 20, 45, 275
    plot_width = width - left - right
    signal_values = [
        read.trace_channels[base][index]
        for base in "ACGT"
        for index in range(sample_start, sample_end + 1)
    ]
    scale = max(1, max(signal_values, default=1))

    def x_position(index: int) -> float:
        fraction = (index - sample_start) / (sample_end - sample_start)
        return width - right - fraction * plot_width if reverse else left + fraction * plot_width

    def polyline(channel: list[int]) -> str:
        return " ".join(
            f"{x_position(index):.1f},"
            f"{top + plot_height - channel[index] / scale * plot_height:.1f}"
            for index in range(sample_start, sample_end + 1)
        )

    paths = []
    for display_base in "ACGT":
        raw_base = COMPLEMENT[display_base] if reverse else display_base
        paths.append(
            f'<polyline points="{polyline(read.trace_channels[raw_base])}" '
            f'stroke="{TRACE_COLORS[display_base]}"/>'
        )
    annotations = []
    for raw_index in range(first, last):
        base = read.sequence[raw_index] if raw_index < len(read.sequence) else "N"
        base = COMPLEMENT.get(base, "N") if reverse else base
        quality = read.qualities[raw_index] if raw_index < len(read.qualities) else 0
        x = x_position(read.peak_locations[raw_index])
        annotations.append(
            f'<text x="{x:.1f}" y="350" text-anchor="middle" font-family="monospace" '
            f'font-size="15" fill="{TRACE_COLORS.get(base, "#555")}">{base}</text>'
            f'<text x="{x:.1f}" y="371" text-anchor="middle" font-family="monospace" '
            f'font-size="11" fill="#5c6770">{quality}</text>'
        )
    center_x = x_position(read.peak_locations[center_base_index])
    title = html.escape(label)
    document = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}"
role="img"><title>{title}</title><rect width="{width}" height="{height}" fill="white"/>
<rect x="{center_x - 8:.1f}" y="{top}" width="16" height="335" fill="#ffe7a8"/>
<g fill="none" stroke-width="1.3" stroke-linejoin="round">{''.join(paths)}</g>
<line x1="{left}" y1="{top + plot_height}" x2="{width - right}"
y2="{top + plot_height}" stroke="#87949d"/>
<line x1="{center_x:.1f}" y1="{top}" x2="{center_x:.1f}" y2="380"
stroke="#b42318" stroke-width="1.5"/>
<text x="{left}" y="25" font-family="system-ui" font-size="15">{title}</text>
{''.join(annotations)}
<text x="{left}" y="405" font-family="system-ui" font-size="12" fill="#5c6770">
Reference-oriented bases; numbers are instrument Phred scores</text></svg>"""
    path.write_text(document, encoding="utf-8")
    return True


def write_qc_dashboard_svg(path: Path, rows: list[dict]) -> bool:
    if not rows:
        return False
    width = 1200
    row_height = 28
    height = 75 + row_height * len(rows)
    label_width, plot_width = 285, 680
    bars = []
    for index, row in enumerate(rows):
        y = 50 + index * row_height
        fraction = float(row.get("q20_fraction") or 0)
        color = "#16883e" if row.get("status") == "PASS" else "#b42318"
        label = html.escape(f"{row.get('sample')} {row.get('direction')}")
        bars.append(
            f'<text x="15" y="{y + 14}" font-family="system-ui" font-size="12">{label}</text>'
            f'<rect x="{label_width}" y="{y}" width="{plot_width}" height="18" fill="#e5e9ec"/>'
            f'<rect x="{label_width}" y="{y}" width="{plot_width * fraction:.1f}" '
            f'height="18" fill="{color}"/>'
            f'<text x="{label_width + plot_width + 12}" y="{y + 14}" '
            f'font-family="monospace" font-size="12">Q20 {fraction:.1%}; '
            f'n={row.get("trimmed_length", 0)}</text>'
        )
    document = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}"
role="img"><title>Batch read QC dashboard</title>
<rect width="{width}" height="{height}" fill="white"/>
<text x="15" y="25" font-family="system-ui" font-size="16">Batch read QC — Q20 fraction</text>
{''.join(bars)}</svg>"""
    path.write_text(document, encoding="utf-8")
    return True


def write_vcf(
    path: Path, reference_name: str, reference_length: int, variants: Iterable[Variant]
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.3\n")
        handle.write("##source=CMHS_SangerFlow_Pipeline\n")
        handle.write(f"##contig=<ID={reference_name},length={reference_length}>\n")
        handle.write('##FILTER=<ID=LowQual,Description="Consensus quality below threshold">\n')
        handle.write('##FILTER=<ID=Ambiguous,Description="Ambiguous IUPAC consensus call">\n')
        handle.write('##ALT=<ID=AMBIG,Description="Ambiguous consensus base">\n')
        handle.write('##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample identifier">\n')
        handle.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">\n')
        handle.write('##INFO=<ID=NOTE,Number=1,Type=String,Description="Additional evidence">\n')
        handle.write('##INFO=<ID=STRANDS,Number=1,Type=String,Description="Supporting reads">\n')
        handle.write('##INFO=<ID=FQ,Number=1,Type=Integer,Description="Forward Phred score">\n')
        handle.write('##INFO=<ID=RQ,Number=1,Type=Integer,Description="Reverse Phred score">\n')
        handle.write(
            '##INFO=<ID=FPR,Number=1,Type=Float,Description="Forward secondary peak ratio">\n'
        )
        handle.write(
            '##INFO=<ID=RPR,Number=1,Type=Float,Description="Reverse secondary peak ratio">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for variant in sorted(variants, key=lambda value: (value.position, value.sample)):
            quality = f"{variant.quality:.1f}"
            note = variant.note.replace(";", ",").replace(" ", "_") or "."
            strands = variant.strand_support or "."
            forward_q = variant.forward_quality if variant.forward_quality is not None else "."
            reverse_q = variant.reverse_quality if variant.reverse_quality is not None else "."
            forward_ratio = (
                variant.forward_peak_ratio if variant.forward_peak_ratio is not None else "."
            )
            reverse_ratio = (
                variant.reverse_peak_ratio if variant.reverse_peak_ratio is not None else "."
            )
            info = (
                f"SAMPLE={safe_name(variant.sample)};TYPE={variant.kind};NOTE={note};"
                f"STRANDS={strands};FQ={forward_q};RQ={reverse_q};"
                f"FPR={forward_ratio};RPR={reverse_ratio}"
            )
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
    quality_files: list[dict] | None = None,
    evidence_files: list[dict] | None = None,
    dashboard_file: str | None = None,
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
    quality_cards = (
        "".join(
            "<article><h3>"
            f"{html.escape(str(item['sample']))} — {html.escape(str(item['direction']))}"
            "</h3>"
            f"<img loading='lazy' src='../traces/{html.escape(str(item['file']))}' "
            f"alt='{html.escape(str(item['sample']))} quality profile'></article>"
            for item in (quality_files or [])
        )
        or "<p>No quality profiles were generated.</p>"
    )
    evidence_cards = (
        "".join(
            f"<article><h3>{html.escape(str(item['title']))}</h3>"
            f"<img loading='lazy' src='../traces/{html.escape(str(item['file']))}' "
            f"alt='{html.escape(str(item['title']))}'></article>"
            for item in (evidence_files or [])
        )
        or "<p>No variant-centred figures were generated.</p>"
    )
    dashboard = (
        f"<img src='../traces/{html.escape(dashboard_file)}' alt='Batch read QC dashboard'>"
        if dashboard_file
        else "<p>No batch dashboard was generated.</p>"
    )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>CMHS SangerFlow Pipeline report</title>
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
</style></head><body><main><p>Human Genomics Solutions</p>
<h1>CMHS SangerFlow Pipeline analysis report</h1>
<div class="card notice"><strong>Research use:</strong> calls require chromatogram review and
independent validation before clinical use.</div>
<section class="card"><h2>Samples</h2>{table(samples)}</section>
<section class="card"><h2>Variants</h2>{table(variants)}</section>
<section class="card"><h2>Batch QC</h2>{dashboard}</section>
<section class="card"><h2>Variant evidence</h2>{evidence_cards}</section>
<section class="card"><h2>Base-quality profiles</h2>{quality_cards}</section>
<section class="card"><h2>Chromatograms</h2>{trace_cards}</section>
<section class="card"><h2>Run metadata</h2><pre>{payload}</pre></section>
</main></body></html>"""
    path.write_text(document, encoding="utf-8")
