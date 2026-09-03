"""End-to-end orchestration for CMHS SangerFlow Pipeline analysis runs."""

from __future__ import annotations

import hashlib
import json
import platform
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import Bio
from Bio import SeqIO
from Bio.Seq import Seq

from . import __version__
from .abi import read_abi
from .alignment import build_consensus, reverse_complement_read
from .evidence import enrich_variant
from .models import ReadData, SampleInput, Variant
from .output import (
    safe_name,
    write_alignment,
    write_fasta,
    write_fastq,
    write_html_report,
    write_multi_fasta,
    write_pair_alignment,
    write_peak_table,
    write_qc_dashboard_svg,
    write_quality_svg,
    write_trace_svg,
    write_trace_window_svg,
    write_tsv,
    write_vcf,
)
from .quality import quality_summary, trim_and_mask
from .variants import analyze_reference


@dataclass(frozen=True, slots=True)
class RunConfig:
    error_cutoff: float = 0.05
    mask_below: int = 20
    min_read_length: int = 80
    quality_delta: int = 8
    min_variant_quality: int = 20
    min_overlap: int = 30
    min_overlap_identity: float = 0.80
    min_reference_identity: float = 0.80
    min_reference_coverage: float = 0.80
    call_mixed_peaks: bool = False
    mixed_peak_ratio: float = 0.33
    mixed_peak_min_signal: int = 100
    auto_orient: bool = True
    orientation_delta: float = 0.05

    def validate(self) -> None:
        if not 0 < self.error_cutoff < 1:
            raise ValueError("error_cutoff must be between 0 and 1")
        if self.mask_below < 0 or self.min_variant_quality < 0:
            raise ValueError("quality thresholds cannot be negative")
        if self.min_read_length < 1:
            raise ValueError("min_read_length must be positive")
        if self.quality_delta < 0:
            raise ValueError("quality_delta cannot be negative")
        if self.min_overlap < 1:
            raise ValueError("min_overlap must be positive")
        if not 0 <= self.min_overlap_identity <= 1:
            raise ValueError("min_overlap_identity must be between 0 and 1")
        if not 0 <= self.min_reference_identity <= 1:
            raise ValueError("min_reference_identity must be between 0 and 1")
        if not 0 <= self.min_reference_coverage <= 1:
            raise ValueError("min_reference_coverage must be between 0 and 1")
        if not 0 <= self.mixed_peak_ratio <= 1:
            raise ValueError("mixed_peak_ratio must be between 0 and 1")
        if self.mixed_peak_min_signal < 0:
            raise ValueError("mixed_peak_min_signal cannot be negative")
        if not 0 <= self.orientation_delta <= 1:
            raise ValueError("orientation_delta must be between 0 and 1")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_reference(path: Path) -> tuple[str, str]:
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise ValueError(f"Reference FASTA must contain exactly one sequence; found {len(records)}")
    sequence = str(records[0].seq).upper().replace("U", "T")
    if not sequence or set(sequence) - set("ACGTN"):
        raise ValueError("Reference must be a non-empty DNA sequence containing only A,C,G,T,N")
    return records[0].id, sequence


def _apply_mixed_peak_calls(read: ReadData, config: RunConfig) -> None:
    if not config.call_mixed_peaks:
        return
    sequence = list(read.sequence)
    for peak in read.peak_evidence:
        if (
            peak.primary_signal >= config.mixed_peak_min_signal
            and peak.secondary_ratio >= config.mixed_peak_ratio
            and peak.quality >= config.mask_below
            and peak.called_base == peak.primary_base
        ):
            sequence[peak.index] = peak.suggested_iupac
    read.sequence = "".join(sequence)


def _prepare_read(path: Path, config: RunConfig, *, reverse: bool = False) -> ReadData:
    read = read_abi(path)
    _apply_mixed_peak_calls(read, config)
    trim_and_mask(
        read,
        error_cutoff=config.error_cutoff,
        mask_below=config.mask_below,
        min_length=config.min_read_length,
    )
    if reverse and read.trimmed_sequence:
        reverse_complement_read(read)
    return read


def _reference_fit(reference: str, sequence: str) -> float:
    try:
        _, alignment = analyze_reference(
            "orientation", reference, sequence, [40] * len(sequence), min_quality=0
        )
    except ValueError:
        return 0.0
    return alignment.identity * alignment.consensus_coverage


def _auto_orient(read: ReadData, reference: str, config: RunConfig) -> None:
    if not config.auto_orient or not read.trimmed_sequence:
        return
    current_score = _reference_fit(reference, read.trimmed_sequence)
    opposite = str(Seq(read.trimmed_sequence).reverse_complement())
    opposite_score = _reference_fit(reference, opposite)
    if opposite_score > current_score + config.orientation_delta:
        reverse_complement_read(read)
        read.orientation_corrected = True


def run_pipeline(
    samples: list[SampleInput],
    reference_path: Path,
    output_dir: Path,
    config: RunConfig,
) -> dict:
    config.validate()
    reference_path = Path(reference_path).resolve()
    output_dir = Path(output_dir).resolve()
    reference_name, reference = _load_reference(reference_path)
    output_names = [safe_name(item.sample) for item in samples]
    if len(output_names) != len(set(output_names)):
        raise ValueError("Sample names are not unique after filesystem-safe normalization")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError(f"Output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    reads_dir = output_dir / "reads"
    consensus_dir = output_dir / "consensus"
    alignments_dir = output_dir / "alignments"
    traces_dir = output_dir / "traces"
    reports_dir = output_dir / "reports"
    for directory in (reads_dir, consensus_dir, alignments_dir, traces_dir, reports_dir):
        directory.mkdir()

    read_rows: list[dict] = []
    sample_rows: list[dict] = []
    variants: list[Variant] = []
    all_reads: list[ReadData] = []
    trace_files: list[dict] = []
    quality_files: list[dict] = []
    evidence_files: list[dict] = []
    consensus_records: list[tuple[str, str]] = []
    manifest = {str(reference_path): _sha256(reference_path)}

    for item in samples:
        sample = safe_name(item.sample)
        row = {
            "sample": item.sample,
            "forward": str(item.forward or ""),
            "reverse": str(item.reverse or ""),
            "status": "FAIL",
            "consensus_length": 0,
            "conflicts": 0,
            "ambiguous_bases": 0,
            "overlap_bases": "",
            "overlap_identity": "",
            "reference_start": "",
            "reference_end": "",
            "reference_identity": "",
            "reference_coverage": "",
            "variants": 0,
            "message": "",
        }
        try:
            forward = _prepare_read(item.forward, config) if item.forward else None
            reverse = _prepare_read(item.reverse, config, reverse=True) if item.reverse else None
            if forward:
                _auto_orient(forward, reference, config)
            if reverse:
                _auto_orient(reverse, reference, config)
            for direction, read in (("forward", forward), ("reverse", reverse)):
                if read is None:
                    continue
                all_reads.append(read)
                manifest[str(read.path)] = _sha256(read.path)
                summary = quality_summary(read)
                summary.update({"sample": item.sample, "direction": direction})
                read_rows.append(summary)
                if read.trimmed_sequence:
                    write_fastq(
                        reads_dir / f"{sample}.{direction}.trimmed.fastq",
                        f"{item.sample}_{direction}",
                        read.trimmed_sequence,
                        read.trimmed_qualities,
                    )
                trace_filename = f"{sample}.{direction}.svg"
                if write_trace_svg(traces_dir / trace_filename, read):
                    trace_files.append(
                        {"sample": item.sample, "direction": direction, "file": trace_filename}
                    )
                quality_filename = f"{sample}.{direction}.quality.svg"
                if write_quality_svg(traces_dir / quality_filename, read):
                    quality_files.append(
                        {"sample": item.sample, "direction": direction, "file": quality_filename}
                    )
            consensus = build_consensus(forward, reverse, quality_delta=config.quality_delta)
            row.update(
                {
                    "consensus_length": len(consensus.sequence),
                    "conflicts": consensus.conflicts,
                    "ambiguous_bases": consensus.ambiguous_bases,
                }
            )
            if forward and reverse:
                row["overlap_bases"] = consensus.overlap_bases
                row["overlap_identity"] = f"{consensus.overlap_identity or 0:.4f}"
                if consensus.overlap_bases < config.min_overlap:
                    raise ValueError(
                        f"paired-read overlap {consensus.overlap_bases} is below "
                        f"minimum {config.min_overlap}"
                    )
                if (consensus.overlap_identity or 0) < config.min_overlap_identity:
                    raise ValueError(
                        f"paired-read overlap identity {consensus.overlap_identity:.3f} is below "
                        f"minimum {config.min_overlap_identity:.3f}"
                    )
            write_fasta(
                consensus_dir / f"{sample}.consensus.fasta", item.sample, consensus.sequence
            )
            consensus_records.append((item.sample, consensus.sequence))
            write_pair_alignment(alignments_dir / f"{sample}.reads.fasta", item.sample, consensus)
            sample_variants, reference_alignment = analyze_reference(
                item.sample,
                reference,
                consensus.sequence,
                consensus.qualities,
                min_quality=config.min_variant_quality,
            )
            write_alignment(
                alignments_dir / f"{sample}.reference.fasta",
                item.sample,
                reference_alignment.aligned_reference,
                reference_alignment.aligned_consensus,
            )
            row.update(
                {
                    "reference_start": reference_alignment.reference_start,
                    "reference_end": reference_alignment.reference_end,
                    "reference_identity": f"{reference_alignment.identity:.4f}",
                    "reference_coverage": f"{reference_alignment.consensus_coverage:.4f}",
                }
            )
            if reference_alignment.identity < config.min_reference_identity:
                raise ValueError(
                    f"reference identity {reference_alignment.identity:.3f} is below minimum "
                    f"{config.min_reference_identity:.3f}; verify sample pairing and reference"
                )
            if reference_alignment.consensus_coverage < config.min_reference_coverage:
                raise ValueError(
                    f"reference alignment covers {reference_alignment.consensus_coverage:.3f} "
                    f"of the consensus, below minimum {config.min_reference_coverage:.3f}; "
                    "verify sample pairing and reference"
                )
            enriched_variants: list[Variant] = []
            for variant_number, variant in enumerate(sample_variants, start=1):
                enriched, centers = enrich_variant(
                    variant, reference_alignment, consensus, forward, reverse
                )
                enriched_variants.append(enriched)
                for direction, center in centers.items():
                    read = forward if direction == "forward" else reverse
                    assert read is not None
                    evidence_filename = (
                        f"{sample}.variant-{variant_number}-{variant.position}."
                        f"{direction}.svg"
                    )
                    title = (
                        f"{item.sample} {variant.kind} {variant.reference}>"
                        f"{variant.alternate} at {reference_name}:{variant.position} "
                        f"({direction})"
                    )
                    if write_trace_window_svg(
                        traces_dir / evidence_filename,
                        read,
                        center,
                        title,
                        reverse=read.is_reverse_complemented,
                    ):
                        evidence_files.append({"title": title, "file": evidence_filename})
            variants.extend(enriched_variants)
            row.update(
                {
                    "status": "PASS",
                    "variants": len(enriched_variants),
                    "message": "",
                }
            )
        except Exception as exc:
            row["message"] = str(exc)
        sample_rows.append(row)

    write_multi_fasta(consensus_dir / "all_consensus.fasta", consensus_records)

    read_fields = [
        "sample",
        "direction",
        "read",
        "file",
        "raw_length",
        "trim_start",
        "trim_end",
        "trimmed_length",
        "mean_q_raw",
        "mean_q_trimmed",
        "continuous_q20_length",
        "q20_fraction",
        "q30_fraction",
        "n_fraction",
        "gc_fraction",
        "median_primary_signal",
        "median_secondary_signal",
        "median_secondary_ratio",
        "reference_orientation",
        "orientation_corrected",
        "status",
    ]
    sample_fields = list(sample_rows[0]) if sample_rows else []
    variant_rows = [asdict(variant) for variant in variants]
    variant_fields = [
        "sample",
        "position",
        "reference",
        "alternate",
        "kind",
        "quality",
        "filter",
        "note",
        "strand_support",
        "forward_quality",
        "reverse_quality",
        "forward_peak_ratio",
        "reverse_peak_ratio",
    ]
    write_tsv(reports_dir / "read_qc.tsv", read_rows, read_fields)
    write_tsv(reports_dir / "samples.tsv", sample_rows, sample_fields)
    write_tsv(reports_dir / "variants.tsv", variant_rows, variant_fields)
    write_peak_table(reports_dir / "peak_evidence.tsv", all_reads)
    write_vcf(reports_dir / "variants.vcf", reference_name, len(reference), variants)
    dashboard_filename = "batch_qc.svg"
    dashboard_created = write_qc_dashboard_svg(traces_dir / dashboard_filename, read_rows)

    metadata = {
        "sangerflow_version": __version__,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "biopython": Bio.__version__,
        "platform": platform.platform(),
        "reference": str(reference_path),
        "reference_name": reference_name,
        "reference_length": len(reference),
        "configuration": asdict(config),
        "sample_count": len(samples),
        "passed_samples": sum(row["status"] == "PASS" for row in sample_rows),
        "failed_samples": sum(row["status"] != "PASS" for row in sample_rows),
        "input_sha256": manifest,
    }
    (reports_dir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_html_report(
        reports_dir / "report.html",
        metadata,
        sample_rows,
        variant_rows,
        trace_files,
        quality_files,
        evidence_files,
        dashboard_filename if dashboard_created else None,
    )
    return metadata
