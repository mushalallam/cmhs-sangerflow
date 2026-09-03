"""Command-line interface for SangerFlow."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .abi import read_abi
from .pairing import auto_pair, read_sample_sheet
from .pipeline import RunConfig, run_pipeline
from .quality import quality_summary, trim_and_mask


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sangerflow",
        description="Auditable Sanger chromatogram QC, consensus, and variant analysis",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser("run", help="analyze one or more Sanger samples")
    run.add_argument("--input", type=Path, required=True, help="directory containing AB1 files")
    run.add_argument("--reference", type=Path, required=True, help="single-record reference FASTA")
    run.add_argument("--output", type=Path, required=True, help="new or empty output directory")
    run.add_argument("--sample-sheet", type=Path, help="CSV with sample,forward,reverse columns")
    run.add_argument(
        "--allow-single", action="store_true", help="allow forward-only or reverse-only samples"
    )
    run.add_argument("--error-cutoff", type=float, default=0.05, help="Mott trimming error cutoff")
    run.add_argument("--mask-below", type=int, default=20, help="mask bases below this Phred score")
    run.add_argument("--min-read-length", type=int, default=80)
    run.add_argument(
        "--quality-delta", type=int, default=8, help="quality difference to resolve read conflicts"
    )
    run.add_argument("--min-variant-quality", type=int, default=20)
    run.add_argument("--min-overlap", type=int, default=30)
    run.add_argument("--min-overlap-identity", type=float, default=0.80)
    run.add_argument("--min-reference-identity", type=float, default=0.80)
    run.add_argument("--min-reference-coverage", type=float, default=0.80)
    run.add_argument(
        "--call-mixed-peaks",
        action="store_true",
        help="experimental: convert strong secondary peaks to IUPAC calls",
    )
    run.add_argument("--mixed-peak-ratio", type=float, default=0.33)
    run.add_argument("--mixed-peak-min-signal", type=int, default=100)

    inspect = subparsers.add_parser("inspect", help="inspect one ABI file without running analysis")
    inspect.add_argument("ab1", type=Path)
    inspect.add_argument("--error-cutoff", type=float, default=0.05)
    inspect.add_argument("--mask-below", type=int, default=20)
    inspect.add_argument("--min-read-length", type=int, default=80)
    return parser


def _run(args: argparse.Namespace) -> int:
    if args.sample_sheet:
        samples = read_sample_sheet(args.sample_sheet, args.input)
        if not args.allow_single:
            incomplete = [
                sample.sample for sample in samples if not sample.forward or not sample.reverse
            ]
            if incomplete:
                raise ValueError(
                    f"Single-read samples require --allow-single: {', '.join(incomplete)}"
                )
    else:
        samples = auto_pair(args.input, allow_single=args.allow_single)
    config = RunConfig(
        error_cutoff=args.error_cutoff,
        mask_below=args.mask_below,
        min_read_length=args.min_read_length,
        quality_delta=args.quality_delta,
        min_variant_quality=args.min_variant_quality,
        min_overlap=args.min_overlap,
        min_overlap_identity=args.min_overlap_identity,
        min_reference_identity=args.min_reference_identity,
        min_reference_coverage=args.min_reference_coverage,
        call_mixed_peaks=args.call_mixed_peaks,
        mixed_peak_ratio=args.mixed_peak_ratio,
        mixed_peak_min_signal=args.mixed_peak_min_signal,
    )
    metadata = run_pipeline(samples, args.reference, args.output, config)
    print(json.dumps(metadata, indent=2, sort_keys=True))
    return 2 if metadata["failed_samples"] else 0


def _inspect(args: argparse.Namespace) -> int:
    read = read_abi(args.ab1)
    trim_and_mask(
        read,
        error_cutoff=args.error_cutoff,
        mask_below=args.mask_below,
        min_length=args.min_read_length,
    )
    result = quality_summary(read)
    result["peak_evidence_count"] = len(read.peak_evidence)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if read.trimmed_sequence else 2


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return _run(args) if args.command == "run" else _inspect(args)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
