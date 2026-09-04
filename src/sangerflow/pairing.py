"""Explicit sample-sheet loading and conservative filename-based pairing."""

from __future__ import annotations

import csv
import re
from pathlib import Path

from .models import SampleInput

PAIR_PATTERNS = (
    re.compile(
        r"^(?P<sample>.+?)[_.-](?P<direction>forward|reverse|f|r)(?:[_.-].*)?\.ab(?:1|i)$",
        re.IGNORECASE,
    ),
    re.compile(r"^(?P<sample>.+?)(?P<direction>F|R)(?:[_.-].*)?\.ab(?:1|i)$"),
)


def _direction(value: str) -> str:
    return "forward" if value.lower() in {"f", "forward"} else "reverse"


def infer_sample_and_direction(path: Path) -> tuple[str, str]:
    for pattern in PAIR_PATTERNS:
        match = pattern.match(path.name)
        if match:
            sample = match.group("sample").rstrip("_.-")
            if sample:
                return sample, _direction(match.group("direction"))
    raise ValueError(
        f"Cannot infer sample/direction from {path.name!r}; use a sample sheet with "
        "sample,forward,reverse columns"
    )


def auto_pair(input_dir: Path, *, allow_single: bool = False) -> list[SampleInput]:
    input_dir = Path(input_dir).resolve()
    files = sorted(
        path for path in input_dir.iterdir() if path.is_file() and path.suffix.lower() == ".ab1"
    )
    if not files:
        raise ValueError(f"No .ab1 files found directly inside {input_dir}")
    grouped: dict[str, dict[str, Path]] = {}
    for path in files:
        sample, direction = infer_sample_and_direction(path)
        slot = grouped.setdefault(sample, {})
        if direction in slot:
            raise ValueError(
                f"Multiple {direction} reads inferred for sample {sample!r}: "
                f"{slot[direction].name}, {path.name}. Use a sample sheet."
            )
        slot[direction] = path.resolve()

    samples: list[SampleInput] = []
    for sample, reads in sorted(grouped.items()):
        if not allow_single and set(reads) != {"forward", "reverse"}:
            missing = "reverse" if "forward" in reads else "forward"
            raise ValueError(f"Sample {sample!r} is missing a {missing} read")
        samples.append(SampleInput(sample, reads.get("forward"), reads.get("reverse")))
    return samples


def read_sample_sheet(path: Path, input_dir: Path | None = None) -> list[SampleInput]:
    path = Path(path).resolve()
    base = Path(input_dir).resolve() if input_dir else path.parent
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"sample", "forward", "reverse"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Sample sheet must contain sample,forward,reverse columns")
        samples: list[SampleInput] = []
        seen: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            sample = (row.get("sample") or "").strip()
            if not sample:
                raise ValueError(f"Missing sample name on sample-sheet line {line_number}")
            if sample in seen:
                raise ValueError(f"Duplicate sample {sample!r} on line {line_number}")
            seen.add(sample)

            def resolve(value: str | None, *, row_number: int = line_number) -> Path | None:
                value = (value or "").strip()
                if not value:
                    return None
                candidate = Path(value)
                result = candidate if candidate.is_absolute() else base / candidate
                if not result.is_file():
                    raise ValueError(f"Input file not found on line {row_number}: {result}")
                return result.resolve()

            forward = resolve(row.get("forward"))
            reverse = resolve(row.get("reverse"))
            if not forward and not reverse:
                raise ValueError(f"Sample {sample!r} has no reads on line {line_number}")
            samples.append(SampleInput(sample, forward, reverse))
    if not samples:
        raise ValueError("Sample sheet contains no samples")
    return samples
