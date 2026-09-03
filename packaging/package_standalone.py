"""Create a standalone release archive and SHA-256 sidecar."""

from __future__ import annotations

import argparse
import hashlib
import tarfile
import zipfile
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--platform-label", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--format", choices=("zip", "tar.gz"), required=True)
    parser.add_argument("--output", type=Path, default=Path("release-assets"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.binary.is_file():
        raise SystemExit(f"Standalone binary not found: {args.binary}")

    args.output.mkdir(parents=True, exist_ok=True)
    stem = f"CMHS-SangerFlow-{args.version}-{args.platform_label}"
    suffix = ".zip" if args.format == "zip" else ".tar.gz"
    archive = args.output / f"{stem}{suffix}"
    quickstart = Path(__file__).with_name("QUICKSTART.md")

    if args.format == "zip":
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
            binary_info = zipfile.ZipInfo(f"{stem}/{args.binary.name}")
            binary_info.external_attr = 0o755 << 16
            bundle.writestr(binary_info, args.binary.read_bytes())
            bundle.write(quickstart, f"{stem}/QUICKSTART.md")
    else:
        with tarfile.open(archive, "w:gz") as bundle:
            bundle.add(args.binary, arcname=f"{stem}/{args.binary.name}")
            bundle.add(quickstart, arcname=f"{stem}/QUICKSTART.md")

    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    checksum = archive.with_name(f"{archive.name}.sha256")
    checksum.write_text(f"{digest}  {archive.name}\n", encoding="utf-8")
    print(archive)
    print(checksum)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
