"""Run release checks against a freshly built standalone executable."""

from __future__ import annotations

import argparse
import json
import subprocess


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary")
    parser.add_argument("--version", required=True)
    args = parser.parse_args()

    version = subprocess.run(
        [args.binary, "--version"], check=True, capture_output=True, text=True
    ).stdout.strip()
    if version != f"sangerflow {args.version}":
        raise SystemExit(f"Unexpected version output: {version!r}")

    completed = subprocess.run(
        [args.binary, "doctor"], check=True, capture_output=True, text=True
    )
    diagnosis = json.loads(completed.stdout)
    if diagnosis.get("status") != "PASS" or diagnosis.get("distribution") != "standalone":
        raise SystemExit(f"Standalone self-check failed: {diagnosis}")
    print(completed.stdout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
