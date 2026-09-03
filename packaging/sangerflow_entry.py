"""PyInstaller entry point for the standalone CMHS SangerFlow executable."""

from sangerflow.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
