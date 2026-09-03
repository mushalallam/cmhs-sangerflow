# CMHS SangerFlow Pipeline — standalone quick start

This package contains a standalone command. Python, Conda, and Biopython do not need to be
installed separately.

## 1. Check the download

Compare the archive's SHA-256 checksum with its matching `.sha256` file on the GitHub release.

- macOS: `shasum -a 256 <archive>`
- Linux: `sha256sum <archive>`
- Windows PowerShell: `Get-FileHash <archive> -Algorithm SHA256`

## 2. Verify the program

- macOS or Linux: open a terminal in this folder and run `./sangerflow doctor`
- Windows: open PowerShell in this folder and run `.\sangerflow.exe doctor`

The result should contain `"status": "PASS"` and `"distribution": "standalone"`.

## 3. Run an analysis

Use `./sangerflow` on macOS/Linux or `.\sangerflow.exe` on Windows:

```text
sangerflow run --input chromatograms --sample-sheet samples.csv --reference reference.fasta --output results
```

Run `sangerflow run --help` to see every option. See the project README for sample-sheet format,
output descriptions, safety notes, and troubleshooting.

## First-launch security messages

These community builds are not yet code-signed. Only continue after checking the SHA-256 value.

- macOS: if Gatekeeper blocks the program, use **System Settings → Privacy & Security → Open
  Anyway**, then try again.
- Windows: Microsoft Defender SmartScreen may show an unrecognized-app warning. Select **More
  info → Run anyway** only after verifying the checksum and GitHub source.

CMHS SangerFlow Pipeline is research-use software and is not a medical device.
