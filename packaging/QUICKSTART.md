# CMHS SangerFlow Pipeline — standalone quick start

This package contains a standalone command. Python, Conda, and Biopython do not need to be
installed separately.

## 1. Check the download

Compare the archive's SHA-256 checksum with its matching `.sha256` file on the GitHub release.
The digest GitHub displays beside the `.sha256` sidecar is the digest of that small sidecar—not
the digest expected for the application archive.

- macOS: `shasum -a 256 <archive>`
- Linux: `sha256sum <archive>`
- Windows PowerShell: `Get-FileHash <archive> -Algorithm SHA256`

## 2. Open the graphical interface

- macOS: double-click `Start CMHS SangerFlow.command`
- Windows: double-click `Start CMHS SangerFlow.bat`
- Linux: run `./Start-CMHS-SangerFlow.sh`

The interface opens in your browser but runs only on this computer. Chromatograms are not
uploaded to an internet service.

## 3. Verify the program

- macOS or Linux: open a terminal in this folder and run `./sangerflow doctor`
- Windows: open PowerShell in this folder and run `.\sangerflow.exe doctor`

The result should contain `"status": "PASS"` and `"distribution": "standalone"`.

## 4. Run an analysis

Use `./sangerflow` on macOS/Linux or `.\sangerflow.exe` on Windows:

```text
sangerflow run --input chromatograms --sample-sheet samples.csv --reference reference.fasta --output results
```

Run `sangerflow run --help` to see every option. See the project README for sample-sheet format,
output descriptions, safety notes, and troubleshooting.

## First-launch security messages

These community builds are not yet code-signed. Only continue after checking the SHA-256 value.

- macOS: if Gatekeeper blocks the program, use **System Settings → Privacy & Security → Open
  Anyway**, then try again. An unsigned first launch can also appear frozen for more than a
  minute while macOS performs online checks. After verifying the downloaded ZIP against its
  published SHA-256 value, stop the stalled launch with **Control-C** and run
  `xattr -dr com.apple.quarantine "/path/to/extracted/CMHS-SangerFlow-folder"` once.
- Windows: Microsoft Defender SmartScreen may show an unrecognized-app warning. Select **More
  info → Run anyway** only after verifying the checksum and GitHub source.
- Linux: if extraction removed executable permissions, run `chmod +x sangerflow
  Start-CMHS-SangerFlow.sh`. The current Linux build requires x86-64 and glibc 2.35 or newer;
  ARM64 and musl/Alpine are not currently packaged.

These workarounds are temporary. Warning-free macOS distribution requires Apple Developer ID
signing and notarization. Windows trust requires signing, and Microsoft Store distribution is
the most reliable way to avoid SmartScreen download warnings.

CMHS SangerFlow Pipeline is research-use software and is not a medical device.
