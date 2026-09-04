# Changelog

## 0.3.0 - 2026-09-04

- Added a local-only graphical interface with native file selection and editable read pairing.
- Added standard paired-read, single-read, mixed-peak, and advanced presets.
- Added live per-sample progress, plain-language completion messages, and report/folder actions.
- Added automatic timestamped result directories and a patient-free synthetic demonstration.
- Added double-click launchers to the standalone macOS and Windows downloads.
- Added a Linux graphical-interface launcher and GUI checks to every native release build.

## 0.2.1 - 2026-09-03

- Added self-contained native downloads for Intel/Apple Silicon macOS, Windows, and Linux.
- Added `sangerflow doctor`, SHA-256 release checksums, and automated native smoke tests.

## 0.2.0 - 2026-09-03

- Adopted the Human Genomics Solutions / CMHS SangerFlow Pipeline branding.
- Added automatic reference-based orientation verification and correction tracking.
- Added variant-centred, reference-oriented chromatogram evidence figures.
- Added per-read Phred profiles and a batch Q20 dashboard.
- Added continuous Q20 length, GC fraction, and primary/secondary signal metrics.
- Added forward/reverse quality, peak-ratio, and strand-support evidence to TSV and VCF.

## 0.1.0 - 2026-09-03

- Initial research release.
- ABI base-call, quality, peak-location, and dye-channel parsing.
- Mott trimming and configurable low-quality masking.
- Conservative sample-sheet and filename pairing.
- Quality-aware paired-read consensus with IUPAC ambiguity calls.
- Reference alignment and SNV/indel reporting in TSV and VCF.
- Reproducibility metadata, SHA-256 input manifest, and human-readable HTML report.
- Four-channel chromatogram SVG rendering and peak-evidence tables.
- Paired-read overlap and reference identity/coverage safety gates.
- Aggregate multi-sample consensus FASTA output.
