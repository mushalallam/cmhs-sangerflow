# Changelog

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
