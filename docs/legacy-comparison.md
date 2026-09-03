# Relationship to the legacy ASAP workflow

SangerFlow was developed as an independent, modern implementation after evaluating laboratory
workflows built around the GPLv3-licensed Automated Sanger Analysis Pipeline (ASAP). It does not
contain, bundle, or invoke ASAP code. If ASAP contributed to a scientific workflow, cite the
original publication:

> Singh A, Bhatia P. Automated Sanger Analysis Pipeline (ASAP): A Tool for Rapidly Analyzing
> Sanger Sequencing Data with Minimum User Interference. *J Biomol Tech.* 2016. PMID: 27790076.

## Review of the local 2021 workflow

The reviewed script provided useful forward/reverse trimming, consensus generation, batch
alignment, and optional exon translation. Several implementation details make it unsuitable as a
reproducible unattended pipeline today:

- tool paths are hard-coded to one computer;
- filenames are concatenated into shell commands without safe argument handling;
- intermediate files are deleted through shell commands;
- samples are paired by their position in two filename lists rather than an explicit manifest;
- wildcard collection can include consensus files left by an earlier run;
- all optional external tools are required during startup, even when a mode does not use them;
- the 2021 batch script advertises several modes but implements only paired `FR` processing;
- output lacks input checksums, parameter capture, structured QC metrics, and machine-readable
  variant calls.

## Migration coverage

SangerFlow replaces the core chromatogram workflow with ABI parsing, quality trimming,
forward/reverse consensus, reference alignment, chromatogram plots, batch aggregation, variant
tables, VCF, and reproducibility metadata. It has no dependency on SEQTK, EMBOSS, ClustalW,
BLAST+, or SeaView.

SangerFlow 0.1 is not a command-line-compatible replacement for ASAP. In particular, automatic
exon extraction, six-frame translation, and protein alignment are intentionally not implemented.
A future protein-consequence feature should use an explicitly versioned annotated reference and
declared transcript rather than selecting a frame heuristically.
