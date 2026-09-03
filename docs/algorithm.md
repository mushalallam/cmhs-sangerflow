# Algorithm and interpretation guide

## Inputs

Each sample may have one forward read, one reverse read, or both. A single-record DNA FASTA is
required as the coordinate reference. A sample sheet is recommended because filename inference
cannot establish biological identity.

## ABI parsing

Biopython reads the called sequence and `PCON` Phred values. When present, SangerFlow also reads
`PLOC`, `FWO_1`, and `DATA9` through `DATA12` to report the strongest and second-strongest dye
signals at each called peak. Missing raw channels do not invalidate instrument base calls, but
the report will contain no peak-level evidence for that file.

## Trimming and masking

SangerFlow uses the Mott maximum-scoring segment. For each Phred quality `Q`, the base score is:

```text
error_cutoff - 10^(-Q/10)
```

The default error cutoff is `0.05`. Within the retained segment, bases below Q20 are converted to
`N`. Reads shorter than 80 bases after trimming fail by default. All thresholds are recorded in
`run_metadata.json`.

## Reverse reads and paired consensus

Reverse reads and their quality arrays are reverse-complemented together. Forward and reverse
reads are globally aligned. Matching calls retain the higher quality. At a disagreement:

1. `N` yields to an informative base.
2. A quality difference of at least eight selects the higher-quality call.
3. Otherwise, compatible calls produce an IUPAC ambiguity code.

Terminal read overhangs are retained. The pair alignment is written for inspection.
By default, a pair fails QC unless it has at least 30 bases of direct overlap with 80% exact
identity. This prevents unrelated or incorrectly assigned reads from silently producing a
consensus.

## Reference alignment and variants

The consensus is globally aligned to the supplied reference with reduced terminal-gap penalties.
Terminal gaps are removed before variant extraction because reference sequence outside the
covered amplicon is noncoverage. Internal mismatches and gaps generate SNV, insertion, or deletion
records. Positions are one-based reference coordinates.

The aligned interval, identity, and fraction of the consensus covered by the reference alignment
are reported per sample. Both identity and consensus coverage must be at least 80% by default.
These are screening safeguards, not evidence that a particular reference accession or sample
assignment is biologically correct.

Ambiguous IUPAC disagreements are represented as `<AMBIG>` in VCF and retain the possible bases in
the `NOTE` field. Calls below the configured minimum consensus quality receive `LowQual`. Indels
are anchored for VCF output but are not currently left-normalized across repetitive sequence.

The combined VCF is intentionally site-like: sample identity is stored in the `SAMPLE` INFO field.
The TSV is the canonical complete call table for multi-sample runs.

## Interpretation limitations

- Automatic filename pairing validates uniqueness, not patient identity.
- Peak ratios are instrument- and chemistry-dependent and require laboratory validation.
- Mixed peaks may represent heterozygosity, mosaicism, contamination, artifacts, or noise.
- The pipeline does not phase variants or assign diploid genotypes.
- Large structural changes and complex repeat-associated indels are outside its intended scope.
- Alignment to an incorrect transcript, strand, build, or paralog produces misleading calls.
- Calls near primer sites or low-quality ends require particular caution.
- Software QC does not replace inspection of electropherograms.

## Validation before clinical use

A laboratory should establish truth sets covering wild type, heterozygous SNVs, homozygous SNVs,
small indels, difficult repeats, low-signal traces, contamination, and failed reactions. Thresholds,
reference accessions, acceptable coverage, manual-review criteria, sensitivity, specificity, and
version-control procedures must be documented independently of this repository.
