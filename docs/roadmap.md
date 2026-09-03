# Development roadmap

CMHS SangerFlow Pipeline prioritizes reviewable evidence and reproducibility over the number of
features. Planned work is grouped by the validation it requires.

## Next priorities

1. Primer-aware amplicon verification and optional primer removal.
2. Versioned GenBank/GFF and transcript input with coding and protein consequences.
3. A printable PDF review report with analyst and reviewer sign-off fields.
4. Reference-guided assembly for tiled amplicons with more than two reads.
5. A redistributable truth set and documented comparison against established Sanger tools.

## Requires dedicated validation

- Heterozygous-indel trace deconvolution.
- Diploid genotype assignment.
- HGVS normalization across repetitive regions and transcript versions.
- Any claim of diagnostic or clinical suitability.

These capabilities will not be inferred from consensus sequence alone. Where a mature external
engine is used, it will be an explicit optional dependency with its version, parameters, and
license recorded in run metadata.
