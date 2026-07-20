# Publication-analysis archive

Files in this directory are historical scripts used during exploratory or
publication analysis. They are not installed, are not part of the supported
StORF-Reporter CLI, and are excluded from release validation.

Some scripts retain hard-coded historical filesystem paths and assumptions
about intermediate CD-HIT files. They are preserved for provenance only and
must not be presented as portable toolkit components. Reproducible equivalents
should be implemented as versioned workflows under `workflows/` with explicit
inputs, environments, and expected outputs before reuse.

Supported analysis code lives under `src/StORF_Reporter`; supported workflow
entry points are documented in the repository README.
