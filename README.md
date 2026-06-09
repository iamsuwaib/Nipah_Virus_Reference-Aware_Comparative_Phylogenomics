# Nipah Virus Reference-Aware Comparative Genomics

This repository contains reproducible analysis code, accession-level inputs, curated metadata, and non-sensitive derived outputs for the study:

**Reference-aware phylogenomics reveals country-lineage structure and glycoprotein variation in public Nipah virus genomes**

## Overview

The study analyzes public complete and near-complete Nipah virus (NiV) genomes using a reference-aware comparative framework. The workflow integrates whole-genome phylogenetic reconstruction, country/genotype-aware SNP profiling, recurrent nonsynonymous variant analysis, statistical testing, and F/G glycoprotein substitution mapping.

The repository is organized to support reproducibility from public GenBank accessions while excluding manuscript drafts, large exploratory intermediates, and private local files.

## Reference-Aware Design

Public NiV genomes are unevenly distributed across countries, hosts, collection years, and outbreaks. A single prototype-reference comparison can misclassify lineage-level divergence as country-specific variation. This workflow therefore applies country/genotype-aware comparisons:

- Bangladesh sequences: `AY988601.1`
- India sequences: `MH396625.1`
- Malaysia/prototype-related sequences: `AF212302.2`
- Hendra virus `NC_001906.3`: rooted phylogenetic outgroup only

This design keeps phylogenetic interpretation, SNP profiling, and F/G glycoprotein analysis aligned with known NiV lineage structure.

## Repository Layout

| Directory | Description |
|---|---|
| `data/accession_lists/` | Accession lists and accession-level context used in the analysis |
| `data/metadata/` | Curated metadata for country, host, year, and lineage annotation |
| `data/reference_manifest/` | Reference and outgroup sequence manifest |
| `scripts/` | Analysis code for the reference-aware comparative workflow |
| `results/summary_tables/` | Non-sensitive derived summary tables |
| `results/figures/` | Derived visual outputs used to support the analysis |
| `results/itol_annotations/` | Tree annotation outputs for phylogenetic visualization |
| `docs/` | Repository notes and file manifest |

## Reproducibility

Create a fresh Python environment and install the required packages:

```bash
pip install -r requirements.txt
```

The Python scripts reproduce the country/reference-aware variant, protein, and statistical analyses from the accession-level inputs and curated metadata provided here. The full phylogenetic component also requires MAFFT and IQ-TREE 2.

Some outputs depend on external alignment and tree reconstruction steps; the provided accession lists, metadata, reference manifest, scripts, and derived summary tables document the analysis state used for the associated manuscript.

## Interpretation Notes

- `AF212302.2` is not treated as a universal reference for all countries.
- `NC_001906.3` is used only as the Hendra virus outgroup for rooted phylogenetic visualization.
- AY988601.1 F-protein `S207L` and `G252D` are treated as known reference-correction positions and excluded from biological interpretation.
- Prototype-relative F/G contrasts are used for broad lineage context, not as country-specific mutation calls when a country/genotype reference is available.

## Supplementary Materials

The master supplementary Excel workbook and supplementary file index are submitted with the manuscript as journal supplementary material. This repository provides the code and non-sensitive derived outputs needed to support reproducible analysis from public accession-level inputs.

## Citation

Please cite the associated manuscript when available:

Suwaib et al. **Reference-aware phylogenomics reveals country-lineage structure and glycoprotein variation in public Nipah virus genomes.**

## License

Code in this repository is released under the MIT License unless otherwise stated.
