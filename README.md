# snRNA-seq of the Rhesus macaque uterus — downstream analysis scripts

Analysis scripts accompanying the master's thesis *Single-nucleus transcriptomics
of the Rhesus macaque uterus: intrauterine infection and IL-1 blockade*
(Samuel Evan Cipolletti, University of Parma, 2026).

The study profiles uterine tissue from *Macaca mulatta* across intra-amniotic
saline controls, animals exposed to *E. coli* and treated with antibiotics, and
animals additionally treated with Anakinra, stratified by preterm labour (PTL)
outcome.

## Scope of this repository

This repository contains **only the downstream analyses that I carried out
myself**: pseudobulk differential expression, functional enrichment, and
cell-cell communication.

The upstream part of the pipeline — Cell Ranger, quality control, doublet
removal, normalisation, integration, clustering and cell-type annotation — was
performed by other members of the project and is **not** included here. The
scripts below start from the already annotated Seurat object.

## Scripts

| Script | What it does |
| --- | --- |
| `pseudobulk_unified.R` | Aggregates counts into one pseudobulk profile per sample, for the whole tissue and for each of the 17 clusters (resolution 0.5), and runs the pairwise group comparisons with DESeq2. |
| `ENRICHMENT_gprofiler2_FIXED_skip.R` | Functional enrichment of the DESeq2 results with g:Profiler (gprofiler2), over GO, KEGG and Reactome, plus the summary plots. |
| `cellchat_procedure1_res05.R` | Builds one CellChat object per experimental group and infers cell-cell communication. |
| `cellchat_procedure2_res05.R` | Pairwise comparison between conditions: differential interactions, pathway ranking, and up/down-regulated ligand-receptor pairs. |
| `cellchat_procedure2_replot.R` | Regenerates the CellChat comparison figures in the more readable form used in the thesis, reusing the objects and tables produced by the previous script. |

Each script starts with a short header describing its inputs, its outputs and
the choices it makes.

## Notes for readers

- **The data are not included.** The Seurat object and the sequencing data are
  not part of this repository. The scripts are provided for inspection, to
  document how the results in the thesis were produced.
- **Paths are specific to the machine they ran on.** Input and output
  directories are hard-coded near the top of each script and point to the HPC
  cluster used for the analysis. They would need to be adapted before running
  the scripts elsewhere.
- **"GSEA" in file and variable names is a misnomer.** The analysis labelled
  that way is an ordered-query enrichment with g:Profiler, not a classical gene
  set enrichment analysis. The names were kept as they were when the analyses
  were run, so that the files match the outputs; the thesis describes the method
  correctly.
- **The LPS group is not analysed in the thesis.** Some comparisons involving it
  are still defined in the scripts, and their outputs were simply not used.
- Comments inside the code are in Italian; the headers are in English.

## Software

Analyses were run in R on a SLURM cluster, using Seurat, DESeq2, gprofiler2 and
CellChat, with ggplot2, patchwork and ComplexHeatmap for the figures.

<!-- TODO: add the exact R and package versions from sessionInfo() -->

## Contact

Samuel Evan Cipolletti — University of Parma.
