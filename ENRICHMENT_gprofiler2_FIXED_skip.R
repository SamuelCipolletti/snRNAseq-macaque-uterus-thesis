#!/usr/bin/env Rscript

################################################################################
# Functional enrichment of pseudobulk DEGs (g:Profiler, gprofiler2)
#
# Reads the DESeq2 result tables produced by pseudobulk_unified.R and, for each
# level and comparison, queries g:Profiler (Macaca mulatta; GO:BP, GO:MF, GO:CC,
# KEGG, Reactome; FDR < 0.05):
#   - GO: up- and down-regulated DEGs (padj < 0.05, |log2FC| >= 1), analysed
#     as two separate ordered queries;
#   - "GSEA": all genes with a valid adjusted p-value, ranked by log2FC, as a
#     single ordered query. Despite the file names, this is g:Profiler
#     ordered-query enrichment, not a classical GSEA.
#
# Related terms are grouped into broader labels for plotting. Outputs: one CSV
# per level, comparison and analysis, aggregated CSVs, and two PDFs of dot
# plots (by comparison across levels, and by level across comparisons).
# Existing enrichment CSVs are read back instead of being recomputed.
################################################################################

library(gprofiler2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
library(grid)

# 0. CONFIGURAZIONE
# ==============================================================================

BASE_DIR <- "/hpc/home/samuelevan.cipolletti/snRNAseq/Pseudobulk/res_05"

CELLTYPES <- c(
  "pseudobulk_total",
  "decidual_epithelial_cluster0",
  "decidual_epithelial_cluster6",
  "decidual_epithelial_cluster9",
  "decidual_epithelial_cluster11",
  "decidual_epithelial_cluster12",
  "decidual_stromal_cluster7",
  "smooth_muscle_cluster10",
  # --- 10 cluster aggiunti ---
  "lymphatic_endothelial_cluster1",
  "vascular_endothelial_cluster2",
  "fibroblasts_cluster3",
  "myeloid_cluster4",
  "macrophages_monocytes_cluster5",
  "lymphoid_cluster8",
  "macrophages_cluster13",
  "vascular_endothelial_cluster14",
  "perivascular_cluster15",
  "fibroblasts_cluster16"
)

CELLTYPES_LABELS <- c(
  "pseudobulk_total"              = "Total\nPseudobulk",
  "decidual_epithelial_cluster0"  = "Epithelial\ncluster 0",
  "decidual_epithelial_cluster6"  = "Epithelial\ncluster 6",
  "decidual_epithelial_cluster9"  = "Epithelial\ncluster 9",
  "decidual_epithelial_cluster11" = "Epithelial\ncluster 11",
  "decidual_epithelial_cluster12" = "Epithelial\ncluster 12",
  "decidual_stromal_cluster7"     = "Decidual\nStromal",
  "smooth_muscle_cluster10"       = "Smooth\nMuscle",
  "lymphatic_endothelial_cluster1"  = "Lymphatic\nendothelial 1",
  "vascular_endothelial_cluster2"   = "Vascular\nendothelial 2",
  "fibroblasts_cluster3"            = "Fibroblasts\ncluster 3",
  "myeloid_cluster4"                = "Myeloid\ncluster 4",
  "macrophages_monocytes_cluster5"  = "Macrophages/\nmonocytes 5",
  "lymphoid_cluster8"               = "Lymphoid\ncluster 8",
  "macrophages_cluster13"           = "Macrophages\ncluster 13",
  "vascular_endothelial_cluster14"  = "Vascular\nendothelial 14",
  "perivascular_cluster15"          = "Perivascular\ncluster 15",
  "fibroblasts_cluster16"           = "Fibroblasts\ncluster 16"
)

# Numero di nuclei per livello — letto DINAMICAMENTE da cell_counts.csv
# (prodotto da pseudobulk_unified.R), non piu' hardcoded.
cellcount_file <- file.path(BASE_DIR, "cell_counts.csv")
if (file.exists(cellcount_file)) {
  cc <- read.csv(cellcount_file, stringsAsFactors = FALSE)
  NUCLEI_PER_CT <- setNames(cc$n_cells, cc$level)[CELLTYPES]
  if (any(is.na(NUCLEI_PER_CT))) {
    warning("Nuclei mancanti in cell_counts.csv per: ",
            paste(CELLTYPES[is.na(NUCLEI_PER_CT)], collapse = ", "))
  }
} else {
  # Non blocco qui: il run-script vuole comunque il file (accanto ai CSV DESeq2),
  # e controlla piu' avanti; il plot-script fa fallback sui nuclei salvati
  # nell'.rds. Metto NA come segnaposto.
  warning("cell_counts.csv non trovato in ", BASE_DIR,
          " — NUCLEI_PER_CT impostato a NA (il run-script lo richiede; ",
          "il plot-script usa i valori dall'.rds).")
  NUCLEI_PER_CT <- setNames(rep(NA_integer_, length(CELLTYPES)), CELLTYPES)
}

COMPARISONS <- list(
  C1  = "DE_C1_Control_vs_EcoliAbx.csv",
  C2  = "DE_C2_Control_vs_EcoliAbx_PTL.csv",
  C3  = "DE_C3_Control_vs_EcoliAbx_noPTL.csv",
  C4  = "DE_C4_EcoliAbx_noPTL_vs_PTL.csv",
  C5  = "DE_C5_EcoliAbx_vs_EcoliAbxAnak.csv",
  C5A = "DE_C5A_EcoliAbx_PTL_vs_EcoliAbxAnak_PTL.csv",
  C5B = "DE_C5B_EcoliAbx_noPTL_vs_EcoliAbxAnak_noPTL.csv",
  C6  = "DE_C6_EcoliAbx_PTL_vs_EcoliAbxAnak_noPTL.csv",
  C7  = "DE_C7_Control_vs_EcoliAbxAnak.csv",
  C8  = "DE_C8_Control_vs_EcoliAbxAnak_PTL.csv",
  C9  = "DE_C9_Control_vs_EcoliAbxAnak_noPTL.csv",
  C10 = "DE_C10_EcoliAbxAnak_PTL_vs_noPTL.csv",
  # --- confronti LPS ---
  C11 = "DE_C11_LPS_vs_Control.csv",
  C13 = "DE_C13_LPS_vs_EcoliAbx_noPTL.csv",
  C14 = "DE_C14_LPS_vs_EcoliAbx_PTL.csv"
)

COMPARISON_LABELS <- list(
  C1  = "C1: Control vs EcoliAbx (PTL+noPTL)",
  C2  = "C2: Control vs EcoliAbx PTL",
  C3  = "C3: Control vs EcoliAbx noPTL",
  C4  = "C4: EcoliAbx noPTL vs PTL",
  C5  = "C5: EcoliAbx (PTL+noPTL) vs EcoliAbxAnak (PTL+noPTL)",
  C5A = "C5A: EcoliAbx PTL vs EcoliAbxAnak PTL",
  C5B = "C5B: EcoliAbx noPTL vs EcoliAbxAnak noPTL",
  C6  = "C6: EcoliAbx PTL vs EcoliAbxAnak noPTL",
  C7  = "C7: Control vs EcoliAbxAnak (PTL+noPTL)",
  C8  = "C8: Control vs EcoliAbxAnak PTL",
  C9  = "C9: Control vs EcoliAbxAnak noPTL",
  C10 = "C10: EcoliAbxAnak PTL vs EcoliAbxAnak noPTL",
  C11 = "C11: LPS vs Control",
  C13 = "C13: LPS vs EcoliAbx noPTL (asimmetrico: 16h vs esito)",
  C14 = "C14: LPS vs EcoliAbx PTL (asimmetrico: 16h vs esito)"
)

# Sigle brevi ma parlanti per l'asse X del PDF per-livello (12 colonne).
# Convenzione: Ctrl=Control, Eco=EcoliAbx, Anak=EcoliAbxAnak; esito tra ().
# Il codice originale resta tra [] cosi' il collegamento a C1..C10 e alla
# legenda estesa e' immediato. \n spezza in due righe corte per non
# sovrapporre le colonne quando ruotate.
COMPARISON_SHORT <- c(
  C1  = "Ctrl vs Eco\n(all) [C1]",
  C2  = "Ctrl vs Eco\n(PTL) [C2]",
  C3  = "Ctrl vs Eco\n(noPTL) [C3]",
  C4  = "Eco noPTL\nvs PTL [C4]",
  C5  = "Eco vs Anak\n(all) [C5]",
  C5A = "Eco vs Anak\n(PTL) [C5A]",
  C5B = "Eco vs Anak\n(noPTL) [C5B]",
  C6  = "Eco PTL vs\nAnak noPTL [C6]",
  C7  = "Ctrl vs Anak\n(all) [C7]",
  C8  = "Ctrl vs Anak\n(PTL) [C8]",
  C9  = "Ctrl vs Anak\n(noPTL) [C9]",
  C10 = "Anak PTL\nvs noPTL [C10]",
  C11 = "LPS vs\nControl [C11]",
  C13 = "LPS vs\nEcoli noPTL [C13]",
  C14 = "LPS vs\nEcoli PTL [C14]"
)

PADJ_THRESH      <- 0.05
LFC_THRESH       <- 1.0
MIN_INTERSECTION <- 3
TOP_N            <- 10

ORGANISM <- "mmulatta"
SOURCES  <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")

# Ogni livello scrive i suoi CSV GO/GSEA in <BASE_DIR>/<livello>/Csv/Enrichment/
# (definito dentro il loop principale, dipende da ct).
#
# Gli output cross-celltype (aggregati + PDF comparativo) vanno in una
# cartella "Summary", affiancata ai livelli, per analisi che vanno oltre
# il singolo celltype.
SOURCE_COLORS <- c(
  "GO:BP" = "#66C2A5",
  "GO:MF" = "#FC8D62",
  "GO:CC" = "#8DA0CB",
  "KEGG"  = "#E78AC3",
  "REAC"  = "#A6D854"
)

# ==============================================================================
# 0b. PERCORSI DERIVATI + oggetto intermedio run->plots
# ==============================================================================

# Cartella Summary (aggregati cross-celltype + PDF Summary)
SUMMARY_DIR <- file.path(BASE_DIR, "Summary")
dir.create(file.path(SUMMARY_DIR, "Csv", "Enrichment"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(SUMMARY_DIR, "Plots", "Enrichment"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. MAPPA DI ACCORPAMENTO TERMINI
# ==============================================================================

TERM_MAP <- c(
  # MOVIMENTO CELLULARE
  "cell migration"                                          = "cell migration /\nmovement",
  "cell motility"                                           = "cell migration /\nmovement",
  "locomotion"                                              = "cell migration /\nmovement",
  "regulation of cell migration"                            = "cell migration /\nmovement",
  "positive regulation of cell migration"                   = "cell migration /\nmovement",
  "negative regulation of cell migration"                   = "cell migration /\nmovement",
  "regulation of cell motility"                             = "cell migration /\nmovement",
  "positive regulation of cell motility"                    = "cell migration /\nmovement",
  "negative regulation of cell motility"                    = "cell migration /\nmovement",
  "regulation of locomotion"                                = "cell migration /\nmovement",
  "positive regulation of locomotion"                       = "cell migration /\nmovement",
  "negative regulation of locomotion"                       = "cell migration /\nmovement",
  "cell chemotaxis"                                         = "cell migration /\nmovement",
  "chemotaxis"                                              = "cell migration /\nmovement",
  "negative regulation of neutrophil migration"             = "cell migration /\nmovement",
  "negative regulation of neutrophil chemotaxis"            = "cell migration /\nmovement",
  "negative regulation of granulocyte chemotaxis"           = "cell migration /\nmovement",
  "cellular response to chemokine"                          = "cell migration /\nmovement",
  "response to chemokine"                                   = "cell migration /\nmovement",
  "chemokine-mediated signaling pathway"                    = "cell migration /\nmovement",
  "taxis"                                                   = "cell migration /\nmovement",
  # PROLIFERAZIONE
  "cell population proliferation"                           = "cell proliferation\n& growth",
  "regulation of cell population proliferation"             = "cell proliferation\n& growth",
  "negative regulation of cell population proliferation"    = "cell proliferation\n& growth",
  "positive regulation of cell population proliferation"    = "cell proliferation\n& growth",
  "smooth muscle cell proliferation"                        = "cell proliferation\n& growth",
  "regulation of smooth muscle cell proliferation"          = "cell proliferation\n& growth",
  "muscle cell proliferation"                               = "cell proliferation\n& growth",
  "epithelial cell proliferation"                           = "cell proliferation\n& growth",
  "negative regulation of epithelial cell proliferation"    = "cell proliferation\n& growth",
  "regulation of epithelial cell proliferation"             = "cell proliferation\n& growth",
  "keratinocyte proliferation"                              = "cell proliferation\n& growth",
  # APOPTOSI
  "apoptotic process"                                       = "apoptosis /\nprogrammed cell death",
  "programmed cell death"                                   = "apoptosis /\nprogrammed cell death",
  "cell death"                                              = "apoptosis /\nprogrammed cell death",
  "regulation of apoptotic process"                         = "apoptosis /\nprogrammed cell death",
  "regulation of programmed cell death"                     = "apoptosis /\nprogrammed cell death",
  "positive regulation of apoptotic process"                = "apoptosis /\nprogrammed cell death",
  "positive regulation of programmed cell death"            = "apoptosis /\nprogrammed cell death",
  "negative regulation of neuron apoptotic process"         = "apoptosis /\nprogrammed cell death",
  # INFIAMMAZIONE
  "inflammatory response"                                   = "inflammatory /\nimmune response",
  "humoral immune response"                                 = "inflammatory /\nimmune response",
  "complement activation"                                   = "inflammatory /\nimmune response",
  "defense response"                                        = "inflammatory /\nimmune response",
  "myeloid leukocyte activation"                            = "inflammatory /\nimmune response",
  "regulation of immune system process"                     = "inflammatory /\nimmune response",
  "response to cytokine"                                    = "inflammatory /\nimmune response",
  "cellular response to cytokine stimulus"                  = "inflammatory /\nimmune response",
  "response to peptide"                                     = "inflammatory /\nimmune response",
  # SEGNALAZIONE
  "intracellular signal transduction"                       = "intracellular\nsignaling",
  "regulation of intracellular signal transduction"         = "intracellular\nsignaling",
  "negative regulation of intracellular signal transduction" = "intracellular\nsignaling",
  "intracellular signaling cassette"                        = "intracellular\nsignaling",
  "cell surface receptor signaling pathway"                 = "intracellular\nsignaling",
  "cell surface receptor protein tyrosine kinase signaling pathway" = "intracellular\nsignaling",
  "enzyme-linked receptor protein signaling pathway"        = "intracellular\nsignaling",
  "signal transduction"                                     = "intracellular\nsignaling",
  "signaling"                                               = "intracellular\nsignaling",
  "cell communication"                                      = "intracellular\nsignaling",
  "regulation of signaling"                                 = "intracellular\nsignaling",
  "regulation of cell communication"                        = "intracellular\nsignaling",
  # SVILUPPO VASCOLARE
  "vasculature development"                                 = "vascular\ndevelopment",
  "blood vessel development"                                = "vascular\ndevelopment",
  "blood vessel morphogenesis"                              = "vascular\ndevelopment",
  "circulatory system development"                          = "vascular\ndevelopment",
  "negative regulation of angiogenesis"                     = "vascular\ndevelopment",
  "negative regulation of vasculature development"          = "vascular\ndevelopment",
  "negative regulation of blood vessel morphogenesis"       = "vascular\ndevelopment",
  "cardiac chamber morphogenesis"                           = "vascular\ndevelopment",
  # METABOLISMO PROTEICO
  "protein metabolic process"                               = "protein metabolism\n& modification",
  "protein modification process"                            = "protein metabolism\n& modification",
  "post-translational protein modification"                 = "protein metabolism\n& modification",
  "protein modification by small protein conjugation"       = "protein metabolism\n& modification",
  "protein modification by small protein removal"           = "protein metabolism\n& modification",
  "macromolecule modification"                              = "protein metabolism\n& modification",
  "macromolecule catabolic process"                         = "protein metabolism\n& modification",
  "catabolic process"                                       = "protein metabolism\n& modification",
  "regulation of primary metabolic process"                 = "protein metabolism\n& modification",
  # TRASPORTO
  "intracellular transport"                                 = "intracellular transport\n& localization",
  "intracellular protein localization"                      = "intracellular transport\n& localization",
  "macromolecule localization"                              = "intracellular transport\n& localization",
  "establishment of localization in cell"                   = "intracellular transport\n& localization",
  "protein localization to organelle"                       = "intracellular transport\n& localization",
  "vesicle-mediated transport"                              = "intracellular transport\n& localization",
  # ORGANELLI
  "organelle organization"                                  = "organelle\norganization",
  "regulation of organelle organization"                    = "organelle\norganization",
  "regulation of organelle assembly"                        = "organelle\norganization",
  # ADESIONE
  "cell adhesion"                                           = "cell adhesion\n& junction",
  "cell-substrate adhesion"                                 = "cell adhesion\n& junction",
  "cell junction organization"                              = "cell adhesion\n& junction",
  "cell junction assembly"                                  = "cell adhesion\n& junction",
  "regulation of cell junction assembly"                    = "cell adhesion\n& junction",
  "cell junction"                                           = "cell adhesion\n& junction",
  # SVILUPPO ANATOMICO
  "anatomical structure development"                        = "anatomical development\n& morphogenesis",
  "anatomical structure morphogenesis"                      = "anatomical development\n& morphogenesis",
  "anatomical structure formation involved in morphogenesis" = "anatomical development\n& morphogenesis",
  "developmental process"                                   = "anatomical development\n& morphogenesis",
  "multicellular organism development"                      = "anatomical development\n& morphogenesis",
  "system development"                                      = "anatomical development\n& morphogenesis",
  "tissue development"                                      = "anatomical development\n& morphogenesis",
  "animal organ development"                                = "anatomical development\n& morphogenesis",
  "animal organ morphogenesis"                              = "anatomical development\n& morphogenesis",
  # STRESS
  "response to stress"                                      = "stress /\nstimulus response",
  "cellular response to stress"                             = "stress /\nstimulus response",
  "response to stimulus"                                    = "stress /\nstimulus response",
  "cellular response to stimulus"                           = "stress /\nstimulus response",
  "response to chemical"                                    = "stress /\nstimulus response",
  "cellular response to chemical stimulus"                  = "stress /\nstimulus response",
  "response to oxygen-containing compound"                  = "stress /\nstimulus response",
  "response to oxidative stress"                            = "stress /\nstimulus response",
  "response to endogenous stimulus"                         = "stress /\nstimulus response",
  # SINAPSI
  "synapse"                                                 = "synapse &\nneurotransmission",
  "synapse organization"                                    = "synapse &\nneurotransmission",
  "synapse assembly"                                        = "synapse &\nneurotransmission",
  "synaptic membrane"                                       = "synapse &\nneurotransmission",
  "postsynapse"                                             = "synapse &\nneurotransmission",
  "postsynaptic density"                                    = "synapse &\nneurotransmission",
  "postsynaptic density membrane"                           = "synapse &\nneurotransmission",
  "postsynaptic specialization"                             = "synapse &\nneurotransmission",
  "postsynaptic specialization membrane"                    = "synapse &\nneurotransmission",
  "postsynaptic membrane"                                   = "synapse &\nneurotransmission",
  "neuron to neuron synapse"                                = "synapse &\nneurotransmission",
  "asymmetric synapse"                                      = "synapse &\nneurotransmission",
  "excitatory synapse"                                      = "synapse &\nneurotransmission",
  "glutamatergic synapse"                                   = "synapse &\nneurotransmission",
  "neurotransmitter receptor complex"                       = "synapse &\nneurotransmission",
  "ionotropic glutamate receptor complex"                   = "synapse &\nneurotransmission",
  "vesicle-mediated transport in synapse"                   = "synapse &\nneurotransmission",
  "negative regulation of synaptic transmission"            = "synapse &\nneurotransmission",
  "regulation of postsynaptic membrane potential"           = "synapse &\nneurotransmission",
  # TRASPORTO IONI
  "monoatomic ion transport"                                = "ion transport &\nchannel activity",
  "monoatomic ion transmembrane transport"                  = "ion transport &\nchannel activity",
  "monoatomic cation transport"                             = "ion transport &\nchannel activity",
  "monoatomic cation transmembrane transport"               = "ion transport &\nchannel activity",
  "potassium ion transport"                                 = "ion transport &\nchannel activity",
  "potassium ion transmembrane transport"                   = "ion transport &\nchannel activity",
  "transmembrane transport"                                 = "ion transport &\nchannel activity",
  "sodium ion homeostasis"                                  = "ion transport &\nchannel activity",
  "monoatomic ion channel activity"                         = "ion transport &\nchannel activity",
  "channel activity"                                        = "ion transport &\nchannel activity",
  "gated channel activity"                                  = "ion transport &\nchannel activity",
  "ligand-gated channel activity"                           = "ion transport &\nchannel activity",
  "monoatomic ion transmembrane transporter activity"       = "ion transport &\nchannel activity",
  "monoatomic cation transmembrane transporter activity"    = "ion transport &\nchannel activity",
  "passive transmembrane transporter activity"              = "ion transport &\nchannel activity",
  "transmembrane transporter activity"                      = "ion transport &\nchannel activity",
  # BINDING NUCLEOTIDI
  "ATP binding"                                             = "nucleotide\nbinding (MF)",
  "adenyl ribonucleotide binding"                           = "nucleotide\nbinding (MF)",
  "adenyl nucleotide binding"                               = "nucleotide\nbinding (MF)",
  "purine ribonucleoside triphosphate binding"              = "nucleotide\nbinding (MF)",
  "purine nucleotide binding"                               = "nucleotide\nbinding (MF)",
  "nucleotide binding"                                      = "nucleotide\nbinding (MF)",
  "nucleoside phosphate binding"                            = "nucleotide\nbinding (MF)",
  "heterocyclic compound binding"                           = "nucleotide\nbinding (MF)",
  "purine ribonucleotide binding"                           = "nucleotide\nbinding (MF)",
  # BINDING PROTEICO
  "enzyme binding"                                          = "protein / enzyme\nbinding (MF)",
  "protein binding"                                         = "protein / enzyme\nbinding (MF)",
  "protein-macromolecule adaptor activity"                  = "protein / enzyme\nbinding (MF)",
  "molecular adaptor activity"                              = "protein / enzyme\nbinding (MF)",
  "identical protein binding"                               = "protein / enzyme\nbinding (MF)",
  "protein-containing complex binding"                      = "protein / enzyme\nbinding (MF)",
  # CITOPLASMA/NUCLEO CC
  "cytoplasm"                                               = "cytoplasm /\nnucleus (CC)",
  "cytosol"                                                 = "cytoplasm /\nnucleus (CC)",
  "nucleoplasm"                                             = "cytoplasm /\nnucleus (CC)",
  "nuclear body"                                            = "cytoplasm /\nnucleus (CC)",
  "cytoskeleton"                                            = "cytoplasm /\nnucleus (CC)",
  "microtubule cytoskeleton"                                = "cytoplasm /\nnucleus (CC)",
  "microtubule organizing center"                           = "cytoplasm /\nnucleus (CC)",
  # EXTRACELLULARE CC
  "extracellular space"                                     = "extracellular\nspace & matrix (CC)",
  "extracellular region"                                    = "extracellular\nspace & matrix (CC)",
  "extracellular matrix"                                    = "extracellular\nspace & matrix (CC)",
  "extracellular exosome"                                   = "extracellular\nspace & matrix (CC)",
  "extracellular vesicle"                                   = "extracellular\nspace & matrix (CC)",
  "extracellular organelle"                                 = "extracellular\nspace & matrix (CC)",
  "extracellular membrane-bounded organelle"                = "extracellular\nspace & matrix (CC)",
  "external encapsulating structure"                        = "extracellular\nspace & matrix (CC)",
  "basement membrane"                                       = "extracellular\nspace & matrix (CC)",
  # MEMBRANA PLASMATICA CC
  "plasma membrane"                                         = "plasma membrane\n(CC)",
  "plasma membrane region"                                  = "plasma membrane\n(CC)",
  "cell periphery"                                          = "plasma membrane\n(CC)",
  "cell surface"                                            = "plasma membrane\n(CC)",
  "basolateral plasma membrane"                             = "plasma membrane\n(CC)",
  "apical plasma membrane"                                  = "plasma membrane\n(CC)",
  "apical part of cell"                                     = "plasma membrane\n(CC)",
  "lateral plasma membrane"                                 = "plasma membrane\n(CC)",
  "membrane raft"                                           = "plasma membrane\n(CC)",
  "membrane microdomain"                                    = "plasma membrane\n(CC)",
  # COMPLESSI INTRACELLULARI CC
  "catalytic complex"                                       = "intracellular\ncomplexes (CC)",
  "transferase complex"                                     = "intracellular\ncomplexes (CC)",
  "intracellular protein-containing complex"                = "intracellular\ncomplexes (CC)",
  "organelle membrane"                                      = "intracellular\ncomplexes (CC)",
  "endomembrane system"                                     = "intracellular\ncomplexes (CC)",
  "Golgi apparatus"                                         = "intracellular\ncomplexes (CC)",
  "endoplasmic reticulum"                                   = "intracellular\ncomplexes (CC)",
  "mitochondrion"                                           = "intracellular\ncomplexes (CC)",
  "ribonucleoprotein granule"                               = "intracellular\ncomplexes (CC)"
)

# ==============================================================================
# 2. HELPER CONDIVISI (fix_parents, accorpamento termini)
# ==============================================================================

fix_parents <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  for (col in colnames(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- sapply(df[[col]], function(x) paste(unlist(x), collapse = ","))
    }
  }
  if ("parents" %in% colnames(df)) df$parents <- as.character(df$parents)
  return(df)
}

# Applica la mappa di accorpamento a un vettore di term_name
# Regole di accorpamento per parola chiave (ordine = priorita).
# Catturano famiglie di termini automaticamente, robusto a termini nuovi.
KEYWORD_RULES <- list(
  c("(T cell|B cell|lymphocyte|leukocyte|immun|inflammat|cytokine|interleukin|interferon|chemokine|antigen|MHC|complement|defense response|natural killer|macrophage|neutrophil|monocyte|granulocyte|mast cell|dendritic cell|innate|adaptive|tumor necrosis|toll-like|NF-kappa|hypersensitivit|histamine|phagocyt)", "immune /\ninflammatory"),
  c("(apopto|cell death|programmed cell death|necropt|pyropto|anoikis|cell killing)", "apoptosis /\ncell death"),
  c("(migrat|motil|chemotax|locomot|taxis|extravasation|cell projection)", "cell migration /\nmovement"),
  c("(transcription|DNA-binding|RNA polymerase|chromatin|gene expression|gene silencing|spliceosom|mRNA|RNA splic|RNA processing|RNA metaboli|histone|nucleosome|transcript)", "transcription /\nRNA processing"),
  c("(proliferat|cell cycle|mitotic|mitosis|meiotic|meiosis|cell division|cell growth|DNA replication|DNA repair|chromosome segregation|centrosome|spindle|kinetochore)", "proliferation /\ncell cycle"),
  c("(vascul|angiogen|blood vessel|circulatory|heart|cardiac|endotheli|aorta|atri|ventric|coronary)", "vascular /\ncardiac development"),
  c("(neuron|synap|axon|dendrit|glial|nervous system|brain|neural|cerebr|cerebell|forebrain|hindbrain|cognition|memory|learning|behavior)", "neuronal /\nsynaptic"),
  c("(muscle|myofibr|sarcomere|contracti|actomyosin|myoblast|myotube|myosin)", "muscle /\ncontraction"),
  c("(adhesion|junction|integrin|cadherin|focal adh|desmosome)", "cell adhesion /\njunction"),
  c("(lipid|cholesterol|fatty acid|sterol|phospholipid|triglycerid|lipoprotein|ceramide|sphingo|glycerolipid|steroid|prostagland|icosanoid)", "lipid /\nsteroid metabolism"),
  c("(ubiquitin|proteasom|protein catabolic|protein modification|protein folding|protein localization|protein transport|protein phosphor|protein metaboli|protein stab|protein matur|peptidyl|protein processing|protein targeting|translation)", "protein metabolism\n& modification"),
  c("(transport|transmembrane|channel|symport|exchang|\\bimport\\b|\\bexport\\b|secretion|exocyt|endocyt|vesicle-mediated|localization to|establishment of.*localization)", "transport /\nlocalization"),
  c("(metaboli|biosynthe|catabolic|glycoly|oxidati|respirat|mitochond|nucleotide|nucleoside|purine|pyrimidine|amino acid|carbohydrate|glucose|glycogen)", "metabolism /\nbiosynthesis"),
  c("(differentiat|development|morphogen|maturation|fate commit|specification|pattern|organ|epitheli|mesenchym|gastrul|embryo|tissue|regenerat)", "development /\nmorphogenesis"),
  c("(kinase|phosphat|MAPK|GTPase|receptor signaling|signal transduction|signaling pathway|signaling|cascade)", "signaling"),
  c("(stress|response to|stimulus|hypoxia|oxidative|heat|wound|detoxif|starvation)", "stress /\nstimulus response"),
  c("(homeostasis|regulation of biological quality|pH reduction|pH elevation|redox)", "homeostasis"),
  c("(extracellular space|extracellular region|extracellular matrix|extracellular exosome|extracellular vesicle|external encapsulating|basement membrane)", "extracellular\nspace & matrix (CC)"),
  c("(plasma membrane|cell periphery|cell surface|apical|basolateral|membrane raft|microvillus|microdomain)", "plasma membrane\n(CC)"),
  c("(cytoplasm|cytosol|nucleoplasm|nuclear|cytoskelet|microtubule|secretory granule|granule|organelle|complex|reticulum|Golgi|ribosom|mitochond|vesicle|lysosom|vacuol|endosom|membrane)", "cellular\ncomponent (CC)"),
  c("(binding|activity|regulator|transporter|inhibitor)$", "binding /\nactivity (MF)"),
  c("(regulation of cellular process|regulation of biological process|regulation of metabolic process|biological regulation|cellular process)", "regulation of\ncellular process")
)

# Classifica un singolo termine con le regole keyword
classify_keyword <- function(term) {
  t <- tolower(term)
  for (rule in KEYWORD_RULES) {
    if (grepl(tolower(rule[1]), t, perl = TRUE)) {
      return(rule[2])
    }
  }
  return(NA_character_)
}

apply_term_map <- function(term_names) {
  # Livello 1: mappa esplicita (controllo preciso sui termini comuni)
  mapped <- TERM_MAP[term_names]

  # Livello 2: regole keyword per i non mappati
  unmapped_idx <- which(is.na(mapped))
  if (length(unmapped_idx) > 0) {
    kw <- vapply(term_names[unmapped_idx], classify_keyword, character(1))
    mapped[unmapped_idx] <- kw
  }

  # Livello 3: troncamento + a capo per quel che resta
  still_idx <- which(is.na(mapped))
  if (length(still_idx) > 0) {
    original <- term_names[still_idx]
    shortened <- ifelse(
      nchar(original) > 35,
      paste0(substr(original, 1, 18), "\n", substr(original, 19, 50)),
      original
    )
    mapped[still_idx] <- shortened
  }

  return(unname(mapped))
}

# 3. LOOP PRINCIPALE
# ==============================================================================

results_by_comp <- list()
deg_counts      <- list()   # conta DEG UP/DOWN per celltype x confronto

for (ct in CELLTYPES) {

  message("\n", strrep("=", 60))
  message("CELLTYPE: ", ct)
  message(strrep("=", 60))

  ct_csv_dir <- file.path(BASE_DIR, ct, "Csv")
  if (!dir.exists(ct_csv_dir)) { warning("Directory non trovata: ", ct_csv_dir); next }

  out_enrich_ct <- file.path(BASE_DIR, ct, "Csv", "Enrichment")
  dir.create(out_enrich_ct, recursive = TRUE, showWarnings = FALSE)

  for (comp_name in names(COMPARISONS)) {

    message("\n  Confronto: ", comp_name)

    csv_file <- file.path(ct_csv_dir, COMPARISONS[[comp_name]])
    if (!file.exists(csv_file)) {
      alt_files <- list.files(ct_csv_dir,
                              pattern = paste0("DE_", comp_name, "_.*\\.csv$"),
                              full.names = TRUE)
      if (length(alt_files) == 0) { message("    CSV non trovato — skip"); next }
      csv_file <- alt_files[1]
    }

    res <- tryCatch(
      read.csv(csv_file, row.names = 1, stringsAsFactors = FALSE),
      error = function(e) { message("    ERRORE lettura: ", e$message); NULL }
    )
    if (is.null(res)) next
    if (!all(c("log2FoldChange", "padj") %in% colnames(res))) { next }

    res$gene <- rownames(res)
    res      <- res %>% filter(!is.na(log2FoldChange), !is.na(padj))

    sig        <- res %>% filter(padj < PADJ_THRESH, abs(log2FoldChange) >= LFC_THRESH)
    genes_UP   <- sig %>% filter(log2FoldChange > 0) %>% arrange(desc(log2FoldChange)) %>% pull(gene)
    genes_DOWN <- sig %>% filter(log2FoldChange < 0) %>% arrange(log2FoldChange) %>% pull(gene)

    message("    DEG UP: ", length(genes_UP), " | DOWN: ", length(genes_DOWN))

    # Salva i conteggi DEG
    if (!comp_name %in% names(deg_counts)) deg_counts[[comp_name]] <- list()
    deg_counts[[comp_name]][[ct]] <- list(up = length(genes_UP), down = length(genes_DOWN))

    # Riduce quello che resta in RAM (results_by_comp) al minimo necessario
    # per Summary + bubble plot per-livello: solo le colonne usate a valle,
    # e solo i top TOP_N termini per source — esattamente quello che finira'
    # comunque in qualsiasi pagina. Il CSV su disco (write.csv sotto) resta
    # SEMPRE completo, non filtrato: questa riduzione riguarda solo l'oggetto
    # accumulato in memoria per tutta la durata dello script (causa dell'OOM
    # con migliaia di termini GSEA grezzi tenuti per 96 combinazioni).
    slim_for_memory <- function(df) {
      if (is.null(df) || nrow(df) == 0) return(df)
      keep_cols <- intersect(
        c("term_name", "p_value", "intersection_size", "source",
          "significant", "direction", "celltype", "comparison"),
        colnames(df)
      )
      df %>%
        select(all_of(keep_cols)) %>%
        group_by(source) %>%
        slice_min(p_value, n = TOP_N, with_ties = FALSE) %>%
        ungroup()
    }

    run_GO <- function(gene_list, direction) {
      go_csv <- file.path(out_enrich_ct, paste0("GO_", comp_name, "_", ct, "_", direction, ".csv"))
      # SKIP: se il CSV esiste gia', lo rileggo invece di richiamare gprofiler2.
      if (file.exists(go_csv)) {
        df <- tryCatch(read.csv(go_csv, stringsAsFactors = FALSE),
                       error = function(e) NULL)
        if (!is.null(df) && nrow(df) > 0) {
          message("    GO ", direction, ": letto da CSV esistente (", nrow(df), " termini) [skip]")
          return(slim_for_memory(df))
        }
        # CSV vuoto/illeggibile: ricalcolo sotto
      }
      if (length(gene_list) < 5) return(NULL)

      # BUG NOTO gprofiler2: con evcodes=TRUE, in alcuni edge case (poche
      # significant terms) l'assegnazione interna di evidence_codes va in
      # errore "replacement has 1 row, data has 0". Fallback automatico:
      # riprova senza evcodes (si perdono solo le evidence codes, non i
      # termini) invece di scartare il confronto/direzione.
      run_gost_go <- function(evcodes) {
        gost(query = gene_list, organism = ORGANISM, ordered_query = TRUE,
             multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
             evcodes = evcodes, user_threshold = PADJ_THRESH,
             correction_method = "fdr", sources = SOURCES)
      }

      gp <- tryCatch(run_gost_go(TRUE), error = function(e) {
        message("    GO ", direction, ": evcodes=TRUE fallito (", e$message,
                "), riprovo senza evcodes")
        tryCatch(run_gost_go(FALSE), error = function(e2) {
          message("    GO ", direction, " ERRORE (anche senza evcodes): ", e2$message)
          NULL
        })
      })

      if (is.null(gp) || is.null(gp$result) || nrow(gp$result) == 0) return(NULL)
      tryCatch({
        df <- fix_parents(gp$result)
        df <- df %>% filter(intersection_size >= MIN_INTERSECTION)
        if (nrow(df) == 0) { message("    GO ", direction, ": 0 termini dopo filtro intersection_size"); return(NULL) }
        df$direction <- direction; df$celltype <- ct; df$comparison <- comp_name
        write.csv(df, go_csv, row.names = FALSE)
        message("    GO ", direction, ": ", nrow(df), " termini")
        slim_for_memory(df)
      }, error = function(e) { message("    GO ", direction, " ERRORE post-processing: ", e$message); NULL })
    }

    go_up   <- run_GO(genes_UP,   "UP")
    go_down <- run_GO(genes_DOWN, "DOWN")

    ranked_genes <- res %>% arrange(desc(log2FoldChange)) %>% pull(gene)
    gsea_res     <- NULL

    gsea_csv <- file.path(out_enrich_ct, paste0("GSEA_", comp_name, "_", ct, ".csv"))
    # SKIP: se il CSV GSEA esiste gia', lo rileggo.
    if (file.exists(gsea_csv)) {
      df <- tryCatch(read.csv(gsea_csv, stringsAsFactors = FALSE),
                     error = function(e) NULL)
      if (!is.null(df) && nrow(df) > 0) {
        message("    GSEA: letto da CSV esistente (", nrow(df), " termini) [skip]")
        gsea_res <- slim_for_memory(df)
      }
    }
    if (is.null(gsea_res) && length(ranked_genes) >= 10) {
      # Stesso fallback evcodes del blocco GO (vedi commento sopra).
      run_gost_gsea <- function(evcodes) {
        gost(query = ranked_genes, organism = ORGANISM, ordered_query = TRUE,
             multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
             evcodes = evcodes, user_threshold = PADJ_THRESH,
             correction_method = "fdr", sources = SOURCES)
      }
      gp_gsea <- tryCatch(run_gost_gsea(TRUE), error = function(e) {
        message("    GSEA: evcodes=TRUE fallito (", e$message, "), riprovo senza evcodes")
        tryCatch(run_gost_gsea(FALSE), error = function(e2) {
          message("    GSEA ERRORE (anche senza evcodes): ", e2$message)
          NULL
        })
      })
      gsea_res <- if (is.null(gp_gsea) || is.null(gp_gsea$result) || nrow(gp_gsea$result) == 0) {
        NULL
      } else {
        tryCatch({
          df <- fix_parents(gp_gsea$result)
          df <- df %>% filter(intersection_size >= MIN_INTERSECTION)
          if (nrow(df) == 0) {
            message("    GSEA: 0 termini dopo filtro intersection_size")
            NULL
          } else {
            df$celltype <- ct; df$comparison <- comp_name
            write.csv(df, gsea_csv, row.names = FALSE)
            message("    GSEA: ", nrow(df), " termini")
            slim_for_memory(df)
          }
        }, error = function(e) { message("    GSEA ERRORE post-processing: ", e$message); NULL })
      }
    }

    if (!comp_name %in% names(results_by_comp)) results_by_comp[[comp_name]] <- list()
    results_by_comp[[comp_name]][[ct]] <- list(go_up = go_up, go_down = go_down, gsea = gsea_res)

  }
}

# ==============================================================================
# 4. CSV AGGREGATI (Summary — cross-celltype)
# ==============================================================================

all_GO_list <- all_GSEA_list <- list()
for (comp_name in names(results_by_comp)) {
  for (ct in names(results_by_comp[[comp_name]])) {
    d <- results_by_comp[[comp_name]][[ct]]
    if (!is.null(d$go_up))   all_GO_list[[paste0(ct,"_",comp_name,"_UP")]]   <- d$go_up
    if (!is.null(d$go_down)) all_GO_list[[paste0(ct,"_",comp_name,"_DOWN")]] <- d$go_down
    if (!is.null(d$gsea))    all_GSEA_list[[paste0(ct,"_",comp_name)]]        <- d$gsea
  }
}
summary_csv_dir <- file.path(SUMMARY_DIR, "Csv", "Enrichment")
if (length(all_GO_list)   > 0) write.csv(bind_rows(all_GO_list),
  file.path(summary_csv_dir, "GO_ALL_celltypes_all_comparisons.csv"), row.names = FALSE)
if (length(all_GSEA_list) > 0) write.csv(bind_rows(all_GSEA_list),
  file.path(summary_csv_dir, "GSEA_ALL_celltypes_all_comparisons.csv"), row.names = FALSE)

# ==============================================================================

# ==============================================================================
# PARAMETRI PLOT  <<< MODIFICA QUI PER CAMBIARE L'ASPETTO DEI GRAFICI >>>
# ==============================================================================
# Ogni voce e' usata piu' sotto nelle funzioni di disegno. Cambiando un numero
# qui cambi tutti i PDF senza cercare valori sparsi nel codice.
PP <- list(

  # ---- Dimensioni pagina PDF (pollici) ----
  summary_pdf_w   = 30,   # PDF Summary (X = livello)
  summary_pdf_h   = 14,   # alzato: piu' spazio verticale per le etichette Y
  perlevel_pdf_w  = 16,   # PDF per-livello (X = confronto) — largo per la legenda a lato
  perlevel_pdf_h  = 15,   # alto: molte righe -> evita sovrapposizione etichette

  # ---- Quanti termini tenere per (gruppo) source ----
  top_n_terms     = TOP_N,             # dal file comune
  min_intersection = MIN_INTERSECTION, # dal file comune

  # ---- Dot plot: dimensione dei punti (N geni) ----
  point_size_range = c(1.5, 7),
  point_alpha      = 0.85,

  # ---- Testo assi / legenda dentro i pannelli ----
  base_size       = 8,
  axis_x_angle    = 30,
  axis_x_size     = 7,
  axis_y_size     = 6.5,
  axis_y_lineheight = 0.85,
  legend_text_size  = 6,
  legend_title_size = 7,
  panel_title_size  = 9,

  # ---- Header di pagina (make_page) ----
  page_title_size   = 14,   # titolo confronto/livello
  page_nuclei_size  = 9,    # riga nuclei
  page_header_size  = 11,   # header pannello (GO UP / DOWN / GSEA)
  page_title_col    = "#1B2631",
  page_nuclei_col   = "#555555",
  # altezze relative delle 4 righe: titolo, nuclei, header, pannello
  page_row_heights  = c(0.06, 0.045, 0.05, 0.845),

  # ---- Colori header per tipo di pagina ----
  col_go_up   = "#148F77",
  col_go_down = "#C0392B",
  col_gsea    = "#185FA5",

  # ---- Disambiguazione punti multicolore ----
  # Se TRUE: quando una stessa etichetta-riga (accorpata) compare con PIU'
  # source, ognuna diventa una riga Y separata col suffisso della source,
  # cosi' due punti di source diversa non cadono piu' nella stessa cella
  # (niente piu' "punti bicolore" sovrapposti). Se una riga ha una sola
  # source resta pulita, senza suffisso.
  split_rows_by_source = TRUE,
  # Suffisso source: TRUE = in linea "[BP]" (etichetta piu' bassa, meno righe
  # fisiche); FALSE = a capo "\n(GO:BP)" (piu' alto, causava sovrapposizioni).
  source_tag_inline    = TRUE,

  # ---- Filtro anti-sparsita' ----
  # Tiene un termine solo se compare in almeno questo numero di colonne
  # (= confronti nel per-livello, = livelli nel Summary). Con 2 spariscono le
  # righe a punto singolo, tipiche delle KEGG disease-pathway (COVID, Malaria,
  # Chagas...) che gonfiavano l'asse Y. Metti 1 per disattivare il filtro.
  min_cols_present     = 2,

  # ---- Legenda-codici confronti (solo PDF per-livello, X = C1..C10) ----
  # Nel per-livello l'asse X sono i codici brevi (C1, C2...): questo riquadro
  # ne elenca il significato esteso. Ora messo A LATO (destra), non sotto:
  # sotto risultava schiacciato e disordinato.
  show_comparison_legend = TRUE,
  comparison_legend_size = 7,     # dimensione testo del riquadro-codici
  comparison_legend_width = 0.26, # frazione di larghezza pagina riservata al riquadro

  # ---- Etichette asse X del per-livello ----
  # Nome INTERO del confronto (senza prefisso Cx:), mandato a capo a questa
  # larghezza in caratteri. Piu' piccolo = piu' righe corte = colonne strette.
  axis_x_wrap = 14
)

# ------------------------------------------------------------------------------
# HELPER — testo della legenda-codici confronti (C1 = ..., C2 = ...).
# Costruito una volta da COMPARISON_LABELS; le etichette gia' contengono il
# prefisso "C1: ...", quindi le uso cosi' come sono, una per riga.
# ------------------------------------------------------------------------------
build_comparison_legend_text <- function() {
  labs <- vapply(names(COMPARISONS), function(cn) {
    lab <- COMPARISON_LABELS[[cn]]
    if (is.null(lab)) cn else lab
  }, character(1))
  n <- length(labs)
  # Due colonne per contenere l'altezza (12 righe -> ~6). La colonna sinistra
  # e' incolonnata a larghezza fissa; monospace nel theme mantiene l'allineamento.
  half <- ceiling(n / 2)
  left  <- labs[seq_len(half)]
  right <- labs[(half + 1):n]
  if (length(right) < length(left)) right <- c(right, rep("", length(left) - length(right)))
  colw <- max(nchar(left)) + 3
  rows <- vapply(seq_along(left), function(i) {
    paste0(formatC(left[i], width = -colw, flag = " "), right[i])
  }, character(1))
  paste(rows, collapse = "\n")
}
COMPARISON_LEGEND_TEXT <- build_comparison_legend_text()

# ------------------------------------------------------------------------------
# HELPER — disambigua le righe che aggregano piu' source.
# Prende il data.frame di plotting (deve avere term_name_grouped + source) e
# restituisce lo stesso df con term_name_grouped modificato SOLO per le
# etichette che compaiono con >1 source: a quelle appende " (GO:BP)" ecc.
# Le etichette con una sola source restano invariate.
# ------------------------------------------------------------------------------
disambiguate_source_rows <- function(df) {
  if (!isTRUE(PP$split_rows_by_source) || nrow(df) == 0) return(df)
  multi <- df %>%
    distinct(term_name_grouped, source) %>%
    count(term_name_grouped, name = "n_src") %>%
    filter(n_src > 1) %>%
    pull(term_name_grouped)
  if (length(multi) == 0) return(df)
  needs_tag <- df$term_name_grouped %in% multi
  # Tag source compatto: "GO:BP" -> "BP", "KEGG" -> "KEGG", "REAC" -> "REAC"
  short_src <- sub("^GO:", "", df$source[needs_tag])
  if (isTRUE(PP$source_tag_inline)) {
    # in linea: "... [BP]" — non aggiunge una riga fisica di testo
    df$term_name_grouped[needs_tag] <- paste0(
      df$term_name_grouped[needs_tag], " [", short_src, "]"
    )
  } else {
    df$term_name_grouped[needs_tag] <- paste0(
      df$term_name_grouped[needs_tag], "\n(", df$source[needs_tag], ")"
    )
  }
  df
}

# ------------------------------------------------------------------------------
# HELPER — filtro anti-sparsita': tiene solo i termini (etichette-riga)
# presenti in >= PP$min_cols_present colonne distinte (col = comparison o
# celltype a seconda del pannello). Elimina le righe a punto singolo.
# ------------------------------------------------------------------------------
filter_sparse_rows <- function(df, col) {
  if (is.null(PP$min_cols_present) || PP$min_cols_present <= 1 || nrow(df) == 0) return(df)
  keep <- df %>%
    distinct(term_name_grouped, .data[[col]]) %>%
    count(term_name_grouped, name = "n_cols") %>%
    filter(n_cols >= PP$min_cols_present) %>%
    pull(term_name_grouped)
  df %>% filter(term_name_grouped %in% keep)
}

# ==============================================================================
# FUNZIONE 1 — pannello dot plot per il PDF Summary (X = livello)
# ==============================================================================

make_panel <- function(panel_data, panel_title) {

  if (is.null(panel_data) || nrow(panel_data) == 0) {
    p_empty <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = "No significant\nterms",
               size = 3.5, color = "grey60", hjust = 0.5, vjust = 0.5) +
      labs(title = panel_title) +
      theme_void() +
      theme(
        plot.title      = element_text(size = 9, face = "bold", hjust = 0.5, color = "grey40"),
        plot.background = element_rect(fill = "grey97", color = "grey80", linewidth = 0.5)
      )
    return(p_empty)
  }

  # Applica accorpamento
  panel_data$term_name_grouped <- apply_term_map(panel_data$term_name)

  # Top N termini per source per celltype, poi prendi i termini unici
  plot_df <- panel_data %>%
    filter(significant == TRUE, intersection_size >= PP$min_intersection) %>%
    group_by(celltype, source) %>%
    slice_min(p_value, n = PP$top_n_terms) %>%
    ungroup() %>%
    mutate(
      neg_log10_p    = -log10(p_value),
      celltype_label = CELLTYPES_LABELS[celltype]
    )

  if (nrow(plot_df) == 0) {
    p_empty <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant\nterms",
               size = 3.5, color = "grey60", hjust = 0.5, vjust = 0.5) +
      labs(title = panel_title) +
      theme_void() +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
    return(p_empty)
  }

  # Separa in righe distinte le etichette che aggregano piu' source
  # (evita i punti bicolore sovrapposti nella stessa cella)
  plot_df <- disambiguate_source_rows(plot_df)

  # Filtro anti-sparsita': via le righe presenti in una sola colonna
  plot_df <- filter_sparse_rows(plot_df, "celltype")
  if (nrow(plot_df) == 0) {
    return(make_panel(NULL, panel_title))
  }

  # Ordine termini: per source poi per p_value
  term_order <- plot_df %>%
    group_by(term_name_grouped) %>%
    summarise(max_p = max(neg_log10_p), source = first(source), .groups = "drop") %>%
    arrange(source, desc(max_p)) %>%
    pull(term_name_grouped) %>%
    unique()

  plot_df$term_name_grouped <- factor(plot_df$term_name_grouped, levels = rev(term_order))
  plot_df$celltype_label    <- factor(plot_df$celltype_label, levels = CELLTYPES_LABELS)

  ggplot(plot_df,
         aes(x     = celltype_label,
             y     = term_name_grouped,
             size  = intersection_size,
             color = source)) +
    geom_point(alpha = PP$point_alpha) +
    scale_size_continuous(range = PP$point_size_range, name = "N genes") +
    scale_color_manual(values = SOURCE_COLORS, name = "Source", drop = FALSE) +
    labs(title = panel_title, x = NULL, y = NULL) +
    theme_bw(base_size = PP$base_size) +
    theme(
      axis.text.x      = element_text(angle = PP$axis_x_angle, hjust = 1, size = PP$axis_x_size),
      axis.text.y      = element_text(size = PP$axis_y_size, lineheight = PP$axis_y_lineheight),
      plot.title       = element_text(size = PP$panel_title_size, face = "bold", hjust = 0.5),
      legend.text      = element_text(size = PP$legend_text_size),
      legend.title     = element_text(size = PP$legend_title_size),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

# ==============================================================================
# FUNZIONE 2 — header di pagina (4 righe: titolo / nuclei / header / pannello)
# ==============================================================================

make_page <- function(title_text, nuclei_text, header_text, header_col, panel) {
  title_grob  <- textGrob(title_text,  gp = gpar(fontsize = PP$page_title_size,  fontface = "bold",   col = PP$page_title_col))
  nuclei_grob <- textGrob(nuclei_text, gp = gpar(fontsize = PP$page_nuclei_size, fontface = "italic", col = PP$page_nuclei_col))
  header_grob <- textGrob(header_text, gp = gpar(fontsize = PP$page_header_size, fontface = "bold",   col = header_col))

  layout_matrix <- rbind(c(1), c(2), c(3), c(4))
  arranged <- arrangeGrob(
    title_grob, nuclei_grob, header_grob, panel,
    layout_matrix = layout_matrix,
    heights = PP$page_row_heights
  )
  grid.newpage()
  grid.draw(arranged)
}

# ==============================================================================
# FUNZIONE 3 — pannello dot plot per il PDF per-livello (X = confronto)
# ==============================================================================

make_panel_by_comparison <- function(panel_data) {
  empty_panel <- function() {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant\nterms",
               size = 3.5, color = "grey60", hjust = 0.5, vjust = 0.5) +
      theme_void() +
      theme(plot.background = element_rect(fill = "grey97", color = "grey80", linewidth = 0.5))
  }

  if (is.null(panel_data) || nrow(panel_data) == 0) return(empty_panel())

  panel_data$term_name_grouped <- apply_term_map(panel_data$term_name)

  plot_df <- panel_data %>%
    filter(significant == TRUE, intersection_size >= PP$min_intersection) %>%
    group_by(comparison, source) %>%
    slice_min(p_value, n = PP$top_n_terms) %>%
    ungroup() %>%
    mutate(neg_log10_p = -log10(p_value))

  if (nrow(plot_df) == 0) return(empty_panel())

  # Separa in righe distinte le etichette che aggregano piu' source
  plot_df <- disambiguate_source_rows(plot_df)

  # Filtro anti-sparsita': via le righe presenti in un solo confronto
  plot_df <- filter_sparse_rows(plot_df, "comparison")
  if (nrow(plot_df) == 0) return(empty_panel())

  term_order <- plot_df %>%
    group_by(term_name_grouped) %>%
    summarise(max_p = max(neg_log10_p), source = first(source), .groups = "drop") %>%
    arrange(source, desc(max_p)) %>%
    pull(term_name_grouped) %>%
    unique()

  plot_df$term_name_grouped <- factor(plot_df$term_name_grouped, levels = rev(term_order))
  plot_df$comparison <- factor(plot_df$comparison, levels = names(COMPARISONS))

  # Etichette asse X: nome INTERO del confronto, mandato a capo automaticamente
  # (str_wrap) per stare nella colonna. Tolgo il prefisso "Cx:" ridondante
  # (c'e' la legenda a lato). Orizzontali e centrate: un nome incolonnato su
  # piu' righe corte ingombra meno di uno lungo ruotato.
  full_names <- vapply(levels(plot_df$comparison), function(cn) {
    lab <- COMPARISON_LABELS[[cn]]
    if (is.null(lab)) return(cn)
    sub("^C[0-9A-Z]+:\\s*", "", lab)   # rimuove "C1: ", "C5A: ", ecc.
  }, character(1))
  x_labels <- vapply(full_names, function(s) {
    paste(strwrap(s, width = PP$axis_x_wrap), collapse = "\n")
  }, character(1))

  p <- ggplot(plot_df,
         aes(x = comparison, y = term_name_grouped,
             size = intersection_size, color = source)) +
    geom_point(alpha = PP$point_alpha) +
    scale_size_continuous(range = PP$point_size_range, name = "N genes") +
    scale_color_manual(values = SOURCE_COLORS, name = "Source", drop = FALSE) +
    scale_x_discrete(labels = x_labels) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = PP$base_size) +
    theme(
      axis.text.x        = element_text(angle = 0, hjust = 0.5, vjust = 1,
                                        size = PP$axis_x_size, lineheight = 0.85),
      axis.text.y         = element_text(size = PP$axis_y_size, lineheight = PP$axis_y_lineheight),
      legend.text         = element_text(size = PP$legend_text_size),
      legend.title        = element_text(size = PP$legend_title_size),
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_blank()
    )

  p
}

# ------------------------------------------------------------------------------
# HELPER — grob della legenda-codici a lato (una voce per riga), da affiancare
# a destra del pannello per-livello. Restituisce NULL se disattivata.
# ------------------------------------------------------------------------------
make_comparison_legend_grob <- function() {
  if (!isTRUE(PP$show_comparison_legend)) return(NULL)
  labs <- vapply(names(COMPARISONS), function(cn) {
    lab <- COMPARISON_LABELS[[cn]]
    if (is.null(lab)) cn else lab
  }, character(1))
  body <- paste(labs, collapse = "\n")
  txt <- paste0("Legenda confronti [codice]:\n\n", body)
  textGrob(
    txt, x = 0.02, y = 0.98, just = c("left", "top"),
    gp = gpar(fontsize = PP$comparison_legend_size, fontfamily = "mono",
              lineheight = 1.15, col = "grey20")
  )
}

# ==============================================================================
# PDF 1 — Summary cross-celltype (3 pagine per confronto: GO UP / DOWN / GSEA)
#         X = livello. Output: Summary/Plots/Enrichment/GO_GSEA_per_comparison_v3.pdf
# ==============================================================================
message("\n", strrep("=", 60))
message("Generazione PDF Summary (3 pagine per confronto: GO UP, GO DOWN, GSEA)...")

dir.create(file.path(SUMMARY_DIR, "Plots", "Enrichment"), recursive = TRUE, showWarnings = FALSE)
out_pdf <- file.path(SUMMARY_DIR, "Plots", "Enrichment", "GO_GSEA_per_comparison_v3.pdf")

pdf(out_pdf, width = PP$summary_pdf_w, height = PP$summary_pdf_h)

for (comp_name in names(COMPARISONS)) {

  message("  Confronto: ", comp_name)
  comp_data <- results_by_comp[[comp_name]]

  if (is.null(comp_data) || length(comp_data) == 0) {
    grid.newpage()
    grid.text(paste0(comp_name, " — nessun dato"), gp = gpar(fontsize = 16, col = "grey50"))
    next
  }

  # Aggrega dati per pannello (su tutti gli 8 livelli)
  go_up_all   <- bind_rows(lapply(comp_data, function(x) x$go_up))
  go_down_all <- bind_rows(lapply(comp_data, function(x) x$go_down))
  gsea_all    <- bind_rows(lapply(comp_data, function(x) x$gsea))

  # Conta DEG totali per confronto (somma su tutti i livelli)
  deg_up_ct <- deg_down_ct <- list()
  for (ct in CELLTYPES) {
    if (!is.null(deg_counts[[comp_name]][[ct]])) {
      deg_up_ct[[ct]]   <- deg_counts[[comp_name]][[ct]]$up
      deg_down_ct[[ct]] <- deg_counts[[comp_name]][[ct]]$down
    }
  }
  total_up   <- sum(unlist(deg_up_ct),   na.rm = TRUE)
  total_down <- sum(unlist(deg_down_ct), na.rm = TRUE)

  title_text  <- COMPARISON_LABELS[[comp_name]]
  nuclei_text <- paste(
    sapply(CELLTYPES, function(ct) {
      paste0(gsub("\n", " ", CELLTYPES_LABELS[ct]), ": ",
             format(NUCLEI_PER_CT[ct], big.mark = ","), " nuclei")
    }),
    collapse = "   |   "
  )

  # ---- Pagina 1: GO UP ----
  p_go_up <- make_panel(if (nrow(go_up_all) > 0) go_up_all else NULL, "")
  make_page(title_text, nuclei_text,
            paste0("GO UP  —  UP: ", total_up),
            PP$col_go_up, p_go_up)

  # ---- Pagina 2: GO DOWN ----
  p_go_down <- make_panel(if (nrow(go_down_all) > 0) go_down_all else NULL, "")
  make_page(title_text, nuclei_text,
            paste0("GO DOWN  —  DOWN: ", total_down),
            PP$col_go_down, p_go_down)

  # ---- Pagina 3: GSEA ----
  p_gsea <- make_panel(if (nrow(gsea_all) > 0) gsea_all else NULL, "")
  make_page(title_text, nuclei_text,
            "GSEA — ranked list (all genes)",
            PP$col_gsea, p_gsea)

  message("    ", comp_name, " OK (3 pagine)")
}

dev.off()
message("\nPDF salvato: ", out_pdf)

# ==============================================================================
# PDF 2 — Bubble plot per-livello (3 pagine: GO UP / DOWN / GSEA)
#         X = confronto (C1..C10, C5A, C5B). Un PDF per livello.
#         Output: <livello>/Plots/Enrichment/Bubble_plots/GO_GSEA_bubble_<livello>.pdf
# ==============================================================================
for (ct in CELLTYPES) {

  message("  Livello: ", ct)

  ct_go_up_all   <- bind_rows(lapply(results_by_comp, function(x) x[[ct]]$go_up))
  ct_go_down_all <- bind_rows(lapply(results_by_comp, function(x) x[[ct]]$go_down))
  ct_gsea_all    <- bind_rows(lapply(results_by_comp, function(x) x[[ct]]$gsea))

  total_up_ct <- sum(sapply(names(COMPARISONS), function(cn) {
    v <- deg_counts[[cn]][[ct]]$up
    if (is.null(v)) 0 else v
  }))
  total_down_ct <- sum(sapply(names(COMPARISONS), function(cn) {
    v <- deg_counts[[cn]][[ct]]$down
    if (is.null(v)) 0 else v
  }))

  bubble_dir <- file.path(BASE_DIR, ct, "Plots", "Enrichment", "Bubble_plots")
  dir.create(bubble_dir, recursive = TRUE, showWarnings = FALSE)
  out_pdf_ct <- file.path(bubble_dir, paste0("GO_GSEA_bubble_", ct, ".pdf"))

  pdf(out_pdf_ct, width = PP$perlevel_pdf_w, height = PP$perlevel_pdf_h)

  title_text_ct  <- paste0(gsub("\n", " ", CELLTYPES_LABELS[ct]), " — tutti i confronti")
  nuclei_text_ct <- paste0(format(NUCLEI_PER_CT[ct], big.mark = ","), " nuclei")

  # Affianca a destra di ogni pannello la legenda-codici confronti.
  panel_with_legend <- function(panel) {
    leg <- make_comparison_legend_grob()
    if (is.null(leg)) return(panel)
    arrangeGrob(panel, leg, ncol = 2,
                widths = c(1 - PP$comparison_legend_width, PP$comparison_legend_width))
  }

  make_page(title_text_ct, nuclei_text_ct,
            paste0("GO UP  —  UP totale (tutti i confronti): ", total_up_ct),
            PP$col_go_up, panel_with_legend(make_panel_by_comparison(ct_go_up_all)))

  make_page(title_text_ct, nuclei_text_ct,
            paste0("GO DOWN  —  DOWN totale (tutti i confronti): ", total_down_ct),
            PP$col_go_down, panel_with_legend(make_panel_by_comparison(ct_go_down_all)))

  make_page(title_text_ct, nuclei_text_ct,
            "GSEA — ranked list (all genes)",
            PP$col_gsea, panel_with_legend(make_panel_by_comparison(ct_gsea_all)))

  dev.off()
  message("    Salvato: ", out_pdf_ct)
}

message("
Script grafici completato.")
