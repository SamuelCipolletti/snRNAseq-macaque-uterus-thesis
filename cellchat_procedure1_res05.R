#!/usr/bin/env Rscript
###############################################################################
# CellChat, Procedure 1: cell-cell communication inference per condition
#
# Builds a separate CellChat object for each experimental group, using the
# 17 clusters (resolution 0.5) as cell groups, the SCT-normalised expression
# data and the full human ligand-receptor database (CellChatDB.human).
# Groups: saline control, the four E. coli + antibiotics groups (with or
# without Anakinra, each split by PTL outcome), and two merged groups in
# which PTL and no-PTL animals are combined. LPS is not included.
#
# Communication probabilities are computed with the trimean method, cell
# groups with fewer than 10 cells are filtered out, and each object is saved
# as an .rds file used as input by Procedure 2.
###############################################################################

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(future)
})

## ---------------------------------------------------------------------------
## Parametri (unica sezione da toccare)
## ---------------------------------------------------------------------------
SEURAT_RDS   <- "/hpc/scratch/samuelevan.cipolletti/res05_annotated.rds"
OUT_DIR      <- "/hpc/home/samuelevan.cipolletti/snRNAseq/CellChat/res_05"
GROUP_BY     <- "integrated_snn_res.0.5"   # 17 cluster numerici, 0-16
CONDITION_COL<- "group"                    # colonna delle condizioni
ASSAY        <- "SCT"                       # SCT@data e' normalizzato e popolato;
                                            # coerente con lo spazio in cui sono
                                            # definiti i cluster 0-16
N_WORKERS    <- 4                           # ridotti da 8: multisession duplica i dati per worker (fix OOM)
MIN_CELLS    <- 10                          # soglia filterCommunication
COMPUTE_TYPE <- "triMean"                   # metodo media per gruppo (default v2)

# Definizione dei gruppi da processare.
# Ciascun elemento: nome_output = vettore di condizioni (dalla colonna 'group').
# Un solo valore = condizione pura; piu' valori = gruppo FUSO (nuclei uniti).
# IA_saline = controllo. LPS escluso.
GROUP_DEFS <- list(
  IA_saline                    = "IA_saline",
  IA_Ecoli_Abx_no_PTL          = "IA_Ecoli_Abx_no_PTL",
  IA_Ecoli_Abx_PTL             = "IA_Ecoli_Abx_PTL",
  IA_Ecoli_Abx_Anakinra_no_PTL = "IA_Ecoli_Abx_Anakinra_no_PTL",
  IA_Ecoli_Abx_Anakinra_PTL    = "IA_Ecoli_Abx_Anakinra_PTL",
  # --- gruppi FUSI (PTL + noPTL insieme), servono per C1, C5, C7 ---
  EcoliAbx_all     = c("IA_Ecoli_Abx_PTL", "IA_Ecoli_Abx_no_PTL"),
  EcoliAbxAnak_all = c("IA_Ecoli_Abx_Anakinra_PTL", "IA_Ecoli_Abx_Anakinra_no_PTL")
)

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# future: alza il limite sulla dimensione dei globals esportati ai worker
# (default 500 MiB troppo basso: i gruppi grandi hanno matrici ~1 GiB) e
# imposta RNG parallel-safe per riproducibilita' del permutation test.
options(future.globals.maxSize = 8 * 1024^3)   # 8 GiB
options(future.rng.onMisuse = "ignore")        # RNG gestito via future.seed nei chiamanti

# Database L-R: human, DATABASE COMPLETO (incluso non-protein signaling).
CellChatDB     <- CellChatDB.human
CellChatDB.use <- CellChatDB            # tutto il database, nessuna esclusione

cat("== Carico l'oggetto Seurat ==\n")
obj <- readRDS(SEURAT_RDS)
cat(sprintf("Nuclei totali: %d\n", ncol(obj)))

# Sanity check colonne
stopifnot(GROUP_BY      %in% colnames(obj@meta.data))
stopifnot(CONDITION_COL %in% colnames(obj@meta.data))
stopifnot(ASSAY         %in% Assays(obj))

# assay di default = quello che useremo, evita ambiguita' nel subset (Seurat v5)
DefaultAssay(obj) <- ASSAY

## ---------------------------------------------------------------------------
## Funzione: costruisci + inferisci CellChat per un gruppo (puro o fuso)
## ---------------------------------------------------------------------------
run_one_group <- function(group_name, conds) {

  out_file <- file.path(OUT_DIR, sprintf("cellchat_%s.rds", group_name))

  # skip-if-exists (idempotenza: non ri-computa run gia' completate)
  if (file.exists(out_file)) {
    cat(sprintf("[SKIP] %s : output gia' presente (%s)\n", group_name, out_file))
    return(invisible(NULL))
  }

  is_fused <- length(conds) > 1
  cat(sprintf("\n===== GRUPPO: %s %s=====\n",
              group_name, if (is_fused) "(FUSO) " else ""))
  cat(sprintf("  Condizioni incluse: %s\n", paste(conds, collapse = " + ")))

  # --- subset ai nuclei di questa/e condizione/i (%in% gestisce anche i fusi) ---
  cells_keep <- rownames(obj@meta.data)[obj@meta.data[[CONDITION_COL]] %in% conds]
  sub <- subset(obj, cells = cells_keep)
  cat(sprintf("  Nuclei in %s: %d\n", group_name, ncol(sub)))

  # --- input per CellChat ---
  # matrice normalizzata (layer 'data' dell'assay scelto) + metadata
  # Seurat v5: si usa 'layer' (slot e' deprecato)
  data.input <- GetAssayData(sub, assay = ASSAY, layer = "data")

  # 'labels' come factor con livelli ordinati numericamente.
  # NB: CellChat rifiuta etichette che sono numeri puri ("Cell labels cannot
  # contain `0`"), quindi anteponiamo il prefisso "C" -> C0, C1, ... C16.
  # Restano univoci; la mappa C<n> -> tipo cellulare si applica a valle.
  # 'labels' come factor con livelli ordinati numericamente.
  # NB: CellChat rifiuta etichette che sono numeri puri ("Cell labels cannot
  # contain `0`"), quindi anteponiamo il prefisso "cluster" -> cluster0..cluster16.
  # Uso "cluster" (non "C") per non confondere con i confronti C1..C10.
  # Restano univoci; la mappa cluster<n> -> tipo cellulare si applica a valle.
  clu_num <- as.integer(as.character(sub@meta.data[[GROUP_BY]]))
  lev     <- paste0("cluster", sort(unique(clu_num)))
  meta <- data.frame(
    labels    = factor(paste0("cluster", clu_num), levels = lev),
    samples   = as.character(sub@meta.data[["sample"]]),  # ID animale = campione biologico
    row.names = colnames(sub)
  )
  meta$samples <- factor(meta$samples)

  cat(sprintf("  Cluster presenti: %s\n", paste(levels(meta$labels), collapse = ", ")))
  cat(sprintf("  Campioni (animali): %s\n", paste(levels(meta$samples), collapse = ", ")))

  # --- crea oggetto CellChat ---
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  cellchat@DB <- CellChatDB.use

  # --- pre-processing ---
  cellchat <- subsetData(cellchat)  # necessario anche usando tutto il DB
  future::plan("multisession", workers = N_WORKERS)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # --- inferenza comunicazione ---
  # raw.use = TRUE: usa i dati grezzi (no smoothing PPI), coerente col default v2
  cellchat <- computeCommunProb(cellchat, type = COMPUTE_TYPE, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  # --- centrality (utile per Procedure 2 e per i signaling role plot) ---
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  # torna a esecuzione sequenziale prima del salvataggio
  future::plan("sequential")

  saveRDS(cellchat, file = out_file)
  cat(sprintf("  [OK] salvato: %s\n", out_file))

  # libera memoria tra un gruppo e l'altro
  rm(sub, data.input, cellchat); gc()
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## Loop principale
## ---------------------------------------------------------------------------
for (gname in names(GROUP_DEFS)) {
  tryCatch(
    run_one_group(gname, GROUP_DEFS[[gname]]),
    error = function(e) {
      cat(sprintf("  [ERRORE] %s : %s\n", gname, conditionMessage(e)))
    }
  )
}

cat("\n== FATTO. Oggetti CellChat in:", OUT_DIR, "==\n")
