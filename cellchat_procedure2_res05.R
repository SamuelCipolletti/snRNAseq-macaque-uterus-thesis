#!/usr/bin/env Rscript
###############################################################################
# CellChat, Procedure 2: pairwise comparison between conditions
#
# Takes the per-group CellChat objects produced by Procedure 1 and runs the
# pairwise comparisons C1-C10 (plus C5A and C5B). In each comparison the two
# objects are merged, with dataset1 as reference and dataset2 as test, so that
# "up" means higher in the test group, consistent with the DESeq2 convention.
#
# Per comparison it produces: overall number and strength of interactions,
# differential heatmap and circle plot by cell-type pair, pathway ranking by
# information flow, and a DEG-based analysis of up- and down-regulated
# ligand-receptor pairs exported as CSV and shown as bubble plots. The merged
# object is saved as .rds and a flag file marks completed comparisons, which
# are skipped on rerun.
###############################################################################

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(ComplexHeatmap)
})

## ---------------------------------------------------------------------------
## Parametri
## ---------------------------------------------------------------------------
IN_DIR  <- "/hpc/home/samuelevan.cipolletti/snRNAseq/CellChat/res_05"          # .rds da Procedure 1
OUT_DIR <- "/hpc/home/samuelevan.cipolletti/snRNAseq/CellChat/res_05/Comparisons"

# soglie DEG per L-R up/down (coerenti col tutorial v2)
THRESH_PC <- 0.1
THRESH_FC <- 0.05
THRESH_P  <- 0.05
LIGAND_LOGFC <- 0.05   # soglia |logFC| ligando per up/down
TOP_N_BUBBLE <- 40     # n. max di coppie L-R nel bubble plot (top per |logFC|); NA = tutte

## ---------------------------------------------------------------------------
## Definizione dei CONFRONTI
## Ogni elemento: list(ref = <gruppo dataset1>, test = <gruppo dataset2>)
## I nomi dei gruppi corrispondono ai file cellchat_<nome>.rds di Procedure 1.
## ---------------------------------------------------------------------------
COMPARISONS <- list(
  C1  = list(ref = "IA_saline",                    test = "EcoliAbx_all"),                  # Control vs EcoliAbx (fuso)
  C2  = list(ref = "IA_saline",                    test = "IA_Ecoli_Abx_PTL"),              # Control vs EcoliAbx PTL
  C3  = list(ref = "IA_saline",                    test = "IA_Ecoli_Abx_no_PTL"),           # Control vs EcoliAbx noPTL
  C4  = list(ref = "IA_Ecoli_Abx_no_PTL",          test = "IA_Ecoli_Abx_PTL"),              # noPTL vs PTL (infezione)
  C5  = list(ref = "EcoliAbx_all",                 test = "EcoliAbxAnak_all"),              # EcoliAbx vs EcoliAbxAnak (fuso)
  C5A = list(ref = "IA_Ecoli_Abx_PTL",             test = "IA_Ecoli_Abx_Anakinra_PTL"),     # EcoliAbx PTL vs +Anak PTL
  C5B = list(ref = "IA_Ecoli_Abx_no_PTL",          test = "IA_Ecoli_Abx_Anakinra_no_PTL"),  # EcoliAbx noPTL vs +Anak noPTL
  C6  = list(ref = "IA_Ecoli_Abx_PTL",             test = "IA_Ecoli_Abx_Anakinra_no_PTL"),  # EcoliAbx PTL vs EcoliAbxAnak noPTL
  C7  = list(ref = "IA_saline",                    test = "EcoliAbxAnak_all"),              # Control vs EcoliAbxAnak (fuso)
  C8  = list(ref = "IA_saline",                    test = "IA_Ecoli_Abx_Anakinra_PTL"),     # Control vs EcoliAbxAnak PTL
  C9  = list(ref = "IA_saline",                    test = "IA_Ecoli_Abx_Anakinra_no_PTL"),  # Control vs EcoliAbxAnak noPTL
  C10 = list(ref = "IA_Ecoli_Abx_Anakinra_no_PTL", test = "IA_Ecoli_Abx_Anakinra_PTL")      # EcoliAbxAnak noPTL vs PTL
)

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

rds_path <- function(group_name) file.path(IN_DIR, sprintf("cellchat_%s.rds", group_name))

## ---------------------------------------------------------------------------
## Funzione: esegui un confronto
## ---------------------------------------------------------------------------
run_comparison <- function(cmp_name, ref_name, test_name) {

  cmp_dir <- file.path(OUT_DIR, cmp_name)
  done_flag <- file.path(cmp_dir, "_DONE")

  if (file.exists(done_flag)) {
    cat(sprintf("[SKIP] %s : gia' completato\n", cmp_name))
    return(invisible(NULL))
  }

  cat(sprintf("\n===== CONFRONTO %s : %s (ref) vs %s (test) =====\n",
              cmp_name, ref_name, test_name))

  f_ref  <- rds_path(ref_name)
  f_test <- rds_path(test_name)
  if (!file.exists(f_ref) || !file.exists(f_test)) {
    cat(sprintf("  [ERRORE] manca un .rds (%s o %s) — Procedure 1 completata?\n",
                basename(f_ref), basename(f_test)))
    return(invisible(NULL))
  }

  dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

  # --- carica e merge (dataset1 = ref, dataset2 = test) ---
  cc.ref  <- readRDS(f_ref)
  cc.test <- readRDS(f_test)
  object.list <- list()
  object.list[[ref_name]]  <- cc.ref
  object.list[[test_name]] <- cc.test
  cellchat <- mergeCellChat(object.list, add.names = names(object.list),
                            cell.prefix = TRUE)

  ## --- (1) confronto globale: n. interazioni e forza ---
  gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2),
                             measure = "weight")
  ggsave(file.path(cmp_dir, "01_compareInteractions.pdf"),
         gg1 + gg2, width = 8, height = 4)

  ## --- (2) differenze per coppia di cell type ---
  # heatmap differenziale (n. interazioni e forza)
  pdf(file.path(cmp_dir, "02_diff_heatmap.pdf"), width = 12, height = 6)
  print(netVisual_heatmap(cellchat) +
        netVisual_heatmap(cellchat, measure = "weight"))
  dev.off()

  # circle plot differenziale + pannello legenda nomi cluster
  # layout: 3 pannelli (n. interazioni | forza | legenda cluster)
  clabs <- levels(cellchat@idents$joint)
  if (is.null(clabs)) clabs <- levels(cellchat@idents[[1]])
  ccols <- tryCatch(scPalette(length(clabs)), error = function(e) rainbow(length(clabs)))

  pdf(file.path(cmp_dir, "02_diff_circle.pdf"), width = 16, height = 6)
  layout(matrix(c(1,2,3), nrow = 1), widths = c(1, 1, 0.5))
  par(xpd = TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
  # pannello legenda: nomi dei cluster con relativo colore
  par(mar = c(1, 0, 2, 0))
  plot.new()
  title("Cluster labels", cex.main = 1.1)
  legend("center", legend = clabs, col = ccols, pch = 19, pt.cex = 1.8,
         bty = "n", ncol = 1, cex = 0.9, y.intersp = 1.1)
  dev.off()

  ## --- (3) ranking pathway per information flow ---
  # stacked = TRUE con test di Wilcoxon: pathway significativamente diversi
  gg_rank <- rankNet(cellchat, mode = "comparison", measure = "weight",
                     stacked = TRUE, do.stat = TRUE)
  ggsave(file.path(cmp_dir, "03_rankNet_pathways.pdf"),
         gg_rank, width = 6, height = 10)

  ## --- (4) analisi DEG-based: L-R up/down-regolati nel dataset2 (test) ---
  # pos.dataset = test => "up" e' piu' alto nel gruppo test (non-reference)
  pos.dataset   <- test_name
  features.name <- paste0(pos.dataset, ".merged")

  cellchat <- tryCatch({
    identifyOverExpressedGenes(
      cellchat, group.dataset = "datasets",
      pos.dataset = pos.dataset, features.name = features.name,
      only.pos = FALSE, thresh.pc = THRESH_PC,
      thresh.fc = THRESH_FC, thresh.p = THRESH_P,
      group.DE.combined = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  [WARN] identifyOverExpressedGenes fallito: %s\n",
                conditionMessage(e)))
    NULL
  })

  if (!is.null(cellchat)) {
    net <- netMappingDEG(cellchat, features.name = features.name,
                         variable.all = TRUE)

    # up nel test (ligando up nel dataset2)
    net.up <- subsetCommunication(cellchat, net = net, datasets = test_name,
                                  ligand.logFC = LIGAND_LOGFC, receptor.logFC = NULL)
    # down nel test = up nel ref (ligando up nel dataset1)
    net.down <- subsetCommunication(cellchat, net = net, datasets = ref_name,
                                    ligand.logFC = -LIGAND_LOGFC, receptor.logFC = NULL)

    write.csv(net.up,   file.path(cmp_dir, "04_LR_up_in_test.csv"),   row.names = FALSE)
    write.csv(net.down, file.path(cmp_dir, "04_LR_down_in_test.csv"), row.names = FALSE)
    cat(sprintf("  L-R up nel test: %d | down nel test: %d\n",
                nrow(net.up), nrow(net.down)))

    # bubble plot up/down: limita alle top-N coppie L-R per |logFC| (leggibilita')
    # net.up/net.down contengono la colonna ligand.logFC dalla mappatura DEG.
    top_pairs <- function(df, n = TOP_N_BUBBLE) {
      # coppia L-R unica ordinata per |ligand.logFC| decrescente
      lfc_col <- intersect(c("ligand.logFC", "receptor.logFC"), colnames(df))[1]
      if (!is.na(n) && !is.null(lfc_col) && lfc_col %in% colnames(df)) {
        ord <- order(abs(df[[lfc_col]]), decreasing = TRUE)
        df  <- df[ord, , drop = FALSE]
      }
      pr <- unique(df[["interaction_name"]])
      if (!is.na(n) && length(pr) > n) pr <- pr[seq_len(n)]
      data.frame(interaction_name = pr, stringsAsFactors = FALSE)
    }

    if (nrow(net.up) > 0) {
      pairLR.up <- top_pairs(net.up)
      gg_up <- tryCatch(
        netVisual_bubble(cellchat, pairLR.use = pairLR.up, comparison = c(1,2),
                         angle.x = 90, remove.isolate = TRUE,
                         title.name = sprintf("Up in %s (top %d L-R)",
                                              test_name, nrow(pairLR.up))),
        error = function(e) NULL)
      if (!is.null(gg_up))
        ggsave(file.path(cmp_dir, "05_bubble_up.pdf"), gg_up, width = 12, height = 10,
               limitsize = FALSE)
    }
    if (nrow(net.down) > 0) {
      pairLR.down <- top_pairs(net.down)
      gg_dn <- tryCatch(
        netVisual_bubble(cellchat, pairLR.use = pairLR.down, comparison = c(1,2),
                         angle.x = 90, remove.isolate = TRUE,
                         title.name = sprintf("Down in %s (top %d L-R)",
                                              test_name, nrow(pairLR.down))),
        error = function(e) NULL)
      if (!is.null(gg_dn))
        ggsave(file.path(cmp_dir, "05_bubble_down.pdf"), gg_dn, width = 12, height = 10,
               limitsize = FALSE)
    }
  }

  # salva l'oggetto merged per esplorazione successiva
  saveRDS(cellchat, file.path(cmp_dir, sprintf("merged_%s.rds", cmp_name)))

  file.create(done_flag)  # flag di completamento (idempotenza)
  cat(sprintf("  [OK] %s completato -> %s\n", cmp_name, cmp_dir))

  rm(cc.ref, cc.test, object.list, cellchat); gc()
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## Loop principale
## ---------------------------------------------------------------------------
for (cmp_name in names(COMPARISONS)) {
  d <- COMPARISONS[[cmp_name]]
  tryCatch(
    run_comparison(cmp_name, d$ref, d$test),
    error = function(e) {
      cat(sprintf("  [ERRORE] %s : %s\n", cmp_name, conditionMessage(e)))
    }
  )
}

cat("\n== FATTO. Confronti in:", OUT_DIR, "==\n")
