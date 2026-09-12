#!/usr/bin/env Rscript

# ================================================================
# Pseudobulk differential expression analysis (DESeq2)
#
# Aggregates RNA counts from the annotated Seurat object into one
# pseudobulk profile per sample, for the whole tissue and for each
# of the 17 clusters (resolution 0.5). For each level, it runs the
# pairwise comparisons between experimental groups with DESeq2
# (C1-C10, C5A, C5B and the LPS comparisons C11, C13, C14).
#
# Outputs per level: DESeq2 result tables (CSV), heatmaps of the
# top DEGs, UP/DOWN barplots, dispersion plots and VST matrices.
# A summary of DEG counts across levels and comparisons is also
# written. Existing outputs are skipped when the script is rerun.
# ================================================================

library(Seurat)
library(Matrix)
library(DESeq2)
library(dplyr)

# --- pacchetti per i plot (guardati: se assenti, i plot vengono saltati) ---
have_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)
have_ggplot   <- requireNamespace("ggplot2",  quietly = TRUE)
if (have_pheatmap) suppressMessages(library(pheatmap)) else
  message("ATTENZIONE: pacchetto 'pheatmap' non disponibile -> heatmap saltate.")
if (have_ggplot)   suppressMessages(library(ggplot2)) else
  message("ATTENZIONE: pacchetto 'ggplot2' non disponibile -> barplot saltati.")
suppressMessages(library(grid))

# ================================================================
# PERCORSI
# ================================================================
seurat_path <- "/hpc/scratch/samuelevan.cipolletti/res05_annotated.rds"
base_outdir <- "/hpc/home/samuelevan.cipolletti/snRNAseq/Pseudobulk/res_05"
dir.create(base_outdir, showWarnings = FALSE, recursive = TRUE)

# ================================================================
# LIVELLI TARGET — "total" = tutte le cellule, nessun filtro cluster
#
#   decidual_epithelial_cluster0  = cluster 0   (possibile sottotipo, es. glandulare/luminale)
#   decidual_epithelial_cluster6  = cluster 6
#   decidual_epithelial_cluster9  = cluster 9
#   decidual_epithelial_cluster11 = cluster 11
#   decidual_epithelial_cluster12 = cluster 12
#   decidual_stromal_cluster7     = cluster 7    (~8.375 nuclei)
#   smooth_muscle_cluster10       = cluster 10   (~4.928 nuclei)
#
# NOTA: la vecchia fusione "decidual_epithelial" (0+6+11) è stata
# sostituita dai 5 livelli separati sopra, per
# indagare eventuali sottotipi epiteliali distinti (glandulare vs
# luminale) prima di aggregare.
# ================================================================
cluster_map <- list(
  "pseudobulk_total"              = NULL,   # NULL = nessun filtro, tutte le cellule
  "decidual_epithelial_cluster0"  = c("0"),
  "decidual_epithelial_cluster6"  = c("6"),
  "decidual_epithelial_cluster9"  = c("9"),
  "decidual_epithelial_cluster11" = c("11"),
  "decidual_epithelial_cluster12" = c("12"),
  "decidual_stromal_cluster7"     = c("7"),
  "smooth_muscle_cluster10"       = c("10"),
  # --- 10 cluster aggiunti (estensione a tutti i 17 livelli) ---
  "lymphatic_endothelial_cluster1"  = c("1"),
  "vascular_endothelial_cluster2"   = c("2"),
  "fibroblasts_cluster3"            = c("3"),
  "myeloid_cluster4"                = c("4"),
  "macrophages_monocytes_cluster5"  = c("5"),
  "lymphoid_cluster8"               = c("8"),
  "macrophages_cluster13"           = c("13"),
  "vascular_endothelial_cluster14"  = c("14"),
  "perivascular_cluster15"          = c("15"),
  "fibroblasts_cluster16"           = c("16")
)

# ================================================================
# CARICAMENTO OGGETTO SEURAT
# ================================================================
message("Caricamento oggetto Seurat...")
obj <- readRDS(seurat_path)
obj@meta.data$cluster_char <- as.character(obj@meta.data$integrated_snn_res.0.5)
message("Cluster presenti: ", paste(sort(unique(obj@meta.data$cluster_char)), collapse = ", "))

# ================================================================
# RICOSTRUZIONE MATRICE COUNTS COMPLETA
# ================================================================
message("Ricostruzione matrice counts...")
layers      <- obj[["RNA"]]@layers
gene_names  <- rownames(obj[["RNA"]])
counts_full <- do.call(cbind, lapply(layers, function(x) x))
rownames(counts_full) <- gene_names
colnames(counts_full) <- colnames(obj)
message(sprintf("Matrice completa: %d geni x %d cellule",
                nrow(counts_full), ncol(counts_full)))

# ================================================================
# METADATA A LIVELLO DI SAMPLE
# ================================================================
meta_full <- obj@meta.data[, c("sample", "group")]
meta_full <- meta_full[!duplicated(meta_full$sample), ]
rownames(meta_full) <- meta_full$sample

meta_full$group7 <- NA
meta_full$group7[grepl("saline",                    meta_full$group)] <- "Control"
meta_full$group7[grepl("Ecoli_Abx_no_PTL",          meta_full$group)] <- "Ecoli_Abx_noPTL"
meta_full$group7[grepl("Ecoli_Abx_PTL",             meta_full$group)] <- "Ecoli_Abx_PTL"
meta_full$group7[grepl("Ecoli_Abx_Anakinra_no_PTL", meta_full$group)] <- "Ecoli_Abx_Anak_noPTL"
meta_full$group7[grepl("Ecoli_Abx_Anakinra_PTL",    meta_full$group)] <- "Ecoli_Abx_Anak_PTL"
meta_full$group7[grepl("LPS",                       meta_full$group)] <- "LPS"
meta_full$group7 <- factor(meta_full$group7)

message("Campioni nei metadata: ", nrow(meta_full))
print(table(meta_full$group7))

# ================================================================
# CONFIG PLOTTING (heatmap + barplot + summary)
# ================================================================
TOP_N_HEATMAP <- 30   # numero massimo di DEG mostrati nell'heatmap per confronto

level_labels <- c(
  pseudobulk_total              = "Tessuto totale (pseudobulk)",
  decidual_epithelial_cluster0  = "Decidual (glandular/luminal) epithelial \u2014 cluster 0",
  decidual_epithelial_cluster6  = "Decidual (glandular/luminal) epithelial \u2014 cluster 6",
  decidual_epithelial_cluster9  = "Decidual (glandular/luminal) epithelial \u2014 cluster 9",
  decidual_epithelial_cluster11 = "Decidual (glandular/luminal) epithelial \u2014 cluster 11",
  decidual_epithelial_cluster12 = "Decidual (glandular/luminal) epithelial \u2014 cluster 12",
  decidual_stromal_cluster7     = "Decidual stromal (and epithelial) \u2014 cluster 7",
  smooth_muscle_cluster10       = "(Myometrial) smooth muscle \u2014 cluster 10",
  lymphatic_endothelial_cluster1  = "Lymphatic endothelial \u2014 cluster 1",
  vascular_endothelial_cluster2   = "Vascular endothelial \u2014 cluster 2",
  fibroblasts_cluster3            = "Fibroblasts \u2014 cluster 3",
  myeloid_cluster4                = "Myeloid cells (mixture) \u2014 cluster 4",
  macrophages_monocytes_cluster5  = "Macrophages/monocytes \u2014 cluster 5",
  lymphoid_cluster8               = "Lymphoid cells (mixture) \u2014 cluster 8",
  macrophages_cluster13           = "Macrophages (mixture myeloid) \u2014 cluster 13",
  vascular_endothelial_cluster14  = "Vascular endothelial \u2014 cluster 14",
  perivascular_cluster15          = "Perivascular \u2014 cluster 15",
  fibroblasts_cluster16           = "Fibroblasts \u2014 cluster 16"
)

# ordine preferito dei gruppi nelle colonne dell'heatmap
group_order <- c("Control", "LPS", "EcoliAbx", "Ecoli_Abx_noPTL", "Ecoli_Abx_PTL",
                 "EcoliAbx_Anakinra", "Ecoli_Abx_Anak_noPTL", "Ecoli_Abx_Anak_PTL")
group_colors <- c(
  Control              = "#2ca25f",
  LPS                  = "#e6b800",
  EcoliAbx             = "#f4a259",
  Ecoli_Abx_noPTL      = "#f4a259",
  Ecoli_Abx_PTL        = "#c0392b",
  EcoliAbx_Anakinra    = "#8e7cc3",
  Ecoli_Abx_Anak_noPTL = "#8e7cc3",
  Ecoli_Abx_Anak_PTL   = "#674ea7"
)
dir_colors  <- c(UP = "#c0392b", DOWN = "#2c5f8a")
comp_order  <- c("C1","C2","C3","C4","C5","C5A","C5B","C6","C7","C8","C9","C10",
                 "C11","C13","C14")

# accumulatore conteggi DEG (per barplot per-livello e summary heatmap)
DEG_COUNTS <- data.frame(level = character(), comparison = character(),
                         n_up = integer(), n_down = integer(), total = integer(),
                         stringsAsFactors = FALSE)

# ================================================================
# FUNZIONI
# ================================================================

save_res <- function(res, name, outdir, skip_existing = TRUE) {
  outfile <- file.path(outdir, paste0("DE_", name, ".csv"))
  if (skip_existing && file.exists(outfile)) {
    message("    SKIP (esiste già): ", basename(outfile))
    return(invisible(NULL))
  }
  write.csv(as.data.frame(res), file = outfile, row.names = TRUE)
  message("    Salvato: ", basename(outfile))
}

# DESeq2 con gruppi base (group7), nessuna fusione
make_dds <- function(pb_counts, meta, groups_to_keep) {
  keep       <- meta$group7 %in% groups_to_keep
  sub_meta   <- meta[keep, , drop = FALSE]
  sub_counts <- pb_counts[, rownames(sub_meta), drop = FALSE]
  sub_meta$group7 <- droplevels(sub_meta$group7)

  dds <- DESeqDataSetFromMatrix(
    countData = round(sub_counts),
    colData   = sub_meta,
    design    = ~ group7
  )
  DESeq(dds)
}

# DESeq2 con gruppi FUSI (per C1, C5, C7)
# fuse_map: lista nome_nuovo = c(gruppo1, gruppo2)
make_dds_fused <- function(pb_counts, meta, groups_to_keep, fuse_map) {
  keep     <- meta$group7 %in% groups_to_keep
  sub_meta <- meta[keep, , drop = FALSE]

  sub_meta$group_fused <- as.character(sub_meta$group7)
  for (new_name in names(fuse_map)) {
    old_names <- fuse_map[[new_name]]
    sub_meta$group_fused[sub_meta$group_fused %in% old_names] <- new_name
  }
  sub_meta$group_fused <- factor(sub_meta$group_fused)

  sub_counts <- pb_counts[, rownames(sub_meta), drop = FALSE]

  dds <- DESeqDataSetFromMatrix(
    countData = round(sub_counts),
    colData   = sub_meta,
    design    = ~ group_fused
  )
  DESeq(dds)
}

# ================================================================
# FUNZIONE: esegue i 12 confronti su un dato pseudobulk + metadata
# ================================================================
prettify_comp <- function(name) gsub("_", " ", name)

# ---- Heatmap Top-N DEG (z-score VST per campione) per un singolo confronto ----
make_deg_heatmap <- function(dds, res, name, level_outdir, group_var,
                             top_n = TOP_N_HEATMAP, skip_existing = TRUE) {
  if (!have_pheatmap) return(invisible(NULL))
  level_dir <- dirname(level_outdir)      # <base>/<livello>
  level_key <- basename(level_dir)        # <livello>
  plots_dir <- file.path(level_dir, "Plots", "Heatmaps")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  outfile <- file.path(plots_dir, paste0("heatmap_", name, ".pdf"))
  if (skip_existing && file.exists(outfile)) {
    message("    SKIP heatmap (esiste): ", basename(outfile)); return(invisible(NULL))
  }

  rd  <- as.data.frame(res); rd <- rd[!is.na(rd$padj), ]
  sig <- rd[rd$padj < 0.05 & abs(rd$log2FoldChange) > 1, , drop = FALSE]
  if (nrow(sig) < 2) {
    message("    Heatmap saltata (", nrow(sig), " DEG significativi): ", name)
    return(invisible(NULL))
  }
  sig <- sig[order(sig$padj), , drop = FALSE]
  top <- head(sig, top_n)
  # UP prima, poi DOWN; dentro ciascun blocco per |log2FC| decrescente
  top <- top[order(top$log2FoldChange < 0, -abs(top$log2FoldChange)), , drop = FALSE]

  vsd <- tryCatch(vst(dds, blind = FALSE),
                  error = function(e) varianceStabilizingTransformation(dds, blind = FALSE))
  mat   <- assay(vsd)[rownames(top), , drop = FALSE]
  mat_z <- t(scale(t(mat))); mat_z[is.na(mat_z)] <- 0

  grp <- as.character(colData(dds)[[group_var]])
  ord <- order(match(grp, group_order))
  mat_z <- mat_z[, ord, drop = FALSE]; grp <- grp[ord]

  ann_col <- data.frame(Group = grp);          rownames(ann_col) <- colnames(mat_z)
  ann_row <- data.frame(Direction = ifelse(top$log2FoldChange > 0, "UP", "DOWN"))
  rownames(ann_row) <- rownames(top)

  grp_cols <- group_colors[intersect(names(group_colors), unique(grp))]
  miss <- setdiff(unique(grp), names(grp_cols))
  if (length(miss)) grp_cols <- c(grp_cols, setNames(rep("#999999", length(miss)), miss))
  ann_colors <- list(Group = grp_cols, Direction = dir_colors)

  lab <- level_labels[level_key]; if (is.na(lab)) lab <- level_key
  ttl <- paste0(lab, "\n", prettify_comp(name),
                "\nTop ", nrow(top), " DEG \u2014 z-score VST per campione")

  pal <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  ph <- pheatmap(mat_z, color = pal,
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 annotation_col = ann_col, annotation_row = ann_row,
                 annotation_colors = ann_colors,
                 show_colnames = TRUE, show_rownames = TRUE,
                 fontsize_row = 8, fontsize_col = 8,
                 breaks = seq(-2, 2, length.out = 101),
                 main = ttl, silent = TRUE)
  h <- max(4, 0.22 * nrow(mat_z) + 2.5)
  w <- max(6, 0.45 * ncol(mat_z) + 5)
  pdf(outfile, width = w, height = h); grid::grid.draw(ph$gtable); dev.off()
  message("    Heatmap salvato: ", basename(outfile))
}

# ---- Calcola res, salva CSV, accumula conteggi, genera heatmap ----
emit <- function(dds, contrast, name, level_outdir) {
  res <- results(dds, contrast = contrast)
  save_res(res, name, level_outdir)
  rd  <- as.data.frame(res); rd <- rd[!is.na(rd$padj), ]
  sig <- rd[rd$padj < 0.05 & abs(rd$log2FoldChange) > 1, , drop = FALSE]
  sigla <- sub("_.*", "", name)
  DEG_COUNTS <<- rbind(DEG_COUNTS, data.frame(
    level = basename(dirname(level_outdir)), comparison = sigla,
    n_up = sum(sig$log2FoldChange > 0), n_down = sum(sig$log2FoldChange < 0),
    total = nrow(sig), stringsAsFactors = FALSE))
  tryCatch(make_deg_heatmap(dds, res, name, level_outdir, group_var = contrast[1]),
           error = function(e) message("    Heatmap ERRORE (", name, "): ", e$message))
}

run_all_comparisons <- function(pb_counts, meta, outdir) {

  # C1: Control vs EcoliAbx (PTL+noPTL fusi)
  tryCatch({
    dds <- make_dds_fused(
      pb_counts, meta,
      groups_to_keep = c("Control", "Ecoli_Abx_noPTL", "Ecoli_Abx_PTL"),
      fuse_map       = list("EcoliAbx" = c("Ecoli_Abx_noPTL", "Ecoli_Abx_PTL"))
    )
    emit(dds, c("group_fused", "EcoliAbx", "Control"),
              "C1_Control_vs_EcoliAbx", outdir)
    message("    C1 completato")
  }, error = function(e) message("    C1 ERRORE: ", e$message))

  # C2: Control vs EcoliAbx_PTL
  tryCatch({
    dds <- make_dds(pb_counts, meta, c("Control", "Ecoli_Abx_PTL"))
    emit(dds, c("group7", "Ecoli_Abx_PTL", "Control"),
              "C2_Control_vs_EcoliAbx_PTL", outdir)
    message("    C2 completato")
  }, error = function(e) message("    C2 ERRORE: ", e$message))

  # C3: Control vs EcoliAbx_noPTL
  tryCatch({
    dds <- make_dds(pb_counts, meta, c("Control", "Ecoli_Abx_noPTL"))
    emit(dds, c("group7", "Ecoli_Abx_noPTL", "Control"),
              "C3_Control_vs_EcoliAbx_noPTL", outdir)
    message("    C3 completato")
  }, error = function(e) message("    C3 ERRORE: ", e$message))

  # C4: EcoliAbx_noPTL vs EcoliAbx_PTL
  tryCatch({
    dds <- make_dds(pb_counts, meta, c("Ecoli_Abx_noPTL", "Ecoli_Abx_PTL"))
    emit(dds, c("group7", "Ecoli_Abx_PTL", "Ecoli_Abx_noPTL"),
              "C4_EcoliAbx_noPTL_vs_PTL", outdir)
    message("    C4 completato")
  }, error = function(e) message("    C4 ERRORE: ", e$message))

  # C5: EcoliAbx vs EcoliAbx_Anakinra (entrambi fusi PTL+noPTL)
  tryCatch({
    dds <- make_dds_fused(
      pb_counts, meta,
      groups_to_keep = c("Ecoli_Abx_noPTL", "Ecoli_Abx_PTL",
                         "Ecoli_Abx_Anak_noPTL", "Ecoli_Abx_Anak_PTL"),
      fuse_map = list(
        "EcoliAbx"          = c("Ecoli_Abx_noPTL", "Ecoli_Abx_PTL"),
        "EcoliAbx_Anakinra" = c("Ecoli_Abx_Anak_noPTL", "Ecoli_Abx_Anak_PTL")
      )
    )
    emit(dds, c("group_fused", "EcoliAbx_Anakinra", "EcoliAbx"),
              "C5_EcoliAbx_vs_EcoliAbxAnak", outdir)
    message("    C5 completato")
  }, error = function(e) message("    C5 ERRORE: ", e$message))

  # C5A: EcoliAbx PTL vs EcoliAbxAnak PTL (diretto, no fusione)
  tryCatch({
    dds <- make_dds(pb_counts, meta, c("Ecoli_Abx_PTL", "Ecoli_Abx_Anak_PTL"))
    emit(dds, c("group7", "Ecoli_Abx_Anak_PTL", "Ecoli_Abx_PTL"),
              "C5A_EcoliAbx_PTL_vs_EcoliAbxAnak_PTL", outdir)
    message("    C5A completato")
  }, error = function(e) message("    C5A ERRORE: ", e$message))

  # C5B: EcoliAbx noPTL vs EcoliAbxAnak noPTL (diretto, no fusione)
  tryCatch({
    dds <- make_dds(pb_counts, meta, c("Ecoli_Abx_noPTL", "Ecoli_Abx_Anak_noPTL"))
    emit(dds, c("group7", "Ecoli_Abx_Anak_noPTL", "Ecoli_Abx_noPTL"),
              "C5B_EcoliAbx_noPTL_vs_EcoliAbxAnak_noPTL", outdir)
    message("    C5B completato")
  }, error = function(e) message("    C5B ERRORE: ", e$message))

  # C6: EcoliAbx_PTL vs EcoliAbx_Anakinra_noPTL (indiretto)
  tryCatch({
    dds <- make_dds(pb_counts, meta, c("Ecoli_Abx_PTL", "Ecoli_Abx_Anak_noPTL"))
    emit(dds, c("group7", "Ecoli_Abx_Anak_noPTL", "Ecoli_Abx_PTL"),
              "C6_EcoliAbx_PTL_vs_EcoliAbxAnak_noPTL", outdir)
    message("    C6 completato")
  }, error = function(e) message("    C6 ERRORE: ", e$message))

  # C7: Control vs EcoliAbx_Anakinra (PTL+noPTL fusi)
  tryCatch({
    dds <- make_dds_fused(
      pb_counts, meta,
      groups_to_keep = c("Control", "Ecoli_Abx_Anak_noPTL", "Ecoli_Abx_Anak_PTL"),
      fuse_map       = list("EcoliAbx_Anakinra" = c("Ecoli_Abx_Anak_noPTL", "Ecoli_Abx_Anak_PTL"))
    )
    emit(dds, c("group_fused", "EcoliAbx_Anakinra", "Control"),
              "C7_Control_vs_EcoliAbxAnak", outdir)
    message("    C7 completato")
  }, error = function(e) message("    C7 ERRORE: ", e$message))

  # C8: Control vs EcoliAbx_Anakinra_PTL
  tryCatch({
    dds <- make_dds(pb_counts, meta, c("Control", "Ecoli_Abx_Anak_PTL"))
    emit(dds, c("group7", "Ecoli_Abx_Anak_PTL", "Control"),
              "C8_Control_vs_EcoliAbxAnak_PTL", outdir)
    message("    C8 completato")
  }, error = function(e) message("    C8 ERRORE: ", e$message))

  # C9: Control vs EcoliAbx_Anakinra_noPTL
  tryCatch({
    dds <- make_dds(pb_counts, meta, c("Control", "Ecoli_Abx_Anak_noPTL"))
    emit(dds, c("group7", "Ecoli_Abx_Anak_noPTL", "Control"),
              "C9_Control_vs_EcoliAbxAnak_noPTL", outdir)
    message("    C9 completato")
  }, error = function(e) message("    C9 ERRORE: ", e$message))

  # C10: EcoliAbxAnak PTL vs EcoliAbxAnak noPTL
  tryCatch({
    dds <- make_dds(pb_counts, meta, c("Ecoli_Abx_Anak_PTL", "Ecoli_Abx_Anak_noPTL"))
    emit(dds, c("group7", "Ecoli_Abx_Anak_PTL", "Ecoli_Abx_Anak_noPTL"),
              "C10_EcoliAbxAnak_PTL_vs_noPTL", outdir)
    message("    C10 completato")
  }, error = function(e) message("    C10 ERRORE: ", e$message))

  # ================================================================
  # CONFRONTI LPS (C11, C13, C14)
  # LPS = riferimento/denominatore -> log2FC positivo = piu' alto
  # nel gruppo di test (Control o E. coli).
  # NOTA: gli LPS sono a endpoint fisso 16h e NON hanno
  # stratificazione PTL; C13/C14 sono quindi asimmetrici.
  # ================================================================

  # C11: Control vs LPS
  tryCatch({
    dds <- make_dds(pb_counts, meta, c("LPS", "Control"))
    emit(dds, c("group7", "Control", "LPS"),
              "C11_LPS_vs_Control", outdir)
    message("    C11 completato")
  }, error = function(e) message("    C11 ERRORE: ", e$message))

  # C13: EcoliAbx_noPTL vs LPS
  tryCatch({
    dds <- make_dds(pb_counts, meta, c("LPS", "Ecoli_Abx_noPTL"))
    emit(dds, c("group7", "Ecoli_Abx_noPTL", "LPS"),
              "C13_LPS_vs_EcoliAbx_noPTL", outdir)
    message("    C13 completato")
  }, error = function(e) message("    C13 ERRORE: ", e$message))

  # C14: EcoliAbx_PTL vs LPS
  tryCatch({
    dds <- make_dds(pb_counts, meta, c("LPS", "Ecoli_Abx_PTL"))
    emit(dds, c("group7", "Ecoli_Abx_PTL", "LPS"),
              "C14_LPS_vs_EcoliAbx_PTL", outdir)
    message("    C14 completato")
  }, error = function(e) message("    C14 ERRORE: ", e$message))
}

# ================================================================
# LOOP SUI LIVELLI (total + per-cluster)
# ================================================================
for (level_name in names(cluster_map)) {

  clusters_target <- cluster_map[[level_name]]

  message("\n================================================================")
  if (is.null(clusters_target)) {
    message("LIVELLO: total (tutte le cellule)")
  } else {
    message("LIVELLO: ", level_name, " (cluster: ", paste(clusters_target, collapse = "+"), ")")
  }
  message("================================================================")

  level_outdir <- file.path(base_outdir, level_name, "Csv")
  dir.create(level_outdir, showWarnings = FALSE, recursive = TRUE)

  # ---- SKIP livello intero se gia' completo -----------------------------
  # save_res() salta solo la scrittura del CSV, ma make_dds()/make_dds_fused()
  # richiamano comunque DESeq(dds) ogni volta (fit pesante, non saltato). Per
  # evitare di rifare inutilmente tutti i fit sui livelli gia' completati in
  # run precedenti, se ci sono gia' tutti e 15 i CSV DE + VST/coldata si salta
  # l'intero livello (pseudobulk, DESeq2, plot compresi).
  existing_de <- list.files(level_outdir, pattern = "^DE_.*\\.csv$")
  vst_dir_chk <- file.path(level_outdir, "VST")
  vst_ok <- file.exists(file.path(vst_dir_chk, paste0("vst_", level_name, ".tsv"))) &&
            file.exists(file.path(vst_dir_chk, paste0("coldata_", level_name, ".tsv")))
  if (length(existing_de) >= 15 && vst_ok) {
    message("  SKIP livello (gia' completo: ", length(existing_de),
            " CSV DE + VST/coldata presenti): ", level_name)
    next
  }

  # ---- Subset nuclei ----
  if (is.null(clusters_target)) {
    cells_lvl <- rownames(obj@meta.data)
  } else {
    cells_lvl <- rownames(obj@meta.data)[obj@meta.data$cluster_char %in% clusters_target]
  }
  message(sprintf("  Nuclei selezionati: %d", length(cells_lvl)))

  if (length(cells_lvl) < 50) {
    message("  SKIP: troppo pochi nuclei (minimo 50)")
    next
  }

  # ---- Pseudobulk ----
  counts_lvl        <- counts_full[, cells_lvl, drop = FALSE]
  cell_samples      <- obj@meta.data[cells_lvl, "sample"]
  sample_factor_lvl <- factor(cell_samples)
  design_mat_lvl     <- model.matrix(~ sample_factor_lvl - 1)
  colnames(design_mat_lvl) <- levels(sample_factor_lvl)
  pb_counts_lvl      <- counts_lvl %*% design_mat_lvl
  message(sprintf("  Pseudobulk: %d geni x %d campioni",
                  nrow(pb_counts_lvl), ncol(pb_counts_lvl)))

  # ---- Allinea metadata ----
  samples_present <- colnames(pb_counts_lvl)
  meta_lvl        <- meta_full[samples_present, , drop = FALSE]

  if (any(is.na(meta_lvl$group7))) {
    message("  Campioni esclusi (LPS/NA): ",
            paste(rownames(meta_lvl)[is.na(meta_lvl$group7)], collapse = ", "))
    meta_lvl       <- meta_lvl[!is.na(meta_lvl$group7), , drop = FALSE]
    pb_counts_lvl  <- pb_counts_lvl[, rownames(meta_lvl), drop = FALSE]
  }

  message("  Campioni nel modello: ", nrow(meta_lvl))
  print(table(meta_lvl$group7))

  message("  Avvio 12 confronti DESeq2...")
  run_all_comparisons(pb_counts_lvl, meta_lvl, level_outdir)

  # ---- Barplot UP/DOWN per confronto (questo livello) ----
  tryCatch({
    if (have_ggplot) {
      dcl <- DEG_COUNTS[DEG_COUNTS$level == level_name, , drop = FALSE]
      if (nrow(dcl) > 0) {
        dcl$comparison <- factor(dcl$comparison, levels = comp_order)
        bar_df <- rbind(
          data.frame(comparison = dcl$comparison, n =  dcl$n_up,   dir = "UP (log2FC>1)"),
          data.frame(comparison = dcl$comparison, n = -dcl$n_down, dir = "DOWN (log2FC<-1)")
        )
        lab <- level_labels[level_name]; if (is.na(lab)) lab <- level_name
        p <- ggplot(bar_df, aes(x = comparison, y = n, fill = dir)) +
          geom_col() + geom_hline(yintercept = 0, color = "black") +
          scale_fill_manual(values = c("UP (log2FC>1)" = "#c0392b",
                                       "DOWN (log2FC<-1)" = "#2c5f8a")) +
          labs(title = lab,
               subtitle = "Significant DEGs per comparison | padj<0.05, |log2FC|>1",
               x = NULL, y = "Number of DEGs", fill = NULL) +
          theme_minimal(base_size = 11) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        bar_dir <- file.path(base_outdir, level_name, "Plots", "Barplot")
        dir.create(bar_dir, recursive = TRUE, showWarnings = FALSE)
        ggsave(file.path(bar_dir, paste0("barplot_", level_name, ".pdf")),
               p, width = 8, height = 5)
        message("  Barplot salvato: ", level_name)
      }
    }
  }, error = function(e) message("  Barplot ERRORE: ", e$message))

  # ---- QC: dispersion plot (fit unico con TUTTI i gruppi di questo livello,
  #      indipendente dai 12 confronti — serve solo a valutare il fit del
  #      modello DESeq2, non produce risultati di DE) ----
  tryCatch({
    dds_qc <- make_dds(pb_counts_lvl, meta_lvl, levels(meta_lvl$group7))
    disp_dir <- file.path(base_outdir, level_name, "Plots", "DESeq2", "Dispersion")
    dir.create(disp_dir, recursive = TRUE, showWarnings = FALSE)
    pdf(file.path(disp_dir, paste0("dispersion_", level_name, ".pdf")), width = 7, height = 6)
    plotDispEsts(dds_qc, main = paste0("Dispersion estimates \u2014 ", level_name))
    dev.off()
    message("  Dispersion plot salvato: ", level_name)
  }, error = function(e) message("  Dispersion plot ERRORE: ", e$message))

  # ---- Export VST + coldata per l'HTML explorer (generate_html.py) ----
  # Matrice VST (tutti i geni, tutti i campioni di questo livello) + coldata
  # sample->condizione, in <livello>/Csv/VST/. Fit DESeq2 indipendente (non
  # riusa dds_qc) per evitare qualsiasi ambiguità se il blocco QC sopra fallisce
  # per questo livello: qui ricalcoliamo da zero, con skip-if-exists come il
  # resto della pipeline.
  tryCatch({
    vst_dir      <- file.path(base_outdir, level_name, "Csv", "VST")
    vst_file     <- file.path(vst_dir, paste0("vst_",     level_name, ".tsv"))
    coldata_file <- file.path(vst_dir, paste0("coldata_", level_name, ".tsv"))
    if (file.exists(vst_file) && file.exists(coldata_file)) {
      message("  SKIP VST/coldata (esistono gia'): ", level_name)
    } else {
      dir.create(vst_dir, recursive = TRUE, showWarnings = FALSE)
      dds_vst  <- make_dds(pb_counts_lvl, meta_lvl, levels(meta_lvl$group7))
      vsd_full <- tryCatch(vst(dds_vst, blind = FALSE),
                           error = function(e) varianceStabilizingTransformation(dds_vst, blind = FALSE))
      vst_mat  <- as.data.frame(assay(vsd_full))
      write.table(vst_mat, vst_file, sep = "\t", quote = FALSE, col.names = NA)

      coldata_df <- data.frame(
        sample    = colnames(vst_mat),
        condition = as.character(colData(dds_vst)$group7),
        stringsAsFactors = FALSE
      )
      write.table(coldata_df, coldata_file, sep = "\t", quote = FALSE, row.names = FALSE)
      message("  VST + coldata salvati: ", level_name,
              " (", nrow(vst_mat), " geni x ", ncol(vst_mat), " campioni)")
    }
  }, error = function(e) message("  VST export ERRORE (", level_name, "): ", e$message))

  # Salva il conteggio reale di nuclei per questo livello, così generate_html.py
  # può mostrarlo in app senza valori hardcoded.
  cellcount_file <- file.path(base_outdir, "cell_counts.csv")
  cc_row <- data.frame(level = level_name, n_cells = length(cells_lvl))
  if (file.exists(cellcount_file)) {
    cc_existing <- read.csv(cellcount_file, stringsAsFactors = FALSE)
    cc_existing <- cc_existing[cc_existing$level != level_name, , drop = FALSE]
    cc_all <- rbind(cc_existing, cc_row)
  } else {
    cc_all <- cc_row
  }
  write.csv(cc_all, cellcount_file, row.names = FALSE)
  message("  Conteggio cellule salvato: ", level_name, " = ", length(cells_lvl))

  message("  === ", level_name, " COMPLETATO ===")
}

# ================================================================
# SUMMARY: heatmap conteggi DEG (livelli x confronti) + barplot immagine 1
# ================================================================
tryCatch({
  summ_dir <- file.path(base_outdir, "Summary")
  dir.create(summ_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(DEG_COUNTS, file.path(summ_dir, "DEG_counts_summary.csv"), row.names = FALSE)

  if (have_pheatmap && nrow(DEG_COUNTS) > 0) {
    m <- xtabs(total ~ level + comparison, data = DEG_COUNTS)
    m <- as.matrix(m)
    cp <- comp_order[comp_order %in% colnames(m)]
    m  <- m[, cp, drop = FALSE]
    rn <- level_labels[rownames(m)]; rn[is.na(rn)] <- rownames(m)[is.na(rn)]
    rownames(m) <- rn
    pal <- colorRampPalette(c("white", "#fddbc7", "#B2182B"))(100)
    ph <- pheatmap(m, color = pal, cluster_rows = (nrow(m) > 2), cluster_cols = FALSE,
                   display_numbers = TRUE, number_format = "%.0f", fontsize_number = 9,
                   main = "Total DEGs per comparison x cell type\n(padj<0.05, |log2FC|>1)",
                   silent = TRUE)
    pdf(file.path(summ_dir, "summary_DEG_heatmap.pdf"),
        width = max(8, 0.7 * ncol(m) + 5), height = max(4, 0.55 * nrow(m) + 3))
    grid::grid.draw(ph$gtable); dev.off()
    message("Summary heatmap salvato in: ", summ_dir)
  }
}, error = function(e) message("Summary heatmap ERRORE: ", e$message))

# ================================================================
# FINE
# ================================================================
message("\n================================================================")
message("PSEUDOBULK UNIFICATO COMPLETATO")
message("Output in: ", base_outdir)
message("  <livello>/Csv/                 -> DE_C1_*.csv ... DE_C14_*.csv (+ C5A, C5B, C11, C13, C14)")
message("  <livello>/Plots/Heatmaps/      -> heatmap_<confronto>.pdf (Top-N DEG, z-score VST)")
message("  <livello>/Plots/Barplot/       -> barplot_<livello>.pdf (UP/DOWN per confronto)")
message("  <livello>/Plots/DESeq2/Dispersion/ -> dispersion_<livello>.pdf (QC)")
message("  <livello>/Csv/VST/            -> vst_<livello>.tsv + coldata_<livello>.tsv (per generate_html.py)")
message("  Summary/                       -> summary_DEG_heatmap.pdf + DEG_counts_summary.csv")
message("  cell_counts.csv                -> nuclei per livello")
message("  --- dettaglio CSV per livello ---")
for (lv in names(cluster_map)) {
  message(sprintf("  %-40s -> DE_C1_*.csv ... DE_C14_*.csv (+ C5A, C5B, C11, C13, C14)",
                  paste0(lv, "/Csv/")))
}
message("================================================================")
