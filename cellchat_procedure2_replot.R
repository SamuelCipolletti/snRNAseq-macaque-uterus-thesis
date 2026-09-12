#!/usr/bin/env Rscript
###############################################################################
# CellChat, Procedure 2 (replot): regenerates the plots only
#
# Reloads the merged object saved by cellchat_procedure2_res05.R for each
# comparison and redraws a more readable set of figures, without repeating the
# merge, the DEG analysis or any other computation. The ligand-receptor pairs
# come from the CSV files already written by that script.
#
# Figures per comparison: differential heatmaps with a cluster-to-cell-type
# legend page; a multipage differential circle plot with one page per cluster,
# shown both as source and as target; and multipage bubble plots of the top
# up- and down-regulated ligand-receptor pairs, one page per sender cluster.
# Plots keep the cluster numbers; cell-type names appear only in titles and
# legends.
###############################################################################

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(ComplexHeatmap)
  library(grid)
  library(scales)   # trans/breaks per la scala colore del bubble plot
})

## ---------------------------------------------------------------------------
## Parametri (coerenti con cellchat_procedure2_res05.R)
## ---------------------------------------------------------------------------
OUT_DIR      <- "/hpc/home/samuelevan.cipolletti/snRNAseq/CellChat/res_05/Comparisons"
TOP_N_BUBBLE <- 40     # n. max coppie L-R nel bubble (top per |logFC|); NA = tutte

# confronti da rigenerare: NULL = tutti quelli presenti in OUT_DIR con merged_*.rds
CMP_ONLY <- NULL       # es. c("C1","C8") per rifarne solo alcuni

## Mappatura cluster -> cell type (Res 0.5, from the consensus annotation;
## see thesis, Methods). I GRAFICI restano coi numeri di cluster;
## questa mappa serve SOLO per la legenda affiancata al circle plot.
CELLTYPE_MAP <- c(
  "0"  = "decidual (glandular/luminal) epithelial cells",
  "1"  = "lymphatic endothelial cells",
  "2"  = "(vascular) endothelial cells",
  "3"  = "decidual stromal (and perivascular) cells",
  "4"  = "macrophages (possible mixture of myeloid cells)",
  "5"  = "macrophages (and monocytes; also Hofbauer cells)",
  "6"  = "decidual (glandular/luminal) epithelial cells",
  "7"  = "decidual stromal (and epithelial cells)",
  "8"  = "mixture of lymphoid cells",
  "9"  = "decidual (luminal/glandular) epithelial cells",
  "10" = "(myometrial) smooth muscle cells",
  "11" = "decidual (glandular/luminal) epithelial cells",
  "12" = "decidual (glandular/luminal) epithelial cells; possible mixture of non-immune cells",
  "13" = "macrophages (possible mixture of myeloid cells)",
  "14" = "(vascular) endothelial cells",
  "15" = "perivascular cells",
  "16" = "decidual stromal cells (and fibroblasts)"
)

## ---------------------------------------------------------------------------
## Helper: seleziona le top-N coppie L-R per |logFC| da un data.frame CSV
## ---------------------------------------------------------------------------
top_pairs <- function(df, n = TOP_N_BUBBLE) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if (!"interaction_name" %in% colnames(df)) return(NULL)
  lfc_col <- intersect(c("ligand.logFC", "receptor.logFC"), colnames(df))
  if (length(lfc_col) >= 1 && !is.na(n)) {
    ord <- order(abs(df[[lfc_col[1]]]), decreasing = TRUE)
    df  <- df[ord, , drop = FALSE]
  }
  pr <- unique(df[["interaction_name"]])
  if (!is.na(n) && length(pr) > n) pr <- pr[seq_len(n)]
  data.frame(interaction_name = pr, stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## Bubble plot custom (leggibilita' migliorata)
##
## Ricostruisce il ggplot da netVisual_bubble(..., return.data = TRUE) per
## avere pieno controllo su etichette, scala colore e font.
## $communication contiene: source.target (asse X "sender -> receiver (dataset)"),
##   interaction_name_2 (asse Y "Ligand - Receptor"), prob (Commun.Prob.), pval.
## ---------------------------------------------------------------------------

# accorcia le etichette X: se il sender e' unico (caso sources.use = i singolo)
# rimuove il "senderN -> " ridondante, lasciando "-> receiver (dataset)".
shorten_xlabels <- function(x) {
  senders <- sub("^(.*?) -> .*$", "\\1", x)
  if (length(unique(senders)) == 1L) {
    return(sub("^.*? -> ", "\u2192 ", x))
  }
  x
}

build_bubble <- function(cellchat, pairLR.use, comparison, sources.use,
                         title.name) {
  bd <- tryCatch(
    netVisual_bubble(cellchat, pairLR.use = pairLR.use, comparison = comparison,
                     sources.use = sources.use, remove.isolate = TRUE,
                     return.data = TRUE),
    error = function(e) NULL)
  if (is.null(bd)) return(NULL)

  df <- bd$communication
  if (is.null(df) || nrow(df) == 0) return(NULL)

  df$source.target      <- factor(df$source.target,
                                  levels = unique(df$source.target))
  df$interaction_name_2 <- factor(df$interaction_name_2,
                                  levels = rev(unique(df$interaction_name_2)))

  xlev   <- levels(df$source.target)
  xshort <- shorten_xlabels(xlev)

  has_psize <- "pval" %in% colnames(df) && length(unique(df$pval)) > 1

  p <- ggplot(df, aes(x = source.target, y = interaction_name_2))
  if (has_psize) {
    p <- p + geom_point(aes(color = prob, size = factor(pval)))
    p <- p + scale_size_manual(
      name   = "p-value",
      values = c("1" = 4, "2" = 2.6, "3" = 1.4),
      labels = c("1" = "< 0.01", "2" = "< 0.05", "3" = "ns"),
      drop   = TRUE)
  } else {
    p <- p + geom_point(aes(color = prob), size = 2.8)
  }

  p +
    scale_x_discrete(labels = xshort) +
    scale_color_viridis_c(
      option = "plasma",
      name   = "Commun.\nProb.",
      trans  = "sqrt",                          # espande i valori bassi
      breaks = scales::pretty_breaks(n = 4),
      guide  = guide_colorbar(barheight = grid::unit(3, "cm"),
                              barwidth  = grid::unit(0.35, "cm"))
    ) +
    labs(title = title.name) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
      axis.text.y  = element_text(size = 9, face = "italic"),  # L-R in corsivo
      axis.title   = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title   = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      legend.key.size = grid::unit(0.4, "cm")
    )
}

## ---------------------------------------------------------------------------
## Rigenera i grafici per un singolo confronto
## ---------------------------------------------------------------------------
replot_one <- function(cmp_name) {
  cmp_dir <- file.path(OUT_DIR, cmp_name)
  rds     <- list.files(cmp_dir, pattern = "^merged_.*\\.rds$", full.names = TRUE)
  if (length(rds) == 0) {
    cat(sprintf("[SKIP] %s : nessun merged_*.rds\n", cmp_name)); return(invisible())
  }
  cat(sprintf("\n===== REPLOT %s =====\n", cmp_name))
  cellchat <- readRDS(rds[1])

  ## --- setup label / colori / nomenclatura cluster = cell type ---
  clabs <- levels(cellchat@idents$joint)
  if (is.null(clabs)) clabs <- levels(cellchat@idents[[1]])
  ncl   <- length(clabs)
  ccols <- tryCatch(scPalette(ncl), error = function(e) rainbow(ncl))

  # estrae il numero del cluster da ciascuna label (es. "cluster3" -> "3")
  cnum  <- gsub("[^0-9]", "", clabs)
  ctype <- CELLTYPE_MAP[cnum]
  ctype[is.na(ctype)] <- "(cell type non mappato)"
  leg_txt <- sprintf("%s = %s", clabs, ctype)

  # pannello riutilizzabile: legenda nomenclatura cluster = cell type
  draw_legend_panel <- function() {
    par(mar = c(1, 0, 2, 0))
    plot.new()
    title("Cluster = cell type (Res 0.5)", cex.main = 1.0, adj = 0)
    legend("center", legend = leg_txt, col = ccols, pch = 19, pt.cex = 1.4,
           bty = "n", ncol = 1, cex = 0.75, y.intersp = 1.15, text.width = 0)
  }

  ## --- (0) heatmap differenziale + pannello nomenclatura ---
  # netVisual_heatmap restituisce oggetti ComplexHeatmap: li affianchiamo e
  # aggiungiamo a destra un pannello di testo con cluster = cell type.
  legend_grob <- function() {
    grid::grid.newpage()
    n <- length(leg_txt)
    grid::grid.text("Cluster = cell type (Res 0.5)", x = 0.02, y = 0.98,
                    just = c("left", "top"), gp = grid::gpar(fontface = "bold", cex = 0.85))
    ys <- seq(0.93, 0.05, length.out = n)
    for (k in seq_len(n)) {
      grid::grid.points(x = grid::unit(0.04, "npc"), y = grid::unit(ys[k], "npc"),
                        pch = 19, size = grid::unit(0.6, "char"),
                        gp = grid::gpar(col = ccols[k]))
      grid::grid.text(leg_txt[k], x = 0.08, y = ys[k], just = c("left", "center"),
                      gp = grid::gpar(cex = 0.62))
    }
  }
  ht_ok <- tryCatch({
    ht1 <- netVisual_heatmap(cellchat)
    ht2 <- netVisual_heatmap(cellchat, measure = "weight")
    pdf(file.path(cmp_dir, "02_diff_heatmap.pdf"), width = 18, height = 8)
    # pagina 1: le due heatmap affiancate
    ComplexHeatmap::draw(ht1 + ht2, ht_gap = grid::unit(1, "cm"))
    # pagina 2: nomenclatura cluster = cell type
    legend_grob()
    dev.off()
    TRUE
  }, error = function(e) { cat(sprintf("  [WARN heatmap] %s\n", conditionMessage(e))); FALSE })
  if (ht_ok) cat("  [OK] 02_diff_heatmap.pdf (heatmap + pagina nomenclatura)\n")

  ## --- (A) circle plot differenziale: PDF MULTIPAGINA, 1 pagina per cluster ---
  # ogni pagina mostra gli archi che coinvolgono il cluster di riferimento:
  #   pannello sx = il cluster come SORGENTE verso tutti (i -> *)
  #   pannello dx = il cluster come BERSAGLIO da tutti  (* -> i)
  # entrambi come differenza del N. di interazioni. Pannello nomenclatura a dx.
  pdf(file.path(cmp_dir, "02_diff_circle.pdf"), width = 16, height = 7)
  for (i in seq_len(ncl)) {
    layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.9))
    par(xpd = TRUE)
    tryCatch(
      netVisual_diffInteraction(cellchat, weight.scale = TRUE,
             sources.use = i,
             title.name = sprintf("%s come SORGENTE (n. interazioni)", clabs[i])),
      error = function(e) { plot.new(); title(sprintf("%s (sorgente): n/d", clabs[i])) })
    tryCatch(
      netVisual_diffInteraction(cellchat, weight.scale = TRUE,
             targets.use = i,
             title.name = sprintf("%s come BERSAGLIO (n. interazioni)", clabs[i])),
      error = function(e) { plot.new(); title(sprintf("%s (bersaglio): n/d", clabs[i])) })
    draw_legend_panel()
  }
  dev.off()
  cat("  [OK] 02_diff_circle.pdf (multipagina, 1 pagina/cluster)\n")

  ## --- bubble plot up/down (top-N per |logFC|) dai CSV ---
  test_name <- tail(names(cellchat@net), 1)   # dataset2 = test
  f_up   <- file.path(cmp_dir, "04_LR_up_in_test.csv")
  f_down <- file.path(cmp_dir, "04_LR_down_in_test.csv")

  ## --- (B) bubble plot: PDF MULTIPAGINA, 1 pagina per cluster-sorgente ---
  # ogni pagina: cluster i come sender verso TUTTI i target; righe = top-N L-R.
  do_bubble <- function(csv, tag, out_pdf) {
    if (!file.exists(csv)) { cat(sprintf("  [skip bubble %s] manca %s\n", tag, basename(csv))); return() }
    df <- read.csv(csv, stringsAsFactors = FALSE)
    pr <- top_pairs(df)
    if (is.null(pr) || nrow(pr) == 0) { cat(sprintf("  [skip bubble %s] nessuna coppia\n", tag)); return() }
    plots <- list()
    for (i in seq_len(ncl)) {
      gg <- build_bubble(
        cellchat, pairLR.use = pr, comparison = c(1, 2), sources.use = i,
        title.name = sprintf("%s in %s | sender: %s = %s (top %d L-R)",
                             tag, test_name, clabs[i], ctype[i], nrow(pr)))
      if (!is.null(gg)) plots[[length(plots) + 1]] <- gg
    }
    if (length(plots) == 0) { cat(sprintf("  [skip bubble %s] nessuna pagina valida\n", tag)); return() }
    pdf(out_pdf, width = 14, height = 10)
    for (p in plots) print(p)
    dev.off()
    cat(sprintf("  [OK] %s (multipagina, %d pagine)\n", basename(out_pdf), length(plots)))
  }

  do_bubble(f_up,   "Up",   file.path(cmp_dir, "05_bubble_up.pdf"))
  do_bubble(f_down, "Down", file.path(cmp_dir, "05_bubble_down.pdf"))

  rm(cellchat); gc()
  invisible()
}

## ---------------------------------------------------------------------------
## Main
## ---------------------------------------------------------------------------
if (is.null(CMP_ONLY)) {
  cmps <- list.dirs(OUT_DIR, recursive = FALSE, full.names = FALSE)
  cmps <- cmps[grepl("^C", cmps)]
} else {
  cmps <- CMP_ONLY
}
cat(sprintf("Confronti da rigenerare: %s\n", paste(cmps, collapse = ", ")))
for (cn in cmps) replot_one(cn)
cat("\n==== REPLOT COMPLETATO ====\n")
