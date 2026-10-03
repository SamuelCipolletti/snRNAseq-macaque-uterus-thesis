#!/usr/bin/env python3
"""
build_demo_explorer.py
 
Builds a demo version of the snRNA-seq DEG explorer from a simulated dataset.
All values are randomly generated (fixed seed, 2026) and none comes from the
thesis data.
 
The script writes a synthetic dataset with the folder layout expected by
generate_html.py: DESeq2-style DE tables, VST matrices with sample metadata,
GO over-representation and GSEA tables, UMAP coordinates, marker dotplot data
and nuclei counts. Gene symbols are real human symbols; their statistics are
random. Cell types (13 generic clusters), conditions (Control and three
treatments, 4 samples each) and comparisons (C1-C6) are generic.
 
generate_html.py is run on a copy in demo_build/_gen, with cell type,
comparison and subtitle labels replaced. The resulting HTML gets a "simulated
data" banner, a noindex tag and a mobile layout, and a QR code pointing to
the public URL is generated.
 
Usage:
    python build_demo_explorer.py <path/to/generate_html.py> [public_URL]
 
Output (demo_build/):
    index.html                    demo page to publish
    qr_demo_explorer.png / .svg   QR code linking to the public URL
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 2026
DEFAULT_URL = "https://samuelcipolletti.github.io/snRNAseq-macaque-uterus-thesis/demo/"

HERE = Path(__file__).resolve().parent
OUT = HERE / "demo_build"
DATA = OUT / "data"
TERM_POOL_DIRS = [Path("/mnt/project"), HERE]   # dove cercare nomi di termini GO (solo ontologia)

rng = np.random.default_rng(SEED)

# ── Universo geni (simboli umani reali, statistiche inventate) ───────────────
def _fam(prefix, items):
    return [f"{prefix}{i}" for i in items]

GENES = sorted(set(
    _fam("RPL", "3 4 5 6 7 7A 8 9 10 10A 11 12 13 13A 14 15 17 18 18A 19 21 22 23 23A 24 26 27 27A 28 29 30 31 32 34 35 35A 36 37 37A 38 39 41".split())
    + _fam("RPS", "2 3 3A 4X 5 6 7 8 9 10 11 12 13 14 15 15A 16 17 18 19 20 21 23 24 25 26 27 27A 28 29".split())
    + "COL1A1 COL1A2 COL3A1 COL4A1 COL4A2 COL5A1 COL5A2 COL6A1 COL6A2 COL6A3 COL12A1 COL14A1 COL15A1 COL18A1".split()
    + _fam("MMP", "1 2 3 7 9 10 11 12 14 19".split())
    + _fam("CXCL", "1 2 3 5 8 10 11 12 14".split())
    + _fam("CCL", "2 3 4 5 8 13 14 19 20 21".split())
    + """IL1A IL1B IL1RN IL1R1 IL6 IL6R IL10 IL15 IL18 IL33 IL7R TNF TNFAIP3 TNFAIP6 TNFRSF1A TNFRSF1B
    NFKBIA NFKB1 RELA RELB KRT5 KRT7 KRT8 KRT18 KRT19 HLA-A HLA-B HLA-C HLA-E HLA-DRA HLA-DRB1 HLA-DPA1
    HLA-DPB1 HLA-DQA1 B2M CD74 PTGS1 PTGS2 PTGER2 PTGER4 PTGES HPGD OXTR GJA1 PGR ESR1 AR NR3C1 SPP1 CD44
    FN1 VIM DCN LUM POSTN SPARC TIMP1 TIMP2 TIMP3 SERPINE1 SERPINF1 THBS1 THBS2 TNC VCAN BGN FBN1 ELN LOX
    IGFBP1 IGFBP2 IGFBP3 IGFBP5 IGFBP7 PRL FOXO1 HAND2 HOXA10 HOXA11 WT1 PDGFRA PDGFRB ACTA2 DES CNN1 MYH11
    TAGLN MYLK LMOD1 PECAM1 VWF CDH5 KDR FLT1 TEK EMCN PROX1 LYVE1 PDPN MMRN1 PTPRC CD68 CD14 CD163 MRC1
    LYZ C1QA C1QB C1QC FCGR3A ITGAM CSF1R CD3E CD3D CD2 CD8A CD4 NKG7 GNLY GZMA GZMB PRF1 KLRD1 MS4A1
    CD79A EPCAM CDH1 MUC1 PAX8 SOX17 FOXJ1 CD24 RGS5 MCAM NOTCH3 KCNJ8 SOCS3 STAT1 STAT3 IRF1 ISG15 IFIT1
    IFIT3 MX1 OAS1 S100A8 S100A9 SAA1 LCN2 CHI3L1 GAPDH ACTB MALAT1 NEAT1 FOS JUN JUNB EGR1 ATF3 HSPA1A
    HSPA1B DUSP1 KLF2 KLF4 ZFP36 IER2 MT1X MT2A MT1E SOD2 VEGFA ANGPT1 ANGPT2 HIF1A EGFR ERBB2 TGFB1 TGFB2
    TGFB3 TGFBR1 TGFBR2 SMAD3 SMAD7 BMP2 BMP4 WNT4 WNT5A CDKN1A MKI67 TOP2A CCNB1 ITGB1 ITGA5 ITGAV ITGB3
    ITGB5 LAMA4 LAMB1 LAMC1 NID1 HSPG2 CCN1 CCN2 CXCR4 CD9 CD63 CD81 ANXA1 ANXA2 LGALS1 LGALS3 CALD1 TPM1
    TPM2 MYL9 FHL2 CRYAB HSPB1 PLAT PLAU PLAUR F3 THBD EDN1 NOS3 APOE APOD CLU CFD C3 C7 PTX3 MGP""".split()
))

# ── Disegno sperimentale GENERICO (non rispecchia lo studio) ─────────────────
# Tipi cellulari, condizioni e confronti sono volutamente neutri: la demo mostra
# come funziona l'explorer, non il dataset della tesi. Per cambiarli basta
# modificare questi tre blocchi.
CELL_TYPES = {   # cluster -> (chiave cartella, etichetta, 5 marker)
    0:  ("epithelial",            "Epithelial cells",            "EPCAM KRT8 KRT18 CDH1 MUC1"),
    1:  ("stromal",               "Stromal cells",               "VIM PDGFRA IGFBP5 IGFBP7 WT1"),
    2:  ("fibroblasts",           "Fibroblasts",                 "DCN LUM COL1A1 COL3A1 POSTN"),
    3:  ("smooth_muscle",         "Smooth muscle cells",         "ACTA2 DES CNN1 MYH11 TAGLN"),
    4:  ("pericytes",             "Pericytes",                   "RGS5 MCAM NOTCH3 KCNJ8 PDGFRB"),
    5:  ("vascular_endothelial",  "Vascular endothelial cells",  "PECAM1 VWF CDH5 KDR EMCN"),
    6:  ("lymphatic_endothelial", "Lymphatic endothelial cells", "PROX1 LYVE1 CCL21 PDPN MMRN1"),
    7:  ("macrophages",           "Macrophages",                 "CD68 CD163 MRC1 C1QA C1QB"),
    8:  ("monocytes",             "Monocytes",                   "CD14 LYZ S100A8 S100A9 FCGR3A"),
    9:  ("dendritic",             "Dendritic cells",             "HLA-DRA HLA-DQA1 CD74 HLA-DPA1 HLA-DPB1"),
    10: ("t_cells",               "T cells",                     "CD3E CD3D CD2 CD4 CD8A"),
    11: ("nk_cells",              "NK cells",                    "NKG7 GNLY GZMB PRF1 KLRD1"),
    12: ("b_cells",               "B cells",                     "MS4A1 CD79A CXCR4 HLA-DRB1 CD24"),
}
CONDITIONS = ["Control", "Treatment_A", "Treatment_B", "Treatment_C"]
SAMPLES_PER_COND = 4
COMPARISON_LABELS = {
    "C1": "C1 — Treatment A vs Control",
    "C2": "C2 — Treatment B vs Control",
    "C3": "C3 — Treatment C vs Control",
    "C4": "C4 — Treatment B vs Treatment A",
    "C5": "C5 — Treatment C vs Treatment A",
    "C6": "C6 — Treatment C vs Treatment B",
}
COMPARISONS = list(COMPARISON_LABELS)
SUBTITLE = "Example dataset — simulated tissue"
GSEA_TERMS = 60   # termini GSEA simulati per confronto e livello (pesa molto sulla dimensione)

N_CLU = len(CELL_TYPES)
CLUSTER_TYPES = {k: v[1] for k, v in CELL_TYPES.items()}
MARKERS = {k: v[2] for k, v in CELL_TYPES.items()}
ANALYSIS_TYPES = {"pseudobulk_total": "Total pseudobulk"}
ANALYSIS_TYPES.update({f"{v[0]}_cluster{k}": f"{v[1]} — cluster {k}" for k, v in CELL_TYPES.items()})
ANALYSES = list(ANALYSIS_TYPES)
MARKERS = {k: v.split() for k, v in MARKERS.items()}


def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1)
    return out


TERM_KEYWORDS = (r"immun|inflamm|matrix|collagen|cytokine|angiogen|muscle|migration|adhesion|leukocyte|"
                 r"chemokine|chemotaxis|interferon|wound|vascul|apoptot|hormone|steroid|"
                 r"prostaglandin|lipopolysaccharide|bacteri|interleukin|tumor necrosis|NF-kappaB|contraction|"
                 r"translation|ribosom|hypoxia|TGF|growth factor|decidua|uter|pregnan|parturition|integrin|"
                 r"endothel|macrophage|\bT cell|defense")


# ── Pool di termini (solo nomi/ID di ontologia, nessun p-value reale) ────────
def term_pool():
    rows = []
    for d in TERM_POOL_DIRS:
        for f in list(d.glob("GO_*.csv")) + list(d.glob("GSEA_*.csv")):
            try:
                rows.append(pd.read_csv(f, usecols=["term_id", "term_name", "source", "term_size"]))
            except Exception:
                pass
    if rows:
        pool = pd.concat(rows).drop_duplicates("term_id").reset_index(drop=True)
        # solo termini plausibili per il tessuto (evita curiosita' tipo "tooth eruption")
        pool = pool[pool["term_name"].str.contains(TERM_KEYWORDS, case=False, regex=True)].reset_index(drop=True)
        if len(pool) > 200:
            return pool
    # ripiego minimo se non ci sono CSV a disposizione
    names = ["inflammatory response", "extracellular matrix organization", "response to lipopolysaccharide",
             "cytokine-mediated signaling pathway", "angiogenesis", "cell adhesion", "wound healing",
             "regulation of cell migration", "response to hypoxia", "leukocyte chemotaxis",
             "collagen fibril organization", "muscle contraction", "translation", "ribosome biogenesis",
             "antigen processing and presentation", "response to interferon-gamma", "apoptotic process",
             "regulation of cell population proliferation", "smooth muscle contraction", "vasculature development"]
    return pd.DataFrame({"term_id": [f"GO:DEMO{i:03d}" for i in range(len(names))],
                         "term_name": names, "source": "GO:BP", "term_size": rng.integers(50, 1500, len(names))})


def enrich_df(pool, n, genes, comp, analysis, direction=None, gsea=False):
    t = pool.sample(n=min(n, len(pool)), random_state=int(rng.integers(1e9))).reset_index(drop=True)
    p = np.array([float(f"{x:.3g}") for x in np.sort(10 ** -rng.uniform(1.4, 9 if not gsea else 14, len(t)))])
    isz = [int(rng.integers(3, max(4, min(len(genes), 25)) + 1)) for _ in range(len(t))]
    df = pd.DataFrame({
        "query": "query_1", "significant": True, "p_value": p, "term_size": t["term_size"],
        "query_size": len(genes), "intersection_size": isz,
        "precision": np.round(np.array(isz) / max(len(genes), 1), 4), "recall": 0.0,
        "term_id": t["term_id"], "source": t["source"], "term_name": t["term_name"],
        "effective_domain_size": 19968, "source_order": 0, "parents": "", "evidence_codes": "",
        "intersection": [",".join(rng.choice(genes, size=min(k, len(genes)), replace=False)) for k in isz],
    })
    if direction:
        df["direction"] = direction
    df["celltype"] = analysis
    df["comparison"] = comp
    return df


# ── Generazione dataset ──────────────────────────────────────────────────────
def build_dataset():
    if DATA.exists():
        shutil.rmtree(DATA)
    DATA.mkdir(parents=True)
    pool = term_pool()
    G = np.array(GENES)
    samples = [f"{c}_{i}" for c in CONDITIONS for i in range(1, SAMPLES_PER_COND + 1)]
    sample_cond = [s.rsplit("_", 1)[0] for s in samples]

    # UMAP: blob gaussiani, uno per cluster
    centers = rng.uniform(-12, 12, size=(N_CLU, 2))
    sizes = rng.integers(150, 900, size=N_CLU)
    um = []
    for k in range(N_CLU):
        xy = centers[k] + rng.normal(0, rng.uniform(0.6, 1.4), size=(sizes[k], 2))
        cond_p = rng.dirichlet(np.ones(len(CONDITIONS)) * 4)
        conds = rng.choice(CONDITIONS, size=sizes[k], p=cond_p)
        samp = [f"{c}_{rng.integers(1, SAMPLES_PER_COND + 1)}" for c in conds]
        um.append(pd.DataFrame({"UMAP_1": xy[:, 0], "UMAP_2": xy[:, 1], "cluster": k,
                                "celltype": CLUSTER_TYPES[k], "condition": conds, "sample": samp}))
    umap = pd.concat(um, ignore_index=True)
    umap.to_csv(DATA / "umap_data.tsv", sep="\t", index=False)

    # Markers dotplot + ranghi
    mgenes = [g for k in range(N_CLU) for g in MARKERS[k]]
    dot = []
    for g in dict.fromkeys(mgenes):
        own = [k for k in range(N_CLU) if g in MARKERS[k]]
        for k in range(N_CLU):
            hi = k in own
            dot.append({"avg.exp": 0, "pct.exp": round(float(rng.uniform(55, 95) if hi else rng.uniform(0, 25)), 2),
                        "features.plot": g, "id": k,
                        "avg.exp.scaled": round(float(rng.uniform(1.2, 2.5) if hi else rng.uniform(-1, 0.4)), 3)})
    pd.DataFrame(dot).to_csv(DATA / "dotplot_data_res05.csv", index=False)
    pd.DataFrame([{"gene": g, "cluster": k, "rank": r + 1}
                  for k in range(N_CLU) for r, g in enumerate(MARKERS[k])]).to_csv(DATA / "top_marker_ranks.csv", index=False)

    # Conteggi nuclei
    counts = {"pseudobulk_total": len(umap)}
    for a in ANALYSES[1:]:
        k = int(a.rsplit("cluster", 1)[1])
        counts[a] = int((umap["cluster"] == k).sum())
    pd.DataFrame({"level": list(counts), "n_cells": list(counts.values())}).to_csv(DATA / "cell_counts.csv", index=False)

    # Per livello: DE, VST, enrichment
    for a in ANALYSES:
        csv = DATA / a / "Csv"; enr = csv / "Enrichment"; vstd = csv / "VST"
        for d in (csv, enr, vstd):
            d.mkdir(parents=True, exist_ok=True)

        base = np.exp(rng.normal(5, 1.6, len(G)))
        for comp in COMPARISONS:
            lfc = rng.normal(0, 0.45, len(G))
            pv = rng.uniform(0, 1, len(G))
            n_de = int(rng.choice([0, 3, 8, 15, 30, 55], p=[.1, .15, .2, .25, .2, .1]))
            if n_de:
                idx = rng.choice(len(G), n_de, replace=False)
                lfc[idx] = rng.choice([-1, 1], n_de) * rng.uniform(1.1, 4, n_de)
                pv[idx] = 10 ** -rng.uniform(4, 12, n_de)
            pv = np.array([float(f"{x:.3g}") for x in pv])        # 3 cifre significative:
            padj = np.array([float(f"{x:.3g}") for x in bh(pv)])  # HTML molto piu' leggero
            de = pd.DataFrame({"baseMean": np.round(base * rng.uniform(.8, 1.2, len(G)), 2),
                               "log2FoldChange": lfc, "lfcSE": np.abs(rng.normal(.3, .1, len(G))),
                               "stat": np.round(lfc / 0.3, 3), "pvalue": pv, "padj": padj}, index=G)
            de.to_csv(csv / f"DE_{comp}_{a}.csv")

            for direc, mask in (("UP", (padj < .05) & (lfc > 1)), ("DOWN", (padj < .05) & (lfc < -1))):
                genes = list(G[mask])
                if len(genes) >= 4:
                    enrich_df(pool, int(rng.integers(5, 30)), genes, comp, a, direc).to_csv(
                        enr / f"GO_{comp}_{a}_{direc}.csv", index=False)
            enrich_df(pool, GSEA_TERMS, list(G), comp, a, gsea=True).to_csv(enr / f"GSEA_{comp}_{a}.csv", index=False)

        # VST: livello base per gene + rumore + piccoli effetti di condizione
        mu = np.log2(base + 1)
        eff = {c: rng.normal(0, .35, len(G)) for c in CONDITIONS}
        mat = np.column_stack([mu + eff[c] + rng.normal(0, .3, len(G)) for c in sample_cond])
        pd.DataFrame(np.round(mat, 3), index=G, columns=samples).to_csv(vstd / f"vst_{a}.tsv", sep="\t")
        pd.DataFrame({"sample": samples, "condition": sample_cond}).to_csv(
            vstd / f"coldata_{a}.tsv", sep="\t", index=False)
    print(f"Dataset simulato: {len(G)} geni, {len(ANALYSES)} livelli, {len(COMPARISONS)} confronti, {len(umap)} nuclei UMAP")


# ── Patch per mobile + banner ────────────────────────────────────────────────
HEAD_INJECT = r"""
<meta name="robots" content="noindex, nofollow">
<style id="demo-mobile">
  #demo-banner { background:var(--warn); color:#fff; font-size:13px; line-height:1.4; padding:7px 14px;
                 text-align:center; font-weight:600; position:fixed; left:0; right:0; bottom:0; z-index:70; }
  body { padding-bottom:34px; }
  #demo-banner span { font-weight:400; opacity:.95; }
  #demo-filters-btn { display:none; }
  @media (max-width: 820px) {
    body { flex-direction:column; padding-bottom:0; }
    body > #demo-banner { position:sticky; top:0; bottom:auto; z-index:60; }
    #sidebar { width:100%; min-width:0; border-right:none; border-bottom:1px solid var(--border2);
               padding:12px 14px; overflow:visible; }
    #sidebar h1 { font-size:19px; padding-right:4px; }
    #sidebar .sub { margin-bottom:10px; }
    #demo-filters-btn { display:block; width:100%; padding:10px; border-radius:8px; font-size:15px;
               font-family:inherit; font-weight:600; background:var(--bg2); color:var(--accent);
               border:1px solid var(--accent); cursor:pointer; }
    #sidebar.demo-collapsed > :not(h1):not(.sub):not(#demo-filters-btn) { display:none !important; }
    select, input[type=text], input[type=number] { font-size:16px; }   /* evita lo zoom automatico su iOS */
    #main { padding:12px; width:100%; min-width:0; }
    .theme-toggle-btn { top:auto; bottom:14px; right:14px; padding:8px 12px; }
    .tabs { overflow-x:auto; flex-wrap:nowrap; -webkit-overflow-scrolling:touch; scrollbar-width:none; }
    .tab { white-space:nowrap; padding:10px 12px; }
    .layout { flex-direction:column; align-items:stretch; }
    .plot-col { width:100%; }
    .side-col { width:100%; min-width:0; position:static; max-height:none; }
    .ov-row-small { flex-direction:column; gap:10px; }
    .vst-grid { grid-template-columns:1fr; }
    .tbl-toolbar, .enrich-toolbar, .heat-controls, .volcano-controls { flex-wrap:wrap; justify-content:flex-start; }
    .heat-controls input[type=range] { width:100%; }
    #deg-table-wrap { overflow-x:auto; }
    #volcano { height:420px !important; }
    #umap-plot { height:620px !important; }
    #markers-plot-wrap { max-height:none !important; }
    .func-tooltip { display:none !important; }
  }
</style>
"""

BODY_INJECT = r"""
<div id="demo-banner">DEMO · simulated data <span>— every value is randomly generated and does not represent study results.</span></div>
"""

SCRIPT_INJECT = r"""
<script>
(function(){
  var sb = document.getElementById('sidebar'); if (!sb) return;
  var btn = document.createElement('button');
  btn.id = 'demo-filters-btn'; btn.type = 'button';
  var sub = sb.querySelector('.sub');
  (sub ? sub.after(btn) : sb.prepend(btn));
  function set(c){ sb.classList.toggle('demo-collapsed', c);
    btn.textContent = c ? 'Show filters' : 'Hide filters';
    window.dispatchEvent(new Event('resize')); }
  btn.onclick = function(){ set(!sb.classList.contains('demo-collapsed')); };
  set(window.matchMedia('(max-width: 820px)').matches);
})();
</script>
"""


PLOTLY_PATCH = r"""
<script>
// Su schermi stretti: legende Plotly orizzontali sotto il grafico, cosi' il
// grafico usa tutta la larghezza invece di cederne meta' alla legenda.
(function(){
  if (!window.Plotly || !window.matchMedia('(max-width: 820px)').matches) return;
  ['react','newPlot'].forEach(function(fn){
    var orig = Plotly[fn];
    Plotly[fn] = function(gd, data, layout, config){
      if (layout && layout.showlegend !== false && !(layout.legend && layout.legend.orientation === 'h')) {
        layout = Object.assign({}, layout, {legend: Object.assign({}, layout.legend || {},
                  {orientation:'h', x:0, xanchor:'left', y:-0.12, yanchor:'top'})});
      }
      var id = (typeof gd === 'string') ? gd : (gd && gd.id);
      if (id === 'enrich-dot' && layout && Array.isArray(data)) {
        // etichette dei termini accorciate: altrimenti occupano tutta la larghezza
        var ys = [];
        data.forEach(function(t){ (t.y || []).forEach(function(v){
          if (typeof v === 'string' && ys.indexOf(v) < 0) ys.push(v); }); });
        if (ys.length) {
          layout = Object.assign({}, layout, {yaxis: Object.assign({}, layout.yaxis || {}, {
            tickmode:'array', tickvals: ys,
            ticktext: ys.map(function(v){ return v.length > 24 ? v.slice(0, 23) + '\u2026' : v; }),
            tickfont: Object.assign({}, (layout.yaxis||{}).tickfont || {}, {size: 10})}),
            margin: Object.assign({}, layout.margin || {}, {l: 150, r: 12})});
        }
      }
      return orig.call(this, gd, data, layout, config);
    };
  });
})();
</script>
"""


def patch_html(html):
    html = html.replace("<title>snRNA-seq Pseudobulk DEG Explorer</title>",
                        "<title>DEMO — snRNA-seq DEG Explorer (simulated data)</title>", 1)
    tag = '<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>'
    html = html.replace(tag, tag + PLOTLY_PATCH, 1)
    html = html.replace("</head>", HEAD_INJECT + "</head>", 1)
    html = html.replace("<body>", "<body>" + BODY_INJECT, 1)
    i = html.rfind("</body>")
    return html[:i] + SCRIPT_INJECT + html[i:]


def _replace_dict(src, name, new):
    """Sostituisce il blocco `NAME = {...}` (fino alla prima riga '}') nella copia."""
    i = src.index(f"{name} = {{")
    j = src.index("\n}\n", i) + 3
    return src[:i] + f"{name} = " + repr(new).replace("', '", "',\n    '") + "\n" + src[j:]


def relabel_copy(path):
    """Mette etichette generiche nella COPIA di lavoro di generate_html.py
    (l'originale resta intatto): tipi cellulari, confronti e sottotitolo."""
    src = path.read_text(encoding="utf-8")
    src = _replace_dict(src, "ANALYSIS_TYPES", ANALYSIS_TYPES)
    src = _replace_dict(src, "COMPARISON_LABELS", COMPARISON_LABELS)
    old_sub = "Macaca mulatta — uterine cell types"
    if old_sub not in src:
        sys.exit("Sottotitolo non trovato in generate_html.py: lo script e' cambiato, aggiorna relabel_copy().")
    src = src.replace(old_sub, SUBTITLE)
    # pulsanti del tab VST che citano i gruppi dello studio
    src = src.replace("'⊕ Fused (ctrl · EcoliAbx · +Anak)'", "'⊕ Merged groups'")
    src = src.replace("'⊕ PTL / noPTL split'", "'⊕ Subgroups'")
    path.write_text(src, encoding="utf-8")


def make_qr(url):
    import qrcode
    import qrcode.image.svg
    qr = qrcode.QRCode(error_correction=qrcode.constants.ERROR_CORRECT_M, box_size=20, border=4)
    qr.add_data(url); qr.make(fit=True)
    qr.make_image(fill_color="black", back_color="white").save(OUT / "qr_demo_explorer.png")
    qr.make_image(image_factory=qrcode.image.svg.SvgPathImage).save(str(OUT / "qr_demo_explorer.svg"))


def main():
    if len(sys.argv) < 2:
        sys.exit("Uso: python3 build_demo_explorer.py <path/generate_html.py> [URL_pubblico]")
    gen = Path(sys.argv[1]).expanduser().resolve()
    url = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_URL
    OUT.mkdir(exist_ok=True)

    build_dataset()

    # copia isolata di generate_html.py (lo script originale non viene toccato)
    work = OUT / "_gen"; work.mkdir(exist_ok=True)
    shutil.copy(gen, work / "generate_html.py")
    relabel_copy(work / "generate_html.py")
    raw = OUT / "_raw.html"
    env = dict(os.environ,
               DEG_ENRICH_DIR=str(DATA), DEG_UMAP_TSV=str(DATA / "umap_data.tsv"),
               DEG_MARKERS_CSV=str(DATA / "dotplot_data_res05.csv"),
               DEG_MARKER_RANKS_CSV=str(DATA / "top_marker_ranks.csv"))
    env.pop("DEG_DATA_DIR", None)
    subprocess.run([sys.executable, str(work / "generate_html.py"), str(DATA), str(raw)], check=True, env=env)

    html = patch_html(raw.read_text(encoding="utf-8"))
    (OUT / "index.html").write_text(html, encoding="utf-8")
    raw.unlink()
    make_qr(url)
    mb = (OUT / "index.html").stat().st_size / 1024 / 1024
    print(f"\nindex.html: {mb:.1f} MB  |  QR -> {url}")


if __name__ == "__main__":
    main()
