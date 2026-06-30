"""
CardioIsoMap — gene-searchable cardiac long-read isoform atlas.

Run locally:  streamlit run app/app.py
Data:         data/processed/*.parquet  (built by code/01_build_tables.py and 02_build_structures.py)

NOTE: abundance is derived from per-sample counts.
"""
from __future__ import annotations

import io
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# --------------------------------------------------------------------------- #
# Config / constants
# --------------------------------------------------------------------------- #
import re

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
SASHIMI_DIR = PROC / "sashimi"


def safe_name(g: str) -> str:
    """Match the filename sanitisation used by code/03_render_sashimi.R."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(g))


def _img_base() -> str:
    """Base URL for remotely-hosted sashimi PNGs (jsDelivr CDN).

    Set via Streamlit secret `IMG_BASE` (deployed) or env var
    CARDIOISOMAP_IMG_BASE (local). Empty -> serve local files in data/processed/sashimi.
    """
    try:
        val = st.secrets.get("IMG_BASE", "")
    except Exception:
        val = ""
    return (val or os.environ.get("CARDIOISOMAP_IMG_BASE", "")).rstrip("/")


IMG_BASE = _img_base()

SAMPLES = ["sample1311", "sample1518", "sample1532", "sample1535", "sample1561", "sample1662"]
SAMPLE_PHENO = {
    "sample1561": "NFH", "sample1662": "NFH",
    "sample1532": "DCM", "sample1535": "DCM",
    "sample1311": "IHD", "sample1518": "IHD",
}
PHENOTYPES = ["NFH", "DCM", "IHD"]
PHENO_LABEL = {"NFH": "NFH (non-failing)", "DCM": "DCM", "IHD": "IHD"}
PHENO_COLOR = {"NFH (non-failing)": "#ef767a", "DCM": "#456990", "IHD": "#49beaa"}
# MAGNET (short-read) etiology groups
MAG_ORDER = ["Non-Failing", "DCM", "HCM", "PPCM"]
MAG_COLOR = {"Non-Failing": "#2196F3", "DCM": "#F44336",
             "HCM": "#FF9800", "PPCM": "#9C27B0"}

st.set_page_config(page_title="CardioIsoMap", page_icon="🫀", layout="wide")


# --------------------------------------------------------------------------- #
# Data loading (cached)
# --------------------------------------------------------------------------- #
@st.cache_data(show_spinner=False)
def load_isoforms() -> pd.DataFrame:
    return pd.read_parquet(PROC / "isoforms.parquet")


@st.cache_data(show_spinner=False)
def load_genes() -> pd.DataFrame:
    return pd.read_parquet(PROC / "genes.parquet")


@st.cache_data(show_spinner=False)
def load_exons() -> pd.DataFrame:
    return pd.read_parquet(PROC / "exons.parquet")


@st.cache_data(show_spinner=False)
def load_magnet_tpm() -> pd.DataFrame:
    """MAGNET short-read isoform TPM (wide: transcript_id + 366 sample cols)."""
    p = PROC / "magnet_isoform_tpm.parquet"
    return pd.read_parquet(p) if p.exists() else pd.DataFrame()


@st.cache_data(show_spinner=False)
def load_magnet_samples() -> pd.DataFrame:
    """MAGNET sample -> etiology (disease group)."""
    p = PROC / "magnet_samples.csv"
    return pd.read_csv(p) if p.exists() else pd.DataFrame(columns=["Run", "etiology"])


@st.cache_data(show_spinner=False)
def load_magnet_enst_tpm() -> pd.DataFrame:
    """MAGNET TPM for reference ENSTs matched by FSM/ISM isoforms (ref_enst + samples)."""
    p = PROC / "magnet_enst_tpm.parquet"
    return pd.read_parquet(p) if p.exists() else pd.DataFrame()


@st.cache_data(show_spinner=False)
def load_isoform_ref_map() -> pd.DataFrame:
    """Known isoform (FSM/ISM) -> matched reference ENST."""
    p = PROC / "isoform_ref_map.csv"
    return pd.read_csv(p) if p.exists() else pd.DataFrame(
        columns=["isoform", "ref_enst", "structural_category"])


@st.cache_data(show_spinner=False)
def load_gencode_exons() -> pd.DataFrame:
    """GENCODE v43 reference exon models (empty frame if not built)."""
    p = PROC / "gencode_exons.parquet"
    if not p.exists():
        return pd.DataFrame(columns=["gene", "transcript_id", "transcript_name",
                                     "chrom", "strand", "start", "end"])
    return pd.read_parquet(p)


@st.cache_data(show_spinner=False)
def gene_list() -> list[str]:
    g = pd.read_csv(PROC / "gene_index.csv")
    return sorted(g["gene"].dropna().astype(str).tolist())


@st.cache_data(show_spinner=False)
def available_sashimi() -> set[str]:
    """Genes that have a rendered coverage/sashimi PNG (from the manifest)."""
    mf = PROC / "sashimi_manifest.csv"
    if not mf.exists():
        return set()
    return set(pd.read_csv(mf)["gene"].astype(str))


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def external_links(gene: str) -> str:
    return (
        f"[Ensembl](https://www.ensembl.org/Multi/Search/Results?q={gene}) · "
        f"[NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/?term={gene}) · "
        f"[GeneCards](https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}) · "
        f"[GTEx](https://gtexportal.org/home/gene/{gene})"
    )


def usage_stack_figure(sub: pd.DataFrame, top_n: int, grouping: str):
    """Stacked horizontal bars of isoform usage composition per phenotype (or sample).

    Each bar sums to 100%: it shows what fraction of the gene's reads each isoform
    accounts for, so an isoform that takes over in DCM/IHD is an isoform switch.
    Built from raw counts within the displayed isoform set.
    """
    if grouping == "Per sample":
        cols = SAMPLES
        labels = {s: f"{s} ({SAMPLE_PHENO[s]})" for s in cols}
    else:
        cols = [f"count_{p}" for p in PHENOTYPES]
        labels = {f"count_{p}": PHENO_LABEL[p] for p in PHENOTYPES}
    cols = [c for c in cols if c in sub.columns]
    if sub.empty or not cols or sub[cols].to_numpy().sum() == 0:
        return None

    d = sub[["isoform"] + cols].copy()
    d["_tot"] = d[cols].sum(axis=1)
    d = d.sort_values("_tot", ascending=False)
    top_ids = d["isoform"].head(top_n).tolist()
    d["grp"] = d["isoform"].where(d["isoform"].isin(top_ids), "Other")
    agg = d.groupby("grp", sort=False)[cols].sum()
    prop = agg.div(agg.sum(axis=0).replace(0, np.nan), axis=1).fillna(0)

    long = (prop.reset_index().melt(id_vars="grp", var_name="col", value_name="usage"))
    long["group"] = long["col"].map(labels)
    iso_order = top_ids + (["Other"] if "Other" in agg.index else [])
    # phenotype/sample order, top of chart = first listed
    grp_axis = [labels[c] for c in cols][::-1]

    palette = (px.colors.qualitative.Set3 + px.colors.qualitative.Pastel)
    cmap = {iso: palette[i % len(palette)] for i, iso in enumerate(top_ids)}
    cmap["Other"] = "#d9d9d9"

    fig = px.bar(long, x="usage", y="group", color="grp", orientation="h",
                 category_orders={"group": grp_axis, "grp": iso_order},
                 color_discrete_map=cmap)
    fig.update_layout(
        barmode="stack", height=110 + 64 * len(cols),
        margin=dict(l=10, r=10, t=10, b=10),
        font=dict(size=16),
        xaxis=dict(tickformat=".0%", title="Isoform usage (share of gene reads)",
                   range=[0, 1], tickfont=dict(size=15), title_font=dict(size=16)),
        yaxis_title=None, yaxis=dict(tickfont=dict(size=15)),
        legend=dict(title_text="Isoform", traceorder="normal", font=dict(size=14)))
    fig.update_traces(hovertemplate="%{y}<br>%{fullData.name}: %{x:.1%}<extra></extra>")
    return fig


def _magnet_rows(sub: pd.DataFrame):
    """Per-isoform MAGNET TPM rows for a gene, from two sources:
       - novel isoforms: their own PB.x.y row in the salmon matrix
       - known isoforms (FSM/ISM): the TPM of their matched GENCODE ENST
    Returns (rows DataFrame indexed by isoform, dict isoform->source-label).
    """
    tpm = load_magnet_tpm()
    if tpm.empty:
        return pd.DataFrame(), {}
    samples = [c for c in tpm.columns if c != "transcript_id"]
    want = set(sub["isoform"])
    parts, src = [], {}

    # novel isoforms — direct PB rows
    nov = tpm[tpm["transcript_id"].isin(want)].set_index("transcript_id")
    for i in nov.index:
        src[i] = "novel (direct)"
    if not nov.empty:
        parts.append(nov[samples])

    # known isoforms — via FSM/ISM -> reference ENST
    rmap = load_isoform_ref_map()
    enst = load_magnet_enst_tpm()
    remaining = want - set(nov.index)
    if not rmap.empty and not enst.empty and remaining:
        rm = rmap[rmap["isoform"].isin(remaining)]
        eidx = enst.set_index("ref_enst")
        rm = rm[rm["ref_enst"].isin(eidx.index)]
        if not rm.empty:
            mapped = eidx.loc[rm["ref_enst"].values, samples].copy()
            mapped.index = rm["isoform"].values
            for iso, cat in zip(rm["isoform"], rm["structural_category"]):
                src[iso] = "FSM (matched ENST)" if cat == "full-splice_match" \
                    else "ISM (matched ENST)"
            parts.append(mapped)

    if not parts:
        return pd.DataFrame(), {}
    rows = pd.concat(parts)
    rows = rows[~rows.index.duplicated()]
    rows.index.name = "isoform"
    return rows, src


def magnet_box_figure(sub: pd.DataFrame, top_n: int, groups: list[str], log_y: bool):
    """Box plot of MAGNET short-read TPM per isoform, split by disease group.

    Returns (figure, n_isoforms_plotted). Covers novel isoforms (own salmon row) plus
    FSM/ISM isoforms mapped to their matched GENCODE transcript's TPM.
    """
    samp = load_magnet_samples()
    rows, src = _magnet_rows(sub)
    if rows.empty or samp.empty:
        return None, 0
    order = rows.mean(axis=1).sort_values(ascending=False).head(top_n).index.tolist()
    long = (rows.loc[order].reset_index()
            .melt(id_vars="isoform", var_name="Run", value_name="TPM"))
    long = long.merge(samp, on="Run", how="left")
    long = long[long["etiology"].isin(groups)]
    if long.empty:
        return None, len(order)
    long["source"] = long["isoform"].map(src)
    fig = px.box(
        long, x="isoform", y="TPM", color="etiology",
        category_orders={"isoform": order,
                         "etiology": [g for g in MAG_ORDER if g in groups]},
        color_discrete_map=MAG_COLOR, points="all",
        hover_data={"source": True})
    # one jittered dot per MAGNET sample; boxes show no separate outliers
    fig.update_traces(jitter=0.35, pointpos=0, boxpoints="all",
                      marker=dict(size=6, opacity=0.6),
                      selector=dict(type="box"))
    fig.update_layout(
        height=max(360, 32 * len(order) + 170),
        boxmode="group", margin=dict(l=10, r=10, t=30, b=90),
        font=dict(size=16),
        xaxis_title=None, yaxis_title="TPM (log)" if log_y else "TPM",
        yaxis=dict(title_font=dict(size=16)),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, title_text="",
                    font=dict(size=14)))
    fig.update_xaxes(tickangle=-40, tickfont=dict(size=14))
    fig.update_yaxes(tickfont=dict(size=14))
    if log_y:
        fig.update_yaxes(type="log")
    return fig, len(order)


# Canonical per-structural-category colors — shared by the per-gene overview charts AND the
# Structure figure so the two always match. (Excluded isoforms are drawn as a lighter tint.)
PIE_CAT_COLORS = {
    "full-splice_match":       "#9b1c31",  # FSM  (dark red)
    "incomplete-splice_match": "#f6d3c9",  # ISM  (pale peach)
    "novel_in_catalog":        "#5fa78f",  # NIC  (medium green)
    "novel_not_in_catalog":    "#0f7a4e",  # NNC  (dark green)
}
PIE_CAT_OTHER = "#7fbf9a"  # Others (fusion / antisense / intergenic / genic …)

# Structure figure reuses the same palette (single source of truth).
CAT_COLORS = PIE_CAT_COLORS
CAT_OTHER = PIE_CAT_OTHER
REF_COLOR = "#9e9e9e"  # GENCODE v43 reference transcripts


def _lighten(hex_color: str, frac: float) -> str:
    """Blend a hex color toward white by `frac` (0=unchanged, 1=white)."""
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    r = int(r + (255 - r) * frac)
    g = int(g + (255 - g) * frac)
    b = int(b + (255 - b) * frac)
    return f"#{r:02x}{g:02x}{b:02x}"


def exon_color(category: str, keep: bool) -> str:
    if category == "GENCODE":
        return REF_COLOR
    base = CAT_COLORS.get(category, CAT_OTHER)
    return base if keep else _lighten(base, 0.55)  # excluded → even lighter


def structure_figure(exons: pd.DataFrame, isoforms: list[str], gene: str,
                     meta: dict | None = None):
    """Draw exon/intron transcript models (genomic coordinates).

    `meta` maps transcript_id -> (structural_category, keep_bool); used to color
    exons by SQANTI category, with excluded isoforms drawn in a lighter shade.
    """
    meta = meta or {}
    ex = exons[exons["transcript_id"].isin(isoforms)].copy()
    if ex.empty:
        return None
    gmin, gmax = ex["start"].min(), ex["end"].max()
    strand = ex["strand"].iloc[0]
    order = isoforms[::-1]  # first isoform on top

    fig, ax = plt.subplots(figsize=(7, max(1.0, 0.28 * len(order) + 0.6)))
    cats_present, has_excluded = set(), False
    for i, tid in enumerate(order):
        te = ex[ex["transcript_id"] == tid].sort_values("start")
        if te.empty:
            continue
        cat, keep = meta.get(tid, (None, True))
        cats_present.add(cat)
        has_excluded = has_excluded or (not keep)
        fill = exon_color(cat, keep)
        tmin, tmax = te["start"].min(), te["end"].max()
        ax.plot([tmin, tmax], [i, i], color="#bbb" if not keep else "#888",
                lw=1, zorder=1)  # intron line
        for _, r in te.iterrows():
            ax.add_patch(plt.Rectangle((r["start"], i - 0.3), r["end"] - r["start"], 0.6,
                                       facecolor=fill, edgecolor="#555", linewidth=0.4,
                                       zorder=2))
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=7)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlim(gmin - (gmax - gmin) * 0.02, gmax + (gmax - gmin) * 0.02)
    ax.set_xlabel(f"{ex['chrom'].iloc[0]} position (bp)   strand: {strand}")
    ax.set_title(f"{gene} — transcript structures")
    ax.spines[["top", "right", "left"]].set_visible(False)

    # Legend: one swatch per category present, plus a note on the lighter shade.
    from matplotlib.patches import Patch
    label = {"full-splice_match": "FSM", "incomplete-splice_match": "ISM",
             "novel_in_catalog": "NIC", "novel_not_in_catalog": "NNC"}
    handles = [Patch(facecolor=CAT_COLORS.get(c, CAT_OTHER), edgecolor="#555",
                     label=label.get(c, str(c)))
               for c in CAT_COLORS if c in cats_present]
    if "GENCODE" in cats_present:
        handles.append(Patch(facecolor=REF_COLOR, edgecolor="#555",
                             label="GENCODE v43"))
    if has_excluded:
        handles.append(Patch(facecolor="#e6e6e6", edgecolor="#555",
                             label="lighter = excluded"))
    if handles:
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0),
                  fontsize=7, frameon=False)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------- #
# Sidebar
# --------------------------------------------------------------------------- #
genes = gene_list()
iso_all = load_isoforms()
genes_tbl = load_genes()

st.sidebar.title("🫀 CardioIsoMap")
st.sidebar.caption("Cardiac long-read isoform atlas")

gene = st.sidebar.selectbox(
    "🔍 Search a gene",
    options=["—"] + genes,
    index=0,
    help="Type a gene symbol, e.g. MYH7, TTN, NPPA",
)

st.sidebar.markdown("**Filters**")
keep_only = st.sidebar.toggle("Curated isoforms only (keep)", value=True)
novel_only = st.sidebar.toggle("Novel isoforms only", value=False)
cats = sorted(iso_all["structural_category"].dropna().unique().tolist())
sel_cats = st.sidebar.multiselect("Structural category", cats, default=[])

st.sidebar.markdown("---")
st.sidebar.markdown(
    "**About** — PacBio long-read isoforms across human heart "
    "(NFH, DCM, IHD). Abundance from per-sample counts."
)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def filter_gene(df: pd.DataFrame) -> pd.DataFrame:
    if keep_only:
        df = df[df["keep"]]
    if novel_only:
        df = df[df["is_novel"]]
    if sel_cats:
        df = df[df["structural_category"].isin(sel_cats)]
    return df


if gene == "—":
    # Landing / dataset-summary page
    st.title("🫀 CardioIsoMap")
    st.subheader("A full-length isoform catalog of the healthy and failing human heart")

    kept = iso_all[iso_all["keep"]]
    n_iso = len(iso_all)
    n_keep = int(iso_all["keep"].sum())
    n_excl = n_iso - n_keep
    # "Novel" = ISM + NIC + NNC, not the curated novel_lr set.
    novel_cat = iso_all["structural_category"].isin(
        ["incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"])
    n_novel_all = int(novel_cat.sum())
    n_novel_keep = int((novel_cat & iso_all["keep"]).sum())
    n_novel_verified = int(iso_all["is_novel"].sum())  # curated/verified subset
    n_genes = iso_all["associated_gene"].nunique()
    n_genes_keep = kept["associated_gene"].nunique()

    st.markdown("### Dataset summary")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Genes (total)", f"{n_genes:,}",
              help="Distinct genes with ≥1 detected isoform")
    c2.metric("Genes after filtering", f"{n_genes_keep:,}",
              delta=f"{n_genes_keep - n_genes:,}",
              help="Genes that retain ≥1 isoform after filtering")
    c3.metric("Isoforms (total)", f"{n_iso:,}")
    c4.metric("Isoforms after filtering", f"{n_keep:,}",
              delta=f"-{n_excl:,} removed",
              help="Isoforms with keep_decision == keep")

    c5, c6, c7, c8 = st.columns(4)
    c5.metric("Excluded isoforms", f"{n_excl:,}",
              help="Artifacts removed by ML filter / rescue rules")
    c6.metric("Novel isoforms (total)", f"{n_novel_all:,}",
              help="Novel transcripts: incomplete-splice_match + novel_in_catalog + "
                   "novel_not_in_catalog (i.e. all isoforms except full-splice_match "
                   "and the other reference categories).")
    c7.metric("Novel isoforms (kept)", f"{n_novel_keep:,}",
              delta=f"-{n_novel_all - n_novel_keep:,} removed",
              help=f"Of which {n_novel_verified:,} are in the curated/verified "
                   "novel set (novel_lr) used for the MAGNET & Novel-evidence views.")
    pct_keep = 100 * n_keep / n_iso if n_iso else 0
    c8.metric("% isoforms kept", f"{pct_keep:.1f}%")

    # Before / after-filtering breakdown by structural category
    st.markdown("### Isoforms by structural category")
    cat_all = iso_all["structural_category"].value_counts()
    cat_keep = kept["structural_category"].value_counts()
    cat_tbl = (pd.DataFrame({"All isoforms": cat_all, "After filtering": cat_keep})
               .fillna(0).astype(int))
    cat_tbl["% kept"] = (100 * cat_tbl["After filtering"] /
                         cat_tbl["All isoforms"].replace(0, np.nan)).round(1)
    cat_tbl = cat_tbl.sort_values("All isoforms", ascending=False)
    cc1, cc2 = st.columns([1.1, 1])
    with cc1:
        styled = (cat_tbl.style
                  .format({"All isoforms": "{:,}", "After filtering": "{:,}",
                           "% kept": "{:.1f}"})
                  .set_properties(**{"font-size": "16px"})
                  .set_table_styles([{"selector": "th",
                                      "props": [("font-size", "16px")]}]))
        st.table(styled)
    with cc2:
        melted = cat_tbl.reset_index()
        melted.columns = ["Category"] + list(cat_tbl.columns)
        order = cat_tbl.index[::-1].tolist()  # smallest at top for h-bars
        solid = [PIE_CAT_COLORS.get(c, PIE_CAT_OTHER) for c in order]

        def _rgba(hexc, a):
            h = hexc.lstrip("#")
            return f"rgba({int(h[0:2],16)},{int(h[2:4],16)},{int(h[4:6],16)},{a})"
        faded = [_rgba(c, 0.35) for c in solid]

        figb = go.Figure()
        figb.add_bar(y=order, x=cat_tbl.loc[order, "All isoforms"],
                     orientation="h", marker_color=faded,
                     marker_line_color=solid, marker_line_width=1,
                     name="All isoforms")
        figb.add_bar(y=order, x=cat_tbl.loc[order, "After filtering"],
                     orientation="h", marker_color=solid,
                     name="After filtering")
        figb.update_layout(
            barmode="overlay", height=340, bargap=0.4,
            margin=dict(l=0, r=10, t=10, b=0),
            font=dict(size=16),
            yaxis_title=None, xaxis_title=None,
            legend=dict(orientation="h", yanchor="bottom", y=1.02,
                        xanchor="left", x=0, title_text="",
                        font=dict(size=15)))
        figb.update_xaxes(tickfont=dict(size=15))
        figb.update_yaxes(tickfont=dict(size=15))
        st.plotly_chart(figb, use_container_width=True)

    st.markdown("### Samples")
    st.markdown(
        "<div style='font-size:16px'>"
        "6 human-heart samples across 3 phenotypes (2 each): "
        "<b>NFH</b> (non-failing control: sample1561, sample1662) · "
        "<b>DCM</b> (sample1532, sample1535) · "
        "<b>IHD</b> (sample1311, sample1518). "
        "Abundance is derived from per-sample read counts."
        "</div>",
        unsafe_allow_html=True,
    )

    st.markdown("---")
    st.markdown("**Try an example gene:** `MYH7` · `TTN` · `NPPA` · `SCN5A` · `MYBPC3`")
    st.info("Use the search box in the sidebar to look up a gene.")
    st.stop()


# ---- Gene selected ----
sub = filter_gene(iso_all[iso_all["associated_gene"] == gene])
grow = genes_tbl[genes_tbl["gene"] == gene]

st.title(f"{gene}")
st.markdown(external_links(gene))

if sub.empty:
    st.warning("No isoforms match the current filters for this gene. Loosen the filters in the sidebar.")
    st.stop()

c1, c2, c3, c4 = st.columns(4)
c1.metric("Isoforms shown", f"{len(sub):,}")
c2.metric("Novel", f"{int(sub['is_novel'].sum()):,}")
dom = sub.sort_values("total_count", ascending=False)["isoform"].iloc[0]
c3.metric("Top isoform", dom)
if "n_DET_sig" in grow.columns and len(grow):
    c4.metric("DET-significant", f"{int(grow['n_DET_sig'].iloc[0]):,}")
if len(grow):
    st.caption(f"{grow['chrom'].iloc[0]} · strand {grow['strand'].iloc[0]} · "
               f"{int(grow['n_isoforms'].iloc[0])} total isoforms in dataset")

tab_over, tab_struct, tab_tbl, tab_expr, tab_novel = st.tabs(
    ["Overview", "Structure", "Isoforms table", "Expression", "Novel evidence"]
)

# ---- Overview ----
with tab_over:
    def _hbar(df, x, y, height, color=None, cmap=None, fmt=None):
        fig = px.bar(df, x=x, y=y, orientation="h",
                     color=color, color_discrete_map=cmap)
        fig.update_layout(
            height=height, margin=dict(l=0, r=10, t=4, b=0),
            xaxis_title=None, yaxis_title=None, showlegend=False,
            bargap=0.4, font=dict(size=16),
            yaxis=dict(categoryorder="total ascending"))
        fig.update_xaxes(tickfont=dict(size=15))
        fig.update_yaxes(tickfont=dict(size=15))
        if not color:
            fig.update_traces(marker_color="#68a691")
        if fmt:
            fig.update_traces(hovertemplate="%{y}: %{x" + fmt + "}<extra></extra>")
        return fig

    cc1, cc2 = st.columns([1, 1])
    with cc1:
        st.markdown("**Isoforms by structural category**")
        ch = (sub["structural_category"].value_counts()
              .rename_axis("category").reset_index(name="isoforms"))
        cat_cmap = {c: PIE_CAT_COLORS.get(c, PIE_CAT_OTHER) for c in ch["category"]}
        st.plotly_chart(_hbar(ch, "isoforms", "category",
                              max(140, 34 * len(ch)),
                              color="category", cmap=cat_cmap),
                        use_container_width=True)
    with cc2:
        st.markdown("**Isoform usage (overall, top 10)**")
        top = (sub.sort_values("usage_overall", ascending=False).head(10)
               [["isoform", "usage_overall", "structural_category"]])
        iso_cmap = {c: PIE_CAT_COLORS.get(c, PIE_CAT_OTHER)
                    for c in top["structural_category"]}
        st.plotly_chart(_hbar(top, "usage_overall", "isoform",
                              max(140, 34 * len(top)), fmt=":.3f",
                              color="structural_category", cmap=iso_cmap),
                        use_container_width=True)

    st.markdown("**Gene-total CPM by phenotype** — long-read cohort (NFH · DCM · IHD)")
    pheno_tot = pd.DataFrame({
        "phenotype": [PHENO_LABEL[p] for p in PHENOTYPES],
        "mean_CPM": [sub[f"cpm_{p}"].sum() for p in PHENOTYPES],
    })
    pc1, _ = st.columns([1, 1])  # keep it narrow, not full-width
    with pc1:
        st.plotly_chart(
            _hbar(pheno_tot, "mean_CPM", "phenotype", 150,
                  color="phenotype", cmap=PHENO_COLOR, fmt=":.0f"),
            use_container_width=True)

# ---- Expression ----
with tab_expr:
    st.markdown("#### Isoform expression — long-read (left) vs MAGNET short-read (right)")
    col_lr, col_sr = st.columns(2)

    # -- LEFT: long-read isoform usage (this study) --
    with col_lr:
        st.markdown("**Long-read isoform usage** — this study (6 hearts)")
        st.caption("Each bar = 100% of the gene's reads in that group; segments are isoforms. "
                   "An isoform that grows from NFH → DCM/IHD is an isoform switch. "
                   "From raw long-read counts within the displayed isoforms.")
        u_group = st.radio("Groups", ["Per phenotype", "Per sample"],
                           horizontal=True, key="lr_group")
        u_top = st.slider("Isoforms shown distinctly (rest → 'Other')", 2, 15,
                          min(8, len(sub)), key="lr_top")
        figu = usage_stack_figure(sub, u_top, u_group)
        if figu is not None:
            st.plotly_chart(figu, use_container_width=True)
        else:
            st.info("No long-read counts to plot for the current selection.")

    # -- RIGHT: MAGNET short-read quantification --
    with col_sr:
        st.markdown("**MAGNET short-read isoform expression (TPM)**")
        st.caption("Independent MAGnet RNA-seq cohort (366 hearts), re-quantified with a "
                   "**Salmon decoy-aware** pipeline against a **merged reference** "
                   "(GENCODE v43 + 18,637 novel long-read isoforms). Box = isoform TPM "
                   "distribution per disease group. Novel isoforms use their own salmon row; "
                   "FSM/ISM use their matched GENCODE transcript's TPM (exact for FSM, "
                   "approximate for ISM). Hover a box for the source.")
        groups = st.multiselect("Disease groups", MAG_ORDER,
                                default=["Non-Failing", "DCM"], key="sr_groups")
        sc1, sc2 = st.columns([1, 1])
        box_n = sc1.slider("Max isoforms (box)", 1, 20, min(6, len(sub)), key="sr_top")
        log_y = sc2.toggle("Log y-axis", value=True, key="sr_log")
        if not groups:
            st.info("Select at least one disease group.")
        else:
            figm, n_box = magnet_box_figure(sub, box_n, groups, log_y)
            if figm is not None:
                st.plotly_chart(figm, use_container_width=True)
            else:
                st.info("No MAGNET-quantified isoforms for this gene.")

# ---- Isoforms table ----
with tab_tbl:
    # Informative SQANTI3 columns (identity / structure / coding / support / filter)
    show_cols = [
        "isoform", "structural_category", "subcategory", "associated_transcript",
        "length", "exons", "ref_exons",
        "coding", "ORF_length", "predicted_NMD",
        "all_canonical", "RTS_stage", "within_CAGE_peak", "polyA_motif_found",
        "pct_junct_supported", "ML_filter", "keep_decision",
        "total_count", "usage_overall",
    ]
    show_cols = [c for c in show_cols if c in sub.columns]
    table = sub.sort_values("total_count", ascending=False)[show_cols].reset_index(drop=True)

    col_cfg = {
        "isoform": st.column_config.TextColumn("isoform", help="PacBio isoform ID (PB.x.y)"),
        "structural_category": st.column_config.TextColumn(
            "structural_category", help="SQANTI3 category vs reference (FSM/ISM/NIC/NNC/…)"),
        "subcategory": st.column_config.TextColumn(
            "subcategory", help="SQANTI3 sub-classification (e.g. reference match, mono-exon)"),
        "associated_transcript": st.column_config.TextColumn(
            "assoc_transcript", help="Matched GENCODE reference transcript (ENST), if any"),
        "length": st.column_config.NumberColumn("length", help="Transcript length (bp)"),
        "exons": st.column_config.NumberColumn("exons", help="Number of exons in this isoform"),
        "ref_exons": st.column_config.NumberColumn(
            "ref_exons", help="Exon count of the matched reference transcript"),
        "coding": st.column_config.TextColumn("coding", help="Coding vs non-coding (SQANTI)"),
        "ORF_length": st.column_config.NumberColumn(
            "ORF_len", help="Predicted ORF length (amino acids)"),
        "predicted_NMD": st.column_config.TextColumn(
            "NMD", help="Predicted nonsense-mediated-decay target"),
        "all_canonical": st.column_config.TextColumn(
            "all_canonical", help="All splice junctions use canonical motifs"),
        "RTS_stage": st.column_config.TextColumn(
            "RTS", help="Flagged as a reverse-transcriptase template-switching artifact"),
        "within_CAGE_peak": st.column_config.TextColumn(
            "CAGE_5p", help="5' end falls within a CAGE peak (independent TSS support)"),
        "polyA_motif_found": st.column_config.TextColumn(
            "polyA_motif", help="PolyA motif found near the 3' end"),
        "pct_junct_supported": st.column_config.NumberColumn(
            "junc_sr_%", format="%.0f%%",
            help="% of splice junctions supported by MAGNET short reads"),
        "ML_filter": st.column_config.TextColumn(
            "ML_filter", help="SQANTI machine-learning filter call (Isoform/Artifact)"),
        "keep_decision": st.column_config.TextColumn(
            "keep", help="Final keep/exclude decision"),
        "total_count": st.column_config.NumberColumn(
            "total_count", help="Total long-read count across the 6 samples"),
        "usage_overall": st.column_config.NumberColumn(
            "usage", format="%.3f", help="Isoform's share of the gene's long-read reads"),
    }
    st.dataframe(table, use_container_width=True, height=440,
                 column_config={k: v for k, v in col_cfg.items() if k in show_cols},
                 hide_index=True)
    st.caption("Columns are SQANTI3 classification fields. Hover a header for its meaning. "
               "Per-sample counts, CPM and short-read/DET evidence are in the full CSV below "
               "and the Novel-evidence tab.")
    st.download_button(
        "⬇ Download this gene's isoforms (CSV)",
        data=sub.sort_values("total_count", ascending=False).to_csv(index=False).encode(),
        file_name=f"CardioIsoMap_{gene}_isoforms.csv",
        mime="text/csv",
    )

# ---- Structure ----
with tab_struct:
    exons = load_exons()
    c_s1, c_s2 = st.columns([1, 1])
    with c_s1:
        max_tx = st.slider("Max isoforms to draw", 2, 40, min(15, len(sub)))
    with c_s2:
        include_ref = st.toggle("Include GENCODE v43 reference transcripts",
                                value=False)

    top = sub.sort_values("total_count", ascending=False).head(max_tx)
    tids = top["isoform"].tolist()
    meta = {r["isoform"]: (r["structural_category"], bool(r["keep"]))
            for _, r in top.iterrows()}
    ex_use = exons

    n_ref = 0
    if include_ref:
        gx = load_gencode_exons()
        gx = gx[gx["gene"] == gene].copy()
        if not gx.empty:
            ref_order = (gx.groupby("transcript_id")["start"].min()
                         .sort_values().index.tolist())
            n_ref = len(ref_order)
            ex_use = pd.concat([exons,
                                gx[["gene", "transcript_id", "chrom",
                                    "strand", "start", "end"]]],
                               ignore_index=True)
            tids = tids + ref_order            # reference models at the bottom
            for t in ref_order:
                meta[t] = ("GENCODE", True)

    fig = structure_figure(ex_use, tids, gene, meta=meta)
    if fig is not None:
        if include_ref and n_ref == 0:
            st.caption("No GENCODE v43 transcript found for this gene "
                       "(novel/PacBio-only locus).")
        st.pyplot(fig)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
        st.download_button("⬇ Download figure (PNG)", data=buf.getvalue(),
                           file_name=f"CardioIsoMap_{gene}_structure.png", mime="image/png")
        plt.close(fig)
    else:
        st.info("No exon structures found for these isoforms.")

# ---- Novel evidence ----
with tab_novel:
    nv = sub[sub["is_novel"]].copy()
    if nv.empty:
        st.info("No novel isoforms for this gene under the current filters.")
    else:
        ncols = [
            "isoform", "novel_feature", "n_junctions", "n_junct_supported",
            "pct_junct_supported", "sr_mean_count", "sr_pct_detected",
            "DET_log2FC", "DET_padj", "DET_significant", "DET_direction",
        ]
        ncols = [c for c in ncols if c in nv.columns]
        st.markdown(f"**{len(nv)} novel isoform(s)** with verification evidence")
        st.dataframe(nv.sort_values("total_count", ascending=False)[ncols].reset_index(drop=True),
                     use_container_width=True, height=380)
        st.caption("Junction support and short-read (sr_*) columns quantify evidence; "
                   "DET_* = differential transcript usage (DCM vs NFH).")

st.markdown("---")
st.caption("CardioIsoMap")
