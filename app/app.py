"""
CardioIsoMap — gene-searchable cardiac long-read isoform atlas.

Run locally:  streamlit run app/app.py
Data:         data/processed/*.parquet  (built by code/01_build_tables.py and 02_build_structures.py)

NOTE: abundance is derived from per-sample counts / total (NOT the SQANTI `FL` column).
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


def expression_figure(sub: pd.DataFrame, value: str, grouping: str, top_n: int):
    """Heatmap of isoform x (sample|phenotype) for the chosen value type."""
    sub = sub.sort_values("total_count", ascending=False).head(top_n)
    if sub.empty:
        return None

    if grouping == "Per sample":
        cols = SAMPLES
        col_labels = [f"{s}\n({SAMPLE_PHENO[s]})" for s in cols]
        if value == "Raw count":
            mat = sub[cols]
        elif value == "CPM":
            mat = sub[[f"cpm_{s}" for s in cols]]
        else:  # Usage proportion
            mat = sub[[f"usage_{s}" for s in cols]]
    else:  # Per phenotype (mean)
        cols = PHENOTYPES
        col_labels = [PHENO_LABEL[p] for p in cols]
        if value == "Raw count":
            mat = sub[[f"count_{p}" for p in cols]]
        elif value == "CPM":
            mat = sub[[f"cpm_{p}" for p in cols]]
        else:
            mat = sub[[f"usage_{p}" for p in cols]]

    mat = mat.copy()
    mat.columns = col_labels
    mat.index = sub["isoform"].values

    fig = px.imshow(
        mat,
        labels=dict(x="", y="Isoform", color=value),
        aspect="auto",
        color_continuous_scale="Viridis",
    )
    fig.update_layout(height=max(300, 22 * len(sub) + 120), margin=dict(l=10, r=10, t=30, b=10))
    return fig


# Per-category exon fill colors (kept isoforms). Excluded isoforms are drawn lighter.
CAT_COLORS = {
    "full-splice_match":     "#efc7c2",  # FSM
    "incomplete-splice_match": "#ffe5d4",  # ISM
    "novel_in_catalog":      "#bfd3c1",  # NIC
    "novel_not_in_catalog":  "#68a691",  # NNC
}
CAT_OTHER = "#cccccc"  # fusion / antisense / intergenic / genic


def _lighten(hex_color: str, frac: float) -> str:
    """Blend a hex color toward white by `frac` (0=unchanged, 1=white)."""
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    r = int(r + (255 - r) * frac)
    g = int(g + (255 - g) * frac)
    b = int(b + (255 - b) * frac)
    return f"#{r:02x}{g:02x}{b:02x}"


def exon_color(category: str, keep: bool) -> str:
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

    fig, ax = plt.subplots(figsize=(10, max(1.5, 0.5 * len(order) + 1)))
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
    ax.set_yticklabels(order, fontsize=8)
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
    "(NFH, DCM, IHD). Abundance from per-sample counts (not SQANTI `FL`)."
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
    st.subheader("A gene-searchable atlas of cardiac transcript isoforms")

    kept = iso_all[iso_all["keep"]]
    n_iso = len(iso_all)
    n_keep = int(iso_all["keep"].sum())
    n_excl = n_iso - n_keep
    n_novel_all = int(iso_all["is_novel"].sum())
    n_novel_keep = int(kept["is_novel"].sum())
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
    c6.metric("Novel isoforms (total)", f"{n_novel_all:,}")
    c7.metric("Novel isoforms (kept)", f"{n_novel_keep:,}")
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
        st.dataframe(cat_tbl, use_container_width=True)
    with cc2:
        st.bar_chart(cat_tbl[["All isoforms", "After filtering"]])

    st.markdown("### Samples")
    st.caption(
        "6 human-heart samples across 3 phenotypes (2 each): "
        "**NFH** (non-failing control: sample1561, sample1662) · "
        "**DCM** (sample1532, sample1535) · "
        "**IHD** (sample1311, sample1518). "
        "Abundance is derived from per-sample read counts and `total` — the SQANTI `FL` column is not used."
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

tab_over, tab_expr, tab_cov, tab_tbl, tab_struct, tab_novel = st.tabs(
    ["Overview", "Expression", "Coverage / Sashimi", "Isoforms table", "Structure", "Novel evidence"]
)

# ---- Overview ----
with tab_over:
    cats_here = sub["structural_category"].value_counts()
    cc1, cc2 = st.columns([1, 1])
    with cc1:
        st.markdown("**Isoforms by structural category**")
        st.bar_chart(cats_here)
    with cc2:
        st.markdown("**Isoform usage (overall, top 10)**")
        top = sub.sort_values("usage_overall", ascending=False).head(10)
        st.bar_chart(top.set_index("isoform")["usage_overall"])
    st.markdown("**Mean CPM by phenotype (gene total)**")
    pheno_tot = pd.DataFrame({
        "phenotype": [PHENO_LABEL[p] for p in PHENOTYPES],
        "mean_CPM": [sub[f"cpm_{p}"].sum() for p in PHENOTYPES],
    })
    st.bar_chart(pheno_tot.set_index("phenotype"))

# ---- Expression ----
with tab_expr:
    cc1, cc2, cc3 = st.columns(3)
    value = cc1.radio("Value", ["CPM", "Raw count", "Usage proportion"], horizontal=False)
    grouping = cc2.radio("Columns", ["Per phenotype", "Per sample"], horizontal=False)
    top_n = cc3.slider("Max isoforms", 5, 60, min(25, len(sub)))
    fig = expression_figure(sub, value, grouping, top_n)
    if fig is not None:
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("Nothing to plot for the current selection.")

# ---- Coverage / Sashimi ----
with tab_cov:
    fname = f"{safe_name(gene)}.png"
    local_png = SASHIMI_DIR / fname
    has_image = (gene in available_sashimi()) or local_png.exists()
    if has_image:
        st.markdown(
            "Per-phenotype read **coverage** with splice-junction (**sashimi**) arcs — "
            "**NFH** (blue), **DCM** (red), **IHD** (green) — over transcript models "
            "(known = red, novel `PB.*` = green)."
        )
        if IMG_BASE:
            url = f"{IMG_BASE}/{fname}"
            st.image(url, use_container_width=True)
            st.link_button("⬇ Open / download full-res (PNG)", url)
        elif local_png.exists():
            st.image(str(local_png), use_container_width=True)
            st.download_button("⬇ Download figure (PNG)", data=local_png.read_bytes(),
                               file_name=f"CardioIsoMap_{gene}_coverage.png", mime="image/png")
        else:
            st.warning("Image is listed in the manifest but not reachable "
                       "(no IMG_BASE set and no local file).")
    else:
        st.info(
            f"No coverage plot for **{gene}** — likely a very large or fusion-span locus "
            "that was skipped during rendering. Coverage/sashimi plots are produced in batch "
            "by `code/03_render_sashimi.R` (Gviz)."
        )

# ---- Isoforms table ----
with tab_tbl:
    show_cols = [
        "isoform", "structural_category", "subcategory", "length", "exons",
        "coding", "is_novel", "keep_decision", "total_count",
        "cpm_NFH", "cpm_DCM", "cpm_IHD", "usage_overall",
    ]
    show_cols = [c for c in show_cols if c in sub.columns]
    table = sub.sort_values("total_count", ascending=False)[show_cols].reset_index(drop=True)
    st.dataframe(table, use_container_width=True, height=420)
    st.download_button(
        "⬇ Download this gene's isoforms (CSV)",
        data=sub.sort_values("total_count", ascending=False).to_csv(index=False).encode(),
        file_name=f"CardioIsoMap_{gene}_isoforms.csv",
        mime="text/csv",
    )

# ---- Structure ----
with tab_struct:
    exons = load_exons()
    max_tx = st.slider("Max isoforms to draw", 2, 40, min(15, len(sub)))
    top = sub.sort_values("total_count", ascending=False).head(max_tx)
    tids = top["isoform"].tolist()
    meta = {r["isoform"]: (r["structural_category"], bool(r["keep"]))
            for _, r in top.iterrows()}
    fig = structure_figure(exons, tids, gene, meta=meta)
    if fig is not None:
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
st.caption("CardioIsoMap · built from PacBio long-read SQANTI3 output · abundance excludes the `FL` column.")
