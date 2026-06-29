#!/usr/bin/env python3
"""
CardioIsoMap — Step 4: rebuild the sashimi manifest from the PNGs on disk.

The batch render (code/03_render_sashimi.R) is sharded; the per-shard manifest
merge did not always complete. The deployed Coverage tab needs a definitive list
of which genes have a rendered image (it can't stat a CDN), so we regenerate the
manifest by scanning data/processed/sashimi/*.png and mapping filenames back to
genes via gene_regions.csv (same safe_name sanitisation as the renderer & app).

Output: data/processed/sashimi_manifest.csv  (gene, png, n_tx)
"""
import re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
SASHIMI = PROC / "sashimi"


def safe_name(g: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(g))


def main() -> None:
    regions = pd.read_csv(PROC / "gene_regions.csv")
    regions["safe"] = regions["gene"].map(safe_name)
    # map safe filename stem -> (gene, n_tx); if two genes sanitise to the same
    # stem, keep the one with more transcripts (deterministic).
    regions = regions.sort_values("n_tx", ascending=False).drop_duplicates("safe")
    by_safe = regions.set_index("safe")[["gene", "n_tx"]].to_dict("index")

    pngs = sorted(SASHIMI.glob("*.png"))
    rows, unmatched = [], []
    for p in pngs:
        stem = p.stem
        info = by_safe.get(stem)
        if info is None:
            unmatched.append(stem)
            continue
        rows.append({"gene": info["gene"], "png": p.name, "n_tx": int(info["n_tx"])})

    manifest = pd.DataFrame(rows).sort_values("gene").reset_index(drop=True)
    out = PROC / "sashimi_manifest.csv"
    manifest.to_csv(out, index=False)

    print(f"PNGs on disk:        {len(pngs):,}")
    print(f"Manifest rows:       {len(manifest):,}")
    print(f"Unmatched filenames: {len(unmatched):,}")
    if unmatched[:5]:
        print("  e.g.", unmatched[:5])
    present = set(manifest["gene"])
    for g in ["MYOM2", "MYH7", "TTN", "NPPA"]:
        print(f"  {g}: {'present' if g in present else 'MISSING'}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
