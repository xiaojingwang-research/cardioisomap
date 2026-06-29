# Deployment — CardioIsoMap

Two repos: a small **app repo** (~11 MB, deployed by Streamlit) and a large **figures repo**
(~1.3 GB of PNGs, served free via jsDelivr). `gh` is not installed here, so the commands below
use plain `git` + the GitHub website to create repos. Replace `<USER>` with your GitHub username.

## A. App repo (code + processed tables)
The `.gitignore` already excludes raw data and the 1.3 GB `data/processed/sashimi/` images.
```bash
cd "/Users/wangx11/Library/CloudStorage/Box-Box/000/LR6_may11/CardioIsoMap"
git init
git add app code docs data/processed/*.parquet \
        data/processed/gene_index.csv data/processed/gene_regions.csv \
        data/processed/sashimi_manifest.csv \
        README.md PROJECT_LOG.md .gitignore
git commit -m "CardioIsoMap app + processed tables"
# create an EMPTY repo named 'cardioisomap' on github.com, then:
git branch -M main
git remote add origin https://github.com/<USER>/cardioisomap.git
git push -u origin main
```

## B. Figures repo (the 1.3 GB of sashimi PNGs)
```bash
cd "/Users/wangx11/Library/CloudStorage/Box-Box/000/LR6_may11/CardioIsoMap/data/processed"
mkdir -p /tmp/cardioisomap-figures && cp -R sashimi /tmp/cardioisomap-figures/
cd /tmp/cardioisomap-figures
git init && git add sashimi && git commit -m "per-gene coverage/sashimi PNGs"
git branch -M main
git remote add origin https://github.com/<USER>/cardioisomap-figures.git
git push -u origin main      # one-time ~1.3 GB upload
```
jsDelivr then serves them at:
`https://cdn.jsdelivr.net/gh/<USER>/cardioisomap-figures@main/sashimi/<gene>.png`

Sanity check after push:
```bash
curl -I "https://cdn.jsdelivr.net/gh/<USER>/cardioisomap-figures@main/sashimi/MYOM2.png"   # expect HTTP/2 200
```

## C. Deploy on Streamlit Community Cloud
1. https://share.streamlit.io → sign in with GitHub → **New app**.
2. Repo `cardioisomap`, branch `main`, main file `app/app.py`.
3. **Advanced settings → Secrets**, paste:
   ```toml
   IMG_BASE = "https://cdn.jsdelivr.net/gh/<USER>/cardioisomap-figures@main/sashimi"
   ```
4. Deploy → public URL `https://<name>.streamlit.app`.

Notes:
- For immutable image URLs, pin a tag/commit instead of `@main` (e.g. `@v1`) in `IMG_BASE`.
- If `IMG_BASE` is unset, the app falls back to local `data/processed/sashimi/` (dev only).
- `gh` CLI alternative to creating repos in the browser: `brew install gh && gh auth login`,
  then `gh repo create <USER>/cardioisomap --public --source=. --push`.
