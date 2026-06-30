# Deployment — CardioIsoMap

Single public repo (~73 MB: app code + processed parquet/CSV tables), deployed on
**Streamlit Community Cloud**. No CDN, no second repo, no secrets — the Coverage/Sashimi
tab (the only thing that needed hosted images) has been removed.

**Live app:** https://github.com/xiaojingwang-research/cardioisomap → deployed via share.streamlit.io

## A. Push the repo to GitHub
`.gitignore` excludes raw data, the (legacy) `data/processed/sashimi/` PNGs, and local tooling;
the `data/processed/*.parquet` + `*.csv` tables ARE committed (the app needs them at runtime).
```bash
cd "/Users/wangx11/Library/CloudStorage/Box-Box/000/LR6_may11/CardioIsoMap"
git add -A
git commit -m "your message"
git push                       # to https://github.com/xiaojingwang-research/cardioisomap
```
Auth over HTTPS uses a **Personal Access Token** (Settings → Developer settings → Tokens
(classic) → scope `repo`) as the password, not your account password. Run once to cache it:
`git config --global credential.helper osxkeychain`.

## B. Deploy on Streamlit Community Cloud
1. https://share.streamlit.io → sign in with GitHub → **Create app** → *Deploy a public app from GitHub*.
2. **Repository** `xiaojingwang-research/cardioisomap` · **Branch** `main` · **Main file** `app/app.py`.
3. (optional) set a custom subdomain, then **Deploy**.
4. Live at `https://<name>.streamlit.app`. No secrets needed.

## C. Updating the live app
Streamlit watches the repo — every push to `main` auto-redeploys within ~1 minute:
```bash
git add -A && git commit -m "..." && git push
```

## Notes
- Free tier ≈ 1 GB RAM. The ~73 MB of tables (mostly the two MAGNET TPM parquets, ~52 MB) load
  fine; if memory ever gets tight, downcast/aggregate the MAGNET tables.
- Rebuild processed tables only if raw inputs change — see the `code/01..07` scripts in `README.md`.
- The Sashimi render scripts (`03_render_sashimi.R`, `04_rebuild_manifest.py`) and their 1.3 GB of
  PNGs remain on disk but are no longer used by the app.
