#!/usr/bin/env python3
"""
Tier 1 — Extended dataset com 8991 proteinas AMRFinderPlus + 43 do painel.

Output: data/hf_dataset/smartsepsis_oph_extended.parquet
  Colunas: variant_id, source, drug_classes, esm2_embedding, esm2_model
  (sem ProtT5/PDB — Tier 1 leve. Para multimodal completo ver split 'panel'.)
"""

import os
import json
import hashlib
import numpy as np
import pandas as pd

from config import REPORTS_DIR, BASE_DIR

ESM_MODEL_NAME = "esm2_t30_150M_UR50D"
OUT_DIR = os.path.join(BASE_DIR, "data", "hf_dataset")
AMR_NPZ = os.path.join(REPORTS_DIR, "amrfinderplus", f"embeddings_{ESM_MODEL_NAME}.npz")
ESM_PANEL_DIR = os.path.join(REPORTS_DIR, "protein_scoring", "embeddings")
PANEL_PARQUET = os.path.join(OUT_DIR, "smartsepsis_oph.parquet")


def main():
    print("=" * 70)
    print("TIER 1 — Build extended dataset")
    print("=" * 70)

    # 1. AMRFinderPlus (8991)
    print(f"\nCarregando {AMR_NPZ}...")
    z = np.load(AMR_NPZ, allow_pickle=True)
    amr_X = z["X"]
    amr_names = list(z["names"])
    amr_labels = [list(l) for l in z["labels"]]
    print(f"  AMRFinderPlus: {len(amr_X)} proteinas, dim={amr_X.shape[1]}")

    rows_amr = []
    for name, emb, drugs in zip(amr_names, amr_X, amr_labels):
        rows_amr.append({
            "variant_id": str(name),
            "source": "AMRFinderPlus",
            "drug_classes": list(drugs),
            "esm2_embedding": emb.tolist(),
            "esm2_model": ESM_MODEL_NAME,
        })

    # 2. Panel (43) — adiciona como source="panel" pra rastreabilidade
    print(f"\nCarregando panel (43 do nosso painel)...")
    panel_df = pd.read_parquet(PANEL_PARQUET)
    rows_panel = []
    for _, p in panel_df.iterrows():
        if p["esm2_embedding"] is None or (isinstance(p["esm2_embedding"], list) and not p["esm2_embedding"]):
            continue
        rows_panel.append({
            "variant_id": p["variant_id"],
            "source": "SmartSepsis-Oph panel",
            "drug_classes": list(p["drug_classes"]) if p["drug_classes"] is not None else [],
            "esm2_embedding": p["esm2_embedding"],
            "esm2_model": ESM_MODEL_NAME,
        })
    print(f"  Panel: {len(rows_panel)} entradas")

    # 3. Concat
    rows = rows_amr + rows_panel
    df = pd.DataFrame(rows)
    print(f"\nTotal extended: {len(df)} entradas")

    out_path = os.path.join(OUT_DIR, "smartsepsis_oph_extended.parquet")
    df.to_parquet(out_path, compression="snappy")
    sz = os.path.getsize(out_path) / 1024 / 1024
    print(f"[Saved] {out_path}  ({sz:.1f} MB)")

    # Update manifest
    mfp = os.path.join(OUT_DIR, "manifest.json")
    if os.path.exists(mfp):
        manifest = json.load(open(mfp))
    else:
        manifest = {"dataset": "smartsepsis-oph", "version": "1.1.0", "files": {}}
    with open(out_path, "rb") as f:
        manifest["files"]["smartsepsis_oph_extended.parquet"] = hashlib.sha256(f.read()).hexdigest()
    manifest["n_rows_panel"] = 45
    manifest["n_rows_extended"] = len(df)
    manifest["version"] = "1.1.0"
    with open(mfp, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[Saved] manifest atualizado")


if __name__ == "__main__":
    main()
