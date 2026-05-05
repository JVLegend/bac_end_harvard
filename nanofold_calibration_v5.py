#!/usr/bin/env python3
"""
Item #5 — PDB calibration via NanoFold structural distribution.

Em vez de fittar pLDDT->RMSD (que requer ground truth experimental, que nao temos
para nossos PDBs), fazemos OUTLIER DETECTION: comparamos os 7 descritores
estruturais dos nossos 43 PDBs vs distribuicao do NanoFold val (1000 cadeias)
matched por length bins. Variantes com z-score > 2 em ≥3 descritores sao flagged
como "structurally atypical" — possivel artefato de predicao ou divergencia real.

Output:
  reports/nanofold_calibration/per_variant_zscores.csv
  reports/nanofold_calibration/outliers.json
  reports/nanofold_calibration/metrics.json
"""

from __future__ import annotations

import csv
import json
import os

import numpy as np
import pandas as pd
from datasets import load_dataset

from config import REPORTS_DIR

OUT_DIR = os.path.join(REPORTS_DIR, "nanofold_calibration")
os.makedirs(OUT_DIR, exist_ok=True)
PANEL_STRUCT_CSV = os.path.join(REPORTS_DIR, "phenotype_probe_v3", "structure_features.csv")


def descriptors_from_ca(coords: np.ndarray) -> dict:
    L = len(coords)
    if L < 10: return None
    centroid = coords.mean(axis=0)
    rg = float(np.sqrt(((coords - centroid) ** 2).sum(axis=1).mean()))
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=2))
    mask = np.abs(np.arange(L)[:, None] - np.arange(L)[None, :]) >= 4
    contacts = (dist < 8.0) & mask
    n_pairs = mask.sum() // 2
    contact_density = float(contacts.sum() / 2 / max(n_pairs, 1))
    centered = coords - centroid
    cov = np.cov(centered.T)
    eigvals = np.sort(np.linalg.eigvalsh(cov))[::-1]
    return {
        "L": L,
        "rg": rg,
        "compactness_ratio": rg / (L ** 0.6),
        "contact_density": contact_density,
        "mean_ca_dist_norm": float(dist[mask].mean() / max(np.sqrt(L), 1)),
        "aspect_ratio": float(eigvals[0] / max(eigvals[2], 1e-6)),
    }


def main():
    print("=" * 70)
    print("ITEM #5 — Calibracao estrutural via NanoFold val")
    print("=" * 70)

    # 1. NanoFold val descriptors
    print("\nCarregando NanoFold val + computando descritores...")
    ds = load_dataset("ChrisHayduk/nanofold-public", split="validation")
    nf_rows = []
    for i, row in enumerate(ds):
        coords = np.array(row["ca_coords"], dtype=np.float32)
        mask = np.array(row["ca_mask"], dtype=bool)
        coords = coords[mask]
        d = descriptors_from_ca(coords)
        if d:
            nf_rows.append(d)
        if (i + 1) % 200 == 0:
            print(f"  [{i+1}/{len(ds)}]")
    nf = pd.DataFrame(nf_rows)
    print(f"NanoFold val: {len(nf)} cadeias com descritores")

    # 2. Z-score per length-bin (40-80, 80-120, 120-180, 180-256)
    bins = [(40, 80), (80, 120), (120, 180), (180, 256)]
    bin_stats = {}
    for lo, hi in bins:
        sub = nf[(nf["L"] >= lo) & (nf["L"] < hi)]
        if len(sub) < 10: continue
        bin_stats[f"{lo}-{hi}"] = {
            "n": int(len(sub)),
            "mean": {c: float(sub[c].mean()) for c in ["rg","compactness_ratio","contact_density","mean_ca_dist_norm","aspect_ratio"]},
            "std":  {c: float(sub[c].std())  for c in ["rg","compactness_ratio","contact_density","mean_ca_dist_norm","aspect_ratio"]},
        }
    print(f"Bins NanoFold construidos: {list(bin_stats.keys())}")

    # 3. Score painel
    if not os.path.exists(PANEL_STRUCT_CSV):
        print(f"ERRO: {PANEL_STRUCT_CSV} nao existe — rode structure_features_v3.py antes")
        return
    panel = pd.read_csv(PANEL_STRUCT_CSV)
    # Note: panel CSV has Z-NORMALIZED values (relativo ao painel). Vou re-extrair raw.
    # Para simplicidade aqui: comparar diretamente usando os valores no CSV (ja relativos
    # ao painel). Para comparacao raw, idealmente teriamos os PDBs aqui — pulo.
    print(f"Panel: {len(panel)} variantes")

    # Outlier detection (z-score robust): >2 sigma vs NanoFold pra mesmo bin
    # Para fazer corretamente, precisamos dos valores raw do panel — re-computamos
    PDB_DIR = "public/pdbs"
    panel_raw = []
    for fname in sorted(os.listdir(PDB_DIR)):
        if not fname.endswith(".pdb"): continue
        name = fname.replace(".pdb", "")
        coords = []
        for line in open(os.path.join(PDB_DIR, fname)):
            if line.startswith("ATOM") and line[13:15].strip() == "CA":
                try:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except ValueError: continue
        coords = np.array(coords)
        d = descriptors_from_ca(coords)
        if d:
            d["variant"] = name
            panel_raw.append(d)
    panel_raw = pd.DataFrame(panel_raw)
    print(f"Panel raw: {len(panel_raw)} variantes (lengths {panel_raw['L'].min()}-{panel_raw['L'].max()})")

    # Aviso: nossas variantes vao 218-668aa, mas NanoFold so vai ate 256aa.
    # Para variantes >256, comparamos com bin 180-256 (extrapolacao).
    rows = []
    for _, p in panel_raw.iterrows():
        L = p["L"]
        if L < 80:   bin_key = "40-80"
        elif L < 120: bin_key = "80-120"
        elif L < 180: bin_key = "120-180"
        else: bin_key = "180-256"  # inclui >256 (extrapolacao)
        if bin_key not in bin_stats: continue
        m = bin_stats[bin_key]["mean"]
        s = bin_stats[bin_key]["std"]
        zscores = {}
        for c in ["rg","compactness_ratio","contact_density","mean_ca_dist_norm","aspect_ratio"]:
            z = (p[c] - m[c]) / max(s[c], 1e-6)
            zscores[f"z_{c}"] = round(float(z), 3)
        n_atypical = sum(1 for z in zscores.values() if abs(z) > 2)
        rows.append({
            "variant": p["variant"],
            "L": int(L),
            "bin_compared": bin_key,
            "extrapolated": L > 256,
            **zscores,
            "n_atypical_descriptors": n_atypical,
            "outlier": n_atypical >= 3,
        })

    df_out = pd.DataFrame(rows).sort_values("n_atypical_descriptors", ascending=False)
    df_out.to_csv(os.path.join(OUT_DIR, "per_variant_zscores.csv"), index=False)

    outliers = df_out[df_out["outlier"]].to_dict(orient="records")
    n_extrap = int(df_out["extrapolated"].sum())
    print(f"\n=== Resultados v5 ===")
    print(f"Painel: {len(df_out)} variantes")
    print(f"Extrapolated (>256aa): {n_extrap} (mecA1/2, mcr-1/5)")
    print(f"Outliers (>=3 atypical descritores): {len(outliers)}")
    if outliers:
        for o in outliers[:5]:
            zs = {k: v for k, v in o.items() if k.startswith("z_")}
            print(f"  • {o['variant']:<15} L={o['L']:>4}  bin={o['bin_compared']:>8}  "
                  f"atypical={o['n_atypical_descriptors']}/5  extrap={o['extrapolated']}")

    metrics = {
        "n_panel_evaluated": int(len(df_out)),
        "n_extrapolated_beyond_256aa": n_extrap,
        "n_structural_outliers": int(len(outliers)),
        "outlier_threshold": "n_atypical_descriptors>=3 (z>2 em >=3 dos 5 descritores)",
        "nanofold_bins": bin_stats,
        "outlier_variants": [o["variant"] for o in outliers],
    }
    with open(os.path.join(OUT_DIR, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2)
    with open(os.path.join(OUT_DIR, "outliers.json"), "w") as f:
        json.dump(outliers, f, indent=2)
    print(f"\n[Saved] {OUT_DIR}/per_variant_zscores.csv")
    print(f"[Saved] {OUT_DIR}/metrics.json")


if __name__ == "__main__":
    main()
