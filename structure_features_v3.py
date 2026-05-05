#!/usr/bin/env python3
"""
Item #3 — Structure-aware ensemble (lite).

Em vez de pretreinar encoder em NanoFold (heavy), usamos diretamente os
descritores estruturais das 41 PDBs do nosso painel (geradas via
ColabFold/ESMFold) + 5 PDBs longos (mecA1/2, mcr-1/5) quando AF Server
voltar.

Para cada PDB do painel, extrai 7 descritores:
  - L                : numero de residuos
  - Rg               : radius of gyration
  - compactness_ratio: Rg / L^0.6
  - contact_density  : fracao Ca-Ca pairs <8 A (i,j com |i-j|>=4)
  - mean_ca_dist_norm: media Ca-Ca distance / sqrt(L)
  - aspect_ratio     : eigvalue1/eigvalue3 da inertia tensor
  - mean_plddt       : pLDDT medio (do B-factor da PDB)

Concatena com ESM-2 (640) + ProtT5 (1024) = 1671d.
Re-treina classifier multi-label e compara com v2 (concat sem struct features).

Output:
  reports/phenotype_probe_v3/metrics.json
  reports/phenotype_probe_v3/classifier.joblib
  reports/phenotype_probe_v3/structure_features.csv  (insumo do site)
"""

from __future__ import annotations

import csv
import json
import os
import warnings

import numpy as np
import pandas as pd
import joblib
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import f1_score
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import MultiLabelBinarizer

from config import REPORTS_DIR

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ESM_MODEL_NAME = os.environ.get("ESM_MODEL", "esm2_t30_150M_UR50D")
OUT_DIR = os.path.join(REPORTS_DIR, "phenotype_probe_v3")
PDB_DIR = "public/pdbs"
ESM_EMB_DIR = os.path.join(REPORTS_DIR, "protein_scoring", "embeddings")
PT5_EMB_DIR = os.path.join(REPORTS_DIR, "prott5_ensemble", "embeddings")

FAMILY_LABELS = {
    "mecA":        (["penam","cephalosporin","methicillin"], "antibiotic target replacement"),
    "mecA1":       (["penam","cephalosporin","methicillin"], "antibiotic target replacement"),
    "mecA2":       (["penam","cephalosporin","methicillin"], "antibiotic target replacement"),
    "blaKPC":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
    "blaNDM":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
    "blaOXA":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
    "blaVIM":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
    "blaIMP":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
    "blaGES":      (["carbapenem","cephalosporin"], "antibiotic inactivation"),
    "blaCTX-M":    (["cephalosporin","penam"], "antibiotic inactivation"),
    "blaCTX-M-15": (["cephalosporin","penam"], "antibiotic inactivation"),
    "vanA":        (["glycopeptide"], "antibiotic target alteration"),
    "mcr-1":       (["polymyxin","peptide"], "antibiotic target alteration"),
    "mcr-5":       (["polymyxin","peptide"], "antibiotic target alteration"),
    "qnrS":        (["fluoroquinolone"], "antibiotic target protection"),
    "armA":        (["aminoglycoside"], "antibiotic target alteration"),
}


def labels_for(name):
    if name in FAMILY_LABELS:
        return FAMILY_LABELS[name]
    for k in FAMILY_LABELS:
        if name.startswith(k):
            return FAMILY_LABELS[k]
    return None


def parse_ca_coords(pdb_path: str):
    """Le coordenadas Ca + B-factors (pLDDT) de um PDB."""
    coords, plddts = [], []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[13:15].strip() != "CA":
                continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                b = float(line[60:66])
                coords.append([x, y, z])
                plddts.append(b)
            except ValueError:
                continue
    return np.array(coords), np.array(plddts)


def structure_descriptors(coords: np.ndarray, plddts: np.ndarray) -> dict:
    """Retorna 7 descritores estruturais compactos."""
    L = len(coords)
    if L < 10:
        return None

    # 1. Radius of gyration
    centroid = coords.mean(axis=0)
    rg = float(np.sqrt(((coords - centroid) ** 2).sum(axis=1).mean()))

    # 2. Compactness (Rg vs scaling N^0.6)
    compactness = rg / (L ** 0.6)

    # 3. Contact density: pares Ca-Ca <8A excluindo |i-j|<4
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=2))
    mask = np.abs(np.arange(L)[:, None] - np.arange(L)[None, :]) >= 4
    contacts = (dist < 8.0) & mask
    n_pairs = mask.sum() // 2
    contact_density = float(contacts.sum() / 2 / max(n_pairs, 1))

    # 4. Mean Ca-Ca distance normalized
    mean_ca_dist_norm = float(dist[mask].mean() / max(np.sqrt(L), 1))

    # 5. Aspect ratio via inertia tensor PCA
    centered = coords - centroid
    cov = np.cov(centered.T)
    eigvals = np.sort(np.linalg.eigvalsh(cov))[::-1]
    aspect_ratio = float(eigvals[0] / max(eigvals[2], 1e-6))

    # 6. Mean pLDDT
    mean_plddt = float(plddts.mean())

    return {
        "L": L,
        "rg": round(rg, 3),
        "compactness_ratio": round(compactness, 3),
        "contact_density": round(contact_density, 4),
        "mean_ca_dist_norm": round(mean_ca_dist_norm, 3),
        "aspect_ratio": round(aspect_ratio, 3),
        "mean_plddt": round(mean_plddt, 2),
    }


def main():
    print("=" * 70)
    print("ITEM #3 — Structure-aware ensemble (PDB descriptors)")
    print("=" * 70)
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.isdir(PDB_DIR):
        print(f"ERRO: {PDB_DIR} nao existe.")
        return

    # 1. Extract structure features dos PDBs
    print(f"\nExtraindo descritores estruturais de {PDB_DIR}/...")
    rows = []
    for fname in sorted(os.listdir(PDB_DIR)):
        if not fname.endswith(".pdb"):
            continue
        name = fname.replace(".pdb", "")
        coords, plddts = parse_ca_coords(os.path.join(PDB_DIR, fname))
        desc = structure_descriptors(coords, plddts)
        if desc is None:
            print(f"  [SKIP] {name}")
            continue
        desc["variant"] = name
        rows.append(desc)
        print(f"  {name:<14} L={desc['L']:>4} Rg={desc['rg']:>5.1f}A "
              f"compact={desc['compactness_ratio']:.2f} "
              f"contacts={desc['contact_density']:.3f} "
              f"plddt={desc['mean_plddt']:.1f}")
    print(f"[Saved] {len(rows)} structure feature vectors")

    df_struct = pd.DataFrame(rows)
    cols = ["L", "rg", "compactness_ratio", "contact_density",
            "mean_ca_dist_norm", "aspect_ratio", "mean_plddt"]

    # Z-score normalize (relativo ao painel)
    df_struct[cols] = (df_struct[cols] - df_struct[cols].mean()) / df_struct[cols].std().replace(0, 1)
    struct_dict = {row["variant"]: np.array([row[c] for c in cols], dtype=np.float32)
                   for _, row in df_struct.iterrows()}

    df_struct.to_csv(os.path.join(OUT_DIR, "structure_features.csv"), index=False)

    # 2. Build training set: ESM-2 + ProtT5 + struct
    print("\nMontando ensemble (ESM-2 + ProtT5 + 7 struct features)...")
    X_list, names_list, drugs_list = [], [], []
    for variant in struct_dict:
        labs = labels_for(variant)
        if labs is None:
            continue
        esm_p = os.path.join(ESM_EMB_DIR, f"{variant}__{ESM_MODEL_NAME}.npy")
        pt5_p = os.path.join(PT5_EMB_DIR, f"{variant}.npy")
        if not (os.path.exists(esm_p) and os.path.exists(pt5_p)):
            continue
        emb = np.concatenate([np.load(esm_p), np.load(pt5_p), struct_dict[variant]])
        X_list.append(emb)
        names_list.append(variant)
        drugs_list.append(labs[0])

    if len(X_list) < 5:
        print(f"ERRO: poucos exemplos com 3 modalidades ({len(X_list)})")
        return

    X = np.vstack(X_list)
    print(f"  N={len(X)}, dim={X.shape[1]} (640 ESM-2 + 1024 ProtT5 + 7 struct)")

    mlb = MultiLabelBinarizer()
    Y = mlb.fit_transform(drugs_list)

    # 3. LOO CV
    loo = LeaveOneOut()
    Y_loo = np.zeros_like(Y)
    for tr, te in loo.split(X):
        c = OneVsRestClassifier(LogisticRegression(max_iter=2000, C=1.0))
        c.fit(X[tr], Y[tr])
        Y_loo[te] = c.predict(X[te])

    hamming = float(np.mean(Y_loo == Y))
    exact = float((Y_loo == Y).all(axis=1).mean())
    f1_per = {}
    for j, cls in enumerate(mlb.classes_):
        tp = int(((Y_loo[:, j] == 1) & (Y[:, j] == 1)).sum())
        fp = int(((Y_loo[:, j] == 1) & (Y[:, j] == 0)).sum())
        fn = int(((Y_loo[:, j] == 0) & (Y[:, j] == 1)).sum())
        p = tp / (tp + fp) if (tp + fp) else 0
        r = tp / (tp + fn) if (tp + fn) else 0
        f1 = 2 * p * r / (p + r) if (p + r) else 0
        f1_per[cls] = round(f1, 3)

    print(f"\n=== Resultados v3 (structure-aware) ===")
    print(f"  N samples:    {len(X)}")
    print(f"  Dim ensemble: {X.shape[1]}")
    print(f"  Hamming LOO:  {hamming:.3f}")
    print(f"  Exact match:  {exact:.3f}")
    print(f"  F1 per class: {f1_per}")

    metrics = {
        "version": "v3_structure_aware",
        "ensemble": "ESM-2 (640) + ProtT5 (1024) + struct (7) = 1671d",
        "n_samples": int(len(X)),
        "n_classes": int(len(mlb.classes_)),
        "hamming_acc_loo": round(hamming, 4),
        "exact_match_loo": round(exact, 4),
        "f1_per_class_loo": f1_per,
        "compared_to_v2": {
            "v2_hamming": 0.9406,
            "v2_exact": 0.6279,
            "delta_hamming_pp": round(hamming - 0.9406, 4),
            "delta_exact_pp": round(exact - 0.6279, 4),
        },
    }
    with open(os.path.join(OUT_DIR, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2)

    clf_full = OneVsRestClassifier(LogisticRegression(max_iter=2000, C=1.0))
    clf_full.fit(X, Y)
    joblib.dump({"clf": clf_full, "mlb": mlb, "model": ESM_MODEL_NAME,
                 "ensemble": "ESM2+ProtT5+struct", "struct_cols": cols},
                os.path.join(OUT_DIR, "classifier.joblib"))

    # Para o site: salvar struct features p/ visualizacao
    site_data = {
        "variants": names_list,
        "features": cols,
        "matrix": [[float(v) for v in X[i, -7:]] for i in range(len(X))],
        "metrics": metrics,
    }
    with open("public/data/structure_features.json", "w") as f:
        json.dump(site_data, f, indent=2)

    print(f"\n[Saved] {OUT_DIR}/metrics.json, classifier.joblib")
    print(f"[Saved] public/data/structure_features.json")


if __name__ == "__main__":
    main()
