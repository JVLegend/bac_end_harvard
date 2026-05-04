#!/usr/bin/env python3
"""
Embed o catalogo AMRFinderPlus inteiro com ESM-2 e treina classificador
multi-label de drug class via LogisticRegression OvR.

- Input: data/amrfinderplus/AMRProt.fa (~10k proteinas) + ReferenceGeneCatalog.txt
- Output:
    reports/amrfinderplus/embeddings.npz   (matriz N x D + labels)
    reports/amrfinderplus/classifier.joblib
    reports/amrfinderplus/metrics.json
    reports/amrfinderplus/predictions_painel.csv (predicoes nas 42 variantes nossas)

Uso (M2 16GB):
    /Users/iaparamedicos/envs/dev/bin/python amrfinderplus_embed.py
    # opcional limitar amostragem:
    AMR_LIMIT=2000 /Users/iaparamedicos/envs/dev/bin/python amrfinderplus_embed.py
    # modelo maior:
    ESM_MODEL=esm2_t33_650M_UR50D /Users/iaparamedicos/envs/dev/bin/python amrfinderplus_embed.py
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch
import esm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score, hamming_loss
from sklearn.model_selection import train_test_split
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import MultiLabelBinarizer

import joblib

from config import BASE_DIR, REPORTS_DIR

ESM_MODEL_NAME = os.environ.get("ESM_MODEL", "esm2_t30_150M_UR50D")
LIMIT = int(os.environ.get("AMR_LIMIT", "0"))  # 0 = sem limite
BATCH = int(os.environ.get("AMR_BATCH", "8"))
ESM_MAX_LEN = 1022

OUT_DIR = os.path.join(REPORTS_DIR, "amrfinderplus")
EMB_PATH = os.path.join(OUT_DIR, f"embeddings_{ESM_MODEL_NAME}.npz")
CLF_PATH = os.path.join(OUT_DIR, f"classifier_{ESM_MODEL_NAME}.joblib")
METRICS_PATH = os.path.join(OUT_DIR, f"metrics_{ESM_MODEL_NAME}.json")
PRED_PATH = os.path.join(OUT_DIR, f"predictions_painel_{ESM_MODEL_NAME}.csv")

PROT_FASTA = os.path.join(BASE_DIR, "data", "amrfinderplus", "AMRProt.fa")
CATALOG_TSV = os.path.join(BASE_DIR, "data", "amrfinderplus", "ReferenceGeneCatalog.txt")

PAINEL_PROTS_DIR = os.path.join(REPORTS_DIR, "protein_scoring", "proteins")

MODEL_LAYERS = {
    "esm2_t6_8M_UR50D": 6,
    "esm2_t12_35M_UR50D": 12,
    "esm2_t30_150M_UR50D": 30,
    "esm2_t33_650M_UR50D": 33,
    "esm2_t36_3B_UR50D": 36,
}


def parse_amrprot(path: str) -> list[tuple[str, str]]:
    """Retorna [(header_full, sequence), ...]."""
    out, hdr, seq = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    out.append((hdr, "".join(seq)))
                hdr = line[1:]
                seq = []
            else:
                seq.append(line)
    if hdr is not None:
        out.append((hdr, "".join(seq)))
    return out


def load_drug_class_map(catalog_path: str) -> dict[str, list[str]]:
    """Mapeia MULTIPLAS chaves -> drug_class (allele, gene_family, accessions)."""
    if not os.path.exists(catalog_path):
        print(f"[WARN] catalogo nao encontrado: {catalog_path}")
        return {}
    out = {}
    with open(catalog_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            klass = (row.get("class") or "").strip().lower()
            subklass = (row.get("subclass") or "").strip().lower()
            labels = []
            if klass:
                labels.extend([s.strip() for s in klass.split("/") if s.strip()])
            if subklass and subklass != klass:
                labels.extend([s.strip() for s in subklass.split("/") if s.strip()])
            if not labels:
                continue
            labels = list(dict.fromkeys(labels))
            # registra MULTIPLAS chaves possiveis (todas apontando pra mesmos labels)
            for col in ("allele", "gene_family", "refseq_protein_accession",
                        "genbank_protein_accession"):
                v = (row.get(col) or "").strip()
                if v:
                    out.setdefault(v, labels)
    return out


def header_candidates(hdr: str) -> list[str]:
    """Retorna varios identificadores candidatos do header FASTA AMRFinderPlus."""
    parts = hdr.split("|")
    cands = []
    for i in (0, 3, 4, 5, 7):
        if i < len(parts) and parts[i].strip():
            cands.append(parts[i].strip())
    return cands


def clean_seq(seq: str) -> str:
    """Remove stop codon `*` e qualquer caractere fora do alfabeto ESM-2."""
    valid = set("ACDEFGHIKLMNPQRSTVWYX")
    return "".join(c for c in seq.upper() if c in valid)


def load_esm():
    print(f"[ESM-2] {ESM_MODEL_NAME}")
    model, alphabet = esm.pretrained.load_model_and_alphabet(ESM_MODEL_NAME)
    model.eval()
    device = torch.device("cpu")
    model.to(device)
    return model, alphabet, MODEL_LAYERS[ESM_MODEL_NAME], device


def embed_batch(seqs: list[tuple[str, str]], model, alphabet, layer, device) -> np.ndarray:
    bc = alphabet.get_batch_converter()
    trimmed = [(n, s[:ESM_MAX_LEN]) for n, s in seqs]
    _, _, tokens = bc(trimmed)
    tokens = tokens.to(device)
    with torch.no_grad():
        out = model(tokens, repr_layers=[layer], return_contacts=False)
    rep = out["representations"][layer]
    embs = []
    for i, (_, s) in enumerate(trimmed):
        L = len(s)
        v = rep[i, 1 : 1 + L].mean(0).cpu().numpy().astype(np.float32)
        embs.append(v)
    return np.vstack(embs)


def main():
    print("=" * 70)
    print(f"AMRFinderPlus + ESM-2 ({ESM_MODEL_NAME})")
    print("=" * 70)
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.exists(PROT_FASTA):
        print(f"ERRO: {PROT_FASTA} nao existe.")
        sys.exit(1)

    drug_map = load_drug_class_map(CATALOG_TSV)
    print(f"Catalog: {len(drug_map)} alelos com drug_class")

    fasta_entries = parse_amrprot(PROT_FASTA)
    print(f"FASTA: {len(fasta_entries)} proteinas")

    # Parear com labels — testa multiplos candidatos de chave por header
    paired = []
    unmatched = 0
    for hdr, seq in fasta_entries:
        labels = None
        chosen_key = None
        for cand in header_candidates(hdr):
            if cand in drug_map:
                labels = drug_map[cand]
                chosen_key = cand
                break
        if not labels:
            unmatched += 1
            continue
        cleaned = clean_seq(seq)
        if 30 <= len(cleaned) <= 1022:
            paired.append((chosen_key, cleaned, labels))
    print(f"Pareadas: {len(paired)} | sem label: {unmatched}")

    if LIMIT > 0:
        rng = np.random.RandomState(42)
        idx = rng.choice(len(paired), size=min(LIMIT, len(paired)), replace=False)
        paired = [paired[i] for i in idx]
        print(f"Sub-amostra (AMR_LIMIT={LIMIT}): {len(paired)}")

    # Cache: se ja existe, pular embedding
    if os.path.exists(EMB_PATH):
        print(f"[CACHE] {EMB_PATH} ja existe, pulando embedding.")
        z = np.load(EMB_PATH, allow_pickle=True)
        X = z["X"]
        names = list(z["names"])
        labels_per = [list(l) for l in z["labels"]]
    else:
        # Embed em batches
        model, alphabet, layer, device = load_esm()
        X_list, names, labels_per = [], [], []
        t0 = time.time()
        for start in range(0, len(paired), BATCH):
            chunk = paired[start : start + BATCH]
            batch = [(allele, seq) for allele, seq, _ in chunk]
            try:
                E = embed_batch(batch, model, alphabet, layer, device)
            except torch.cuda.OutOfMemoryError:
                print("OOM no batch — cai pra batch=1")
                E_list = []
                for b in batch:
                    E_list.append(embed_batch([b], model, alphabet, layer, device))
                E = np.vstack(E_list)
            X_list.append(E)
            for allele, _, labs in chunk:
                names.append(allele)
                labels_per.append(labs)
            done = start + len(chunk)
            elapsed = time.time() - t0
            eta = elapsed / max(done, 1) * (len(paired) - done)
            print(
                f"  [{done:>5}/{len(paired)}] elapsed={elapsed:6.1f}s  ETA={eta/60:.1f}min",
                flush=True,
            )
        X = np.vstack(X_list)
        np.savez_compressed(
            EMB_PATH,
            X=X,
            names=np.array(names, dtype=object),
            labels=np.array(labels_per, dtype=object),
        )
        print(f"[Saved] {EMB_PATH}  shape={X.shape}")

    # Treino classificador
    mlb = MultiLabelBinarizer()
    Y = mlb.fit_transform(labels_per)
    print(f"Labels: N={len(Y)}  classes={len(mlb.classes_)}")
    print(f"  classes: {list(mlb.classes_)[:20]}{'...' if len(mlb.classes_) > 20 else ''}")

    X_train, X_test, Y_train, Y_test = train_test_split(
        X, Y, test_size=0.2, random_state=42
    )
    clf = OneVsRestClassifier(LogisticRegression(max_iter=2000, C=1.0, n_jobs=1), n_jobs=-1)
    print(f"Treinando OvR LR em N_train={len(X_train)} ...")
    t0 = time.time()
    clf.fit(X_train, Y_train)
    print(f"  treino: {time.time() - t0:.1f}s")

    Y_pred = clf.predict(X_test)
    f1_macro = f1_score(Y_test, Y_pred, average="macro", zero_division=0)
    f1_micro = f1_score(Y_test, Y_pred, average="micro", zero_division=0)
    hloss = hamming_loss(Y_test, Y_pred)
    exact = float(np.mean((Y_test == Y_pred).all(axis=1)))
    print(f"\nHold-out (20%):")
    print(f"  F1 macro    = {f1_macro:.3f}")
    print(f"  F1 micro    = {f1_micro:.3f}")
    print(f"  Hamming acc = {1 - hloss:.3f}")
    print(f"  Exact match = {exact:.3f}")

    # Re-treinar em tudo pra inferencia final
    clf_full = OneVsRestClassifier(LogisticRegression(max_iter=2000, C=1.0))
    clf_full.fit(X, Y)
    joblib.dump({"clf": clf_full, "mlb": mlb, "model": ESM_MODEL_NAME}, CLF_PATH)
    print(f"[Saved] {CLF_PATH}")

    # Aplicar nas 42 variantes do nosso painel (proteinas em reports/protein_scoring/proteins/)
    panel_preds = []
    if os.path.isdir(PAINEL_PROTS_DIR):
        print(f"\nAplicando classificador nas variantes do painel...")
        # Carregar embeddings cacheados do painel (ja gerados em protein_scoring.py)
        emb_dir = os.path.join(REPORTS_DIR, "protein_scoring", "embeddings")
        for fa in sorted(os.listdir(PAINEL_PROTS_DIR)):
            name = fa.replace(".fasta", "")
            cache = os.path.join(emb_dir, f"{name}__{ESM_MODEL_NAME}.npy")
            if not os.path.exists(cache):
                continue
            v = np.load(cache).reshape(1, -1)
            pred = clf_full.predict(v)[0]
            try:
                proba = clf_full.predict_proba(v)[0]
                top_idx = np.argsort(-proba)[:5]
                top = [(mlb.classes_[k], round(float(proba[k]), 3)) for k in top_idx]
            except Exception:
                top = []
            classes_pred = [mlb.classes_[k] for k in np.where(pred == 1)[0]]
            panel_preds.append(
                {
                    "variant": name,
                    "predicted_classes": ";".join(classes_pred) or "(none)",
                    "top5_with_proba": str(top),
                }
            )

        if panel_preds:
            with open(PRED_PATH, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(panel_preds[0].keys()))
                w.writeheader()
                w.writerows(panel_preds)
            print(f"[Saved] {PRED_PATH}  ({len(panel_preds)} variantes)")

    metrics = {
        "model": ESM_MODEL_NAME,
        "n_total": int(len(Y)),
        "n_classes": int(len(mlb.classes_)),
        "classes": list(mlb.classes_),
        "f1_macro_holdout": round(f1_macro, 4),
        "f1_micro_holdout": round(f1_micro, 4),
        "hamming_acc_holdout": round(1 - hloss, 4),
        "exact_match_holdout": round(exact, 4),
        "panel_variants_predicted": len(panel_preds),
    }
    with open(METRICS_PATH, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"[Saved] {METRICS_PATH}")
    print("\nFeito.")


if __name__ == "__main__":
    main()
