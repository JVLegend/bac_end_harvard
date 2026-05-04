#!/usr/bin/env python3
"""
Fase 5 — ProtT5 ensemble.

Gera embeddings ProtT5 (Rostlab/prot_t5_xl_uniref50) para as variantes do painel
e concatena com ESM-2 pra criar um phenotype probe ensemble (1024 + 640 = 1664 dim).

Hardware: ProtT5 XL ~2.5GB. Em M2 16GB CPU roda ~5-15s por proteina.
Uso:
    /Users/iaparamedicos/envs/dev/bin/python prott5_ensemble.py
"""

from __future__ import annotations

import csv
import json
import os
import time
import warnings

import numpy as np
import pandas as pd
import torch
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import MultiLabelBinarizer
from transformers import T5EncoderModel, T5Tokenizer

from config import REPORTS_DIR

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning)

PROTEINS_DIR = os.path.join(REPORTS_DIR, "protein_scoring", "proteins")
ESM_EMB_DIR = os.path.join(REPORTS_DIR, "protein_scoring", "embeddings")
OUT_DIR = os.path.join(REPORTS_DIR, "prott5_ensemble")
PROT_T5_EMB_DIR = os.path.join(OUT_DIR, "embeddings")

ESM_MODEL_NAME = os.environ.get("ESM_MODEL", "esm2_t30_150M_UR50D")
PT5_MODEL = "Rostlab/prot_t5_xl_uniref50"


def parse_fasta(path: str) -> str:
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            out.append(line.strip())
    return "".join(out)


def load_t5():
    print(f"[ProtT5] carregando {PT5_MODEL} (~2.5GB)...")
    tokenizer = T5Tokenizer.from_pretrained(PT5_MODEL, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(PT5_MODEL)
    model.eval()
    return tokenizer, model


def embed_t5(seq: str, tokenizer, model) -> np.ndarray:
    # ProtT5 espera sequencia com espacos entre residuos e qq letra fora ACDEFGHIKLMNPQRSTVWY -> X
    valid = set("ACDEFGHIKLMNPQRSTVWY")
    cleaned = "".join(c if c in valid else "X" for c in seq.upper())
    s = " ".join(list(cleaned))
    enc = tokenizer(s, add_special_tokens=True, padding=True, return_tensors="pt")
    with torch.no_grad():
        out = model(input_ids=enc["input_ids"], attention_mask=enc["attention_mask"])
    # ultima dim: hidden_state. shape: [1, T, D]. T = len(cleaned) + 1 (EOS no final)
    rep = out.last_hidden_state[0, : len(cleaned)]
    return rep.mean(0).cpu().numpy().astype(np.float32)


# Reaproveita labels do phenotype_probe original
FAMILY_LABELS = {
    "mecA":       (["penam", "cephalosporin", "methicillin"], "antibiotic target replacement"),
    "mecA1":      (["penam", "cephalosporin", "methicillin"], "antibiotic target replacement"),
    "mecA2":      (["penam", "cephalosporin", "methicillin"], "antibiotic target replacement"),
    "blaKPC":     (["carbapenem", "cephalosporin", "penam"], "antibiotic inactivation"),
    "blaNDM":     (["carbapenem", "cephalosporin", "penam"], "antibiotic inactivation"),
    "blaOXA":     (["carbapenem", "cephalosporin", "penam"], "antibiotic inactivation"),
    "blaVIM":     (["carbapenem", "cephalosporin", "penam"], "antibiotic inactivation"),
    "blaIMP":     (["carbapenem", "cephalosporin", "penam"], "antibiotic inactivation"),
    "blaGES":     (["carbapenem", "cephalosporin"],          "antibiotic inactivation"),
    "blaCTX-M":   (["cephalosporin", "penam"],               "antibiotic inactivation"),
    "blaCTX-M-15":(["cephalosporin", "penam"],               "antibiotic inactivation"),
    "vanA":       (["glycopeptide"],                         "antibiotic target alteration"),
    "mcr-1":      (["polymyxin", "peptide"],                 "antibiotic target alteration"),
    "mcr-5":      (["polymyxin", "peptide"],                 "antibiotic target alteration"),
    "qnrS":       (["fluoroquinolone"],                      "antibiotic target protection"),
    "armA":       (["aminoglycoside"],                       "antibiotic target alteration"),
}


def labels_for(name: str):
    base = name.split("-")[0] if name.startswith("bla") else name
    if name in FAMILY_LABELS:
        return FAMILY_LABELS[name]
    for k in FAMILY_LABELS:
        if name.startswith(k):
            return FAMILY_LABELS[k]
    return None


def main():
    print("=" * 70)
    print("FASE 5 — ProtT5 ensemble (ProtT5 + ESM-2 concat)")
    print("=" * 70)
    os.makedirs(PROT_T5_EMB_DIR, exist_ok=True)

    if not os.path.isdir(PROTEINS_DIR):
        print("ERRO: rode protein_scoring.py antes.")
        return

    tokenizer, model = load_t5()
    files = sorted(os.listdir(PROTEINS_DIR))

    print(f"\nGerando embeddings ProtT5 para {len(files)} variantes...")
    t0 = time.time()
    for i, fa in enumerate(files, 1):
        name = fa.replace(".fasta", "")
        cache = os.path.join(PROT_T5_EMB_DIR, f"{name}.npy")
        if os.path.exists(cache):
            continue
        seq = parse_fasta(os.path.join(PROTEINS_DIR, fa))
        if len(seq) < 30 or len(seq) > 1500:
            continue
        e = embed_t5(seq, tokenizer, model)
        np.save(cache, e)
        print(f"  [{i:>2}/{len(files)}] {name:<18} L={len(seq):>4}aa  dim={e.shape[0]}")
    print(f"[ProtT5] {time.time() - t0:.1f}s")

    # Treino ensemble
    X, names, drug_y, mech_y = [], [], [], []
    for fa in files:
        name = fa.replace(".fasta", "")
        labs = labels_for(name)
        if labs is None:
            continue
        esm_path = os.path.join(ESM_EMB_DIR, f"{name}__{ESM_MODEL_NAME}.npy")
        t5_path = os.path.join(PROT_T5_EMB_DIR, f"{name}.npy")
        if not (os.path.exists(esm_path) and os.path.exists(t5_path)):
            continue
        e = np.concatenate([np.load(esm_path), np.load(t5_path)])
        X.append(e)
        names.append(name)
        drug_y.append(labs[0])
        mech_y.append(labs[1])

    if len(X) < 5:
        print(f"Amostras insuficientes ({len(X)})")
        return
    X = np.vstack(X)
    print(f"\nEnsemble dim = {X.shape[1]}  N={len(X)}")

    mlb = MultiLabelBinarizer()
    Y = mlb.fit_transform(drug_y)

    loo = LeaveOneOut()
    Y_loo = np.zeros_like(Y)
    for tr, te in loo.split(X):
        c = OneVsRestClassifier(LogisticRegression(max_iter=2000, C=1.0))
        c.fit(X[tr], Y[tr])
        Y_loo[te] = c.predict(X[te])
    hamming = float(np.mean(Y_loo == Y))
    exact = float(np.mean((Y_loo == Y).all(axis=1)))

    f1_per = {}
    for j, cls in enumerate(mlb.classes_):
        tp = int(((Y_loo[:, j] == 1) & (Y[:, j] == 1)).sum())
        fp = int(((Y_loo[:, j] == 1) & (Y[:, j] == 0)).sum())
        fn = int(((Y_loo[:, j] == 0) & (Y[:, j] == 1)).sum())
        p = tp / (tp + fp) if (tp + fp) else 0
        r = tp / (tp + fn) if (tp + fn) else 0
        f1 = 2 * p * r / (p + r) if (p + r) else 0
        f1_per[cls] = round(f1, 3)

    metrics = {
        "ensemble": "ESM-2 + ProtT5 (concat)",
        "esm_model": ESM_MODEL_NAME,
        "prott5_model": PT5_MODEL,
        "n": int(len(X)),
        "dim": int(X.shape[1]),
        "drug_hamming_acc_loo": round(hamming, 4),
        "drug_exact_match_loo": round(exact, 4),
        "drug_per_class_f1_loo": f1_per,
    }
    out_json = os.path.join(OUT_DIR, "metrics.json")
    with open(out_json, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\n[Saved] {out_json}")
    print(f"  hamming LOO = {hamming:.3f}")
    print(f"  exact match = {exact:.3f}")


if __name__ == "__main__":
    main()
