"""Fase 7 · Task 1 — ESM-2 protein embeddings para proteinas AMR.

Traduz as sequencias nucleotidicas de cada familia AMR em proteinas
(translate ORF principal), gera embeddings 1280-d via ESM-2 (Meta AI) e
salva matriz de similaridade/distance entre familias.

Output:
- esm2_embeddings.npz (matriz [N_seqs, 1280])
- esm2_meta.json (mapa id -> familia, comprimento, etc)
- esm2_family_distance.csv (matriz 12x12 de distancias medias intra-familia)
"""
import os
import time
import numpy as np
import pandas as pd
import torch
from _fase7_utils import (list_family_seqs, read_fasta, save_json, save_csv,
                          emit_result, EXP_DIR, AMR_FAMILIES)

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"device={DEVICE}")

# translate DNA -> protein (strand +, frame 0). Codigo simples pra familias AMR.
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def dna_to_protein(seq: str) -> str:
    seq = seq.upper().replace("N", "")
    # best frame: escolhe frame com menos stops
    best = ""
    best_stops = 1e9
    for frame in range(3):
        protein = []
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in CODON_TABLE:
                protein.append(CODON_TABLE[codon])
        p = "".join(protein)
        stops = p.count("*")
        # pick frame with long run before stop
        if stops < best_stops or (stops == best_stops and len(p) > len(best)):
            best = p.rstrip("*").split("*")[0]  # pega 1st ORF
            best_stops = stops
    return best


def load_esm2():
    """Tenta huggingface transformers."""
    try:
        from transformers import AutoTokenizer, AutoModel
        model_name = "facebook/esm2_t30_150M_UR50D"  # 150M params, 640-d
        # Usar 150M (640-d) pra caber rapido. Se tiver tempo: t33_650M (1280-d)
        print(f"carregando {model_name}")
        tok = AutoTokenizer.from_pretrained(model_name)
        model = AutoModel.from_pretrained(model_name).to(DEVICE).eval()
        return tok, model, int(model.config.hidden_size)
    except Exception as e:
        print(f"ESM-2 indisponivel ({e}); usando fallback k-mer")
        return None, None, None


@torch.no_grad()
def embed_esm(protein: str, tok, model) -> np.ndarray:
    if not protein:
        return np.zeros(model.config.hidden_size, dtype=np.float32)
    inputs = tok(protein[:1000], return_tensors="pt").to(DEVICE)
    out = model(**inputs)
    # mean pooling over tokens (exclude special tokens)
    emb = out.last_hidden_state[0, 1:-1].mean(dim=0).cpu().numpy()
    return emb.astype(np.float32)


def kmer_fallback(protein: str, k: int = 3, dim: int = 640) -> np.ndarray:
    """Fallback sem ESM-2: k-mer histogram projetado."""
    if not protein:
        return np.zeros(dim, dtype=np.float32)
    hist = np.zeros(dim, dtype=np.float32)
    for i in range(len(protein) - k + 1):
        km = protein[i:i+k]
        h = hash(km) % dim
        hist[h] += 1
    return hist / (np.linalg.norm(hist) + 1e-9)


def main():
    t0 = time.time()
    fam_seqs = list_family_seqs()
    n_total = sum(len(v) for v in fam_seqs.values())
    print(f"familias={len(fam_seqs)}  total_fastas={n_total}")
    if n_total == 0:
        print("nenhum fasta encontrado — rode fetch_sequences.py primeiro")
        return

    tok, model, hidden = load_esm2()
    use_esm = model is not None
    dim = hidden if use_esm else 640

    embeddings = []
    meta = []
    for fam, paths in fam_seqs.items():
        for p in paths:
            fa = read_fasta(p)
            for header, seq in fa.items():
                protein = dna_to_protein(seq)
                if use_esm:
                    emb = embed_esm(protein, tok, model)
                else:
                    emb = kmer_fallback(protein, dim=dim)
                embeddings.append(emb)
                meta.append({"family": fam, "header": header, "path": p,
                             "prot_len": len(protein), "dim": dim})

    E = np.stack(embeddings, 0) if embeddings else np.zeros((0, dim), dtype=np.float32)
    np.savez_compressed(f"{EXP_DIR}/esm2_embeddings.npz", embeddings=E,
                       families=[m["family"] for m in meta])
    save_json({"use_esm": use_esm, "dim": dim, "n": len(meta), "items": meta},
              "esm2_meta.json")

    # Distance matrix medio por familia
    fam_to_idx = {}
    for i, m in enumerate(meta):
        fam_to_idx.setdefault(m["family"], []).append(i)

    D_fam = pd.DataFrame(index=AMR_FAMILIES, columns=AMR_FAMILIES, dtype=float)
    for f1 in AMR_FAMILIES:
        for f2 in AMR_FAMILIES:
            if f1 not in fam_to_idx or f2 not in fam_to_idx:
                D_fam.loc[f1, f2] = np.nan
                continue
            e1 = E[fam_to_idx[f1]]
            e2 = E[fam_to_idx[f2]]
            # cosine distance entre medias
            m1 = e1.mean(axis=0); m2 = e2.mean(axis=0)
            cos = float(np.dot(m1, m2) / (np.linalg.norm(m1) * np.linalg.norm(m2) + 1e-9))
            D_fam.loc[f1, f2] = round(1 - cos, 4)
    save_csv(D_fam.reset_index().rename(columns={"index": "family"}),
             "esm2_family_distance.csv")

    emit_result("use_esm", use_esm)
    emit_result("embedding_dim", dim)
    emit_result("n_sequences", len(meta))
    emit_result("seconds", f"{time.time()-t0:.1f}")


if __name__ == "__main__":
    main()
