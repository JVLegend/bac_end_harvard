#!/usr/bin/env python3
"""
ESMFold — predicao de estrutura 3D para as 42 variantes AMR do painel.
Gera arquivos .pdb para cada variante. Renderizaveis com NGL.js no dashboard.

Hardware: M2 16GB OK pra proteinas ate ~400aa (CPU). >500aa pode dar OOM.
Tempo: ~30-90s por proteina em CPU M2.

Uso:
    /Users/iaparamedicos/envs/dev/bin/python protein_structure.py
    PROTEINS=blaKPC-3,blaNDM-1 /Users/iaparamedicos/envs/dev/bin/python protein_structure.py
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import torch
import esm

from config import REPORTS_DIR

PROTEINS_DIR = os.path.join(REPORTS_DIR, "protein_scoring", "proteins")
PDB_DIR = os.path.join(REPORTS_DIR, "protein_structure", "pdbs")

MAX_LEN = int(os.environ.get("ESMFOLD_MAX_LEN", "400"))
ONLY = os.environ.get("PROTEINS", "").split(",") if os.environ.get("PROTEINS") else None


def parse_fasta(path: str) -> str:
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq)


def main():
    print("=" * 70)
    print("ESMFold — predicao 3D para variantes AMR (M2 CPU)")
    print("=" * 70)
    os.makedirs(PDB_DIR, exist_ok=True)

    if not os.path.isdir(PROTEINS_DIR):
        print(f"ERRO: rode protein_scoring.py antes (precisa de {PROTEINS_DIR}).")
        sys.exit(1)

    print("[ESMFold] carregando modelo (download na 1a vez ~3GB)...")
    model = esm.pretrained.esmfold_v1()
    model.eval()
    # M2 16GB: usar CPU; MPS tem bugs com ESMFold
    device = torch.device("cpu")
    model.to(device)
    # Otimizacoes pra CPU
    model.set_chunk_size(64)
    print(f"[ESMFold] pronto. device={device}, chunk=64")

    files = sorted(os.listdir(PROTEINS_DIR))
    if ONLY:
        wanted = set(ONLY)
        files = [f for f in files if f.replace(".fasta", "") in wanted]

    for fa in files:
        name = fa.replace(".fasta", "")
        out = os.path.join(PDB_DIR, f"{name}.pdb")
        if os.path.exists(out):
            print(f"  [SKIP] {name} ja tem PDB")
            continue
        seq = parse_fasta(os.path.join(PROTEINS_DIR, fa))
        if len(seq) > MAX_LEN:
            print(f"  [SKIP] {name} L={len(seq)} > MAX_LEN={MAX_LEN} (M2 OOM risk)")
            continue
        if len(seq) < 30:
            continue

        print(f"  [{name}] L={len(seq)}aa  predicting...", end=" ", flush=True)
        t0 = time.time()
        with torch.no_grad():
            output = model.infer_pdb(seq)
        with open(out, "w") as f:
            f.write(output)
        # plDDT mean (qualidade da predicao)
        plddt_lines = [l for l in output.split("\n") if l.startswith("ATOM") and l[13:15] == "CA"]
        if plddt_lines:
            try:
                avg_plddt = sum(float(l[60:66]) for l in plddt_lines) / len(plddt_lines)
                print(f"{time.time() - t0:.1f}s  pLDDT={avg_plddt:.1f}")
            except Exception:
                print(f"{time.time() - t0:.1f}s")
        else:
            print(f"{time.time() - t0:.1f}s")

    print(f"\n[Saved] PDBs em {PDB_DIR}")


if __name__ == "__main__":
    main()
