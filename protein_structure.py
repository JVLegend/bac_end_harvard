#!/usr/bin/env python3
"""
ESMFold via HuggingFace transformers — predicao de estrutura 3D para as 42 variantes AMR.
Gera arquivos .pdb para cada variante, renderizaveis com NGL.js no dashboard.

Hardware: M2 16GB OK pra proteinas ate ~400aa em CPU (lento, ~2-5min/proteina).
Modelo ~2.5GB no HuggingFace cache.

Uso:
    /Users/iaparamedicos/envs/dev/bin/python protein_structure.py
    PROTEINS=blaKPC-3,blaNDM-1 /Users/iaparamedicos/envs/dev/bin/python protein_structure.py
    ESMFOLD_MAX_LEN=300 ... (proteinas menores rodam mais rapido)
"""

from __future__ import annotations

import os
import sys
import time

import torch
from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

from config import REPORTS_DIR

PROTEINS_DIR = os.path.join(REPORTS_DIR, "protein_scoring", "proteins")
PDB_DIR = os.path.join(REPORTS_DIR, "protein_structure", "pdbs")

MAX_LEN = int(os.environ.get("ESMFOLD_MAX_LEN", "350"))
ONLY = os.environ.get("PROTEINS", "").split(",") if os.environ.get("PROTEINS") else None


def parse_fasta(path: str) -> str:
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq)


def convert_outputs_to_pdb(outputs):
    """Converte saida do ESMFold (HuggingFace) em string PDB."""
    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
    outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
    final_atom_positions = final_atom_positions.cpu().numpy()
    final_atom_mask = outputs["atom37_atom_exists"]
    pdbs = []
    for i in range(outputs["aatype"].shape[0]):
        aa = outputs["aatype"][i]
        pred_pos = final_atom_positions[i]
        mask = final_atom_mask[i]
        resid = outputs["residue_index"][i] + 1
        pred = OFProtein(
            aatype=aa,
            atom_positions=pred_pos,
            atom_mask=mask,
            residue_index=resid,
            b_factors=outputs["plddt"][i],
            chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
        )
        pdbs.append(to_pdb(pred))
    return pdbs


def main():
    print("=" * 70)
    print("ESMFold (HuggingFace transformers) — 3D para variantes AMR")
    print("=" * 70)
    os.makedirs(PDB_DIR, exist_ok=True)

    if not os.path.isdir(PROTEINS_DIR):
        print(f"ERRO: rode protein_scoring.py antes (precisa de {PROTEINS_DIR}).")
        sys.exit(1)

    print("[ESMFold] carregando facebook/esmfold_v1 (download ~2.5GB na 1a vez)...")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1", low_cpu_mem_usage=True, torch_dtype=torch.float32
    )
    model.eval()
    device = torch.device("cpu")
    model.to(device)
    # Reduz memoria pra M2 16GB
    model.esm = model.esm.float()
    model.trunk.set_chunk_size(64)
    print(f"[ESMFold] pronto. device={device}, chunk=64")

    files = sorted(os.listdir(PROTEINS_DIR))
    if ONLY:
        wanted = set(ONLY)
        files = [f for f in files if f.replace(".fasta", "") in wanted]

    print(f"\nProcessando {len(files)} variantes (MAX_LEN={MAX_LEN}aa)\n")
    done = 0
    for fa in files:
        name = fa.replace(".fasta", "")
        out = os.path.join(PDB_DIR, f"{name}.pdb")
        if os.path.exists(out):
            print(f"  [SKIP] {name} ja tem PDB")
            continue
        seq = parse_fasta(os.path.join(PROTEINS_DIR, fa))
        if len(seq) > MAX_LEN:
            print(f"  [SKIP] {name} L={len(seq)} > MAX_LEN={MAX_LEN}")
            continue
        if len(seq) < 30:
            continue

        print(f"  [{name}] L={len(seq)}aa  predicting...", end=" ", flush=True)
        t0 = time.time()
        try:
            inputs = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
            inputs = {k: v.to(device) for k, v in inputs.items()}
            with torch.no_grad():
                outputs = model(**inputs)
            pdbs = convert_outputs_to_pdb(outputs)
            with open(out, "w") as f:
                f.write(pdbs[0])
            avg_plddt = float(outputs["plddt"].mean().item())
            print(f"{time.time() - t0:.1f}s  pLDDT={avg_plddt:.1f}")
            done += 1
        except Exception as e:
            print(f"FAIL: {type(e).__name__}: {str(e)[:80]}")

    print(f"\n[Saved] {done} PDBs em {PDB_DIR}")


if __name__ == "__main__":
    main()
