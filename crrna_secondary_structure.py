#!/usr/bin/env python3
"""
Fase 6 — Estrutura secundaria do crRNA via ViennaRNA (RNAfold).

Avalia hairpins/free-energy de cada guide Cas12a candidato:
  - dG da estrutura MFE
  - probabilidade de hairpin na seed region (1-8nt)
  - score de penalidade pra integrar em covariance_probes.py

Requer ViennaRNA instalado: brew install viennarna (binario RNAfold no PATH)

Uso:
    /Users/iaparamedicos/envs/dev/bin/python crrna_secondary_structure.py
"""

from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
import sys

from config import REPORTS_DIR, GUIDES_DIR

OUT_DIR = os.path.join(REPORTS_DIR, "crrna_structure")
LB_DR = "AATTTCTACTAAGTGTAGAT"  # LbCas12a direct repeat (scaffold)


def run_rnafold(seq: str) -> tuple[str, float]:
    """Retorna (notacao_dot-bracket, dG_mfe). Espera RNAfold no PATH."""
    p = subprocess.run(
        ["RNAfold", "--noPS", "-d2"],
        input=seq.replace("T", "U"),
        capture_output=True,
        text=True,
        timeout=20,
    )
    out = p.stdout.strip().split("\n")
    # linha 2: ".(((...)))..  ( -3.40)"
    m = re.match(r"^([().]+)\s+\(\s*(-?\d+\.\d+)\)", out[1]) if len(out) >= 2 else None
    if not m:
        return "", 0.0
    return m.group(1), float(m.group(2))


def seed_paired_fraction(notation: str, seed_start: int, seed_len: int = 8) -> float:
    """Fracao de bases pareadas na seed region."""
    seed = notation[seed_start : seed_start + seed_len]
    if not seed:
        return 0.0
    return sum(1 for c in seed if c != ".") / len(seed)


def main():
    print("=" * 70)
    print("FASE 6 — Estrutura secundaria crRNA (ViennaRNA)")
    print("=" * 70)

    if shutil.which("RNAfold") is None:
        print("ERRO: RNAfold nao encontrado.")
        print("       Instale: brew install viennarna")
        sys.exit(1)

    if not os.path.isdir(GUIDES_DIR):
        print(f"ERRO: {GUIDES_DIR} nao existe. Rode design_guides.py antes.")
        sys.exit(1)

    os.makedirs(OUT_DIR, exist_ok=True)

    rows = []
    for fname in sorted(os.listdir(GUIDES_DIR)):
        if not fname.endswith(".tsv"):
            continue
        path = os.path.join(GUIDES_DIR, fname)
        gene = fname.replace("_guides.tsv", "").replace(".tsv", "")

        with open(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                spacer = row.get("spacer") or row.get("guide_seq") or row.get("sequence")
                if not spacer:
                    continue
                # Estrutura completa (DR + spacer)
                full = LB_DR + spacer
                notation_full, dG_full = run_rnafold(full)
                # Estrutura so do spacer (auto-fold)
                notation_sp, dG_sp = run_rnafold(spacer)

                # Seed = pos 1-8 do spacer (apos DR no full)
                seed_pair_full = seed_paired_fraction(notation_full, len(LB_DR), 8)
                seed_pair_sp = seed_paired_fraction(notation_sp, 0, 8)

                # Score de penalidade: -dG forte + seed pareada = ruim
                penalty = max(0.0, (-dG_full - 5)) + 5 * seed_pair_full
                rows.append(
                    {
                        "gene": gene,
                        "spacer": spacer,
                        "dG_full_kcal": round(dG_full, 2),
                        "dG_spacer_only": round(dG_sp, 2),
                        "seed_paired_full": round(seed_pair_full, 2),
                        "seed_paired_spacer": round(seed_pair_sp, 2),
                        "structure_penalty": round(penalty, 2),
                        "notation_full": notation_full,
                    }
                )
                print(
                    f"  {gene:<14} dG={dG_full:>5.1f}  seed_paired={seed_pair_full:.2f}  "
                    f"penalty={penalty:.2f}"
                )

    if not rows:
        print("Nenhum guide encontrado.")
        return

    out_csv = os.path.join(OUT_DIR, "crrna_structure_scores.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n[Saved] {out_csv}  ({len(rows)} guides)")

    summary = {
        "n_guides": len(rows),
        "mean_dG_full": round(sum(r["dG_full_kcal"] for r in rows) / len(rows), 2),
        "high_penalty_count": sum(1 for r in rows if r["structure_penalty"] > 5),
    }
    out_json = os.path.join(OUT_DIR, "summary.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[Saved] {out_json}")


if __name__ == "__main__":
    main()
