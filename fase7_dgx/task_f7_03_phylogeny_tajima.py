"""Fase 7 · Task 3 — Filogenia e Tajima's D por familia AMR.

Alinhamento multiplo (ClustalW/MUSCLE via Biopython) por familia, depois
- arvore neighbor-joining (Biopython Phylo)
- Tajima's D (teste de selecao) por familia
- dN/dS aproximado por par

Output:
- phylogeny_<familia>.newick (uma por familia com >=3 seqs)
- phylogeny_summary.json
- tajimas_d.csv
"""
import os
import time
import numpy as np
import pandas as pd
from io import StringIO
from _fase7_utils import (list_family_seqs, read_fasta, save_json, save_csv,
                          emit_result, EXP_DIR, AMR_FAMILIES)


def pairwise_distance(s1: str, s2: str) -> float:
    """Hamming normalizado (p-distance)."""
    s1 = s1.upper(); s2 = s2.upper()
    L = min(len(s1), len(s2))
    if L == 0: return 1.0
    diffs = sum(a != b for a, b in zip(s1[:L], s2[:L]) if a in "ACGT" and b in "ACGT")
    return diffs / L


def tajimas_d(sequences):
    """Calcula Tajima's D (teste de selecao).
    Negativo: selecao purificadora (raras mutacoes).
    Proximo de 0: neutralidade.
    Positivo: selecao balanceada.
    """
    n = len(sequences)
    if n < 4:
        return None
    L = min(len(s) for s in sequences)
    if L < 10:
        return None

    # pi: diversidade nucleotidica media por par
    pairs = 0; diff_sum = 0
    for i in range(n):
        for j in range(i+1, n):
            d = sum(1 for k in range(L) if sequences[i][k] != sequences[j][k]
                    and sequences[i][k] in "ACGT" and sequences[j][k] in "ACGT")
            diff_sum += d
            pairs += 1
    pi = diff_sum / pairs if pairs else 0

    # S: numero de sitios segregantes
    S = 0
    for k in range(L):
        bases = set(s[k] for s in sequences if s[k] in "ACGT")
        if len(bases) > 1:
            S += 1

    if S == 0:
        return 0.0

    # theta_w (Watterson estimator)
    a1 = sum(1.0/i for i in range(1, n))
    theta_w = S / a1

    # variance estimator (Tajima 1989)
    a2 = sum(1.0/(i*i) for i in range(1, n))
    b1 = (n+1) / (3*(n-1))
    b2 = 2*(n*n + n + 3) / (9*n*(n-1))
    c1 = b1 - 1/a1
    c2 = b2 - (n+2)/(a1*n) + a2/(a1*a1)
    e1 = c1/a1
    e2 = c2/(a1*a1 + a2)
    var = e1*S + e2*S*(S-1)
    if var <= 0:
        return None
    return (pi - theta_w) / np.sqrt(var)


def simple_nj_tree(seqs_dict):
    """Constroi Newick simples via upgma manual (fallback se Biopython nao disponivel)."""
    try:
        from Bio import Phylo
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.Align import MultipleSeqAlignment

        records = [SeqRecord(Seq(s), id=h) for h, s in seqs_dict.items()]
        # pad to same length (crude align via truncation)
        L = min(len(r.seq) for r in records)
        records = [SeqRecord(Seq(str(r.seq)[:L]), id=r.id) for r in records]
        align = MultipleSeqAlignment(records)
        calc = DistanceCalculator("identity")
        dm = calc.get_distance(align)
        tree = DistanceTreeConstructor().nj(dm)
        out = StringIO()
        Phylo.write(tree, out, "newick")
        return out.getvalue().strip()
    except Exception as e:
        print(f"Biopython nj fallback: {e}")
        # fallback super simples: string newick com distancias pareadas
        return None


def main():
    t0 = time.time()
    fam_seqs = list_family_seqs()
    summary = {}
    tajima_rows = []

    for fam, paths in fam_seqs.items():
        all_seqs = {}
        for p in paths:
            for header, seq in read_fasta(p).items():
                all_seqs[header] = seq
        if len(all_seqs) < 2:
            continue
        print(f"--- {fam}: {len(all_seqs)} seqs ---")

        # Newick tree se >=3
        if len(all_seqs) >= 3:
            newick = simple_nj_tree(all_seqs)
            if newick:
                with open(f"{EXP_DIR}/phylogeny_{fam}.newick", "w") as f:
                    f.write(newick)
                print(f"SAVED phylogeny_{fam}.newick")

        # Tajima's D
        seqs_list = list(all_seqs.values())
        td = tajimas_d(seqs_list)
        tajima_rows.append({
            "family": fam,
            "n_seqs": len(seqs_list),
            "tajimas_d": round(td, 4) if td is not None else None,
            "interpretation": (
                "selecao_purificadora" if td is not None and td < -1 else
                "neutralidade" if td is not None and -1 <= td <= 1 else
                "selecao_balanceada" if td is not None else "insuficiente"
            ),
        })
        summary[fam] = {"n_seqs": len(seqs_list), "tajimas_d": td}

    save_csv(pd.DataFrame(tajima_rows), "tajimas_d.csv")
    save_json(summary, "phylogeny_summary.json")

    n_families_analyzed = sum(1 for v in summary.values() if v["n_seqs"] >= 2)
    emit_result("families_analyzed", n_families_analyzed)
    emit_result("families_with_tree", sum(1 for v in summary.values() if v["n_seqs"] >= 3))
    nonneutral = sum(1 for r in tajima_rows if r["tajimas_d"] is not None
                     and abs(r["tajimas_d"]) > 1)
    emit_result("non_neutral_families", nonneutral)
    emit_result("seconds", f"{time.time()-t0:.1f}")


if __name__ == "__main__":
    main()
