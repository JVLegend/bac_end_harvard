#!/usr/bin/env python3
"""
Script 3: Design de primers RPA para amplificação isotérmica.
Primers flanqueiam o site do crRNA selecionado.
"""

import os
import csv

from config import TARGETS, SEQUENCES_DIR, GUIDES_DIR, PRIMERS_DIR, RPA
from utils import parse_fasta, gc_content, tm_basic, max_homopolymer, self_complementarity_score


def load_best_guide(gene_name: str) -> dict | None:
    """Carrega o melhor guide do TSV."""
    filepath = os.path.join(GUIDES_DIR, f"{gene_name}_cas12a_guides.tsv")
    if not os.path.exists(filepath):
        return None

    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            return row  # primeiro = melhor
    return None


def score_primer(seq: str) -> dict:
    """Pontua um primer candidato."""
    gc = gc_content(seq)
    tm = tm_basic(seq)
    homo = max_homopolymer(seq)
    self_comp = self_complementarity_score(seq)

    score = 100.0

    # GC
    if RPA["gc_min"] <= gc <= RPA["gc_max"]:
        score += 10
    else:
        score -= 20

    # Tm
    if RPA["tm_min"] <= tm <= RPA["tm_max"]:
        score += 15
    else:
        score -= abs(tm - 60) * 2

    # Homopolímero
    if homo >= 5:
        score -= 30
    elif homo >= 4:
        score -= 10

    # Auto-complementaridade
    score -= self_comp * 3

    # Evitar G/C no 3' end excessivo (clamp)
    if seq[-3:].count("G") + seq[-3:].count("C") >= 3:
        score -= 10

    return {
        "score": round(score, 1),
        "gc": round(gc, 3),
        "tm": round(tm, 1),
        "homopolymer": homo,
        "self_comp": self_comp,
    }


def design_primers_for_target(
    gene_name: str, sequence: str, guide_position: int
) -> list[dict]:
    """Desenha primers RPA flanqueando o site do guide."""
    seq = sequence.upper()
    gene_len = len(seq)
    primer_len = RPA["primer_length_optimal"]

    # Definir janela de amplificação em torno do guide
    amplicon_target = RPA["amplicon_optimal"]
    flank = amplicon_target // 2

    # Região upstream para primer forward
    fwd_region_start = max(0, guide_position - flank - primer_len)
    fwd_region_end = max(0, guide_position - 10)  # pelo menos 10bp antes do guide

    # Região downstream para primer reverse
    rev_region_start = min(gene_len, guide_position + 24 + 10)  # 10bp depois do spacer
    rev_region_end = min(gene_len, guide_position + 24 + flank + primer_len)

    candidates = []

    # Gerar primers forward
    for fwd_start in range(fwd_region_start, fwd_region_end - primer_len + 1, 3):
        fwd_seq = seq[fwd_start : fwd_start + primer_len]
        if len(fwd_seq) < primer_len or "N" in fwd_seq:
            continue

        fwd_scores = score_primer(fwd_seq)

        # Gerar primers reverse (reverse complement da região downstream)
        for rev_end in range(rev_region_start + primer_len, rev_region_end + 1, 3):
            rev_start = rev_end - primer_len
            if rev_start < 0 or rev_end > gene_len:
                continue

            # Primer reverse = reverse complement
            rev_template = seq[rev_start:rev_end]
            if len(rev_template) < primer_len or "N" in rev_template:
                continue

            rev_seq = rev_template[::-1].translate(str.maketrans("ATCG", "TAGC"))
            rev_scores = score_primer(rev_seq)

            # Calcular amplicon
            amplicon_size = rev_end - fwd_start
            if amplicon_size < RPA["amplicon_min"] or amplicon_size > RPA["amplicon_max"]:
                continue

            # Score combinado
            pair_score = (fwd_scores["score"] + rev_scores["score"]) / 2
            # Bonus por Tm match
            tm_diff = abs(fwd_scores["tm"] - rev_scores["tm"])
            pair_score -= tm_diff * 2
            # Bonus por amplicon próximo do ideal
            pair_score -= abs(amplicon_size - RPA["amplicon_optimal"]) * 0.5

            candidates.append(
                {
                    "gene": gene_name,
                    "fwd_seq": fwd_seq,
                    "fwd_start": fwd_start,
                    "fwd_tm": fwd_scores["tm"],
                    "fwd_gc": fwd_scores["gc"],
                    "rev_seq": rev_seq,
                    "rev_start": rev_start,
                    "rev_tm": rev_scores["tm"],
                    "rev_gc": rev_scores["gc"],
                    "amplicon_size": amplicon_size,
                    "tm_diff": round(tm_diff, 1),
                    "pair_score": round(pair_score, 1),
                }
            )

    candidates.sort(key=lambda x: x["pair_score"], reverse=True)
    return candidates[:5]


def save_primers(gene_name: str, primers: list[dict]):
    """Salva primers em TSV."""
    filepath = os.path.join(PRIMERS_DIR, f"{gene_name}_rpa_primers.tsv")
    fieldnames = [
        "rank",
        "gene",
        "fwd_seq",
        "fwd_start",
        "fwd_tm",
        "fwd_gc",
        "rev_seq",
        "rev_start",
        "rev_tm",
        "rev_gc",
        "amplicon_size",
        "tm_diff",
        "pair_score",
    ]

    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for i, p in enumerate(primers, 1):
            p["rank"] = i
            writer.writerow({k: p.get(k, "") for k in fieldnames})

    print(f"  Salvo em: {filepath}")


def main():
    print("=" * 60)
    print("DESIGN PRIMERS RPA - Pipeline CRISPR-Cas12a")
    print("=" * 60)

    os.makedirs(PRIMERS_DIR, exist_ok=True)

    all_primers = {}

    for gene_name in TARGETS:
        print(f"\n{'─' * 40}")
        print(f"[{gene_name}]")

        # Carregar melhor guide
        guide = load_best_guide(gene_name)
        if not guide:
            print(f"  ✗ Nenhum guide encontrado. Execute design_guides.py primeiro!")
            continue

        guide_pos = int(guide["position"])
        print(f"  Guide posição: {guide_pos} (spacer: {guide['spacer_seq']})")

        # Carregar sequência
        fasta_path = os.path.join(SEQUENCES_DIR, f"{gene_name}.fasta")
        if not os.path.exists(fasta_path):
            print(f"  ✗ Sequência não encontrada: {fasta_path}")
            continue

        seqs = parse_fasta(fasta_path)
        sequence = list(seqs.values())[0]

        # Desenhar primers
        primers = design_primers_for_target(gene_name, sequence, guide_pos)

        if primers:
            save_primers(gene_name, primers)
            all_primers[gene_name] = primers

            # Mostrar melhor par
            best = primers[0]
            print(f"\n  Melhor par (score: {best['pair_score']}):")
            print(f"  FWD: 5'-{best['fwd_seq']}-3' (Tm={best['fwd_tm']}°C, GC={best['fwd_gc']*100:.0f}%)")
            print(f"  REV: 5'-{best['rev_seq']}-3' (Tm={best['rev_tm']}°C, GC={best['rev_gc']*100:.0f}%)")
            print(f"  Amplicon: {best['amplicon_size']} bp")
        else:
            print(f"  ✗ Nenhum par de primers válido encontrado")

    return all_primers


if __name__ == "__main__":
    main()
