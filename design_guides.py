#!/usr/bin/env python3
"""
Script 2: Design de crRNAs para Cas12a (DETECTR).
Escaneia sequências gênicas por sites PAM TTTV e gera guides ranqueados.
"""

import os
import csv

from config import TARGETS, SEQUENCES_DIR, GUIDES_DIR, CAS12A, CAS12A_DIRECT_REPEAT
from utils import (
    parse_fasta,
    find_pam_sites,
    extract_spacer,
    score_guide,
    gc_content,
    reverse_complement,
)


def design_guides_for_gene(gene_name: str, sequence: str) -> list[dict]:
    """Desenha e ranqueia guides Cas12a para um gene."""
    spacer_len = CAS12A["spacer_length_optimal"]
    gene_len = len(sequence)

    print(f"  Sequência: {gene_len} bp")
    print(f"  Escaneando sites PAM TTTV...")

    # Encontrar todos os sites PAM
    pam_sites = find_pam_sites(sequence, "TTTV")
    print(f"  Sites PAM encontrados: {len(pam_sites)}")

    if not pam_sites:
        print(f"  ✗ Nenhum site PAM TTTV encontrado!")
        return []

    # Extrair e pontuar cada spacer candidato
    candidates = []
    for site in pam_sites:
        spacer = extract_spacer(sequence, site, spacer_len)

        # Pular spacers com bases ambíguas
        if not spacer or len(spacer) < spacer_len:
            continue
        if any(b not in "ATCG" for b in spacer.upper()):
            continue

        scores = score_guide(spacer, site["position"], gene_len)

        # Filtrar por critérios mínimos
        if scores["gc"] < CAS12A["gc_min"] or scores["gc"] > CAS12A["gc_max"]:
            continue

        candidates.append(
            {
                "gene": gene_name,
                "spacer_seq": spacer,
                "pam_seq": site["pam_seq"],
                "position": site["position"],
                "strand": site["strand"],
                "spacer_length": len(spacer),
                **scores,
                # crRNA completo = direct repeat + spacer
                "crRNA_LbCas12a": CAS12A_DIRECT_REPEAT["LbCas12a"] + spacer,
                "crRNA_AsCas12a": CAS12A_DIRECT_REPEAT["AsCas12a"] + spacer,
            }
        )

    # Ranquear por score
    candidates.sort(key=lambda x: x["score"], reverse=True)

    # Top N guides
    top_n = CAS12A["top_guides"]
    top = candidates[:top_n]

    print(f"  Candidatos válidos: {len(candidates)}")
    print(f"  Top {top_n} selecionados")

    return top


def save_guides(gene_name: str, guides: list[dict]):
    """Salva guides em TSV."""
    filepath = os.path.join(GUIDES_DIR, f"{gene_name}_cas12a_guides.tsv")
    fieldnames = [
        "rank",
        "gene",
        "spacer_seq",
        "pam_seq",
        "position",
        "strand",
        "spacer_length",
        "score",
        "gc",
        "homopolymer",
        "self_comp",
        "poly_t",
        "rel_position",
        "crRNA_LbCas12a",
        "crRNA_AsCas12a",
    ]

    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for i, guide in enumerate(guides, 1):
            guide["rank"] = i
            writer.writerow({k: guide.get(k, "") for k in fieldnames})

    print(f"  Salvo em: {filepath}")


def print_guide_summary(guides: list[dict]):
    """Imprime resumo dos guides."""
    for i, g in enumerate(guides, 1):
        print(f"\n  --- Guide #{i} (score: {g['score']}) ---")
        print(f"  Spacer:   5'-{g['spacer_seq']}-3'")
        print(f"  PAM:      {g['pam_seq']} ({g['strand']} strand)")
        print(f"  Posição:  {g['position']} ({g['rel_position']*100:.1f}% do gene)")
        print(f"  GC:       {g['gc']*100:.1f}%")
        print(f"  Homopoly: {g['homopolymer']}")
        print(f"  crRNA:    5'-{g['crRNA_LbCas12a']}-3' (LbCas12a)")


def main():
    print("=" * 60)
    print("DESIGN GUIDES - Pipeline CRISPR-Cas12a")
    print("=" * 60)

    os.makedirs(GUIDES_DIR, exist_ok=True)

    all_guides = {}

    for gene_name in TARGETS:
        fasta_path = os.path.join(SEQUENCES_DIR, f"{gene_name}.fasta")

        print(f"\n{'─' * 40}")
        print(f"[{gene_name}]")

        if not os.path.exists(fasta_path):
            print(f"  ✗ Arquivo não encontrado: {fasta_path}")
            print(f"  → Execute fetch_sequences.py primeiro!")
            continue

        seqs = parse_fasta(fasta_path)
        if not seqs:
            print(f"  ✗ Nenhuma sequência no arquivo")
            continue

        # Usar primeira sequência
        seq_id = list(seqs.keys())[0]
        sequence = seqs[seq_id]

        guides = design_guides_for_gene(gene_name, sequence)
        if guides:
            save_guides(gene_name, guides)
            print_guide_summary(guides)
            all_guides[gene_name] = guides
        else:
            print(f"  ✗ Nenhum guide válido encontrado")

    # Resumo final
    print(f"\n{'=' * 60}")
    print("RESUMO FINAL")
    print("=" * 60)
    for gene_name, guides in all_guides.items():
        best = guides[0]
        print(f"  {gene_name}: {len(guides)} guides | melhor score={best['score']} | spacer={best['spacer_seq']}")

    return all_guides


if __name__ == "__main__":
    main()
