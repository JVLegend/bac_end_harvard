#!/usr/bin/env python3
"""
Script 4: Verificação de especificidade via NCBI BLAST API.
Valida que guides e primers são específicos para os organismos-alvo.
"""

import os
import csv
import time
import requests

from config import TARGETS, GUIDES_DIR, PRIMERS_DIR, REPORTS_DIR, NCBI_EMAIL


BLAST_API = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"


def submit_blast(sequence: str, program: str = "blastn", database: str = "nt") -> str | None:
    """Submete BLAST ao NCBI e retorna o RID (Request ID)."""
    params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": sequence,
        "FORMAT_TYPE": "JSON2",
        "WORD_SIZE": "7",  # menor para sequências curtas
        "EXPECT": "10",
        "EMAIL": NCBI_EMAIL,
    }

    try:
        response = requests.post(BLAST_API, data=params, timeout=30)
        response.raise_for_status()
        text = response.text

        # Extrair RID
        for line in text.split("\n"):
            if line.strip().startswith("RID ="):
                return line.split("=")[1].strip()
        return None
    except requests.RequestException as e:
        print(f"    Erro ao submeter BLAST: {e}")
        return None


def check_blast_status(rid: str) -> str:
    """Verifica status do BLAST job. Retorna 'WAITING', 'READY', ou 'FAILED'."""
    params = {"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}
    try:
        response = requests.get(BLAST_API, params=params, timeout=30)
        text = response.text
        if "Status=WAITING" in text:
            return "WAITING"
        elif "Status=FAILED" in text:
            return "FAILED"
        elif "Status=READY" in text:
            return "READY"
        return "UNKNOWN"
    except requests.RequestException:
        return "UNKNOWN"


def get_blast_results(rid: str) -> dict | None:
    """Recupera resultados do BLAST."""
    params = {
        "CMD": "Get",
        "FORMAT_TYPE": "JSON2",
        "RID": rid,
    }
    try:
        response = requests.get(BLAST_API, params=params, timeout=60)
        response.raise_for_status()
        return response.json()
    except (requests.RequestException, ValueError) as e:
        print(f"    Erro ao recuperar resultados: {e}")
        return None


def run_blast_check(sequence: str, seq_name: str, max_wait: int = 120) -> list[dict]:
    """Executa BLAST completo e retorna hits."""
    print(f"    Submetendo {seq_name} ao BLAST...")
    rid = submit_blast(sequence)
    if not rid:
        print(f"    ✗ Falha ao submeter BLAST")
        return []

    print(f"    RID: {rid} - aguardando resultados...")

    # Poll até resultado pronto
    elapsed = 0
    while elapsed < max_wait:
        time.sleep(15)
        elapsed += 15
        status = check_blast_status(rid)
        if status == "READY":
            break
        elif status == "FAILED":
            print(f"    ✗ BLAST falhou")
            return []
        print(f"    ... aguardando ({elapsed}s)")

    if elapsed >= max_wait:
        print(f"    ✗ Timeout ({max_wait}s) - tente novamente depois")
        return []

    results = get_blast_results(rid)
    if not results:
        return []

    # Parse hits
    hits = []
    try:
        search = results["BlastOutput2"][0]["report"]["results"]["search"]
        for hit in search.get("hits", [])[:10]:  # top 10 hits
            desc = hit["description"][0]
            hsps = hit["hsps"][0]
            hits.append(
                {
                    "accession": desc.get("accession", ""),
                    "title": desc.get("title", "")[:80],
                    "identity": hsps.get("identity", 0),
                    "align_len": hsps.get("align_len", 0),
                    "mismatches": hsps.get("align_len", 0) - hsps.get("identity", 0),
                    "evalue": hsps.get("evalue", 999),
                    "bit_score": hsps.get("bit_score", 0),
                }
            )
    except (KeyError, IndexError) as e:
        print(f"    Erro ao parsear resultados: {e}")

    return hits


def local_cross_reactivity_check(guides: dict[str, str]) -> list[dict]:
    """
    Verificação local de cross-reactivity entre os guides do painel.
    Compara cada guide contra os outros para garantir que não há cross-match.
    """
    results = []
    gene_names = list(guides.keys())

    for i, gene1 in enumerate(gene_names):
        for gene2 in gene_names[i + 1 :]:
            seq1 = guides[gene1].upper()
            seq2 = guides[gene2].upper()

            # Contar mismatches
            min_len = min(len(seq1), len(seq2))
            matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
            identity = matches / min_len if min_len > 0 else 0

            results.append(
                {
                    "gene1": gene1,
                    "gene2": gene2,
                    "identity": round(identity, 3),
                    "matches": matches,
                    "length": min_len,
                    "cross_reactive": identity > 0.75,  # >75% = risco
                }
            )

    return results


def load_guides() -> dict[str, str]:
    """Carrega melhor spacer de cada gene."""
    guides = {}
    for gene_name in TARGETS:
        filepath = os.path.join(GUIDES_DIR, f"{gene_name}_cas12a_guides.tsv")
        if not os.path.exists(filepath):
            continue
        with open(filepath) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                guides[gene_name] = row["spacer_seq"]
                break
    return guides


def main():
    print("=" * 60)
    print("SPECIFICITY CHECK - Pipeline CRISPR-Cas12a")
    print("=" * 60)

    os.makedirs(REPORTS_DIR, exist_ok=True)

    guides = load_guides()
    if not guides:
        print("✗ Nenhum guide encontrado. Execute design_guides.py primeiro!")
        return

    # 1. Cross-reactivity local
    print("\n--- Cross-reactivity entre guides do painel ---")
    cross_results = local_cross_reactivity_check(guides)
    for cr in cross_results:
        status = "⚠ RISCO" if cr["cross_reactive"] else "✓ OK"
        print(
            f"  {cr['gene1']} vs {cr['gene2']}: "
            f"{cr['identity']*100:.1f}% identidade ({cr['matches']}/{cr['length']} matches) {status}"
        )

    # 2. BLAST check (pode ser lento - internet necessária)
    print("\n--- BLAST check (NCBI) ---")
    print("  NOTA: Requer internet. Pode levar 1-2 min por sequência.")

    blast_results = {}
    for gene_name, spacer in guides.items():
        print(f"\n  [{gene_name}] Spacer: {spacer}")
        hits = run_blast_check(spacer, f"{gene_name}_spacer")
        blast_results[gene_name] = hits

        if hits:
            print(f"    Top hits:")
            for h in hits[:3]:
                print(f"      {h['accession']}: {h['title']} (identity={h['identity']}, e={h['evalue']:.2e})")
        else:
            print(f"    Nenhum hit ou timeout")

        time.sleep(1)  # rate limit

    # 3. Salvar relatório
    report_path = os.path.join(REPORTS_DIR, "specificity_report.tsv")
    with open(report_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["type", "gene", "detail", "result"])

        for cr in cross_results:
            writer.writerow([
                "cross_reactivity",
                f"{cr['gene1']}_vs_{cr['gene2']}",
                f"identity={cr['identity']*100:.1f}%",
                "RISK" if cr["cross_reactive"] else "OK",
            ])

        for gene, hits in blast_results.items():
            for h in hits:
                writer.writerow([
                    "blast",
                    gene,
                    f"{h['accession']}|{h['title']}|identity={h['identity']}",
                    f"e={h['evalue']:.2e}",
                ])

    print(f"\n  Relatório salvo em: {report_path}")

    return {"cross_reactivity": cross_results, "blast": blast_results}


if __name__ == "__main__":
    main()
