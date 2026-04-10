#!/usr/bin/env python3
"""
Script 1: Busca sequências DNA dos genes-alvo no NCBI via Entrez API.
Baixa as sequências nucleotídicas (não proteicas) necessárias para design de crRNA e primers.
"""

import os
import time
import requests

from config import TARGETS, SEQUENCES_DIR, NCBI_BASE_URL, NCBI_EMAIL


def fetch_gene_sequence(accession: str, gene_name: str) -> str | None:
    """Busca sequência nucleotídica do NCBI Entrez."""
    url = f"{NCBI_BASE_URL}/efetch.fcgi"
    params = {
        "db": "nucleotide",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "email": NCBI_EMAIL,
    }

    print(f"  Buscando {gene_name} ({accession}) no NCBI...")
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        content = response.text.strip()
        if content.startswith(">"):
            print(f"  ✓ {gene_name}: sequência obtida com sucesso")
            return content
        else:
            print(f"  ✗ {gene_name}: resposta inesperada do NCBI")
            print(f"    Primeiros 200 chars: {content[:200]}")
            return None
    except requests.RequestException as e:
        print(f"  ✗ {gene_name}: erro na requisição - {e}")
        return None


def save_fasta(content: str, filepath: str):
    """Salva conteúdo FASTA em arquivo."""
    with open(filepath, "w") as f:
        f.write(content + "\n")
    print(f"  Salvo em: {filepath}")


def extract_sequence_from_fasta(fasta_content: str) -> tuple[str, str]:
    """Extrai header e sequência de conteúdo FASTA."""
    lines = fasta_content.strip().split("\n")
    header = lines[0][1:]  # remove >
    sequence = "".join(line.strip() for line in lines[1:] if not line.startswith(">"))
    return header, sequence


def main():
    print("=" * 60)
    print("FETCH SEQUENCES - Pipeline CRISPR-Cas12a")
    print("=" * 60)

    os.makedirs(SEQUENCES_DIR, exist_ok=True)

    results = {}

    for gene_name, target in TARGETS.items():
        accession = target["gene_accession"]
        fasta_path = os.path.join(SEQUENCES_DIR, f"{gene_name}.fasta")

        # Verificar se já existe
        if os.path.exists(fasta_path):
            print(f"\n[{gene_name}] Arquivo já existe: {fasta_path}")
            with open(fasta_path) as f:
                content = f.read()
            header, seq = extract_sequence_from_fasta(content)
            results[gene_name] = {
                "accession": accession,
                "length": len(seq),
                "file": fasta_path,
            }
            print(f"  Comprimento: {len(seq)} bp")
            continue

        # Buscar do NCBI
        print(f"\n[{gene_name}] Buscando sequência...")
        content = fetch_gene_sequence(accession, gene_name)

        if content:
            save_fasta(content, fasta_path)
            header, seq = extract_sequence_from_fasta(content)
            results[gene_name] = {
                "accession": accession,
                "length": len(seq),
                "file": fasta_path,
            }
            print(f"  Comprimento: {len(seq)} bp")
        else:
            results[gene_name] = {"accession": accession, "error": "fetch failed"}

        # Rate limit NCBI (max 3 requests/sec sem API key)
        time.sleep(0.5)

    # Resumo
    print("\n" + "=" * 60)
    print("RESUMO")
    print("=" * 60)
    for gene_name, info in results.items():
        if "error" in info:
            print(f"  ✗ {gene_name}: {info['error']}")
        else:
            print(f"  ✓ {gene_name}: {info['length']} bp → {info['file']}")

    return results


if __name__ == "__main__":
    main()
