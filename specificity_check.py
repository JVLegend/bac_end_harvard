#!/usr/bin/env python3
"""
Script 4: Verificação de especificidade via NCBI BLAST API.
Valida que guides e primers são específicos para os organismos-alvo.
Usa amplicons completos (100-200bp) em vez de spacers curtos (20nt) para resultados confiáveis.
"""

import os
import csv
import time
import requests

from config import TARGETS, SEQUENCES_DIR, GUIDES_DIR, PRIMERS_DIR, REPORTS_DIR, NCBI_EMAIL
from utils import parse_fasta


BLAST_API = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"


def extract_amplicon(gene_name: str) -> dict | None:
    """Extrai o amplicon completo (região entre primers) da sequência do gene."""
    # Carregar sequência
    fasta_path = os.path.join(SEQUENCES_DIR, f"{gene_name}.fasta")
    if not os.path.exists(fasta_path):
        return None
    seqs = parse_fasta(fasta_path)
    sequence = list(seqs.values())[0].upper()

    # Carregar melhor primer pair
    primer_path = os.path.join(PRIMERS_DIR, f"{gene_name}_rpa_primers.tsv")
    if not os.path.exists(primer_path):
        return None
    with open(primer_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        primer = next(reader, None)
    if not primer:
        return None

    # Carregar melhor guide
    guide_path = os.path.join(GUIDES_DIR, f"{gene_name}_cas12a_guides.tsv")
    if not os.path.exists(guide_path):
        return None
    with open(guide_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        guide = next(reader, None)

    fwd_start = int(primer["fwd_start"])
    amplicon_size = int(primer["amplicon_size"])
    amplicon_seq = sequence[fwd_start : fwd_start + amplicon_size]

    return {
        "gene": gene_name,
        "amplicon_seq": amplicon_seq,
        "amplicon_size": len(amplicon_seq),
        "fwd_seq": primer["fwd_seq"],
        "rev_seq": primer["rev_seq"],
        "fwd_start": fwd_start,
        "spacer_seq": guide["spacer_seq"] if guide else "",
    }


def submit_blast(sequence: str, program: str = "blastn", database: str = "nt") -> str | None:
    """Submete BLAST ao NCBI e retorna o RID (Request ID)."""
    params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": sequence,
        "FORMAT_TYPE": "JSON2",
        "WORD_SIZE": "11",
        "EXPECT": "0.05",
        "EMAIL": NCBI_EMAIL,
    }

    try:
        response = requests.post(BLAST_API, data=params, timeout=30)
        response.raise_for_status()
        text = response.text

        for line in text.split("\n"):
            if line.strip().startswith("RID ="):
                return line.split("=")[1].strip()
        return None
    except requests.RequestException as e:
        print(f"    Erro ao submeter BLAST: {e}")
        return None


def check_blast_status(rid: str) -> str:
    """Verifica status do BLAST job."""
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
    """Recupera resultados do BLAST em JSON."""
    params = {
        "CMD": "Get",
        "FORMAT_TYPE": "JSON2",
        "RID": rid,
    }
    try:
        response = requests.get(BLAST_API, params=params, timeout=120)
        response.raise_for_status()
        text = response.text.strip()
        if text.startswith("<!DOCTYPE") or text.startswith("<html") or text.startswith("<p"):
            print(f"    NCBI retornou HTML — tentando formato texto...")
            return get_blast_results_text(rid)
        return response.json()
    except (requests.RequestException, ValueError) as e:
        print(f"    Erro JSON, tentando formato texto...")
        return get_blast_results_text(rid)


def get_blast_results_text(rid: str) -> dict | None:
    """Fallback: recupera resultados em formato Text e parseia a tabela de hits."""
    params = {
        "CMD": "Get",
        "FORMAT_TYPE": "Text",
        "RID": rid,
    }
    try:
        response = requests.get(BLAST_API, params=params, timeout=120)
        response.raise_for_status()
        text = response.text

        # Parsear a seção "Sequences producing significant alignments"
        hits = []
        in_table = False
        for line in text.split("\n"):
            stripped = line.strip()

            if "Sequences producing significant alignments:" in stripped:
                in_table = True
                continue

            if in_table and stripped == "":
                # Linha vazia pode ser separador entre header e dados, ou fim da tabela
                if hits:
                    break
                continue

            if in_table and stripped.startswith(">"):
                break  # Início dos detalhes de alinhamento

            if in_table and stripped:
                # Formato: accession description...  score  evalue  ident
                # O score, evalue e ident estão nas últimas 3 colunas
                parts = stripped.rsplit(None, 3)
                if len(parts) >= 4:
                    desc_part = parts[0]
                    try:
                        bit_score = float(parts[1])
                        evalue_str = parts[2]
                        ident_str = parts[3].replace("%", "")

                        # Extrair accession (primeira palavra)
                        accession = desc_part.split()[0] if desc_part else ""
                        title = desc_part[:100]

                        evalue = float(evalue_str)
                        identity_pct = float(ident_str)

                        hits.append({
                            "accession": accession,
                            "title": title,
                            "identity_pct": identity_pct,
                            "align_len": 0,
                            "mismatches": 0,
                            "evalue": evalue,
                            "bit_score": bit_score,
                        })
                    except (ValueError, IndexError):
                        continue

        if hits:
            return {"_tabular_hits": hits}
        return None
    except (requests.RequestException, ValueError) as e:
        print(f"    Erro no fallback texto: {e}")
        return None


def run_blast_check(sequence: str, seq_name: str, max_wait: int = 300) -> list[dict]:
    """Executa BLAST completo e retorna hits."""
    print(f"    Submetendo {seq_name} ({len(sequence)}bp) ao BLAST...")
    rid = submit_blast(sequence)
    if not rid:
        print(f"    ✗ Falha ao submeter BLAST")
        return []

    print(f"    RID: {rid} — aguardando resultados...")

    elapsed = 0
    while elapsed < max_wait:
        time.sleep(20)
        elapsed += 20
        status = check_blast_status(rid)
        if status == "READY":
            print(f"    ✓ Resultado pronto ({elapsed}s)")
            break
        elif status == "FAILED":
            print(f"    ✗ BLAST falhou")
            return []
        print(f"    ... aguardando ({elapsed}s)")

    if elapsed >= max_wait:
        print(f"    ✗ Timeout ({max_wait}s) — tente novamente depois")
        return []

    results = get_blast_results(rid)
    if not results:
        return []

    # Parse hits — formato JSON2
    hits = []
    if "_tabular_hits" in results:
        # Veio do fallback tabular
        for h in results["_tabular_hits"][:10]:
            hits.append({
                "accession": h["accession"],
                "title": h["accession"],  # tabular não tem título
                "identity": h.get("align_len", 0),
                "identity_pct": h.get("identity_pct", 0),
                "align_len": h.get("align_len", 0),
                "mismatches": h.get("mismatches", 0),
                "evalue": h.get("evalue", 999),
                "bit_score": h.get("bit_score", 0),
            })
    else:
        try:
            search = results["BlastOutput2"][0]["report"]["results"]["search"]
            for hit in search.get("hits", [])[:10]:
                desc = hit["description"][0]
                hsps = hit["hsps"][0]
                align_len = hsps.get("align_len", 0)
                identity = hsps.get("identity", 0)
                identity_pct = (identity / align_len * 100) if align_len > 0 else 0
                hits.append({
                    "accession": desc.get("accession", ""),
                    "title": desc.get("title", "")[:100],
                    "identity": identity,
                    "identity_pct": round(identity_pct, 1),
                    "align_len": align_len,
                    "mismatches": align_len - identity,
                    "evalue": hsps.get("evalue", 999),
                    "bit_score": hsps.get("bit_score", 0),
                })
        except (KeyError, IndexError) as e:
            print(f"    Erro ao parsear JSON: {e}")

    return hits


def local_cross_reactivity_check(guides: dict[str, str]) -> list[dict]:
    """Verificação local de cross-reactivity entre guides do painel."""
    results = []
    gene_names = list(guides.keys())

    for i, gene1 in enumerate(gene_names):
        for gene2 in gene_names[i + 1:]:
            seq1 = guides[gene1].upper()
            seq2 = guides[gene2].upper()

            min_len = min(len(seq1), len(seq2))
            matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
            identity = matches / min_len if min_len > 0 else 0

            results.append({
                "gene1": gene1,
                "gene2": gene2,
                "identity": round(identity, 3),
                "matches": matches,
                "length": min_len,
                "cross_reactive": identity > 0.75,
            })

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
    print("Usando amplicons completos para BLAST confiável")
    print("=" * 60)

    os.makedirs(REPORTS_DIR, exist_ok=True)

    guides = load_guides()
    if not guides:
        print("✗ Nenhum guide encontrado. Execute design_guides.py primeiro!")
        return

    # 1. Cross-reactivity local (spacers)
    print("\n--- Cross-reactivity entre guides do painel ---")
    cross_results = local_cross_reactivity_check(guides)
    for cr in cross_results:
        status = "⚠ RISCO" if cr["cross_reactive"] else "✓ OK"
        print(
            f"  {cr['gene1']} vs {cr['gene2']}: "
            f"{cr['identity']*100:.1f}% identidade ({cr['matches']}/{cr['length']} matches) {status}"
        )

    # 2. Extrair amplicons
    print("\n--- Extraindo amplicons para BLAST ---")
    amplicons = {}
    for gene_name in TARGETS:
        amp = extract_amplicon(gene_name)
        if amp:
            amplicons[gene_name] = amp
            print(f"  ✓ {gene_name}: amplicon {amp['amplicon_size']}bp extraído")
        else:
            print(f"  ✗ {gene_name}: não foi possível extrair amplicon")

    # 3. BLAST com amplicons completos
    print("\n--- BLAST check com amplicons (NCBI) ---")
    print("  NOTA: Amplicons ~150bp dão resultados muito melhores que spacers 20nt.")
    print("  Pode levar 2-5 min por sequência.\n")

    blast_results = {}
    for gene_name, amp in amplicons.items():
        print(f"  [{gene_name}] Amplicon: {amp['amplicon_size']}bp")
        print(f"    Spacer dentro: {amp['spacer_seq']}")
        hits = run_blast_check(amp["amplicon_seq"], f"{gene_name}_amplicon")
        blast_results[gene_name] = hits

        if hits:
            print(f"\n    Top hits:")
            for i, h in enumerate(hits[:5], 1):
                pct = h.get("identity_pct", 0)
                print(
                    f"      {i}. {h['accession']}: {h.get('title', 'N/A')}"
                    f"\n         identity={pct:.1f}% | align={h['align_len']}bp | "
                    f"mismatches={h['mismatches']} | e={h['evalue']:.2e} | score={h['bit_score']}"
                )

            # Análise de especificidade
            top_hit = hits[0]
            top_pct = top_hit.get("identity_pct", 0)
            if top_pct >= 95:
                print(f"\n    ✓ ESPECÍFICO: top hit {top_pct:.1f}% identity — provável match no gene-alvo")
            elif top_pct >= 80:
                print(f"\n    ⚠ ATENÇÃO: top hit {top_pct:.1f}% identity — verificar se é o organismo correto")
            else:
                print(f"\n    ✗ BAIXA IDENTIDADE: top hit {top_pct:.1f}% — pode indicar problema no design")
        else:
            print(f"    Nenhum hit encontrado ou timeout")

        time.sleep(2)  # rate limit entre genes

    # 4. Salvar relatório
    report_path = os.path.join(REPORTS_DIR, "specificity_report.tsv")
    with open(report_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "type", "gene", "query_type", "query_length",
            "hit_accession", "hit_title", "identity_pct",
            "align_len", "mismatches", "evalue", "bit_score", "assessment"
        ])

        for cr in cross_results:
            writer.writerow([
                "cross_reactivity",
                f"{cr['gene1']}_vs_{cr['gene2']}",
                "spacer", "20",
                "", "",
                f"{cr['identity']*100:.1f}",
                cr['length'], cr['length'] - cr['matches'],
                "", "",
                "RISK" if cr["cross_reactive"] else "OK",
            ])

        for gene, hits in blast_results.items():
            amp = amplicons.get(gene, {})
            for h in hits:
                pct = h.get("identity_pct", 0)
                assessment = "SPECIFIC" if pct >= 95 else ("CHECK" if pct >= 80 else "LOW")
                writer.writerow([
                    "blast_amplicon",
                    gene,
                    "amplicon",
                    amp.get("amplicon_size", ""),
                    h["accession"],
                    h.get("title", ""),
                    f"{pct:.1f}",
                    h["align_len"],
                    h["mismatches"],
                    f"{h['evalue']:.2e}",
                    h["bit_score"],
                    assessment,
                ])

    print(f"\n{'=' * 60}")
    print(f"Relatório salvo em: {report_path}")
    print("=" * 60)

    return {"cross_reactivity": cross_results, "blast": blast_results}


if __name__ == "__main__":
    main()
