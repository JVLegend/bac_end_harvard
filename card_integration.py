#!/usr/bin/env python3
"""
Integracao com CARD (Comprehensive Antibiotic Resistance Database).
Baixa dados do CARD, mapeia para alvos existentes, descobre novas variantes
e enriquece CSVs com metadados (mecanismo de resistencia, ontologia ARO, drug class).

Referencia: https://card.mcmaster.ca/
"""

import os
import csv
import json
import time
import zipfile
import tempfile
import requests

from config import BASE_DIR, REPORTS_DIR

# === CARD URLs ===
CARD_DATA_URL = "https://card.mcmaster.ca/latest/data"
CARD_CACHE_DIR = os.path.join(BASE_DIR, "card_data")
CARD_JSON_PATH = os.path.join(CARD_CACHE_DIR, "card.json")
CARD_ENRICHED_CSV = os.path.join(BASE_DIR, "targets_brazil_card.csv")
CARD_DISCOVERY_CSV = os.path.join(REPORTS_DIR, "card_new_variants.csv")

# Familias genicas que rastreamos no pipeline
TRACKED_FAMILIES = {
    "mecA": ["mecA", "mecA1", "mecA2"],
    "blaKPC": ["KPC"],
    "blaNDM": ["NDM"],
    "vanA": ["vanA"],
    "mcr": ["mcr-1", "mcr-5", "MCR"],
    "blaCTX-M": ["CTX-M"],
    "blaOXA-48": ["OXA-48", "OXA-181", "OXA-232"],
    "blaVIM": ["VIM"],
    "blaIMP": ["IMP"],
    "blaGES": ["GES"],
    "qnrS": ["qnrS", "QnrS"],
    "armA": ["armA", "ArmA"],
}

# Mapeamento de drug class CARD -> nomes simplificados PT-BR
DRUG_CLASS_MAP = {
    "carbapenem": "Carbapenems",
    "cephalosporin": "Cefalosporinas",
    "penam": "Beta-lactamicos",
    "methicillin": "Meticilina",
    "vancomycin": "Vancomicina",
    "colistin": "Colistina",
    "polymyxin": "Polimixinas",
    "fluoroquinolone": "Fluoroquinolonas",
    "aminoglycoside": "Aminoglicosideos",
    "cephamycin": "Cefamicinas",
    "monobactam": "Monobactamicos",
}


def download_card_data(force: bool = False) -> str:
    """Baixa o arquivo card.json do CARD. Retorna caminho do JSON."""
    os.makedirs(CARD_CACHE_DIR, exist_ok=True)

    if os.path.exists(CARD_JSON_PATH) and not force:
        age_hours = (time.time() - os.path.getmtime(CARD_JSON_PATH)) / 3600
        if age_hours < 168:  # 7 dias
            print(f"[CARD] Usando cache ({age_hours:.0f}h): {CARD_JSON_PATH}")
            return CARD_JSON_PATH
        print(f"[CARD] Cache expirado ({age_hours:.0f}h), re-baixando...")

    print(f"[CARD] Baixando dados de {CARD_DATA_URL}...")
    try:
        resp = requests.get(CARD_DATA_URL, timeout=120, stream=True)
        resp.raise_for_status()

        # CARD retorna um ZIP contendo card.json
        with tempfile.NamedTemporaryFile(suffix=".tar.bz2", delete=False) as tmp:
            tmp_path = tmp.name
            for chunk in resp.iter_content(chunk_size=8192):
                tmp.write(chunk)

        # Tentar extrair como ZIP primeiro (formato mais comum do download)
        try:
            with zipfile.ZipFile(tmp_path, "r") as zf:
                # Procurar card.json dentro do ZIP
                json_files = [n for n in zf.namelist() if n.endswith("card.json")]
                if json_files:
                    zf.extract(json_files[0], CARD_CACHE_DIR)
                    extracted = os.path.join(CARD_CACHE_DIR, json_files[0])
                    if extracted != CARD_JSON_PATH:
                        os.replace(extracted, CARD_JSON_PATH)
                    print(f"[CARD] Extraido: {CARD_JSON_PATH}")
                else:
                    print("[CARD] card.json nao encontrado no ZIP. Tentando tar.bz2...")
                    _extract_tar_bz2(tmp_path)
        except zipfile.BadZipFile:
            # Pode ser tar.bz2
            _extract_tar_bz2(tmp_path)

        os.unlink(tmp_path)

    except Exception as e:
        print(f"[CARD] Erro no download: {e}")
        if os.path.exists(CARD_JSON_PATH):
            print("[CARD] Usando cache anterior.")
        else:
            print("[CARD] Sem cache disponivel. Usando modo offline.")
            return None

    return CARD_JSON_PATH


def _extract_tar_bz2(tmp_path: str):
    """Extrai card.json de um tar.bz2."""
    import tarfile
    with tarfile.open(tmp_path, "r:bz2") as tar:
        members = [m for m in tar.getmembers() if m.name.endswith("card.json")]
        if members:
            tar.extract(members[0], CARD_CACHE_DIR)
            extracted = os.path.join(CARD_CACHE_DIR, members[0].name)
            if extracted != CARD_JSON_PATH:
                os.replace(extracted, CARD_JSON_PATH)
            print(f"[CARD] Extraido (tar.bz2): {CARD_JSON_PATH}")
        else:
            print("[CARD] card.json nao encontrado no arquivo.")


def parse_card_json(json_path: str) -> list[dict]:
    """
    Faz parse do card.json e extrai entradas relevantes para nossos alvos.
    Retorna lista de dicts com metadados enriquecidos.
    """
    print(f"[CARD] Parsing {json_path}...")
    with open(json_path, encoding="utf-8") as f:
        data = json.load(f)

    entries = []
    skipped = 0

    for key, entry in data.items():
        # Pular metadados do CARD (chaves nao-numericas)
        if not key.isdigit():
            continue

        model_name = entry.get("model_name", "")
        model_type = entry.get("model_type", "")
        aro_accession = entry.get("ARO_accession", "")
        aro_name = entry.get("ARO_name", "")
        aro_description = entry.get("ARO_description", "")

        # Extrair categorias ARO
        categories = entry.get("ARO_category", {})
        drug_classes = []
        resistance_mechanisms = []
        gene_family_names = []

        for cat_key, cat_val in categories.items():
            cat_name = cat_val.get("category_aro_name", "")
            cat_class = cat_val.get("category_aro_class_name", "")

            if cat_class == "Drug Class":
                drug_classes.append(cat_name)
            elif cat_class == "Resistance Mechanism":
                resistance_mechanisms.append(cat_name)
            elif cat_class == "AMR Gene Family":
                gene_family_names.append(cat_name)

        # Extrair sequencias de DNA (se disponiveis)
        dna_seqs = {}
        model_sequences = entry.get("model_sequences", {})
        for seq_key, seq_data in model_sequences.items():
            for acc_key, acc_data in seq_data.items():
                dna_seq = acc_data.get("dna_sequence", {})
                if dna_seq:
                    dna_seqs[dna_seq.get("accession", "")] = {
                        "sequence": dna_seq.get("sequence", ""),
                        "fmin": dna_seq.get("fmin", ""),
                        "fmax": dna_seq.get("fmax", ""),
                    }

        # Extrair NCBI accessions dos parametros do modelo
        model_params = entry.get("model_param", {})
        snp_info = {}
        for param_key, param_data in model_params.items():
            param_type = param_data.get("param_type", "")
            if param_type == "snp":
                snp_val = param_data.get("param_value", {})
                for snp_key, snp_detail in snp_val.items():
                    snp_info[snp_key] = snp_detail

        entries.append({
            "card_model_id": key,
            "model_name": model_name,
            "model_type": model_type,
            "aro_accession": aro_accession,
            "aro_name": aro_name,
            "aro_description": aro_description,
            "drug_classes": "; ".join(sorted(drug_classes)),
            "resistance_mechanisms": "; ".join(sorted(resistance_mechanisms)),
            "gene_families": "; ".join(sorted(gene_family_names)),
            "dna_accessions": "; ".join(dna_seqs.keys()),
            "has_sequence": len(dna_seqs) > 0,
            "snp_count": len(snp_info),
        })

    print(f"[CARD] Total de modelos parseados: {len(entries)} (ignorados: {skipped})")
    return entries


def match_card_to_pipeline(card_entries: list[dict]) -> dict:
    """
    Mapeia entradas CARD para as familias genicas do pipeline.
    Retorna dict {familia: [entradas CARD matchadas]}.
    """
    matched = {family: [] for family in TRACKED_FAMILIES}
    unmatched_relevant = []

    for entry in card_entries:
        name = entry["model_name"]
        aro_name = entry["aro_name"]

        found_family = None
        for family, keywords in TRACKED_FAMILIES.items():
            for kw in keywords:
                if kw.lower() in name.lower() or kw.lower() in aro_name.lower():
                    found_family = family
                    break
            if found_family:
                break

        if found_family:
            matched[found_family].append(entry)
        else:
            # Verificar se e AMR relevante (beta-lactamase, carbapenemase, etc.)
            mechanisms = entry.get("resistance_mechanisms", "").lower()
            drugs = entry.get("drug_classes", "").lower()
            if any(term in mechanisms for term in ["inactivation", "target alteration", "efflux"]):
                if any(term in drugs for term in ["carbapenem", "cephalosporin", "methicillin"]):
                    unmatched_relevant.append(entry)

    # Estatisticas
    print(f"\n[CARD] Mapeamento para familias do pipeline:")
    print(f"{'Familia':<15} {'Entradas CARD':<15} {'Exemplos'}")
    print("-" * 60)
    total_matched = 0
    for family, entries in sorted(matched.items()):
        examples = ", ".join(e["model_name"] for e in entries[:3])
        if len(entries) > 3:
            examples += f" (+{len(entries)-3})"
        print(f"{family:<15} {len(entries):<15} {examples}")
        total_matched += len(entries)

    print(f"\nTotal mapeado: {total_matched}")
    print(f"AMR relevante nao-mapeado: {len(unmatched_relevant)}")

    return matched, unmatched_relevant


def load_existing_variants(csv_path: str = None) -> set:
    """Carrega nomes de variantes ja existentes no pipeline."""
    path = csv_path or os.path.join(BASE_DIR, "targets_brazil_variants.csv")
    existing = set()
    if os.path.exists(path):
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                existing.add(row["name"].strip())
    return existing


def simplify_drug_class(drug_str: str) -> str:
    """Converte drug class CARD para nome simplificado PT-BR."""
    drug_lower = drug_str.lower()
    for eng, ptbr in DRUG_CLASS_MAP.items():
        if eng in drug_lower:
            return ptbr
    return drug_str.split(";")[0].strip() if drug_str else "Outros"


def generate_enriched_csv(matched: dict, existing_variants: set):
    """
    Gera CSV enriquecido combinando dados existentes com metadados CARD.
    """
    os.makedirs(REPORTS_DIR, exist_ok=True)

    # Carregar dados existentes dos dois CSVs
    existing_data = {}
    for csv_name in ["targets_brazil.csv", "targets_brazil_variants.csv"]:
        csv_path = os.path.join(BASE_DIR, csv_name)
        if os.path.exists(csv_path):
            with open(csv_path) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    existing_data[row["name"]] = dict(row)

    # Enriquecer com dados CARD
    enriched = []
    for family, card_entries in matched.items():
        # Dados CARD agregados para a familia
        all_drugs = set()
        all_mechanisms = set()
        all_gene_fams = set()
        aro_ids = set()

        for entry in card_entries:
            if entry["drug_classes"]:
                all_drugs.update(d.strip() for d in entry["drug_classes"].split(";"))
            if entry["resistance_mechanisms"]:
                all_mechanisms.update(m.strip() for m in entry["resistance_mechanisms"].split(";"))
            if entry["gene_families"]:
                all_gene_fams.update(g.strip() for g in entry["gene_families"].split(";"))
            if entry["aro_accession"]:
                aro_ids.add(entry["aro_accession"])

        card_meta = {
            "card_drug_classes": "; ".join(sorted(all_drugs)),
            "card_resistance_mechanisms": "; ".join(sorted(all_mechanisms)),
            "card_gene_families": "; ".join(sorted(all_gene_fams)),
            "card_aro_ids": "; ".join(sorted(aro_ids)),
            "card_variant_count": len(card_entries),
        }

        # Enriquecer variantes existentes desta familia
        for name, data in existing_data.items():
            gene_fam = data.get("gene_family", data.get("name", ""))
            is_in_family = False

            # Verificar se pertence a esta familia
            for kw in TRACKED_FAMILIES.get(family, []):
                if kw.lower() in name.lower() or kw.lower() in gene_fam.lower():
                    is_in_family = True
                    break
            # Match direto pelo nome da familia
            if family.lower() in name.lower() or family.lower() in gene_fam.lower():
                is_in_family = True

            if is_in_family:
                row = {**data, **card_meta, "card_family": family}
                enriched.append(row)

    # Escrever CSV enriquecido
    if enriched:
        fieldnames = list(enriched[0].keys())
        with open(CARD_ENRICHED_CSV, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in enriched:
                writer.writerow(row)
        print(f"\n[CARD] CSV enriquecido salvo: {CARD_ENRICHED_CSV} ({len(enriched)} entradas)")
    else:
        print("[CARD] Nenhuma entrada enriquecida gerada.")

    return enriched


def discover_new_variants(matched: dict, existing_variants: set):
    """
    Identifica variantes no CARD que nao estao no pipeline.
    Gera CSV de descobertas para revisao.
    """
    os.makedirs(REPORTS_DIR, exist_ok=True)
    discoveries = []

    for family, card_entries in matched.items():
        for entry in card_entries:
            name = entry["model_name"]
            # Verificar se ja existe (por nome exato ou parcial)
            already_tracked = any(
                name.lower() in existing.lower() or existing.lower() in name.lower()
                for existing in existing_variants
            )
            if not already_tracked:
                discoveries.append({
                    "card_model_name": name,
                    "pipeline_family": family,
                    "aro_accession": entry["aro_accession"],
                    "aro_name": entry["aro_name"],
                    "drug_classes": entry["drug_classes"],
                    "resistance_mechanisms": entry["resistance_mechanisms"],
                    "has_sequence": entry["has_sequence"],
                    "description": entry["aro_description"][:200],
                })

    if discoveries:
        fieldnames = list(discoveries[0].keys())
        with open(CARD_DISCOVERY_CSV, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in discoveries:
                writer.writerow(row)
        print(f"\n[CARD] Novas variantes descobertas: {len(discoveries)}")
        print(f"[CARD] Salvo em: {CARD_DISCOVERY_CSV}")

        # Resumo por familia
        print(f"\n{'Familia':<15} {'Novas variantes':<18} {'Exemplos'}")
        print("-" * 65)
        by_family = {}
        for d in discoveries:
            fam = d["pipeline_family"]
            if fam not in by_family:
                by_family[fam] = []
            by_family[fam].append(d["card_model_name"])

        for fam, names in sorted(by_family.items()):
            examples = ", ".join(names[:3])
            if len(names) > 3:
                examples += f" (+{len(names)-3})"
            print(f"{fam:<15} {len(names):<18} {examples}")
    else:
        print("\n[CARD] Nenhuma variante nova descoberta (pipeline ja cobre tudo).")

    return discoveries


def generate_card_report(matched: dict, discoveries: list, enriched: list):
    """Gera relatorio textual da integracao CARD."""
    report_path = os.path.join(REPORTS_DIR, "card_integration_report.txt")

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("RELATORIO DE INTEGRACAO CARD\n")
        f.write(f"SmartLab BacEnd - {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write("=" * 70 + "\n\n")

        # Resumo geral
        total_card = sum(len(entries) for entries in matched.values())
        f.write(f"Total de entradas CARD mapeadas: {total_card}\n")
        f.write(f"Variantes enriquecidas: {len(enriched)}\n")
        f.write(f"Novas variantes descobertas: {len(discoveries)}\n\n")

        # Detalhes por familia
        f.write("-" * 70 + "\n")
        f.write("DETALHES POR FAMILIA GENICA\n")
        f.write("-" * 70 + "\n\n")

        for family, entries in sorted(matched.items()):
            f.write(f"\n[{family}] - {len(entries)} entradas no CARD\n")

            if entries:
                # Drug classes
                drugs = set()
                mechs = set()
                for e in entries:
                    if e["drug_classes"]:
                        drugs.update(d.strip() for d in e["drug_classes"].split(";"))
                    if e["resistance_mechanisms"]:
                        mechs.update(m.strip() for m in e["resistance_mechanisms"].split(";"))

                f.write(f"  Drug classes: {', '.join(sorted(drugs))}\n")
                f.write(f"  Mecanismos: {', '.join(sorted(mechs))}\n")
                f.write(f"  Variantes CARD:\n")
                for e in entries[:10]:
                    f.write(f"    - {e['model_name']} (ARO:{e['aro_accession']})\n")
                if len(entries) > 10:
                    f.write(f"    ... e mais {len(entries)-10}\n")

        # Descobertas
        if discoveries:
            f.write(f"\n{'='*70}\n")
            f.write("VARIANTES NOVAS (nao rastreadas no pipeline)\n")
            f.write(f"{'='*70}\n\n")
            for d in discoveries:
                f.write(f"  [{d['pipeline_family']}] {d['card_model_name']}\n")
                f.write(f"    ARO: {d['aro_accession']} | Drugs: {d['drug_classes']}\n")
                f.write(f"    Mecanismo: {d['resistance_mechanisms']}\n")
                f.write(f"    Sequencia disponivel: {'Sim' if d['has_sequence'] else 'Nao'}\n\n")

    print(f"\n[CARD] Relatorio salvo: {report_path}")
    return report_path


def run_card_integration(force_download: bool = False):
    """Executa pipeline completo de integracao CARD."""
    print("=" * 70)
    print("INTEGRACAO CARD - Comprehensive Antibiotic Resistance Database")
    print("SmartLab BacEnd | IA para Medicos")
    print("=" * 70)

    # 1. Download
    json_path = download_card_data(force=force_download)
    if not json_path:
        print("[CARD] Abortando: sem dados disponiveis.")
        return None

    # 2. Parse
    card_entries = parse_card_json(json_path)
    if not card_entries:
        print("[CARD] Nenhuma entrada encontrada no JSON.")
        return None

    # 3. Mapear para pipeline
    matched, unmatched = match_card_to_pipeline(card_entries)

    # 4. Carregar variantes existentes
    existing = load_existing_variants()
    print(f"\n[CARD] Variantes ja rastreadas no pipeline: {len(existing)}")

    # 5. Enriquecer CSV
    enriched = generate_enriched_csv(matched, existing)

    # 6. Descobrir novas variantes
    discoveries = discover_new_variants(matched, existing)

    # 7. Gerar relatorio
    report_path = generate_card_report(matched, discoveries, enriched)

    # Resumo final
    print(f"\n{'='*70}")
    print("INTEGRACAO CARD CONCLUIDA")
    print(f"{'='*70}")
    print(f"  Entradas CARD analisadas:  {len(card_entries)}")
    print(f"  Mapeadas para pipeline:    {sum(len(v) for v in matched.values())}")
    print(f"  Variantes enriquecidas:    {len(enriched)}")
    print(f"  Novas descobertas:         {len(discoveries)}")
    print(f"  Relatorio:                 {report_path}")
    print(f"  CSV enriquecido:           {CARD_ENRICHED_CSV}")
    if discoveries:
        print(f"  CSV descobertas:           {CARD_DISCOVERY_CSV}")

    return {
        "matched": matched,
        "enriched": enriched,
        "discoveries": discoveries,
        "report": report_path,
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Integracao CARD para SmartLab BacEnd")
    parser.add_argument("--force", action="store_true", help="Forcar re-download dos dados CARD")
    args = parser.parse_args()
    run_card_integration(force_download=args.force)
