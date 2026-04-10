#!/usr/bin/env python3
"""
Script 5: Montagem do painel multiplex final.
Gera layout do dispositivo, folha de encomenda de oligos e relatório.
"""

import os
import csv
from datetime import datetime

from config import TARGETS, GUIDES_DIR, PRIMERS_DIR, REPORTS_DIR, CAS12A_DIRECT_REPEAT, REPORTER, CONTROLS


def load_best_guide(gene_name: str) -> dict | None:
    filepath = os.path.join(GUIDES_DIR, f"{gene_name}_cas12a_guides.tsv")
    if not os.path.exists(filepath):
        return None
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            return row
    return None


def load_best_primers(gene_name: str) -> dict | None:
    filepath = os.path.join(PRIMERS_DIR, f"{gene_name}_rpa_primers.tsv")
    if not os.path.exists(filepath):
        return None
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            return row
    return None


def generate_oligo_order(guides: dict, primers: dict) -> list[dict]:
    """Gera folha de encomenda de oligos."""
    oligos = []

    for gene_name in TARGETS:
        guide = guides.get(gene_name)
        primer = primers.get(gene_name)

        if guide:
            # crRNA (com direct repeat)
            spacer = guide["spacer_seq"]
            crRNA = CAS12A_DIRECT_REPEAT["LbCas12a"] + spacer
            oligos.append(
                {
                    "name": f"crRNA_{gene_name}_LbCas12a",
                    "sequence": crRNA,
                    "length": len(crRNA),
                    "type": "crRNA",
                    "modification_5": "none",
                    "modification_3": "none",
                    "scale": "10 nmol",
                    "purification": "PAGE",
                    "notes": f"Guide RNA for Cas12a targeting {gene_name}",
                }
            )

        if primer:
            # Forward primer
            oligos.append(
                {
                    "name": f"RPA_FWD_{gene_name}",
                    "sequence": primer["fwd_seq"],
                    "length": len(primer["fwd_seq"]),
                    "type": "RPA primer",
                    "modification_5": "none",
                    "modification_3": "none",
                    "scale": "25 nmol",
                    "purification": "Standard desalting",
                    "notes": f"RPA forward primer for {gene_name}",
                }
            )
            # Reverse primer
            oligos.append(
                {
                    "name": f"RPA_REV_{gene_name}",
                    "sequence": primer["rev_seq"],
                    "length": len(primer["rev_seq"]),
                    "type": "RPA primer",
                    "modification_5": "none",
                    "modification_3": "none",
                    "scale": "25 nmol",
                    "purification": "Standard desalting",
                    "notes": f"RPA reverse primer for {gene_name}",
                }
            )

    # Reporter fluorescente
    oligos.append(
        {
            "name": "Reporter_ssDNA_FQ",
            "sequence": "TTATTATT",
            "length": 8,
            "type": "Reporter",
            "modification_5": "6-FAM",
            "modification_3": "BHQ-1",
            "scale": "100 nmol",
            "purification": "HPLC",
            "notes": "Trans-cleavage reporter for Cas12a detection",
        }
    )

    return oligos


def estimate_cost(oligos: list[dict]) -> float:
    """Estima custo aproximado dos oligos."""
    total = 0.0
    for oligo in oligos:
        base_cost = len(oligo["sequence"]) * 0.20  # ~$0.20/base
        if oligo["type"] == "crRNA":
            base_cost *= 2.5  # RNA synthesis premium
        if "FAM" in oligo.get("modification_5", ""):
            base_cost += 40  # fluorophore
        if "BHQ" in oligo.get("modification_3", ""):
            base_cost += 30  # quencher
        if oligo["purification"] == "PAGE":
            base_cost += 20
        elif oligo["purification"] == "HPLC":
            base_cost += 30
        total += base_cost
    return round(total, 2)


def generate_report(guides: dict, primers: dict, oligos: list[dict], cost: float):
    """Gera relatório final em texto."""
    report_path = os.path.join(REPORTS_DIR, "panel_report.txt")

    with open(report_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("CRISPR-Cas12a PAPER-BASED DIAGNOSTIC PANEL\n")
        f.write(f"Hackathon - Hospital Fecal Residue Detection\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write("=" * 70 + "\n\n")

        # Layout do dispositivo
        f.write("DEVICE LAYOUT (4 spots)\n")
        f.write("─" * 30 + "\n")
        f.write("+--------+--------+\n")
        f.write("|   P    |  mecA  |\n")
        f.write("|  (+)   |  MRSA  |\n")
        f.write("+--------+--------+\n")
        f.write("|   N    |  KPC   |\n")
        f.write("|  (-)   |  CRE   |\n")
        f.write("+--------+--------+\n\n")

        # Detalhes por alvo
        for gene_name, target in TARGETS.items():
            f.write(f"\n{'─' * 50}\n")
            f.write(f"TARGET: {gene_name} ({target['pathogen']})\n")
            f.write(f"Organism: {target['organism']}\n")
            f.write(f"Gene accession: {target['gene_accession']}\n")
            f.write(f"Clinical: {target['clinical_relevance']}\n")

            guide = guides.get(gene_name)
            if guide:
                f.write(f"\n  crRNA Guide:\n")
                f.write(f"    Spacer:  5'-{guide['spacer_seq']}-3'\n")
                f.write(f"    PAM:     {guide['pam_seq']} ({guide['strand']} strand)\n")
                f.write(f"    Score:   {guide['score']}\n")
                f.write(f"    GC:      {float(guide['gc'])*100:.1f}%\n")
                f.write(f"    crRNA:   5'-{guide['crRNA_LbCas12a']}-3'\n")

            primer = primers.get(gene_name)
            if primer:
                f.write(f"\n  RPA Primers:\n")
                f.write(f"    FWD: 5'-{primer['fwd_seq']}-3'\n")
                f.write(f"         Tm={primer['fwd_tm']}°C, GC={float(primer['fwd_gc'])*100:.0f}%\n")
                f.write(f"    REV: 5'-{primer['rev_seq']}-3'\n")
                f.write(f"         Tm={primer['rev_tm']}°C, GC={float(primer['rev_gc'])*100:.0f}%\n")
                f.write(f"    Amplicon: {primer['amplicon_size']} bp\n")

        # Controles
        f.write(f"\n{'─' * 50}\n")
        f.write("CONTROLS\n")
        for ctrl_type, ctrl in CONTROLS.items():
            f.write(f"  [{ctrl['spot']}] {ctrl['name']}: {ctrl['description']}\n")

        # Reporter
        f.write(f"\n{'─' * 50}\n")
        f.write("REPORTER\n")
        f.write(f"  Type: {REPORTER['type']}\n")
        f.write(f"  {REPORTER['fluorophore']}/{REPORTER['quencher']}\n")
        f.write(f"  Sequence: {REPORTER['sequence']}\n")

        # Protocolo
        f.write(f"\n{'─' * 50}\n")
        f.write("PROTOCOL SUMMARY\n")
        f.write("  1. Collect sample from hospital fecal residue\n")
        f.write("  2. Apply to glass fiber pad (chemical lysis)\n")
        f.write("  3. Fold paper device (lateral flow extraction)\n")
        f.write("  4. Add elution buffer to reaction spots\n")
        f.write("  5. Incubate at 37°C for 20 min (RPA amplification)\n")
        f.write("  6. Cas12a trans-cleavage activates reporter\n")
        f.write("  7. Read fluorescence under UV or smartphone\n")
        f.write("  Time to result: ~30 minutes\n")

        # Oligos
        f.write(f"\n{'─' * 50}\n")
        f.write(f"OLIGO ORDER SHEET ({len(oligos)} oligos)\n")
        f.write(f"Estimated cost: ${cost:.2f}\n\n")
        for oligo in oligos:
            f.write(f"  {oligo['name']}:\n")
            f.write(f"    5'-{oligo['sequence']}-3'\n")
            f.write(f"    Length: {oligo['length']}nt | Scale: {oligo['scale']} | {oligo['purification']}\n")
            if oligo["modification_5"] != "none":
                f.write(f"    5' mod: {oligo['modification_5']}\n")
            if oligo["modification_3"] != "none":
                f.write(f"    3' mod: {oligo['modification_3']}\n")

        f.write(f"\n{'=' * 70}\n")

    print(f"  Relatório salvo em: {report_path}")
    return report_path


def save_oligo_order(oligos: list[dict]):
    """Salva folha de encomenda em TSV."""
    filepath = os.path.join(REPORTS_DIR, "oligo_order.tsv")
    fieldnames = [
        "name",
        "sequence",
        "length",
        "type",
        "modification_5",
        "modification_3",
        "scale",
        "purification",
        "notes",
    ]

    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for oligo in oligos:
            writer.writerow(oligo)

    print(f"  Oligo order salvo em: {filepath}")


def main():
    print("=" * 60)
    print("MULTIPLEX PANEL - Pipeline CRISPR-Cas12a")
    print("=" * 60)

    os.makedirs(REPORTS_DIR, exist_ok=True)

    # Carregar guides e primers
    guides = {}
    primers = {}
    for gene_name in TARGETS:
        g = load_best_guide(gene_name)
        if g:
            guides[gene_name] = g
        p = load_best_primers(gene_name)
        if p:
            primers[gene_name] = p

    if not guides:
        print("✗ Nenhum guide encontrado. Execute o pipeline na ordem correta!")
        return

    print(f"\n  Guides carregados: {list(guides.keys())}")
    print(f"  Primers carregados: {list(primers.keys())}")

    # Gerar folha de oligos
    oligos = generate_oligo_order(guides, primers)
    cost = estimate_cost(oligos)

    print(f"\n  Total de oligos: {len(oligos)}")
    print(f"  Custo estimado: ${cost:.2f}")

    # Salvar
    save_oligo_order(oligos)
    report_path = generate_report(guides, primers, oligos, cost)

    # Imprimir relatório na tela
    print(f"\n{'=' * 60}")
    with open(report_path) as f:
        print(f.read())

    return {"oligos": oligos, "cost": cost, "report": report_path}


if __name__ == "__main__":
    main()
