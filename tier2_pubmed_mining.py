#!/usr/bin/env python3
"""
Tier 2 — PubMed mining para variantes AMR oftalmologicamente reportadas.

Busca abstracts de endoftalmite/queratite com mencao a AMR, depois usa Claude API
para extrair: variantes mencionadas (ex. mecA, blaKPC-3), patogeno, contexto clinico.

Output: reports/tier2_pubmed/
  - hits_raw.json     (PMID + abstract por query)
  - extractions.json  (variantes estruturadas, normalizadas)
  - summary.csv       (gene_family, n_papers, drugs, ocular_contexts)
"""

import csv
import json
import os
import time
import xml.etree.ElementTree as ET

import requests

from config import REPORTS_DIR

OUT_DIR = os.path.join(REPORTS_DIR, "tier2_pubmed")
os.makedirs(OUT_DIR, exist_ok=True)

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = "smartsepsis@iaparamedicos.com.br"

# Queries dirigidas
QUERIES = [
    '("endophthalmitis"[Title/Abstract] OR "endophthalmitis"[MeSH]) AND ("antibiotic resistance" OR "MRSA" OR "carbapenem" OR "ESBL")',
    '("keratitis"[Title/Abstract]) AND ("antibiotic resistance" OR "Pseudomonas aeruginosa") AND ("resistant" OR "MIC")',
    '("microbial keratitis") AND ("multidrug resistant" OR "fluoroquinolone resistance")',
    '("post-cataract endophthalmitis" OR "postoperative endophthalmitis") AND ("Staphylococcus" OR "Pseudomonas" OR "Klebsiella")',
    '("ocular infection" OR "intraocular infection") AND ("mecA" OR "blaKPC" OR "blaNDM" OR "vanA" OR "mcr")',
    '("conjunctivitis" OR "neonatal ophthalmia") AND ("ESBL" OR "MRSA" OR "carbapenem-resistant")',
    '("pediatric endophthalmitis") AND ("antibiotic resistance")',
]

MAX_PER_QUERY = 30  # cap pra nao explodir


def esearch(query: str, max_n: int) -> list[str]:
    """Retorna PMIDs do esearch."""
    r = requests.get(f"{NCBI_BASE}/esearch.fcgi", params={
        "db": "pubmed", "term": query, "retmax": max_n,
        "retmode": "json", "email": NCBI_EMAIL,
    }, timeout=30)
    r.raise_for_status()
    return r.json().get("esearchresult", {}).get("idlist", [])


def efetch_abstracts(pmids: list[str]) -> list[dict]:
    """Retorna [{pmid, title, abstract, year, journal}, ...]."""
    if not pmids: return []
    r = requests.get(f"{NCBI_BASE}/efetch.fcgi", params={
        "db": "pubmed", "id": ",".join(pmids),
        "rettype": "abstract", "retmode": "xml", "email": NCBI_EMAIL,
    }, timeout=60)
    r.raise_for_status()
    out = []
    root = ET.fromstring(r.text)
    for art in root.findall(".//PubmedArticle"):
        pmid = art.findtext(".//PMID") or ""
        title = art.findtext(".//ArticleTitle") or ""
        abstract_parts = [t.text for t in art.findall(".//Abstract/AbstractText") if t.text]
        abstract = " ".join(abstract_parts).strip()
        year = art.findtext(".//PubDate/Year") or art.findtext(".//PubDate/MedlineDate", "")
        journal = art.findtext(".//Journal/Title") or ""
        if abstract:
            out.append({"pmid": pmid, "title": title, "abstract": abstract,
                        "year": year, "journal": journal})
    return out


def main():
    print("=" * 70)
    print("TIER 2 — PubMed mining (sem Claude ainda, so coleta)")
    print("=" * 70)

    all_hits = []
    seen_pmids = set()
    for q in QUERIES:
        print(f"\n[{q[:80]}...]")
        pmids = esearch(q, MAX_PER_QUERY)
        new = [p for p in pmids if p not in seen_pmids]
        seen_pmids.update(new)
        print(f"  {len(pmids)} PMIDs, {len(new)} novos")
        if not new: continue
        articles = efetch_abstracts(new)
        for a in articles:
            a["query"] = q
        all_hits.extend(articles)
        time.sleep(1)  # NCBI rate limit

    # Salvar raw
    out_raw = os.path.join(OUT_DIR, "hits_raw.json")
    with open(out_raw, "w") as f:
        json.dump(all_hits, f, indent=2)
    print(f"\n[Saved] {out_raw}: {len(all_hits)} abstracts unicos")

    # Quick keyword extraction (proxy ate Claude API)
    KEY_GENES = ["mecA", "blaKPC", "blaNDM", "blaOXA-48", "blaOXA-23", "blaVIM", "blaIMP",
                 "blaCTX-M", "vanA", "vanB", "mcr-1", "mcr-5", "qnrS", "armA",
                 "blaTEM", "blaSHV", "cfrA", "ermA", "ermB", "ermC", "tetM"]
    KEY_PATHO = ["Staphylococcus aureus", "Staphylococcus epidermidis", "Pseudomonas aeruginosa",
                 "Klebsiella pneumoniae", "Escherichia coli", "Enterococcus", "Acinetobacter"]
    KEY_OCULAR = ["endophthalmitis", "keratitis", "conjunctivitis", "microbial",
                  "post-cataract", "intraocular", "vitreous", "corneal"]

    extractions = []
    for h in all_hits:
        text = (h["title"] + " " + h["abstract"]).lower()
        genes = [g for g in KEY_GENES if g.lower() in text]
        pathos = [p for p in KEY_PATHO if p.lower() in text]
        contexts = [c for c in KEY_OCULAR if c.lower() in text]
        if genes or pathos:
            extractions.append({
                "pmid": h["pmid"],
                "title": h["title"][:200],
                "year": h["year"],
                "journal": h["journal"],
                "genes_mentioned": genes,
                "pathogens_mentioned": pathos,
                "ocular_contexts": contexts,
            })

    out_ext = os.path.join(OUT_DIR, "extractions.json")
    with open(out_ext, "w") as f:
        json.dump(extractions, f, indent=2)
    print(f"[Saved] {out_ext}: {len(extractions)} papers com gene+patogeno")

    # Summary CSV
    gene_count = {}
    for e in extractions:
        for g in e["genes_mentioned"]:
            gene_count.setdefault(g, []).append(e["pmid"])
    rows = []
    for g, pmids in sorted(gene_count.items(), key=lambda x: -len(x[1])):
        rows.append({"gene": g, "n_papers": len(pmids), "pmids": ",".join(pmids[:5])})
    out_csv = os.path.join(OUT_DIR, "summary.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene","n_papers","pmids"])
        w.writeheader()
        w.writerows(rows)

    print(f"\nTop genes em context oftalmologico:")
    for r in rows[:10]:
        print(f"  {r['gene']:<14} {r['n_papers']:>3} papers")
    print(f"\n[Saved] {out_csv}")


if __name__ == "__main__":
    main()
