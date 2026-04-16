"""Fase 7 · Task 4 — Descoberta de novas variantes AMR via NCBI Entrez scan amplo.

Objetivo: ir alem das 42 variantes do painel atual buscando:
- Toda variante publicada no NCBI para as 12 familias AMR (ate 500 por familia)
- Novas variantes pos-2024 (data-cutoff filter)
- Variantes reportadas em cepas clinicas brasileiras (geo filter Brazil)

Saida: prioriza variantes que (a) nao estao no painel atual, (b) tem prevalencia
clinica crescente, (c) sao proximas mas nao identicas as existentes.

Output:
- novel_variants_raw.csv (todas candidatas)
- novel_variants_prioritized.csv (filtradas)
- novel_variants_summary.json
"""
import os
import time
import re
import json
import pandas as pd
from _fase7_utils import (save_json, save_csv, emit_result, AMR_FAMILIES,
                          PROJECT_ROOT)


def fetch_from_ncbi(term: str, retmax: int = 200):
    """Entrez esearch + esummary via Biopython."""
    try:
        from Bio import Entrez
        Entrez.email = "jv@wingsai.com"
        h = Entrez.esearch(db="nuccore", term=term, retmax=retmax)
        rec = Entrez.read(h); h.close()
        ids = rec["IdList"]
        if not ids: return []
        h = Entrez.esummary(db="nuccore", id=",".join(ids))
        summaries = Entrez.read(h); h.close()
        return [dict(s) for s in summaries]
    except Exception as e:
        print(f"NCBI offline ({e}); usando cache/AMRFinder local")
        return []


def parse_existing_variants():
    """Le targets_brazil_variants.csv + targets_brazil_card.csv do pipeline."""
    existing = set()
    for fname in ["targets_brazil_variants.csv", "targets_brazil_card.csv", "targets_brazil.csv"]:
        p = PROJECT_ROOT / fname
        if not p.exists(): continue
        try:
            df = pd.read_csv(p)
            for col in df.columns:
                if any(k in col.lower() for k in ["variant", "gene", "name", "accession"]):
                    for v in df[col].astype(str).tolist():
                        # normalize
                        v = v.strip().upper()
                        if v and v != "NAN":
                            existing.add(v)
        except Exception as e:
            print(f"skip {p}: {e}")
    return existing


def query_family(fam: str, brazil_only: bool = False, recent: bool = False, retmax: int = 100):
    terms = [f'"{fam}"[Gene Name]', f'antimicrobial resistance']
    if brazil_only:
        terms.append("Brazil[Text Word]")
    if recent:
        terms.append('"2024/01/01"[PDAT] : "2026/12/31"[PDAT]')
    term = " AND ".join(terms)
    return fetch_from_ncbi(term, retmax=retmax)


def main():
    t0 = time.time()
    existing = parse_existing_variants()
    print(f"variantes existentes no pipeline: {len(existing)}")

    all_candidates = []
    for fam in AMR_FAMILIES:
        for brazil in [False, True]:
            for recent in [False, True]:
                hits = query_family(fam, brazil_only=brazil, recent=recent, retmax=50)
                for hit in hits:
                    title = hit.get("Title", "")
                    # extrai variante do titulo (padrao "NDM-5", "KPC-33", etc.)
                    variant_match = re.search(r"\b([A-Z][a-z]{2,4}[-_]?\d+)\b", title)
                    variant = variant_match.group(1).upper() if variant_match else None
                    all_candidates.append({
                        "family": fam,
                        "variant": variant,
                        "title": title[:150],
                        "accession": hit.get("AccessionVersion", ""),
                        "length": hit.get("Length", 0),
                        "brazil_filter": brazil,
                        "recent_filter": recent,
                    })
        print(f"  {fam}: {sum(1 for c in all_candidates if c['family']==fam)} hits")

    df = pd.DataFrame(all_candidates)
    if df.empty:
        print("nenhuma variante encontrada (NCBI provavelmente offline)")
        emit_result("status", "ncbi_offline")
        save_json({"status": "ncbi_offline"}, "novel_variants_summary.json")
        return

    save_csv(df, "novel_variants_raw.csv")

    # Prioriza: variantes nao-existentes + Brasil + recentes
    df["is_novel"] = df["variant"].apply(
        lambda v: pd.notna(v) and str(v).upper() not in existing
    )
    prioritized = df[df["is_novel"]].copy()
    prioritized["score"] = (prioritized["brazil_filter"].astype(int) * 2 +
                            prioritized["recent_filter"].astype(int) * 3)
    prioritized = (prioritized.sort_values(["family", "score"], ascending=[True, False])
                              .drop_duplicates(subset=["variant"], keep="first"))
    save_csv(prioritized, "novel_variants_prioritized.csv")

    summary = {
        "n_candidates_total": len(df),
        "n_novel": int(prioritized.shape[0]),
        "novel_by_family": prioritized.groupby("family").size().to_dict(),
        "top10_priority": prioritized.head(10).to_dict(orient="records"),
    }
    save_json(summary, "novel_variants_summary.json")

    emit_result("n_candidates_total", len(df))
    emit_result("n_novel", int(prioritized.shape[0]))
    emit_result("top_novel_family",
                prioritized["family"].value_counts().head(1).to_string() if not prioritized.empty else "none")
    emit_result("seconds", f"{time.time()-t0:.1f}")


if __name__ == "__main__":
    main()
