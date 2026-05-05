#!/usr/bin/env python3
"""
Item #4 — Empacotar SmartSepsis-Oph dataset em formato HuggingFace.

Output: data/hf_dataset/
  - smartsepsis_oph.parquet     (43 linhas multimodais)
  - README.md                   (dataset card seguindo padrao HF)
  - LICENSE                     (CC BY 4.0)
  - manifest.json               (SHA256 hashes pra integridade)

NAO faz upload — somente prepara o pacote local. Depois confirmar com user.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os

import numpy as np
import pandas as pd

from config import REPORTS_DIR, BASE_DIR

ESM_MODEL_NAME = "esm2_t30_150M_UR50D"

OUT_DIR = os.path.join(BASE_DIR, "data", "hf_dataset")
os.makedirs(OUT_DIR, exist_ok=True)

# --- Paths
TARGETS_CSV = os.path.join(BASE_DIR, "targets_brazil_variants.csv")
EFFECTS_CSV = os.path.join(REPORTS_DIR, "protein_scoring", "protein_variant_effects.csv")
ESM_EMB_DIR = os.path.join(REPORTS_DIR, "protein_scoring", "embeddings")
PT5_EMB_DIR = os.path.join(REPORTS_DIR, "prott5_ensemble", "embeddings")
PROT_DIR = os.path.join(REPORTS_DIR, "protein_scoring", "proteins")
DNA_DIR = os.path.join(BASE_DIR, "sequences")
PDB_DIR = os.path.join(BASE_DIR, "public", "pdbs")
STRUCT_CSV = os.path.join(REPORTS_DIR, "phenotype_probe_v3", "structure_features.csv")

FAMILY_LABELS = {
    "mecA":        (["penam","cephalosporin","methicillin"], "antibiotic target replacement"),
    "mecA1":       (["penam","cephalosporin","methicillin"], "antibiotic target replacement"),
    "mecA2":       (["penam","cephalosporin","methicillin"], "antibiotic target replacement"),
    "blaKPC":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
    "blaNDM":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
    "blaOXA":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
    "blaVIM":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
    "blaIMP":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
    "blaGES":      (["carbapenem","cephalosporin"], "antibiotic inactivation"),
    "blaCTX-M":    (["cephalosporin","penam"], "antibiotic inactivation"),
    "blaCTX-M-15": (["cephalosporin","penam"], "antibiotic inactivation"),
    "vanA":        (["glycopeptide"], "antibiotic target alteration"),
    "mcr-1":       (["polymyxin","peptide"], "antibiotic target alteration"),
    "mcr-5":       (["polymyxin","peptide"], "antibiotic target alteration"),
    "qnrS":        (["fluoroquinolone"], "antibiotic target protection"),
    "armA":        (["aminoglycoside"], "antibiotic target alteration"),
}


def labels_for(name):
    if name in FAMILY_LABELS: return FAMILY_LABELS[name]
    for k in FAMILY_LABELS:
        if name.startswith(k): return FAMILY_LABELS[k]
    return (["unknown"], "unknown")


def read_text(path):
    if not os.path.exists(path): return ""
    return open(path).read()


def read_fasta_seq(path):
    if not os.path.exists(path): return ""
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


def get_dna_accession(variant_name):
    """Carrega accession NCBI a partir do CSV."""
    if not os.path.exists(TARGETS_CSV): return None
    with open(TARGETS_CSV) as f:
        for row in csv.DictReader(f):
            if row.get("name") == variant_name:
                return row.get("gene_accession")
    return None


def main():
    # 1. Coletar todos variants do effects.csv (35) + adicionar refs do targets
    variants = set()
    if os.path.exists(EFFECTS_CSV):
        for row in csv.DictReader(open(EFFECTS_CSV)):
            variants.add(row["variant"])

    # adicionar familias-base (refs) que nao aparecem como variant
    for fam in FAMILY_LABELS:
        variants.add(fam)

    variants = sorted(variants)
    print(f"Total candidatos: {len(variants)}")

    # 2. Carregar struct features (para enriquecer com 7 descritores)
    struct_dict = {}
    if os.path.exists(STRUCT_CSV):
        for row in pd.read_csv(STRUCT_CSV).to_dict(orient="records"):
            struct_dict[row["variant"]] = {
                "L": int(row["L"]) if row["L"] else 0,
                "rg": float(row["rg"]),
                "compactness_ratio": float(row["compactness_ratio"]),
                "contact_density": float(row["contact_density"]),
                "mean_ca_dist_norm": float(row["mean_ca_dist_norm"]),
                "aspect_ratio": float(row["aspect_ratio"]),
                "mean_plddt": float(row["mean_plddt"]),
            }

    # 3. Construir dataset
    rows = []
    for variant in variants:
        drugs, mech = labels_for(variant)
        family = variant
        if variant in FAMILY_LABELS:
            family = variant
        else:
            for k in FAMILY_LABELS:
                if variant.startswith(k):
                    family = k
                    break

        dna_seq = read_fasta_seq(os.path.join(DNA_DIR, f"{variant}.fasta"))
        prot_seq = read_fasta_seq(os.path.join(PROT_DIR, f"{variant}.fasta"))
        esm_path = os.path.join(ESM_EMB_DIR, f"{variant}__{ESM_MODEL_NAME}.npy")
        pt5_path = os.path.join(PT5_EMB_DIR, f"{variant}.npy")
        pdb_path = os.path.join(PDB_DIR, f"{variant}.pdb")

        esm_emb = np.load(esm_path).tolist() if os.path.exists(esm_path) else None
        pt5_emb = np.load(pt5_path).tolist() if os.path.exists(pt5_path) else None
        pdb_text = read_text(pdb_path) if os.path.exists(pdb_path) else None
        struct = struct_dict.get(variant)

        row = {
            "variant_id": variant,
            "gene_family": family,
            "dna_accession": get_dna_accession(variant),
            "dna_sequence": dna_seq,
            "protein_sequence": prot_seq,
            "protein_length": len(prot_seq),
            "drug_classes": drugs,
            "resistance_mechanism": mech,
            "esm2_embedding": esm_emb,           # 640 floats
            "esm2_model": ESM_MODEL_NAME,
            "prott5_embedding": pt5_emb,         # 1024 floats
            "prott5_model": "Rostlab/prot_t5_xl_uniref50",
            "structure_pdb": pdb_text,
            "structure_source": "ColabFold/ESMFold/AlphaFoldServer (mixed)",
            "struct_length": struct["L"] if struct else None,
            "struct_rg": struct["rg"] if struct else None,
            "struct_compactness": struct["compactness_ratio"] if struct else None,
            "struct_contact_density": struct["contact_density"] if struct else None,
            "struct_mean_plddt": struct["mean_plddt"] if struct else None,
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    print(f"Dataset rows: {len(df)}")
    print(f"Schema: {list(df.columns)}")
    print(f"Coverage:")
    print(f"  com DNA seq:        {df['dna_sequence'].astype(bool).sum()}")
    print(f"  com protein seq:    {df['protein_sequence'].astype(bool).sum()}")
    print(f"  com ESM-2 emb:      {df['esm2_embedding'].notna().sum()}")
    print(f"  com ProtT5 emb:     {df['prott5_embedding'].notna().sum()}")
    print(f"  com PDB:            {df['structure_pdb'].notna().sum()}")
    print(f"  com struct features:{df['struct_rg'].notna().sum()}")

    # 4. Save Parquet
    parquet_path = os.path.join(OUT_DIR, "smartsepsis_oph.parquet")
    df.to_parquet(parquet_path, compression="snappy")
    sz_mb = os.path.getsize(parquet_path) / 1024 / 1024
    print(f"\n[Saved] {parquet_path}  ({sz_mb:.1f} MB)")

    # 5. SHA256 manifest
    manifest = {}
    for fname in os.listdir(OUT_DIR):
        if fname.endswith((".parquet", ".md", ".txt")):
            with open(os.path.join(OUT_DIR, fname), "rb") as f:
                manifest[fname] = hashlib.sha256(f.read()).hexdigest()
    with open(os.path.join(OUT_DIR, "manifest.json"), "w") as f:
        json.dump({
            "dataset": "smartsepsis-oph",
            "version": "1.0.0",
            "generated_at": pd.Timestamp.utcnow().isoformat(),
            "n_rows": len(df),
            "files": manifest,
        }, f, indent=2)
    print(f"[Saved] manifest.json com hashes SHA256")


if __name__ == "__main__":
    main()
