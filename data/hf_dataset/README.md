---
license: cc-by-4.0
language:
- en
tags:
- biology
- protein
- AMR
- antimicrobial-resistance
- ophthalmology
- endophthalmitis
- crispr-cas12a
- protein-language-models
- esm2
- prott5
- alphafold
- multimodal
size_categories:
- n<1K
task_categories:
- text-classification
- feature-extraction
- token-classification
pretty_name: SmartSepsis-Oph — multimodal AMR variants for ophthalmology
configs:
- config_name: panel
  data_files:
  - split: train
    path: smartsepsis_oph.parquet
- config_name: extended
  data_files:
  - split: train
    path: smartsepsis_oph_extended.parquet
---

# SmartSepsis-Oph: Multimodal AMR variants for ophthalmology research

Curated multimodal dataset of 43 antimicrobial resistance (AMR) gene variants relevant
to ocular bacterial pathogens (endophthalmitis, keratitis, perioperative prophylaxis),
each annotated with DNA sequence, protein sequence, ESM-2 + ProtT5 embeddings,
predicted 3D structure (PDB), structural descriptors, drug class labels and resistance
mechanism.

Built as the public dataset accompanying the SmartSepsis-Oph research line at
HC-FMUSP × Mass Eye and Ear (Harvard Medical School), led by Dr. Gustavo Sakuno
(oculomics & multi-omics biomarkers).

## Why this dataset

Existing AMR resources (CARD, AMRFinderPlus, NDARO) provide reference sequences but
no aligned multimodal annotations. To train and benchmark protein-language-model
classifiers, structural property predictors, and AI-driven CRISPR diagnostics targeting
**ocular pathogens specifically**, we curated a compact set of 43 clinically prevalent
variants and pre-computed every modality used in our pipeline. Inspired by the
[NanoFold-public](https://huggingface.co/datasets/ChrisHayduk/nanofold-public)
distribution pattern (Hayduk 2026).

## Coverage

| Family | Variants | Drug class | Source |
|---|---|---|---|
| mecA | mecA1, mecA2 | penam, cephalosporin, methicillin (MRSA) | RefSeq |
| blaKPC | KPC-2, 3, 4, 5, 11, 30, 31 | carbapenem, cephalosporin, penam | RefSeq |
| blaNDM | NDM-1, 2, 5, 7 | carbapenem, cephalosporin, penam | RefSeq |
| blaOXA-48 | OXA-48, 181, 232 | carbapenem, cephalosporin, penam | RefSeq |
| blaVIM | VIM-1, 2, 4 | carbapenem, cephalosporin, penam | RefSeq |
| blaIMP | IMP-1, 6 | carbapenem, cephalosporin, penam | RefSeq |
| blaGES | GES-1, 5 | carbapenem, cephalosporin | RefSeq |
| blaCTX-M | CTX-M-2, 8, 9, 14, 27 | cephalosporin, penam (ESBL) | RefSeq |
| vanA | vanA | glycopeptide (vancomycin) | RefSeq |
| mcr | mcr-1, mcr-1.1, mcr-5 | polymyxin, peptide | RefSeq |
| qnrS | qnrS1, qnrS2 | fluoroquinolone | RefSeq |
| armA | armA | aminoglycoside | RefSeq |

## Schema

```python
{
    "variant_id": str,                 # "blaKPC-3"
    "gene_family": str,                # "blaKPC"
    "dna_accession": str,              # "NG_049257.1"
    "dna_sequence": str,
    "protein_sequence": str,           # longest ORF, table 11 (bacterial)
    "protein_length": int,             # AA
    "drug_classes": list[str],         # ["carbapenem", "cephalosporin", "penam"]
    "resistance_mechanism": str,       # "antibiotic inactivation"
    "esm2_embedding": list[float],     # 640d, mean-pooled
    "esm2_model": str,                 # "esm2_t30_150M_UR50D"
    "prott5_embedding": list[float],   # 1024d, mean-pooled
    "prott5_model": str,               # "Rostlab/prot_t5_xl_uniref50"
    "structure_pdb": str,              # full PDB text
    "structure_source": str,           # "ColabFold/ESMFold/AlphaFoldServer"
    "struct_length": int,              # CA atoms
    "struct_rg": float,                # radius of gyration (A)
    "struct_compactness": float,       # Rg / L^0.6
    "struct_contact_density": float,   # fraction of CA-CA pairs <8 A
    "struct_mean_plddt": float         # 0-100, prediction confidence
}
```

## Two configs

| Config | Rows | Size | What's inside |
|---|---|---|---|
| `panel` | 45 (43 com tudo) | 3.2 MB | Multimodal completo: DNA + protein + ESM-2 + ProtT5 + PDB + struct features + drug class |
| `extended` | **9.034** | 34 MB | AMRFinderPlus catalog (8.991) + panel (43): variant_id + source + drug_classes + ESM-2 embedding |

## How to use

```python
from datasets import load_dataset
import numpy as np

# Multimodal panel (curated, full pipeline)
ds_panel = load_dataset("jvlegend/smartsepsis-oph", "panel", split="train")
row = ds_panel[0]
print(row["variant_id"], row["gene_family"], row["protein_length"])
emb = np.array(row["esm2_embedding"])           # 640d
ensemble = np.concatenate([emb, np.array(row["prott5_embedding"])])  # 1664d
print("ensemble shape:", ensemble.shape)

# Extended (9034 entries from AMRFinderPlus + panel, ESM-2 + drug class)
ds_ext = load_dataset("jvlegend/smartsepsis-oph", "extended", split="train")
print(f"Extended: {len(ds_ext)} entries, {len(set(ds_ext['source']))} sources")
# -> Extended: 9034 entries, 2 sources
```

## Dataset construction

1. **Source acquisition** — variants pulled from NCBI RefSeq via Entrez API (NG_*
   accessions curated against AMRFinderPlus / CARD ontology).
2. **Translation** — longest ORF using bacterial genetic code (table 11) via Biopython.
3. **ESM-2 embeddings** — `esm2_t30_150M_UR50D`, mean-pooled across residues.
4. **ProtT5 embeddings** — `Rostlab/prot_t5_xl_uniref50`, mean-pooled.
5. **3D structure prediction** — combination of:
   - ESMFold via HuggingFace transformers (proteins ≤400 aa)
   - ColabFold AF2 (mcr-1, mcr-5 ~540 aa)
   - AlphaFold Server AF3 (mecA1, mecA2 ~665 aa)
   PDB rank_1 selected per variant.
6. **Structure descriptors** — Rg, compactness ratio, contact density, mean Cα-Cα,
   aspect ratio, mean pLDDT computed from CA coordinates.
7. **Drug class / mechanism labels** — derived from CARD ontology (CC BY 4.0)
   harmonized to clinically meaningful classes.

## Companion code

Pipeline source at <https://github.com/JVLegend/smartsepsis>. Includes:
- CRISPR-Cas12a guide design (`design_guides.py`)
- Multi-label OvR classifier with NanoFold-augmented negative calibration
  (`phenotype_probe_v2.py`)
- Structure-aware ensemble (`structure_features_v3.py`)
- Pangenome of 21 K. pneumoniae + E. coli isolates (`pangenome.sh`)

## Personal & sensitive information

**None.** All data is derived from public NCBI RefSeq sequences and predicted
structures. No patient data, biological samples, or PII. Embeddings are
deterministic functions of the public sequences.

## Considerations for use

- **Predicted structures** carry inherent uncertainty (mean pLDDT ~85-95 across the
  panel, but local regions can be lower). Use the per-variant `struct_mean_plddt`
  for downstream weighting.
- **AlphaFold Server (AF3) terms** apply to the mecA1/mecA2 PDBs — research use only,
  no commercial redistribution. For commercial use, regenerate via ColabFold AF2 or
  AlphaFold 2 OpenFold.
- **Class imbalance** — drug class distribution mirrors clinical relevance (heavy on
  β-lactams, lighter on glycopeptide/polymyxin). For balanced training, augment with
  AMRFinderPlus catalog (8,991 sequences) referenced in the companion paper.

## Citation

```bibtex
@dataset{smartsepsis_oph_2026,
  title  = {{SmartSepsis-Oph}: Multimodal AMR variants for ophthalmology research},
  author = {Dias, Jo{\~a}o Victor and Sakuno, Gustavo and Primo, Raul},
  year   = {2026},
  url    = {https://huggingface.co/datasets/JVLegend/smartsepsis-oph},
  doi    = {tbd},
  note   = {Dataset accompanying the SmartSepsis-Oph research line, HC-FMUSP × Mass Eye and Ear (Harvard).}
}
```

## Authors & affiliations

- **João Victor Dias** — CTO & AI Architect, IA para Médicos; PhD candidate HC-FMUSP
- **Dr. Gustavo Sakuno** — Clinical & Scientific Lead, postdoc Harvard Medical School / Mass Eye and Ear; PhD USP — Ophthalmology & Oculomics
- **Raul Primo** — Software Engineer, IA para Médicos

## Conflict of interest

Authors declare affiliation with **IA para Médicos** (project sponsor). No financial
conflicts.

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

You are free to share and adapt for any purpose, including commercial — provided
you give appropriate credit (cite as above) and indicate if changes were made.

## Contact

Issues, corrections, additions: open an issue on
<https://github.com/JVLegend/smartsepsis> or reach out via
<https://www.iaparamedicos.com.br/>.
