# SmartSepsis-Oph

**AI-augmented multimodal diagnostics for ophthalmic infection.**

[![Status](https://img.shields.io/badge/status-Phase%200%3A%20computational%20scaffolding-orange)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Stage](https://img.shields.io/badge/stage-wet--lab%20partnership%20in%20formation-green)]()

> **Important.** SmartSepsis-Oph is a research program in computational
> scaffolding phase. The repository builds the AI layer for ocular-infection
> molecular diagnostics — a multimodal reference dataset, structural and
> foundation-model classifiers for AMR genes, candidate CRISPR-Cas12a guide
> libraries, and a planned interpretation pipeline. **No claims of clinical
> performance are made.** All results in this repository are *in silico*.
> Wet-lab validation is being organized with a partner laboratory at Mass
> Eye and Ear / Harvard Medical School.

---

## What this project is

Ocular infection (endophthalmitis, microbial keratitis, post-injection
infections) is the hardest sample in molecular medicine. Volumes are sub-
microliter. The vitreous is an immune sanctuary — normally sterile —
which means the diagnostic question is **"which organism is here?"** not
"is this organism resistant?". Culture fails in 30% or more of cases.
Existing molecular panels were not designed for this matrix.

SmartSepsis-Oph attacks this gap as an AI-first scaffolding program with
two complementary tracks:

### Track A — AI-augmented metagenomic identification (primary)

A pipeline that takes Nanopore metagenomic reads from ocular samples and
returns: organism identification, resistance-gene catalog, predicted drug-
class susceptibility, and structural annotation of AMR variants. Foundation-
model embeddings (ESM-2, ProtT5) and AlphaFold structural features carry
most of the interpretive load on noisy, low-coverage reads. **This is the
primary scientific track**, developed in collaboration with a renowned
ocular-infection laboratory at Mass Eye and Ear / Harvard Medical School
(senior advisor and host lab — full disclosure of names will follow public
announcement).

### Track B — CRISPR-Cas12a resistance triage (adjunct)

A paper-strip Cas12a + RPA panel for resistance-gene detection. Designed
for time-critical empirical guidance in the OR, or for settings without
nanopore capacity. Not intended as a primary identification tool — its
clinical use is conditional on suspected organism. Track B inherits the
213 computationally designed guide candidates across 12 resistance-gene
families produced in Phase 0.

---

## What this project is *not*

- Not a medical device.
- Not validated experimentally.
- Not a substitute for culture, antibiogram, or clonal surveillance.
- Not approved for any clinical use under ANVISA, FDA, or any regulatory
  body.

Performance metrics (sensitivity, specificity, limit of detection,
time-to-result) will only be reported after experimental validation with
the partner laboratory.

---

## Phased roadmap

| Phase | Scope | Status |
|---|---|---|
| **0 — Computational scaffolding** | Multimodal dataset, foundation-model classifiers, candidate guide library, scoping review | **Current** |
| **1 — Scoping review** | PRISMA-ScR mapping CRISPR-Dx and Nanopore metagenomic Dx for ocular infection | In drafting |
| **2 — Pipeline development** | AI-augmented Nanopore taxonomic + AMR pipeline; benchmarks against curated reads | Planned, partner laboratory engaged |
| **3 — Clinical validation** | Prospective sampling under IRB, performance against culture | Planned |
| **4 — Regulatory** | IVD pathway (RDC 830/2023, equivalents) | Planned |

---

## Team and advisors

- **João Victor Pacheco Dias** — Founder & AI Lead. Doctoral candidate at
  HC-FMUSP (Medical AI). CTO WingsAI. ITU/WHO member, AI for Health.
  Technical advisor CBO.
- **Dr. Gustavo Sakuno** — Co-author. Postdoctoral fellow at Mass Eye and
  Ear / Harvard Medical School. PhD USP, Ophthalmology and Oculomics.
  Clinical and translational lead.
- **Senior scientific advisor** — *Renowned ocular-infection investigator
  at Mass Eye and Ear / Harvard Medical School. Name to be announced
  following formal confirmation. Will provide mentorship and access to her
  laboratory's wet-lab and Nanopore metagenomics capability.*

**We continue to seek** additional collaborators in CRISPR-Dx, ophthalmic
microbiology, and AI-for-genomics. Co-authorship on forthcoming papers
remains on offer for early contributors.

Contact: **iaparamedicos@gmail.com**

---

## What is in this repository

```
.
├── README.md                  # this file
├── LICENSE                    # MIT
├── CONTRIBUTING.md            # how to engage with the project
├── CITATION.cff               # academic citation metadata
├── CHANGELOG.md               # versioned changes
├── MIGRATION.md               # planned src/ migration
├── public/                    # project website (PT-BR / EN)
├── docs/                      # MkDocs Material site (also deployed to gh-pages)
├── preprint/                  # framework preprint, scoping review drafts
├── scripts/                   # operational scripts (stress test, etc.)
├── rtx5090/                   # compute pack for the GPU machine
├── alphafold3_submission/     # AlphaFold 3 ternary-complex input pack
├── src/                       # pipeline modules (migration in progress)
├── data/                      # output artifacts, including HuggingFace dataset
└── *.py                       # design pipeline scripts (see below)
```

Key Python entry points (still at repo root pending [MIGRATION.md](MIGRATION.md)):

- `fetch_sequences.py` — NCBI reference pulls
- `design_guides.py` — Cas12a guide RNA design (Track B)
- `covariance_probes.py` — biophysical re-scoring of guides
- `design_primers.py` — RPA primer design
- `specificity_check.py` — BLAST-based specificity (Track B)
- `evo2_scoring.py` — functional variant impact scoring
- `card_integration.py` — CARD enrichment
- `protein_structure.py` / `prott5_ensemble.py` / `phenotype_probe*.py` — structural and PLM-based AMR characterization (Track A backbone)
- `multiplex_panel.py` — panel optimization
- `run_batch.py` — batch orchestration

## Published artifacts

- **Multimodal dataset** on HuggingFace:
  https://huggingface.co/datasets/JVLegend/smartsepsis-oph
  - `panel` (45 curated variants, ESM-2 + ProtT5 + AlphaFold embeddings + structural descriptors)
  - `extended` (9,034 rows)
- **Candidate guide library (Track B)**:
  `fase7_dgx_results/fase7_results/rnafold_guides.csv` — 213 ranked
  candidate guides across 42 target entries (12 reference gene families).
- **Documentation site**: published via MkDocs Material; build/deploy
  workflow in `.github/workflows/docs.yml`.

## How to engage

If you are a researcher, clinician, or funder who can move this from
Phase 0 forward:

- See [CONTRIBUTING.md](CONTRIBUTING.md).
- Write to **iaparamedicos@gmail.com**.
- Open an issue with the `partnership`, `advisor`, or `funding` label.

## License

MIT — see [LICENSE](LICENSE). Open by design.

## Disclaimer

SmartSepsis-Oph is a research program by **IA para Médicos**. The code and
designs in this repository are for research use only. No claim of clinical
performance is made. No medical device approval has been sought or granted.
Do not use any output of this pipeline to make clinical decisions.
