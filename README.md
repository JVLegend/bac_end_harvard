# SmartSepsis-Oph

**AI-first molecular diagnostics design for ophthalmic infections.**

[![Status](https://img.shields.io/badge/status-Phase%200%3A%20computational%20design-orange)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Stage](https://img.shields.io/badge/stage-seeking%20wet--lab%20partners-red)]()

> **Important — read this first.** SmartSepsis-Oph is a research program in
> computational design phase. This repository contains an AI-driven pipeline
> for designing CRISPR-Cas12a guide-RNA libraries, RPA primers, and a
> paper-strip assay architecture targeted at ophthalmic infections. **No
> claims of clinical performance are made.** All results in this repository
> are *in silico*. Experimental validation is pending and is what we are
> actively seeking collaborators for.

---

## What this project is

Bacterial endophthalmitis can blind an eye in 24–48 hours. Cultures take
48–72. Existing molecular diagnostics are built for sample volumes orders of
magnitude larger than what ophthalmology yields (sub-microliter vitreous tap,
corneal scrape, ocular swab).

SmartSepsis-Oph attacks this gap with an AI-first design pipeline:

- **Guide-RNA library design** (CRISPR-Cas12a) for 12 resistance-gene
  families relevant to ocular infection, prioritized for the Brazilian
  epidemiological context.
- **In silico specificity modeling** against public reference repositories.
- **Isothermal RPA primer design** compatible with thermocycler-free,
  point-of-care workflows.
- **Functional impact scoring** of variants (inspired by Evo 2 / EVEE).
- **Paper-strip assay architecture** under design.

First declared target: ***Staphylococcus aureus* / mecA** — the most common
and best-validated combination in the CRISPR-Dx literature.

## What this project is *not*

- Not a medical device.
- Not validated experimentally.
- Not a substitute for culture, antibiogram, or clonal surveillance.
- Not approved for any clinical use under ANVISA, FDA, or any regulatory
  body.

Performance metrics (sensitivity, specificity, limit of detection,
time-to-result) will only be reported after experimental validation with a
wet-lab partner.

---

## Phased roadmap

| Phase | Scope | Status |
|---|---|---|
| **0 — Computational design** | Guide and primer libraries, in silico specificity, scoring, structural analysis. | **Current** |
| **1 — Wet-lab feasibility** | Single-organism, single-gene validation (S. aureus / mecA) in cultured isolates. | Seeking wet-lab partner |
| **2 — Panel expansion** | Multiplex of priority resistance markers, real-amplicon validation. | Planned |
| **3 — Clinical validation** | Prospective sampling under IRB, performance against culture and PCR. | Planned |
| **4 — Regulatory submission** | IVD (RDC 830/2023) pathway and equivalent. | Planned |

---

## Team

- **João Victor Pacheco Dias** — Founder & AI lead. Doctoral candidate
  HC-FMUSP (Medical AI). CTO WingsAI. ITU/WHO member, AI for Health.
  Technical advisor CBO.
- **Dr. Gustavo Sakuno** — Clinical advisor. Postdoc, Mass Eye and Ear /
  Harvard Medical School. PhD USP, Ophthalmology & Oculomics.

**We are actively seeking:**
- Senior scientific advisors in CRISPR-Dx, ocular microbiology, and
  molecular pathology
- Wet-lab collaborators with CRISPR-Dx or isothermal amplification capability
- Clinical partners with endophthalmitis case volume
- Seed-stage funders

Contact: **iaparamedicos@gmail.com**

Co-authorship on the forthcoming framework preprint is offered to early
collaborators.

---

## Repository layout

```
.
├── README.md                  # this file
├── LICENSE                    # MIT
├── CONTRIBUTING.md            # how to engage with the project
├── CITATION.cff               # academic citation metadata
├── public/                    # project website (PT-BR / EN)
├── src/                       # pipeline modules (work in progress)
├── data/                      # output artifacts (json, csv)
├── card_data/                 # CARD database integration cache
├── *.py                       # design pipeline scripts (see below)
└── docs/                      # design notes and rationale
```

### Pipeline scripts at root

These scripts implement the design pipeline. Stable entry points:

- `fetch_sequences.py` — pulls reference sequences from NCBI
- `design_guides.py` — CRISPR-Cas12a guide RNA design
- `covariance_probes.py` — re-scoring with 18 biophysical features
- `design_primers.py` — RPA primer design (isothermal, 37 °C)
- `specificity_check.py` — BLAST-based specificity analysis
- `evo2_scoring.py` — functional variant impact scoring
- `card_integration.py` — CARD database enrichment
- `clinical_interpreter.py` — natural-language interpretation
- `multiplex_panel.py` — multiplex panel optimization
- `run_batch.py` — batch orchestration

> Note: this layout reflects the project's research-prototype origin. We
> intend to migrate everything under `src/` and provide a stable
> `pyproject.toml` ahead of the framework preprint. Until then, expect rough
> edges.

---

## Reproducing the design pipeline

> *To be expanded.* We are preparing reproducibility documentation alongside
> the framework preprint. If you want to reproduce a specific result now,
> open an issue or write to us.

A minimal flow looks like:

```bash
python fetch_sequences.py
python design_guides.py
python covariance_probes.py
python design_primers.py
python specificity_check.py
python run_batch.py 4
```

---

## How to engage

If you are a researcher, clinician, or funder who can move this from Phase
0 to Phase 1:

- See [CONTRIBUTING.md](CONTRIBUTING.md).
- Write to **iaparamedicos@gmail.com**.
- Open an issue with the `partnership` or `advisor` label.

If you are using ideas from this work in your own research, see
[CITATION.cff](CITATION.cff).

---

## License

MIT — see [LICENSE](LICENSE). We chose a permissive license because we
believe open science accelerates this field. Derivative works that lead to
clinical products should comply with all applicable regulatory frameworks.

---

## Disclaimer

SmartSepsis-Oph is a research program by **IA para Médicos**. The code and
designs in this repository are for research use only. No claim of clinical
performance is made. No medical device approval has been sought or granted.
Do not use any output of this pipeline to make clinical decisions.
