# CRISPR-Cas12a Paper-Based Diagnostic for Hospital Fecal Residue Detection

Hackathon project: a computational pipeline for designing a paper-based CRISPR-Cas12a (DETECTR) diagnostic device that detects antimicrobial-resistant bacteria in hospital fecal residues.

## Overview

A hospital worker collects a sample, applies it to a paper device, and gets a colorimetric/fluorescent result in ~30 minutes. The device detects:

| Spot | Target | Gene | Pathogen |
|------|--------|------|----------|
| 1 | MRSA | mecA | *Staphylococcus aureus* |
| 2 | CRE/CRAB | blaKPC | *Acinetobacter baumannii* |
| P | Positive control | 16S rRNA | Universal bacterial |
| N | Negative control | NTC | No template |

## How it works

```
Sample → Lysis (glass fiber) → Extraction (paper lateral flow)
→ RPA amplification (37°C, 20 min) → Cas12a trans-cleavage → Fluorescence readout
```

## Pipeline

Run scripts in order:

```bash
pip install -r requirements.txt

# 1. Fetch gene DNA sequences from NCBI
python fetch_sequences.py

# 2. Design Cas12a guide RNAs (crRNAs)
python design_guides.py

# 3. Design RPA primers
python design_primers.py

# 4. Check specificity via BLAST (requires internet)
python specificity_check.py

# 5. Generate final panel report + oligo order sheet
python multiplex_panel.py
```

Outputs are generated in `sequences/`, `guides/`, `primers/`, and `reports/` directories.

## Technology

- **CRISPR enzyme**: Cas12a (LbCas12a) — PAM: TTTV
- **Amplification**: RPA (Recombinase Polymerase Amplification) — 37°C, 20 min
- **Detection**: ssDNA-FQ reporter (FAM/BHQ-1) cleaved by Cas12a trans-activity
- **Readout**: UV fluorescence or smartphone camera

## Project structure

```
├── config.py              # Target definitions, Cas12a/RPA parameters
├── utils.py               # Shared functions (reverse complement, GC%, Tm, PAM scan)
├── fetch_sequences.py     # Download gene DNA from NCBI Entrez
├── design_guides.py       # crRNA guide design (PAM scan + scoring)
├── design_primers.py      # RPA primer design
├── specificity_check.py   # BLAST specificity validation
├── multiplex_panel.py     # Final panel assembly + oligo order
├── requirements.txt       # Python dependencies
├── WP_057521704.1/        # NCBI protein data (MecA, S. aureus)
└── WP_063860633.1.zip     # NCBI protein data (KPC-10, A. baumannii)
```

## Data sources

- [WP_057521704.1](https://www.ncbi.nlm.nih.gov/protein/WP_057521704.1) — MecA, *S. aureus*
- [WP_063860633.1](https://www.ncbi.nlm.nih.gov/protein/WP_063860633.1) — KPC-10, *A. baumannii*
- Gene sequences: NCBI RefSeq (NG_047945.1, NG_049243.1)
