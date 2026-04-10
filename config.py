"""
Configuração central do pipeline de diagnóstico CRISPR-Cas12a paper-based.
Hackathon: Detecção de resíduos fecais hospitalares.
"""

import os

# === Diretórios ===
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SEQUENCES_DIR = os.path.join(BASE_DIR, "sequences")
GUIDES_DIR = os.path.join(BASE_DIR, "guides")
PRIMERS_DIR = os.path.join(BASE_DIR, "primers")
REPORTS_DIR = os.path.join(BASE_DIR, "reports")

# === Alvos do painel ===
TARGETS = {
    "mecA": {
        "name": "mecA",
        "description": "Ceftaroline-resistant PBP2a peptidoglycan transpeptidase",
        "organism": "Staphylococcus aureus",
        "protein_accession": "WP_057521704.1",
        "gene_accession": "NG_047945.1",
        "pathogen": "MRSA",
        "clinical_relevance": "Methicillin-resistant S. aureus - major HAI pathogen",
    },
    "blaKPC": {
        "name": "blaKPC",
        "description": "KPC family class A beta-lactamase",
        "organism": "Acinetobacter baumannii",
        "protein_accession": "WP_063860633.1",
        "gene_accession": "NG_049243.1",
        "pathogen": "CRE/CRAB",
        "clinical_relevance": "Carbapenem-resistant - WHO critical priority",
    },
}

# === Parâmetros Cas12a (LbCas12a / AsCas12a) ===
CAS12A = {
    "pam": "TTTV",  # V = A, C, G (não T)
    "pam_strand": "non-target",  # PAM reconhecido na fita não-alvo, upstream
    "spacer_length_min": 20,
    "spacer_length_max": 24,
    "spacer_length_optimal": 20,
    "gc_min": 0.30,
    "gc_max": 0.70,
    "gc_optimal_min": 0.40,
    "gc_optimal_max": 0.60,
    "max_homopolymer": 4,  # máximo de bases idênticas consecutivas
    "top_guides": 5,  # número de guides a reportar por gene
}

# === Parâmetros RPA (Recombinase Polymerase Amplification) ===
RPA = {
    "primer_length_min": 30,
    "primer_length_max": 35,
    "primer_length_optimal": 32,
    "amplicon_min": 100,
    "amplicon_max": 200,
    "amplicon_optimal": 150,
    "tm_min": 54.0,
    "tm_max": 67.0,
    "gc_min": 0.30,
    "gc_max": 0.70,
    "temperature": 37,  # °C
    "time_minutes": 20,
}

# === Reporter Cas12a ===
REPORTER = {
    "type": "ssDNA-FQ",
    "fluorophore": "FAM",
    "quencher": "BHQ-1",
    "sequence": "/56-FAM/TTATT/3BHQ_1/",  # 5nt ssDNA reporter typical
    "description": "ssDNA reporter clivado por atividade trans-cleavage de Cas12a",
}

# === Controles ===
CONTROLS = {
    "positive": {
        "name": "16S rRNA",
        "description": "Universal bacterial 16S rRNA - confirma presença de DNA bacteriano",
        "spot": "P",
    },
    "negative": {
        "name": "NTC",
        "description": "No-template control - confirma ausência de contaminação",
        "spot": "N",
    },
}

# === NCBI Entrez ===
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = "hackathon@example.com"  # required by NCBI policy

# === Cas12a Direct Repeat (scaffold para crRNA) ===
CAS12A_DIRECT_REPEAT = {
    "LbCas12a": "AATTTCTACTAAGTGTAGAT",
    "AsCas12a": "AATTTCTACTCTTGTAGAT",
}
