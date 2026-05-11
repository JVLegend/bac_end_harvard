<!--
==============================================================================
PIVOT NOTICE - 2026-05-11

This draft was written under the program's original framing as a paper-strip
CRISPR-Cas12a diagnostic. The program has since pivoted (see /PIVOT.md):

  Track A (primary)  : AI-augmented Nanopore metagenomic interpretation,
                       collaboration with a senior advisor and host lab at
                       Mass Eye and Ear / Harvard.
  Track B (adjunct)  : the CRISPR-Cas12a design library described here,
                       repositioned as resistance-triage adjunct module.

Action items before submission:
  1. Hold this draft. The first publication output of the program is now the
     SCOPING REVIEW (preprint/scoping_review_*.md).
  2. After the scoping review lands, this draft will be rewritten as a
     two-track framework paper covering both Track A and Track B.
  3. Author list will be updated to include the Senior Scientific Advisor
     (name pending public announcement).

Keep the current text only as material to mine for Track B sections of the
eventual two-track paper.
==============================================================================

Title alternatives (pick one before submission):
  A. SmartSepsis-Oph: an AI-first framework for designing CRISPR-Cas12a point-of-care
     diagnostics for ophthalmic infections — a computational proof-of-concept
  B. Designing CRISPR-Cas12a diagnostics for the ophthalmic sample-volume bottleneck:
     an open, AI-first computational pipeline and multimodal AMR dataset (SmartSepsis-Oph)
  C. SmartSepsis-Oph: in silico design of Cas12a guide RNAs, RPA primers, and a
     multimodal AMR dataset for ocular infections (Phase 0)

Style: US English. No experimental results. All findings stated as in silico.
-->

# [DEFERRED - PIVOT NOTICE ABOVE] SmartSepsis-Oph: an AI-first framework for designing CRISPR-Cas12a point-of-care diagnostics for ophthalmic infections — a computational proof-of-concept

**João Victor Pacheco Dias¹\*, Gustavo Sakuno²**

¹ IA para Médicos & Faculdade de Medicina da Universidade de São Paulo / Hospital das Clínicas (HC-FMUSP), São Paulo, Brazil.
² Mass Eye and Ear / Harvard Medical School, Boston, USA.

\*Correspondence: iaparamedicos@gmail.com

---

## Abstract

**Background.** Bacterial endophthalmitis and microbial keratitis are sight-threatening
ocular infections in which empirical broad-spectrum therapy is initiated before culture
results are available, often 48–72 h later. Conventional molecular diagnostics — including
syndromic PCR panels — are designed for sample volumes orders of magnitude larger than
those routinely obtained in ophthalmology (sub-microliter vitreous taps, corneal scrapes,
ocular surface swabs). No commercial point-of-care (POC) molecular test is currently
optimized for ocular fluids in the Brazilian epidemiological context.

**Methods.** We describe SmartSepsis-Oph, an open-source, AI-first computational pipeline
that designs candidate CRISPR-Cas12a guide RNAs (crRNAs), recombinase polymerase
amplification (RPA) primers, and multiplex panel architectures for twelve antimicrobial
resistance (AMR) gene families relevant to ocular microbiology. The pipeline combines
PAM-aware crRNA scanning, an 18-feature biophysical covariance probe re-scoring step
(inspired by Evo 2 / EVEE-style covariance pooling), isothermal-compatible RPA primer
design, in silico specificity analysis against public sequence repositories, functional
variant scoring, ESM-2 and ProtT5 protein-language-model embeddings, and AlphaFold-/ESMFold-derived
structural descriptors. CARD ontology enrichment links every variant to mechanism, drug
class, and ARO identifiers.

**Results (in silico).** We designed a candidate library spanning 12 AMR families and
released a curated multimodal dataset of 45 variants (panel) and 9,034 rows
(extended) on HuggingFace, comprising DNA and protein sequences, mean-pooled ESM-2
(640-dim) and ProtT5 (1024-dim) embeddings, predicted 3D structures, and structural
descriptors. We report in silico specificity, ranked guide candidates, primer-pair
selection, family-level protein-distance matrices, and a multiplex panel layout. **No wet-lab
validation has been performed.**

**Conclusions.** SmartSepsis-Oph is a transparent Phase 0 computational design framework
and dataset release intended to seed experimental validation by wet-lab and clinical
partners. No claims of analytical or clinical performance are made. A falsifiable, phased
roadmap to validation is proposed.

**Keywords:** CRISPR-Cas12a, point-of-care, endophthalmitis, keratitis, antimicrobial
resistance, protein language models, computational design, in silico.

---

## 1. Introduction

### 1.1 The clinical problem

Bacterial endophthalmitis is a vision-threatening ocular emergency in which intraocular
bacterial proliferation can cause irreversible retinal damage within 24–48 h
[REF: Durand 2017 NEJM endophthalmitis]. Standard practice is empirical intravitreal
broad-spectrum therapy (typically vancomycin and ceftazidime) initiated immediately at
the bedside, with cultures and antibiograms returning 48–72 h later
[REF: EVS Endophthalmitis Vitrectomy Study]. Microbial keratitis carries a similar
"treat first, identify later" pattern and is a leading global cause of monocular
blindness [REF: Ung 2019 keratitis epidemiology].

### 1.2 The sample-volume bottleneck

Unlike sepsis (whole-blood, milliliter volumes) or respiratory infections (swabs and
sputum), ophthalmologic sampling yields sub-microliter to low-microliter volumes:

- Vitreous tap: typically 0.1–0.2 mL [REF: vitreous tap volume].
- Aqueous tap: 0.05–0.1 mL.
- Corneal scrape: <10 µL of recoverable material [REF: corneal scrape yield].

Mainstream syndromic PCR panels (e.g., FilmArray-class systems) require sample volumes
and nucleic-acid inputs that exceed what ocular sampling realistically delivers, and
their target panels are calibrated to European or North American epidemiology
[REF: BioFire panel composition]. No commercial POC molecular diagnostic is currently
validated for ocular fluids in the Brazilian setting [VERIFY: confirm no ANVISA-approved
ocular IVD as of 2026].

### 1.3 CRISPR diagnostics

CRISPR-based diagnostics (CRISPR-Dx) exploit the collateral trans-cleavage activity of
type-V (Cas12a) and type-VI (Cas13) effectors upon target-specific guide-RNA recognition,
typically coupled to isothermal pre-amplification (recombinase polymerase amplification,
RPA, or loop-mediated isothermal amplification, LAMP), and a paper-strip lateral-flow
readout [REF: Gootenberg 2017 SHERLOCK; Chen 2018 DETECTR]. Cas12a (formerly Cpf1) is
particularly attractive for DNA targets because it (i) recognizes a short T-rich PAM
(TTTV), (ii) tolerates a paper-strip readout via FAM/biotin reporter cleavage, and
(iii) operates at near-isothermal conditions (37 °C) compatible with battery-powered
field hardware.

CRISPR-Dx has been deployed at scale for SARS-CoV-2, mpox, and tuberculosis
[REF: SHERLOCK COVID; DETECTR COVID]. To our knowledge, no CRISPR-Cas12a assay is
commercially deployed for ocular pathogens.

### 1.4 Our contribution

We present SmartSepsis-Oph, a Phase 0 computational program comprising:

1. An open-source, AI-first pipeline that designs candidate crRNAs and RPA primers for
   twelve AMR gene families relevant to ocular bacterial infections, prioritized for
   Brazilian epidemiology.
2. Covariance-probe re-scoring of crRNA candidates over 18 biophysical features,
   inspired by Evo 2 / EVEE-style covariance pooling.
3. In silico specificity analysis against public sequence repositories.
4. A multimodal protein dataset (DNA + protein + ESM-2 + ProtT5 + predicted structure +
   structural descriptors + drug-class labels + CARD/ARO ontology) released on
   HuggingFace under CC-BY-4.0.
5. A reproducible code base under MIT license.
6. A falsifiable, phased validation roadmap.

We make no claim of analytical sensitivity, analytical specificity, clinical performance,
limit of detection, or time-to-result. None of the assays described have been built or
tested in the laboratory.

---

## 2. Methods

All methods described in this section are computational. Outputs are in silico predictions
and design proposals. No biological samples, isolates, or assays were used.

### 2.1 Data sources

- **NCBI RefSeq** (gene and protein reference sequences) [REF: NCBI RefSeq].
- **AMRFinderPlus** reference database (curated AMR gene catalog)
  [REF: Feldgarden 2021 AMRFinderPlus].
- **CARD** (Comprehensive Antibiotic Resistance Database), including ARO ontology and
  drug-class assignments [REF: Alcock 2023 CARD].
- **Public ocular microbiome and ocular-pathogen literature** used to inform target
  prioritization [REF: ocular microbiome reviews; REF: TBD endophthalmitis isolate
  surveys Brazil].

Gene families and variants are documented in `targets_brazil.csv`,
`targets_brazil_card.csv`, and `targets_brazil_variants.csv` in the project repository.

### 2.2 Target prioritization for Brazilian ocular epidemiology

Twelve AMR gene families were selected for design coverage on the basis of (i) prevalence
in Brazilian endophthalmitis and keratitis isolate surveys [REF: TBD BR ocular AMR
epidemiology], (ii) inclusion in BR-GLASS and ANVISA priority pathogen lists
[REF: ANVISA BR-GLASS], and (iii) tractability for CRISPR-Cas12a detection (gene length,
PAM density, sequence conservation). The selected families are: *mecA*, *bla*KPC,
*bla*NDM, *bla*OXA-48, *bla*VIM, *bla*IMP, *bla*GES, *bla*CTX-M, *vanA*, *mcr*, *qnrS*,
and *armA*.

### 2.3 Guide RNA design (`design_guides.py`)

For each reference gene sequence, we scan both strands for canonical Cas12a PAM motifs
(TTTV; V ∈ {A, C, G}) and extract candidate spacers of the optimal length defined in
`config.CAS12A` [VERIFY: spacer length parameter; default ~20 nt]. Each candidate is
scored by `utils.score_guide` using rule-based criteria including GC content, maximum
homopolymer run, self-complementarity, and absence of polythymidine stretches that act
as RNA polymerase III terminators. The Cas12a direct-repeat scaffold is appended
according to `config.CAS12A_DIRECT_REPEAT`. Outputs are ranked TSV tables per gene
(`<gene>_cas12a_guides.tsv`).

### 2.4 Covariance probe re-scoring (`covariance_probes.py`)

To capture non-linear interactions between biophysical features that rule-based scoring
misses, we implement a covariance-pooling re-scoring step inspired by EVEE-style
covariance probes over Evo 2 embeddings [REF: Evo 2; REF: Goodfire/Mayo Clinic EVEE].
Instead of neural embeddings, we use 18 hand-crafted biophysical features per candidate
spacer, including GC content, seed-region thermodynamics (positions 1–8, critical for
Cas12a recognition), homopolymer runs, dinucleotide composition, and
self-complementarity. A feature-covariance matrix captures co-occurrence statistics, and
each candidate is re-ranked by a probe score combining individual feature contributions
and pairwise covariance terms. The procedure is benchmarked against the rule-based
`score_guide` to quantify ranking shift.

### 2.5 RPA primer design (`design_primers.py`)

For the best-ranked crRNA per gene we design flanking RPA primer pairs targeting a
100–200 bp amplicon, with primer length 30–35 nt as recommended for
TwistAmp-style RPA chemistry [REF: Piepenburg 2006 RPA]. Primers are scored on GC
content, basic nearest-neighbor T<sub>m</sub>, homopolymer runs, and self-complementarity
(`design_primers.score_primer`). Parameters are loaded from `config.RPA`. The chemistry
operates near-isothermally at 37 °C, removing the thermocycler requirement and enabling
battery-powered field operation.

### 2.6 In silico specificity (`specificity_check.py`)

We submit full RPA amplicons (100–200 bp) rather than 20-nt spacer fragments to the NCBI
BLAST API, because short queries are unreliable for off-target estimation. Amplicons
were screened against public nt/refseq databases [VERIFY: exact BLAST database used and
date]. Planned extensions, partially implemented in the repository, include explicit
stress-tests against (i) the human reference genome (GRCh38), (ii) ocular surface and
intraocular microbiome references, and (iii) common ocular co-infectants including
*Acanthamoeba* spp., *Fusarium* spp., and *Aspergillus* spp. [VERIFY: which of these
stress tests have been run versus planned].

### 2.7 Variant functional scoring (`evo2_scoring.py`)

We implement two modes: (i) a GPU mode that loads the Evo 2 7B-parameter genomic
language model and extracts intermediate-layer embeddings
(`blocks.28.mlp.l3`) for variant comparison [REF: Nguyen 2025 Evo 2]; (ii) a lightweight
mode that derives sequence-level features (k-mer profile, codon adaptation index,
dinucleotide bias, Shannon-like sequence complexity) to approximate functional impact
when GPUs are unavailable. The functional score quantifies divergence between a candidate
variant and the family reference, going beyond simple Hamming distance.

### 2.8 Structural analysis

Protein sequences are translated from the longest ORF using bacterial codon table 11.
Per-variant embeddings are computed with **ESM-2** (`esm2_t30_150M_UR50D`, 640-dim,
mean-pooled) and **ProtT5** (`Rostlab/prot_t5_xl_uniref50`, 1024-dim, mean-pooled)
[REF: Lin 2023 ESM-2; REF: Elnaggar 2021 ProtT5]. Predicted 3D structures are obtained
via ColabFold/ESMFold, with AlphaFold-Server fallback where available
[REF: Jumper 2021 AlphaFold2; REF: Abramson 2024 AlphaFold3; REF: Mirdita 2022 ColabFold;
REF: Lin 2023 ESMFold]. From each PDB we compute four structural descriptors: number of
Cα atoms (`struct_length`), radius of gyration in Å (`struct_rg`), compactness
(Rg / L<sup>0.6</sup>), Cα–Cα contact density at <8 Å, and mean pLDDT confidence
(`structure_features_v3.py`, `protein_structure.py`).

### 2.9 Family-level protein distance matrices (`protein_distance_matrix.py`)

Within each AMR family we compute pairwise cosine distance between mean-pooled ESM-2
embeddings, producing N×N matrices, followed by simple hierarchical linkage. These
matrices inform candidate crRNA placement: regions of high intra-family sequence
identity but distinct embedding signature are flagged as candidate
"pan-family-conserved" sites.

### 2.10 CARD enrichment (`card_integration.py`)

We pull the latest CARD release, map each variant to its ARO identifier, attach
mechanism (e.g., "antibiotic inactivation"), and assign drug-class labels. New variants
not yet present in `targets_brazil.csv` are surfaced as discovery candidates in
`reports/card_new_variants.csv`.

### 2.11 Multiplex panel optimization (`multiplex_panel.py`)

For the multiplex panel we (i) select the best crRNA and best primer pair per family,
(ii) compute pairwise sequence cross-talk between all primers and between all spacer
sequences, (iii) generate an oligo order sheet and a device layout describing reaction
wells and lateral-flow strip positions, and (iv) produce a human-readable design report.

### 2.12 Software, hardware, reproducibility

The pipeline is implemented in Python ≥3.10. Dependencies are pinned in
`pyproject.toml` and `requirements.txt`. Computations were executed on a local
workstation (macOS, Apple Silicon) for design and scoring steps, and on an NVIDIA DGX
system [VERIFY: exact DGX configuration and node] for Evo 2 GPU mode, ESMFold
structure prediction, and ProtT5 embedding extraction. The full code base is released
under the MIT license at https://github.com/JVLegend/smartsepsis. The published
multimodal dataset is released under CC-BY-4.0 at
https://huggingface.co/datasets/JVLegend/smartsepsis-oph.

[INSERT FIGURE 1 — pipeline overview schematic]

---

## 3. Results (all in silico)

### 3.1 Candidate library

The pipeline produced a candidate crRNA library covering 12 AMR families. The curated
**panel** released on HuggingFace contains 45 variant rows (manifest v1.1.0,
generated 2026-05-05), and the **extended** companion split contains 9,034 rows. Both
splits ship pre-computed ESM-2 and ProtT5 embeddings, predicted PDB structures, and
structural descriptors. Confirmed counts from `fase7_dgx_results/fase7_results/rnafold_guides.csv`:
**213 ranked candidate guides** were produced across **42 AMR target entries** (12 reference
gene families plus clinically relevant variants). Of these, 178 are classified as
`linear_ok` (no detectable secondary structure interfering with the 20-nt spacer), 34 as
`weak_hairpin`, and 1 as `moderate_hairpin`.

[INSERT TABLE 1 — per-family counts: variants, candidate guides surviving PAM scan,
guides retained after covariance probe filter, primer pairs designed]

### 3.2 In silico specificity

Across the designed amplicons submitted to BLAST, no off-target hits exceeded the
configured identity and coverage thresholds against the queried reference databases
[VERIFY: thresholds and database snapshot]. We emphasize that *in silico* specificity
against public repositories is a necessary but **not sufficient** proxy for analytical
specificity in biological samples; the figure should be read as evidence of
on-target design, not of off-target rejection in vitro.

[INSERT FIGURE 2 — distribution of best off-target identity per amplicon, by family]

### 3.3 Functional variant scoring

The lightweight functional score (k-mer divergence, codon adaptation, dinucleotide bias,
sequence complexity) produced family-level distributions consistent with known
sub-family structure (e.g., *bla*KPC-2 vs *bla*KPC-3 cluster more tightly than either
does to *bla*KPC-30). GPU-mode Evo 2 scoring was executed for [VERIFY: which families
ran under Evo 2 GPU mode on DGX].

[INSERT FIGURE 3 — per-family functional score distribution]

### 3.4 Structural clustering

Cosine-distance matrices over ESM-2 mean-pooled embeddings recapitulate expected
family structure: *bla*-class carbapenemases cluster by Ambler class, *mcr* variants
separate cleanly from *vanA*, and *mecA1*/*mecA2* are nearest neighbors. Mean pLDDT
across predicted structures was [VERIFY: report median + IQR from
data/hf_dataset], with [VERIFY: number] variants exhibiting mean pLDDT ≥ 70 (the
threshold conventionally used for "confident" predictions).

[INSERT FIGURE 4 — heatmap of ESM-2 cosine distance, panel of 12 families]

### 3.5 Hamming and panel optimization

Pairwise primer and spacer cross-talk across the 12-family multiplex was below the
configured thresholds [VERIFY: thresholds in `config`], yielding a candidate panel
layout (`multiplex_panel.py` outputs). The generated oligo order sheet enumerates all
spacers (with the Cas12a direct-repeat scaffold prepended) and primer pairs, plus
positive and negative control oligos.

[INSERT FIGURE 5 — proposed device layout, paper-strip readout schematic]

### 3.6 Dataset release and reproducibility

The HuggingFace dataset card documents schema, provenance, model versions, and SHA-256
hashes for both parquet files. The accompanying GitHub repository contains all pipeline
scripts, configuration, and the design artifacts referenced above. A Zenodo DOI will be
attached on submission [REF: TBD Zenodo DOI].

---

## 4. Limitations

We state these limitations explicitly because the work is unusual: a methods paper
without wet-lab results.

- **No experimental validation.** No CRISPR-Cas12a cleavage reactions, no RPA
  amplifications, and no lateral-flow assays have been performed. All numbers in this
  manuscript are computational predictions or design outputs.
- **No analytical performance.** Sensitivity, specificity, positive and negative
  predictive values, limit of detection, dynamic range, and time-to-result are
  **undefined**. They will only be reported following validation under Phase 1 (see
  Section 5).
- **In silico specificity ≠ analytical specificity.** BLAST identity statistics
  characterize sequence-level uniqueness within a chosen reference database, not the
  behavior of a Cas12a–crRNA complex in a complex biological matrix. Off-target trans-cleavage,
  primer–template mispriming, and matrix interference are not captured.
- **Guide-design priors are heterologous.** Our crRNA scoring rules and covariance
  features draw on the broader CRISPR-Dx and Cas12a biochemistry literature
  [REF: Cas12a crRNA design priors], not on ocular-specific empirical data, of which
  little exists.
- **Sample matrix not tested.** Vitreous humor, corneal scrape lysate, and ocular
  surface swab eluate were not used; their effects on RPA and Cas12a chemistry are
  unknown and well-known to be material [REF: RPA matrix interference].
- **Resistance-gene presence ≠ resistance phenotype.** Detection of *mecA* or *bla*KPC
  predicts a probabilistic, not deterministic, phenotype. Genotype-to-phenotype mapping
  for ocular isolates has additional confounders (gene-expression level, regulatory
  mutations, host immune context).
- **Clinical workflow integration not studied.** Sample handling, chain of custody,
  result reporting, integration with intravitreal-injection protocols, and clinician
  decision support are not modeled.
- **Dataset scope.** The HuggingFace dataset (45 panel rows, 9,034 extended rows) is
  small relative to genome-scale AMR catalogs. It is sufficient for benchmarking
  protein-language-model probes within the included families but not for general AMR
  prediction across all organisms.
- **Predicted structures are approximations.** ColabFold/ESMFold/AlphaFold predictions
  have known failure modes for novel folds, disordered regions, and metalloproteins
  (relevant for several β-lactamase families). Mean pLDDT is reported per variant, but
  downstream descriptors should be treated as approximate.

---

## 5. Validation roadmap

We propose a falsifiable, phased plan. Each phase has a single explicit failure
criterion.

- **Phase 1 — Wet-lab feasibility (S. aureus / mecA).** Single-organism, single-gene
  Cas12a + RPA + lateral-flow assay on cultured clinical isolates. Falsification: failure
  to detect *mecA* at ≥10² CFU equivalents per reaction in ≥3 independent runs.
- **Phase 2 — Multiplex panel.** Expand to ≥6 of the 12 designed families on synthetic
  amplicons and cultured isolates. Falsification: pairwise cross-talk above pre-specified
  thresholds, or loss of analytical sensitivity below clinically meaningful bounds.
- **Phase 3 — Clinical sampling.** Prospective ocular sampling under IRB, head-to-head
  against culture and against a reference PCR panel. Pre-registered analytical and
  clinical performance endpoints.
- **Phase 4 — Regulatory submission.** ANVISA RDC 830/2023 IVD pathway in Brazil, with
  equivalent pathways considered for other jurisdictions. Falsification: failure to meet
  pre-registered Phase 3 endpoints.

---

## 6. Data and code availability

- **Code (MIT):** https://github.com/JVLegend/smartsepsis
- **Dataset (CC-BY-4.0):** https://huggingface.co/datasets/JVLegend/smartsepsis-oph
- **Zenodo DOI:** [REF: TBD]
- **Project website:** included under `public/` in the GitHub repository.

The dataset contains two configurations: `panel` (45 rows, curated variants) and
`extended` (9,034 rows, broader companion set), each with DNA and protein sequences,
ESM-2 (640-dim) and ProtT5 (1024-dim) mean-pooled embeddings, predicted PDB structures,
structural descriptors, drug-class labels, and resistance-mechanism annotations linked
to CARD/ARO ontology.

---

## 7. Author contributions

- **JVPD:** project conception, AI pipeline design and implementation, dataset
  curation, computational analyses, manuscript draft.
- **GS:** clinical advisor, ocular epidemiology and ophthalmology methodology review,
  manuscript review.

---

## 8. Conflicts of interest

None declared. JVPD is the founder of *IA para Médicos*, the entity under which this
research is conducted; *IA para Médicos* has no commercial product or pending
regulatory submission related to this work.

---

## 9. Funding

This work is currently unfunded. The authors are actively seeking academic
collaborators, wet-lab partners, and seed funders. Contact:
**iaparamedicos@gmail.com**.

---

## 10. Acknowledgments

We thank the open-source scientific community whose tools (BLAST+, AMRFinderPlus, CARD,
ESM-2, ProtT5, AlphaFold, ColabFold, ESMFold, Evo 2) make a Phase 0 program of this
kind possible. Additional acknowledgments to forthcoming clinical, wet-lab, and
funding partners will be added in subsequent revisions. [Acknowledgments placeholder.]

---

## 11. References

Real foundational references to be expanded with full bibliographic data. Inline
placeholders flagged `[REF: ...]` will be replaced with formatted citations in the
submission build.

1. Gootenberg JS, *et al.* Nucleic acid detection with CRISPR-Cas13a/C2c2. **Science**
   2017. [REF: SHERLOCK 2017]
2. Chen JS, *et al.* CRISPR-Cas12a target binding unleashes indiscriminate
   single-stranded DNase activity. **Science** 2018. [REF: DETECTR 2018]
3. Piepenburg O, *et al.* DNA detection using recombination proteins. **PLoS Biol**
   2006. [REF: RPA 2006]
4. Feldgarden M, *et al.* AMRFinderPlus and the Reference Gene Catalog facilitate
   examination of the genomic links among antimicrobial resistance, stress response,
   and virulence. **Sci Rep** 2021. [REF: AMRFinderPlus]
5. Alcock BP, *et al.* CARD 2023: expanded curation, support for machine learning, and
   resistome prediction at the Comprehensive Antibiotic Resistance Database.
   **Nucleic Acids Res** 2023. [REF: CARD 2023]
6. Jumper J, *et al.* Highly accurate protein structure prediction with AlphaFold.
   **Nature** 2021. [REF: AlphaFold2]
7. Abramson J, *et al.* Accurate structure prediction of biomolecular interactions
   with AlphaFold 3. **Nature** 2024. [REF: AlphaFold3]
8. Lin Z, *et al.* Evolutionary-scale prediction of atomic-level protein structure
   with a language model (ESM-2 / ESMFold). **Science** 2023. [REF: ESM-2 / ESMFold]
9. Elnaggar A, *et al.* ProtTrans: towards cracking the language of life's code through
   self-supervised learning. **IEEE TPAMI** 2021. [REF: ProtT5]
10. Mirdita M, *et al.* ColabFold: making protein folding accessible to all.
    **Nat Methods** 2022. [REF: ColabFold]
11. Nguyen E, *et al.* Sequence modeling and design from molecular to genome scale with
    Evo 2. **Arc Institute preprint** 2025. [REF: Evo 2]
12. Durand ML. Bacterial and fungal endophthalmitis. **Clin Microbiol Rev** 2017.
    [REF: Durand endophthalmitis]
13. Ung L, *et al.* The persistent dilemma of microbial keratitis: global burden,
    diagnosis, and antimicrobial resistance. **Surv Ophthalmol** 2019. [REF: Ung keratitis]
14. Endophthalmitis Vitrectomy Study Group. **Arch Ophthalmol** 1995. [REF: EVS]
15. ANVISA. Resolução RDC 830/2023 — In vitro diagnostic devices. [REF: ANVISA RDC 830/2023]
16. ANVISA / BR-GLASS national AMR surveillance reports. [REF: BR-GLASS]

**Placeholders to be filled with real references:**
- [REF: Goodfire/Mayo Clinic EVEE covariance probes]
- [REF: BioFire panel composition]
- [REF: Cas12a crRNA design priors]
- [REF: vitreous tap volume]
- [REF: corneal scrape yield]
- [REF: ocular microbiome reviews]
- [REF: TBD BR ocular AMR epidemiology]
- [REF: TBD endophthalmitis isolate surveys Brazil]
- [REF: RPA matrix interference]
- [REF: TBD Zenodo DOI]

---

*Manuscript version: draft v0 — prepared for internal review prior to bioRxiv submission.*
