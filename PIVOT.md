# Program pivot — 2026-05-11

## Summary

SmartSepsis-Oph pivoted from "AI-first CRISPR paper-strip diagnostic" to
**"AI-augmented multimodal diagnostics for ophthalmic infection"**, with
Nanopore metagenomics as the primary track and CRISPR-Cas12a as an
adjunct module for resistance triage.

## What drove the pivot

Dr. Gustavo Sakuno (Mass Eye and Ear / Harvard, co-author) reviewed the
project and surfaced two correct scientific objections:

1. **Wrong clinical question.** Vitreous is an immune sanctuary —
   normally sterile. The question is *which organism is here?*, not *is
   this organism resistant?*. A CRISPR resistance panel presupposes an
   organism identity that has not been established.
2. **Kitome problem in low-biomass samples.** Pan-bacterial 16S
   amplification on a paper strip is dominated by residual bacterial DNA
   in reagents (well-documented in low-biomass microbiome literature —
   Salter 2014; Doan 2016). Signal-to-noise on sub-microliter sterile-by-
   default samples is unfavorable for paper-strip CRISPR as primary
   identification.

Subsequent conversation revealed that Dr. Sakuno's senior advisor at
Mass Eye and Ear / Harvard has expertise in **Nanopore metagenomic
diagnostics for ocular infection** and is willing to participate as
senior scientific advisor and host the wet-lab work in her laboratory.

## What the program now is

Two complementary tracks:

### Track A (primary) — AI-augmented Nanopore metagenomic interpretation
- Take Nanopore reads from ocular samples
- Apply foundation-model embeddings (ESM-2, ProtT5) and AlphaFold-derived
  structural features to interpret noisy, low-coverage AMR reads
- Output: taxonomic identification + AMR gene catalog + predicted drug
  class + structural annotation
- Wet-lab partner: senior advisor's laboratory at Mass Eye and Ear

### Track B (adjunct) — CRISPR-Cas12a resistance triage
- Paper-strip Cas12a + RPA for rapid resistance-gene detection
- Conditional on suspected organism (not primary ID)
- Use cases: time-critical empirical guidance in OR; low-resource
  settings without Nanopore capacity
- Inherits the 213 candidate guides across 12 resistance-gene families
  produced in earlier Phase 0 work

## What stays vs what changes

### Stays (preserved value)
- HuggingFace multimodal dataset (`JVLegend/smartsepsis-oph`) — becomes
  the backbone of Track A AMR interpretation, more central than before.
- AlphaFold structures for 43 AMR variants — more useful for Track A
  (structural annotation of fragmented Nanopore reads) than for Track B.
- 213 designed CRISPR guides — become Track B adjunct module.
- AlphaFold 3 ternary complex submission pack (`alphafold3_submission/`)
  — still relevant for Track B structural feasibility.
- RTX 5090 compute pack (`rtx5090/`) — still useful; Evo 2 embeddings
  apply to AMR characterization regardless of track.
- BLAST stress-test scripts — still useful for Track B specificity.
- Honest Phase 0 framing throughout.

### Changes (primary scientific direction)
- Primary technology: paper-strip CRISPR → Nanopore metagenomics.
- Primary diagnostic question: resistance detection → pathogen ID +
  AMR interpretation.
- Computational pipeline becomes interpretation-centered (reads →
  taxonomy + AMR), not design-centered (gene → guide RNA library).
- Dr. Sakuno repositioned from "Clinical Advisor" to **co-author**.
- Added "Senior Scientific Advisor" slot (renowned investigator at Mass
  Eye and Ear; name to be announced after formal confirmation).
- Added "Partner Laboratory" slot at Mass Eye and Ear.
- First publication output reordered: **scoping review (Paper 01) →
  before** the framework preprint.

## Files touched by this pivot

- `README.md` — full rewrite with two-track framing
- `CITATION.cff` — Sakuno as proper co-author, advisor placeholder, new
  abstract and keywords, version 0.2.0-phase0
- `public/index.html` — hero, stats, dor.solution, dor.why, diff section,
  team (Sakuno co-author + Senior Advisor placeholder + Partner Lab
  placeholder), footer
- `public/lang.js` — same keys, PT and EN
- `public/paper.html` — Paper 01 reframed as scoping review (P0 in
  drafting); TOC entry updated
- `preprint/draft_v0.md` — header note marking pivot; framework preprint
  deferred until scoping review lands
- `CONTRIBUTING.md` — what we still seek updated

## Files NOT touched and why

- All Python pipeline scripts at repo root — code still valid; some
  scripts become Track A backbone, others Track B; no rename needed yet.
- `alphafold3_submission/` — still useful as Track B structural pack.
- `rtx5090/` — still useful; Evo 2 + stress-test jobs apply regardless.
- `scripts/stress_test_specificity.py` — still useful for Track B.
- `docs/` MkDocs site — content needs revision in next pass; not blocked.

## Next steps after the pivot

1. Finalize and submit scoping review (already scaffolded under
   `preprint/scoping_review_*.md`).
2. Formal confirmation of Senior Scientific Advisor and partner
   laboratory → public announcement → update `CITATION.cff`, `README.md`,
   `public/index.html` team section, and `preprint/` author lists.
3. Rewrite `preprint/draft_v0.md` as a Track A + Track B framework paper
   after scoping review provides the rigorous prior-art baseline.
4. Update `docs/` MkDocs pages with two-track framing.
5. Plan Track A pipeline experiments (Nanopore reads → AMR
   interpretation benchmarks) with partner laboratory.

## Honesty rules carried over

- No clinical performance claims for either track until experimental
  validation.
- Senior Advisor name not published before formal confirmation.
- Track B guide library remains "computationally designed, experimentally
  unvalidated".
- This pivot is documented openly here rather than hidden behind a
  rewrite — anyone reading the repo can trace the program's evolution.
