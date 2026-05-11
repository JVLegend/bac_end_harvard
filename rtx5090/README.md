# SmartSepsis-Oph — RTX 5090 Compute Pack

**Audience:** an LLM coding agent (Claude / Codex / Aider) operating on a
Linux box with **NVIDIA RTX 5090 (32 GB VRAM)**, working directory = this
folder.

**Goal:** run 3 GPU-bound workloads that the user's Mac cannot run locally,
then pack the outputs so they can be pulled back to the main repo.

---

## What you (the agent) need to do

Execute these 3 jobs **in order**. Each step has a `scripts/NN_*.sh` driver
that does the actual work. **Read the script before running**; if anything
looks wrong for this machine, stop and ask the human.

| # | Job | Driver | Approx time | Output |
|---|---|---|---|---|
| 0 | Verify environment (CUDA, disk, deps) | `scripts/00_check_environment.sh` | 1 min | console only |
| 1 | BLAST specificity stress-test (GRCh38 + ocular commensals + Acanthamoeba/fungi) | `scripts/01_run_stress_test.sh` | 2-6 h (download + BLAST) | `outputs/stress_test/` |
| 2 | Boltz-1 ternary complex prediction (LbCas12a + crRNA + dsDNA), 3 top mecA guides | `scripts/02_run_boltz_ternary.sh` | 30-90 min | `outputs/boltz/` |
| 3 | Evo 2 7B embeddings for 12 AMR gene-family reference DNA sequences | `scripts/03_run_evo2_embeddings.sh` | 30-90 min | `outputs/evo2/` |
| ⨠ | Pack everything for upload | `scripts/99_pack_results.sh` | 1 min | `outputs/smartsepsis_5090_<timestamp>.tar.gz` |

---

## Pre-flight checks the agent should make

Before kicking off the heavy jobs, confirm:

1. `nvidia-smi` shows an RTX 5090 with **at least 30 GB VRAM free**.
2. `df -h .` shows **at least 25 GB free** on the disk holding this folder
   (GRCh38 + Boltz model weights + Evo 2 weights).
3. Python 3.11 available (`python3.11 --version`).
4. `git`, `curl`, `wget`, `tar` available.
5. CUDA drivers compatible with PyTorch 2.4+ (Boltz and Evo 2 dependency).
6. Internet access for first run (downloads weights + databases).

If any check fails, abort and tell the human exactly what to fix.

---

## How to run

```bash
# from this directory (rtx5090/)
bash scripts/00_check_environment.sh
bash scripts/01_run_stress_test.sh
bash scripts/02_run_boltz_ternary.sh
bash scripts/03_run_evo2_embeddings.sh
bash scripts/99_pack_results.sh
```

You can run 1, 2, 3 in parallel if VRAM/CPU permits (`tmux` recommended).
Stress-test is mostly CPU + disk; Boltz and Evo 2 are GPU. Boltz and Evo 2
should NOT run at the same time on a single GPU (each peaks ~20 GB).

---

## Inputs (already in this folder)

`inputs/` contains everything the scripts need:

- `cas12a_lba_protein.fasta` — LbCas12a, 1228 aa
- `top5_mecA_crRNAs.fasta` — 5 crRNAs (scaffold + 20-nt spacer)
- `target_dsDNA_mecA-{0,1,2}.fasta` — top + bottom strands of target DNA
- `targets_SUMMARY.tsv` — PAM/position metadata
- `all_guides.csv` — full 213-guide library (input for stress-test)
- `targets_brazil_variants.csv` — 12 gene-family reference catalog

---

## Outputs that need to come back

After `99_pack_results.sh`, upload `outputs/smartsepsis_5090_<timestamp>.tar.gz`
back to the user (Google Drive, scp, GitHub release, whatever they prefer).

Inside the tarball, the user expects:

```
outputs/
├── stress_test/
│   ├── stress_test_<utc>.tsv          # raw hit-level results
│   ├── stress_test_summary.tsv        # per-guide pass/fail per DB
│   └── stress_test.log
├── boltz/
│   ├── <guide_id>/
│   │   ├── boltz_results_*/           # Boltz native output
│   │   ├── plddt.json
│   │   └── render.png                 # optional, if pymol/chimerax available
│   └── boltz_summary.tsv              # mean pLDDT + ipTM per complex
├── evo2/
│   ├── embeddings.npz                 # variant_id -> vector
│   ├── distance_matrix.csv            # cosine distance vs family reference
│   └── evo2.log
└── MANIFEST.json                      # describes everything
```

The user will integrate these into the preprint (`preprint/draft_v0.md`) as
Figures 2, 3, and a refreshed Section 3.3.

---

## Failure modes and what to do

- **Boltz install fails on CUDA mismatch** — try `pip install boltz==1.*`
  and pin torch to the CUDA version that matches `nvidia-smi`. If still
  failing, skip Boltz and run only stress-test + Evo 2.
- **Evo 2 model download fails** — Evo 2 weights are on HuggingFace
  (`arcinstitute/evo2_7b` or similar). Make sure `HF_TOKEN` env var is set
  if the model requires gated access.
- **BLAST OOM** — split the guide CSV into chunks of 50 and run sequentially.
- **Disk fills up** — delete the uncompressed GRCh38 FASTA once `makeblastdb`
  finishes; only the `.nhr/.nin/.nsq` index files are needed.

If something is genuinely broken, **do not invent numbers or fake outputs**.
Leave a `FAILED.txt` in the relevant `outputs/<job>/` subfolder explaining
what blew up. The user prefers honest failure to fake results.

---

## Honesty rules (carry over to all artifacts)

Everything produced here is **in silico**. Do not write log lines or
README snippets that imply experimental validation. Section 3 of the
preprint (`preprint/draft_v0.md`) makes this explicit; preserve that tone.
