#!/usr/bin/env bash
# Evo 2 7B embeddings + cosine-distance matrix for AMR reference sequences.
# Replaces the k-mer proxy currently in evo2_scoring.py with real embeddings.
# GPU bound (~14 GB VRAM with 4-bit quant; ~28 GB FP16).
set -euo pipefail
cd "$(dirname "$0")/.."
source .venv/bin/activate

OUT=outputs/evo2
mkdir -p "$OUT"
LOG="$OUT/evo2.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== Evo 2 start: $(date -u +%FT%TZ) ==="

# Try install evo2 if missing
if ! python -c "import evo2" 2>/dev/null; then
  echo "Installing Evo 2 from PyPI..."
  pip install evo2 2>&1 | tail -5 || {
    echo "PyPI install failed; trying from source..."
    pip install git+https://github.com/ArcInstitute/evo2.git 2>&1 | tail -10 || {
      echo "Evo 2 install failed. Writing FAILED.txt." > "$OUT/FAILED.txt"
      echo "Install steps: see https://github.com/ArcInstitute/evo2" >> "$OUT/FAILED.txt"
      exit 1
    }
  }
fi

REPO_ROOT="$(git rev-parse --show-toplevel)"

python3 <<'PY'
import os, json, csv, numpy as np, pandas as pd
import torch
from pathlib import Path

OUT = Path("outputs/evo2")
OUT.mkdir(parents=True, exist_ok=True)

# Load AMR reference sequences from the HF dataset (panel split, 45 rows)
df = pd.read_parquet("../data/hf_dataset/smartsepsis_oph.parquet")
print(f"Loaded {len(df)} reference variants")
seqs = df[["variant_id","gene_family","dna_sequence"]].dropna()
print(f"With DNA sequence: {len(seqs)}")

# Load Evo 2
print("Loading Evo 2 7B (4-bit if low VRAM)...")
free_gb = torch.cuda.mem_get_info()[0] / 1e9
quant = free_gb < 24
print(f"  Free VRAM: {free_gb:.1f} GB  | quantized: {quant}")
try:
    from evo2 import Evo2
    model = Evo2(model_name="evo2_7b") if not quant else Evo2(model_name="evo2_7b", load_in_4bit=True)
except Exception as e:
    print("Evo 2 load failed:", e)
    (OUT / "FAILED.txt").write_text(f"Evo 2 load: {e}")
    raise

# Embed each sequence (mean-pool last hidden state)
embeddings = {}
for _, row in seqs.iterrows():
    vid = row["variant_id"]; sq = str(row["dna_sequence"])[:8192]  # truncate to ctx
    try:
        emb = model.embed(sq, layer="last", pool="mean")
        embeddings[vid] = emb.cpu().numpy().astype("float32")
        print(f"  {vid}: {embeddings[vid].shape}")
    except Exception as e:
        print(f"  {vid}: FAILED ({e})")

np.savez(OUT / "embeddings.npz", **embeddings)
print(f"Saved {len(embeddings)} embeddings -> {OUT}/embeddings.npz")

# Cosine distance matrix per gene_family vs family reference
fam_groups = seqs.groupby("gene_family")
rows = []
for fam, grp in fam_groups:
    ref_id = grp.iloc[0]["variant_id"]   # first variant as family reference
    if ref_id not in embeddings: continue
    ref = embeddings[ref_id]
    ref = ref / (np.linalg.norm(ref) + 1e-9)
    for _, r in grp.iterrows():
        vid = r["variant_id"]
        if vid not in embeddings: continue
        v = embeddings[vid]; v = v / (np.linalg.norm(v) + 1e-9)
        cos = float(np.dot(ref, v))
        rows.append({"gene_family": fam, "variant_id": vid, "reference_variant": ref_id,
                     "cosine_similarity": cos, "cosine_distance": 1.0 - cos})

pd.DataFrame(rows).to_csv(OUT / "distance_matrix.csv", index=False)
print(f"Wrote distance_matrix.csv ({len(rows)} rows)")
PY

echo "=== Evo 2 done: $(date -u +%FT%TZ) ==="
ls -la "$OUT"
