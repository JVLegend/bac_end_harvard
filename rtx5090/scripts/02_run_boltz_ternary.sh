#!/usr/bin/env bash
# Boltz-1 ternary complex prediction: LbCas12a + crRNA + dsDNA target
# for the top 3 mecA guides. GPU bound (~20 GB VRAM peak).
set -euo pipefail
cd "$(dirname "$0")/.."
source .venv/bin/activate

OUT=outputs/boltz
mkdir -p "$OUT"
LOG="$OUT/boltz.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== Boltz ternary start: $(date -u +%FT%TZ) ==="

# Verify boltz installed
if ! python -c "import boltz" 2>/dev/null; then
  echo "Installing boltz..."
  pip install boltz==1.* 2>&1 | tail -5
fi

# Read Cas12a sequence (single-line)
CAS12A=$(awk '!/^>/' inputs/cas12a_lba_protein.fasta | tr -d '\n')

# Build Boltz YAML input for each of 3 mecA guides
build_yaml() {
  local guide_id="$1"
  local crrna_seq="$2"
  local target_top="$3"
  local target_bot="$4"
  local out_yaml="$OUT/${guide_id}/input.yaml"
  mkdir -p "$OUT/${guide_id}"
  cat > "$out_yaml" <<EOF
version: 1
sequences:
  - protein:
      id: A
      sequence: "${CAS12A}"
  - rna:
      id: B
      sequence: "${crrna_seq}"
  - dna:
      id: C
      sequence: "${target_top}"
  - dna:
      id: D
      sequence: "${target_bot}"
EOF
  echo "  wrote $out_yaml"
}

# Extract crRNA + target sequences for each guide
declare -A CRRNAS
while read -r line; do
  if [[ $line == ">"* ]]; then
    cur="${line#>crRNA_}"; cur="${cur%% *}"
  else
    CRRNAS[$cur]="$line"
  fi
done < inputs/top5_mecA_crRNAs.fasta

for gid in mecA-0 mecA-1 mecA-2; do
  TOP=$(awk -v g="target_top_${gid}" '$0==">"g{f=1;next} /^>/{f=0} f' "inputs/target_dsDNA_${gid}.fasta")
  BOT=$(awk -v g="target_bot_${gid}" '$0==">"g{f=1;next} /^>/{f=0} f' "inputs/target_dsDNA_${gid}.fasta")
  build_yaml "$gid" "${CRRNAS[$gid]}" "$TOP" "$BOT"
done

# Run Boltz on each YAML
for gid in mecA-0 mecA-1 mecA-2; do
  echo ""
  echo "[Boltz] predicting $gid ..."
  boltz predict "$OUT/${gid}/input.yaml" --out_dir "$OUT/${gid}" --use_msa_server 2>&1 | tail -20 || {
    echo "  FAILED for $gid" > "$OUT/${gid}/FAILED.txt"
    continue
  }
done

# Summarize
python3 - <<'PY'
import json, os, glob, csv
out_dir = "outputs/boltz"
rows = [["guide_id","mean_plddt","ipTM","status"]]
for gid in ["mecA-0","mecA-1","mecA-2"]:
    pred = glob.glob(f"{out_dir}/{gid}/boltz_results*/predictions/*/confidence*.json")
    if not pred:
        rows.append([gid,"","","NO_OUTPUT"]); continue
    with open(pred[0]) as f: d = json.load(f)
    rows.append([gid, d.get("complex_plddt",""), d.get("iptm",""), "OK"])
with open(f"{out_dir}/boltz_summary.tsv","w") as f:
    csv.writer(f, delimiter="\t").writerows(rows)
print(open(f"{out_dir}/boltz_summary.tsv").read())
PY

echo "=== Boltz done: $(date -u +%FT%TZ) ==="
