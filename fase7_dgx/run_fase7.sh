#!/bin/bash
# Fase 7 BacEnd — Orchestrator na DGX Spark
# Tasks: ESM-2 embeddings, RNAfold, Phylogenia+Tajima, Novel variants NCBI
set -u

PROJECT="/home/oftalmousp/jv-teste/harvard_bacend"
FASE7_DIR="$PROJECT/fase7_dgx"
RESULTS="$PROJECT/fase7_results"
LOGS="$PROJECT/fase7_logs"
TIMEOUT=1800

mkdir -p "$RESULTS" "$LOGS"
cd "$FASE7_DIR" || { echo "ERR: $FASE7_DIR nao existe"; exit 1; }

# Venv (reaproveita do MIQA)
VENV="/home/oftalmousp/jv-teste/miqa_backend/venv"
if [ -f "$VENV/bin/activate" ]; then
    source "$VENV/bin/activate"
    PY="$VENV/bin/python"
else
    PY=python3
fi

# Deps
$PY -m pip install -q biopython pandas numpy torch transformers 2>&1 | tail -3
$PY -m pip install -q ViennaRNA 2>&1 | tail -2 || echo "(ViennaRNA Python pkg falhou — fallback ativo)"

export FASE7_DIR="$RESULTS"

SUMMARY="$LOGS/fase7_summary.csv"
[ ! -f "$SUMMARY" ] && echo "task,status,seconds,result_line" > "$SUMMARY"

run_task() {
    local name="$1"; local script="$2"
    local log="$LOGS/${name}.log"
    echo ""
    echo "################################################################"
    echo "# $name — $(date +%H:%M:%S)"
    echo "################################################################"
    local t0=$(date +%s)
    if timeout "$TIMEOUT" "$PY" "$script" > "$log" 2>&1; then
        local dt=$(( $(date +%s) - t0 ))
        local rl=$(grep '^RESULT:' "$log" | head -3 | tr '\n' ';' | sed 's/,/;/g' | sed 's/"//g')
        echo "OK  $name (${dt}s) — $rl"
        echo "$name,ok,$dt,\"$rl\"" >> "$SUMMARY"
        tail -5 "$log"
    else
        local dt=$(( $(date +%s) - t0 ))
        echo "FAIL $name (${dt}s) — last 20 lines:"
        tail -20 "$log"
        echo "$name,fail,$dt," >> "$SUMMARY"
    fi
}

echo "========================================"
echo "Fase 7 BacEnd — DGX run  $(date +%Y-%m-%d\ %H:%M)"
echo "python: $PY"
echo "results: $RESULTS"
echo "cuda: $($PY -c 'import torch; print(torch.cuda.is_available(), torch.cuda.device_count())' 2>/dev/null)"
echo "========================================"

run_task "f7_02_rnafold"         task_f7_02_rnafold.py
run_task "f7_03_phylo_tajima"    task_f7_03_phylogeny_tajima.py
run_task "f7_01_esm2"            task_f7_01_esm2_embeddings.py
run_task "f7_04_novel_variants"  task_f7_04_novel_variants.py

echo ""
echo "========================================"
echo "FASE 7 DONE — summary: $SUMMARY"
cat "$SUMMARY"
echo "Files produced:"
ls "$RESULTS" | sort
