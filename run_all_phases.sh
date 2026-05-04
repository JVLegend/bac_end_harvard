#!/usr/bin/env bash
# ====================================================================
# SmartSepsis — Pipeline Fases 4-7 com checkpoint/resume
#
# Uso:    bash run_all_phases.sh
# Resume: rode o mesmo comando — etapas concluidas sao puladas.
# Reset:  rm -rf .run_state    (forca tudo de novo)
#
# Marca cada fase como concluida em .run_state/<fase>.done.
# Scripts internos ja cacheiam embeddings (.npy/.npz) e PDBs, entao
# mesmo dentro de uma fase parcialmente rodada, so o que faltar processa.
# ====================================================================

set -e
set -o pipefail
cd "$(dirname "$0")"

PY="/Users/iaparamedicos/envs/dev/bin/python"
LOG="run_all.log"
STATE_DIR=".run_state"
mkdir -p "$STATE_DIR"

banner() {
    local msg="$1"
    echo ""
    echo "============================================================" | tee -a "$LOG"
    echo " $msg  [$(date +%H:%M:%S)]" | tee -a "$LOG"
    echo "============================================================" | tee -a "$LOG"
}

mark_done() { touch "$STATE_DIR/$1.done"; }
is_done()   { [ -f "$STATE_DIR/$1.done" ]; }

run_phase() {
    local key="$1"; shift
    local label="$1"; shift
    if is_done "$key"; then
        echo "[SKIP] $label  (ja concluido em $STATE_DIR/$key.done)" | tee -a "$LOG"
        return 0
    fi
    banner "$label"
    if "$@" >> "$LOG" 2>&1; then
        mark_done "$key"
        echo "[OK] $label" | tee -a "$LOG"
    else
        local rc=$?
        echo "[FAIL rc=$rc] $label — veja $LOG, rode novamente para retomar" | tee -a "$LOG"
        exit 1
    fi
}

# ====================================================================
# Pre-requisitos (idempotente)
# ====================================================================
if ! command -v RNAfold &> /dev/null; then
    banner "Instalando ViennaRNA"
    brew install viennarna 2>&1 | tee -a "$LOG"
fi

export PATH="/Users/iaparamedicos/envs/dev/bin:$PATH"
if ! command -v prodigal &> /dev/null; then
    echo "[INFO] prodigal nao instalado — Fase 7 sera pulada."
    echo "       Instale: brew install prodigal"
fi

# ====================================================================
# FASE 4 — ESMFold (3D structures)
# ====================================================================
run_phase "fase4_esmfold" "FASE 4 — ESMFold (3D structures, ate 400aa)" \
    bash -c "ESMFOLD_MAX_LEN=400 PROTEINS=blaKPC-3,blaNDM-1,mecA1,vanA,blaCTX-M-2,blaCTX-M-14,blaKPC-30 \"$PY\" protein_structure.py"

# Copiar PDBs pro public/ (idempotente)
mkdir -p public/pdbs
if ls reports/protein_structure/pdbs/*.pdb 1> /dev/null 2>&1; then
    cp -u reports/protein_structure/pdbs/*.pdb public/pdbs/ 2>/dev/null || true
fi

# ====================================================================
# FASE 5a — AMRFinderPlus (CATALOGO COMPLETO ~10k proteinas)
# ====================================================================
run_phase "fase5a_amrfinder_full" "FASE 5a — AMRFinderPlus catalogo COMPLETO (~10k, 30-40min)" \
    "$PY" amrfinderplus_embed.py

# ====================================================================
# FASE 5b — ProtT5 ensemble
# ====================================================================
run_phase "fase5b_prott5" "FASE 5b — ProtT5 ensemble" \
    "$PY" prott5_ensemble.py

# ====================================================================
# FASE 6 — ViennaRNA crRNA structure
# ====================================================================
if ls guides/*.tsv 1> /dev/null 2>&1; then
    run_phase "fase6_viennarna" "FASE 6 — Estrutura secundaria crRNA (ViennaRNA)" \
        "$PY" crrna_secondary_structure.py
else
    echo "[SKIP] Fase 6 — guides/*.tsv nao encontrados (rode design_guides.py antes)" | tee -a "$LOG"
fi

# ====================================================================
# FASE 7 — Pangenome bacteriano (HGT detection)
# ====================================================================
if command -v prodigal &> /dev/null; then
    run_phase "fase7_pangenome" "FASE 7 — Pangenome (download + prodigal + panaroo)" \
        bash pangenome.sh
else
    echo "[SKIP] Fase 7 — prodigal ausente." | tee -a "$LOG"
fi

# ====================================================================
# Deploy
# ====================================================================
if command -v vercel &> /dev/null; then
    banner "Deploy Vercel"
    vercel --prod --yes 2>&1 | tail -5 | tee -a "$LOG"
fi

# ====================================================================
banner "FEITO — checkpoints em $STATE_DIR/"
echo "Resultados:"
echo "  - reports/protein_structure/pdbs/*.pdb     (Fase 4 — 3D)"
echo "  - reports/amrfinderplus/                    (Fase 5a — classificador)"
echo "  - reports/prott5_ensemble/metrics.json      (Fase 5b — ensemble)"
echo "  - reports/crrna_structure/                  (Fase 6 — RNA fold)"
echo "  - reports/pangenome/                        (Fase 7 — HGT/core)"
echo ""
echo "Site: https://smartsepsis.vercel.app/structure.html"
echo ""
echo "Para rerodar uma fase especifica: rm $STATE_DIR/<fase>.done"
