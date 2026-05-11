#!/usr/bin/env bash
# SmartSepsis-Oph rtx5090 — environment preflight
set -uo pipefail
echo "=== SmartSepsis-Oph RTX 5090 preflight ==="
FAIL=0

check() { local label="$1"; shift; if "$@" >/dev/null 2>&1; then echo "  [ OK ] $label"; else echo "  [FAIL] $label"; FAIL=$((FAIL+1)); fi; }

echo "[CUDA / GPU]"
if command -v nvidia-smi >/dev/null; then
  nvidia-smi --query-gpu=name,memory.total,memory.free,driver_version --format=csv,noheader
  FREE_MB=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  if [ "${FREE_MB:-0}" -lt 28000 ]; then echo "  [WARN] GPU free <28 GB ($FREE_MB MB)"; fi
else echo "  [FAIL] nvidia-smi not found"; FAIL=$((FAIL+1)); fi

echo "[Disk]"
df -h . | tail -1
AVAIL_GB=$(df -BG . | tail -1 | awk '{print $4}' | tr -d 'G')
if [ "${AVAIL_GB:-0}" -lt 25 ]; then echo "  [FAIL] need >=25 GB free, have ${AVAIL_GB} GB"; FAIL=$((FAIL+1)); fi

echo "[Tools]"
check "python3.11"  python3.11 --version
check "git"         git --version
check "curl"        curl --version
check "wget"        wget --version
check "tar"         tar --version
check "blastn"      blastn -version
check "makeblastdb" makeblastdb -version

echo "[Python venv]"
if [ ! -d .venv ]; then
  echo "  Creating .venv..."
  python3.11 -m venv .venv && echo "  [ OK ] .venv created"
fi
# shellcheck source=/dev/null
source .venv/bin/activate
pip install --quiet --upgrade pip
echo "[Installing Python deps]"
pip install --quiet -r requirements.txt 2>&1 | tail -5

echo ""
if [ "$FAIL" -gt 0 ]; then
  echo "=== PREFLIGHT FAILED ($FAIL issues) — fix and re-run ==="
  exit 1
fi
echo "=== PREFLIGHT OK — ready to run scripts/01..03 ==="
