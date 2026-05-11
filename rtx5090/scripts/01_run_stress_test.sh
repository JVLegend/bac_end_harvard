#!/usr/bin/env bash
# Stress-test specificity: BLAST guides against GRCh38, ocular commensals,
# and Acanthamoeba/fungal pathogens. CPU + disk heavy. No GPU needed.
set -euo pipefail
cd "$(dirname "$0")/.."
source .venv/bin/activate

OUT=outputs/stress_test
mkdir -p "$OUT" blast_dbs
LOG="$OUT/stress_test.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== Stress-test start: $(date -u +%FT%TZ) ==="

# 1. Build BLAST databases (idempotent: skips if already built)
cd blast_dbs

if [ ! -f human_genome.nhr ]; then
  echo "[1/3] Downloading GRCh38 (~900 MB compressed, ~3 GB uncompressed)..."
  curl -L -o human.fna.gz \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
  gunzip human.fna.gz
  makeblastdb -in human.fna -dbtype nucl -out human_genome -title "GRCh38_p14"
  rm -f human.fna   # free ~3 GB now that index exists
fi

if [ ! -f ocular_microbiome.nhr ]; then
  echo "[2/3] Downloading ocular commensals..."
  : > ocular_micro.fna
  for ACC in NC_004461.1 NC_006085.1 NZ_CP013128.1; do
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=text" >> ocular_micro.fna
    sleep 1
  done
  makeblastdb -in ocular_micro.fna -dbtype nucl -out ocular_microbiome -title "ocular_commensals"
fi

if [ ! -f ocular_fungi_acanthamoeba.nhr ]; then
  echo "[3/3] Downloading Acanthamoeba + ocular fungi..."
  : > fungi_acanth.fna
  for ACC in NW_004457519.1 NC_007194.1 NC_032089.1 NC_037070.1; do
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=text" >> fungi_acanth.fna
    sleep 1
  done
  makeblastdb -in fungi_acanth.fna -dbtype nucl -out ocular_fungi_acanthamoeba -title "ocular_eukaryotes"
fi
cd ..

# 2. Run the stress-test python script (lives in repo's scripts/, two levels up)
REPO_ROOT="$(git rev-parse --show-toplevel)"
GUIDES_CSV="inputs/all_guides.csv"

echo ""
echo "Running stress_test_specificity.py..."
python3 "$REPO_ROOT/scripts/stress_test_specificity.py" \
  --guides "$GUIDES_CSV" \
  --blast-dir blast_dbs \
  --output-dir "$OUT"

echo "=== Stress-test done: $(date -u +%FT%TZ) ==="
ls -la "$OUT"
