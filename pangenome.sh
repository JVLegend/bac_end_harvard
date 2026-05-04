#!/usr/bin/env bash
# ====================================================================
# Fase 7 — Pangenome bacteriano (HGT detection)
#
# Pipeline:
#   1. Baixa ~50 genomas Klebsiella pneumoniae + E. coli do RefSeq
#   2. Anota cada um com prokka (gera .gff)
#   3. Roda panaroo no conjunto (core/accessory genome, HGT signal)
#
# Resume: cada genome eh skip se ja anotado. Panaroo eh skip se outdir existe.
# ====================================================================

set -e
export PATH="/Users/iaparamedicos/envs/dev/bin:$PATH"
cd "$(dirname "$0")"

PY="/Users/iaparamedicos/envs/dev/bin/python"
REPO_DIR="$(pwd)"
PG_DIR="data/pangenome"
GENOMES_DIR="$PG_DIR/genomes"
PROKKA_DIR="$PG_DIR/prokka"
PANAROO_OUT="reports/pangenome/panaroo_out"

mkdir -p "$GENOMES_DIR" "$PROKKA_DIR" "reports/pangenome"

# === Lista de assemblies (RefSeq accessions, prevalentes em UTI BR) ===
# 25 K. pneumoniae + 25 E. coli (mix de ST/serotipo, com AMR conhecido)
ACCESSIONS=(
    # K. pneumoniae representativos
    GCF_000240185.1 GCF_001596075.1 GCF_002240885.1 GCF_002902925.1 GCF_003014415.1
    GCF_003812245.1 GCF_004035915.1 GCF_004151375.1 GCF_004151875.1 GCF_004152545.1
    GCF_009025895.1 GCF_010120655.1 GCF_011045645.1 GCF_012226325.1 GCF_013396265.1
    GCF_013798765.1 GCF_014193465.1 GCF_015097635.1 GCF_017654635.1 GCF_018399305.1
    GCF_019388325.1 GCF_020874975.1 GCF_022377535.1 GCF_023719555.1 GCF_024740105.1
    # E. coli representativos
    GCF_000005845.2 GCF_000008865.2 GCF_000010385.1 GCF_000010745.1 GCF_000017745.1
    GCF_000019425.1 GCF_000026325.2 GCF_000026545.1 GCF_001941055.1 GCF_002007705.1
    GCF_002896145.1 GCF_003952565.1 GCF_004337995.1 GCF_004803135.1 GCF_005886055.1
    GCF_007990505.1 GCF_009832145.1 GCF_011880525.1 GCF_013030255.1 GCF_014167395.1
    GCF_015852495.1 GCF_017639175.1 GCF_019836925.1 GCF_021391275.1 GCF_023650915.1
)

# ====================================================================
# Etapa 1 — Download (resume-friendly)
# ====================================================================
echo "[Pangenome] Etapa 1/3 — download de ${#ACCESSIONS[@]} genomas"
total=${#ACCESSIONS[@]}
i=0
for acc in "${ACCESSIONS[@]}"; do
    i=$((i+1))
    fa="$GENOMES_DIR/${acc}.fna"
    if [ -f "$fa" ] && [ -s "$fa" ]; then
        echo "  [$i/$total] $acc  (cached)"
        continue
    fi

    # Resolver caminho FTP via assembly_summary (ou tentar ambos os bancos)
    base="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${acc:4:3}/${acc:7:3}/${acc:10:3}"
    # Tenta listar com retry (NCBI tras vezes da timeout)
    listing=""
    for try in 1 2 3; do
        listing=$(curl -sL --max-time 30 "$base/" 2>/dev/null | grep -oE "${acc}_[A-Za-z0-9._]+" | head -1)
        [ -n "$listing" ] && break
        sleep 2
    done
    if [ -z "$listing" ]; then
        echo "  [$i/$total] $acc  [SKIP — listing falhou apos 3 tentativas]"
        continue
    fi

    url="$base/$listing/${listing}_genomic.fna.gz"
    echo "  [$i/$total] $acc  baixando $listing..."
    ok=0
    for try in 1 2 3; do
        if curl -fsSL --max-time 120 -o "$fa.gz" "$url"; then
            ok=1; break
        fi
        sleep 3
    done
    if [ "$ok" = "1" ]; then
        gunzip -f "$fa.gz"
    else
        echo "  [$i/$total] $acc  [FAIL apos 3 tentativas de download]"
        rm -f "$fa.gz"
    fi
done

n_genomes=$(ls "$GENOMES_DIR"/*.fna 2>/dev/null | wc -l | tr -d ' ')
echo "[Pangenome] $n_genomes genomas em disco"

# ====================================================================
# Etapa 2 — Anotacao via prodigal (rapido, sem dep tbl2asn)
# ====================================================================
echo ""
echo "[Pangenome] Etapa 2/3 — anotacao prodigal -> GFF panaroo-compativel"
i=0
for fa in "$GENOMES_DIR"/*.fna; do
    name=$(basename "$fa" .fna)
    out="$PROKKA_DIR/$name"
    gff="$out/$name.gff"
    i=$((i+1))
    if [ -f "$gff" ]; then
        echo "  [$i/$n_genomes] $name  (cached gff)"
        continue
    fi
    mkdir -p "$out"
    echo "  [$i/$n_genomes] $name  prodigal rodando..."
    # Prodigal gera GFFv3 + proteinas .faa. Panaroo aceita GFF com FASTA no final.
    if prodigal -i "$fa" -f gff -o "$out/${name}_prodigal.gff" \
                -a "$out/$name.faa" -d "$out/$name.ffn" -q -p meta 2>/dev/null; then
        # Converte GFF prodigal -> formato prokka que panaroo espera, anexa FASTA
        "$PY" "$REPO_DIR/prodigal_to_prokka_gff.py" \
              "$out/${name}_prodigal.gff" "$fa" "$name" "$gff"
    else
        echo "    [WARN] prodigal falhou em $name"
        rm -rf "$out"
    fi
done

# ====================================================================
# Etapa 3 — Panaroo (pangenome inferencia)
# ====================================================================
echo ""
echo "[Pangenome] Etapa 3/3 — panaroo"
if [ -d "$PANAROO_OUT" ] && [ -f "$PANAROO_OUT/gene_presence_absence.csv" ]; then
    echo "  [SKIP] panaroo ja rodou (output em $PANAROO_OUT)"
else
    rm -rf "$PANAROO_OUT"
    # Excluir os *_prodigal.gff intermediarios (sem ##FASTA, panaroo nao aceita)
    gff_files=()
    for g in "$PROKKA_DIR"/*/*.gff; do
        case "$g" in *_prodigal.gff) continue;; esac
        gff_files+=("$g")
    done
    if [ ${#gff_files[@]} -lt 5 ]; then
        echo "  [FAIL] poucos GFFs (${#gff_files[@]}). Precisa pelo menos 5."
        exit 1
    fi
    echo "  Rodando panaroo em ${#gff_files[@]} genomas (~30-90min)..."
    if command -v panaroo &> /dev/null; then
        panaroo -i "${gff_files[@]}" -o "$PANAROO_OUT" --clean-mode sensitive --remove-invalid-genes -t 1
    else
        echo "  [INFO] panaroo nao instalado. Instale com:"
        echo "         $PY -m pip install panaroo"
        exit 1
    fi
fi

# ====================================================================
# Sumario
# ====================================================================
echo ""
echo "[Pangenome] DONE"
if [ -f "$PANAROO_OUT/summary_statistics.txt" ]; then
    cat "$PANAROO_OUT/summary_statistics.txt"
fi
