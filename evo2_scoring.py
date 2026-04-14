#!/usr/bin/env python3
"""
Score funcional de variantes AMR usando Evo 2 (Arc Institute).
Calcula distancia funcional entre variantes usando embeddings genomicos,
em vez de simples contagem de mismatches.

Modos de operacao:
  - gpu: Usa modelo Evo 2 completo (requer CUDA + evo2 instalado)
  - lightweight: Scoring baseado em features de sequencia (sem GPU)

Inspirado em: EVEE (Goodfire/Mayo Clinic) - covariance probes sobre Evo 2

Referencia: https://github.com/ArcInstitute/evo2
"""

import os
import csv
import json
import math
import time
from datetime import datetime
from collections import Counter

from config import BASE_DIR, REPORTS_DIR, SEQUENCES_DIR, NCBI_BASE_URL, NCBI_EMAIL
from utils import gc_content, reverse_complement

# === Configuracao ===
EVO2_REPORTS_DIR = os.path.join(REPORTS_DIR, "evo2_scoring")
EVO2_CACHE_DIR = os.path.join(BASE_DIR, "evo2_cache")
EVO2_MODEL_NAME = "evo2_7b"  # Modelo mais leve, nao requer FP8
EVO2_LAYER = "blocks.28.mlp.l3"  # Embeddings intermediarios (recomendado pelo paper)


# =====================================================================
# MODO GPU - Evo 2 completo
# =====================================================================

def check_evo2_available() -> bool:
    """Verifica se Evo 2 esta instalado e GPU disponivel."""
    try:
        import torch
        if not torch.cuda.is_available():
            print("[EVO2] CUDA nao disponivel.")
            return False
        import evo2
        print(f"[EVO2] Evo 2 disponivel. GPU: {torch.cuda.get_device_name(0)}")
        return True
    except ImportError:
        print("[EVO2] Pacote evo2 nao instalado.")
        return False


def load_evo2_model():
    """Carrega modelo Evo 2."""
    from evo2 import Evo2
    import torch
    print(f"[EVO2] Carregando modelo {EVO2_MODEL_NAME}...")
    model = Evo2(EVO2_MODEL_NAME)
    print(f"[EVO2] Modelo carregado.")
    return model


def get_evo2_logits(model, sequence: str):
    """Obtem logits (probabilidades por posicao) do Evo 2."""
    import torch
    input_ids = torch.tensor(
        model.tokenizer.tokenize(sequence),
        dtype=torch.int,
    ).unsqueeze(0).to("cuda:0")

    outputs, _ = model(input_ids)
    logits = outputs[0]  # (seq_len, vocab_size)
    return logits.detach().cpu()


def get_evo2_embeddings(model, sequence: str, layer: str = None):
    """Obtem embeddings intermediarios do Evo 2."""
    import torch
    layer = layer or EVO2_LAYER
    input_ids = torch.tensor(
        model.tokenizer.tokenize(sequence),
        dtype=torch.int,
    ).unsqueeze(0).to("cuda:0")

    outputs, embeddings = model(
        input_ids,
        return_embeddings=True,
        layer_names=[layer],
    )
    return embeddings[layer].detach().cpu()


def compute_log_likelihood_gpu(model, sequence: str) -> float:
    """Calcula log-likelihood da sequencia segundo Evo 2."""
    import torch
    import torch.nn.functional as F

    logits = get_evo2_logits(model, sequence)
    probs = F.softmax(logits, dim=-1)

    # Mapear bases para indices do tokenizer
    base_to_idx = {}
    for base in "ACGT":
        tokens = model.tokenizer.tokenize(base)
        if tokens:
            base_to_idx[base] = tokens[0]

    total_ll = 0.0
    count = 0
    for i, base in enumerate(sequence[1:], 1):  # Skip first (no context)
        if base in base_to_idx:
            idx = base_to_idx[base]
            prob = probs[i - 1, idx].item()
            if prob > 0:
                total_ll += math.log(prob)
                count += 1

    return total_ll / count if count > 0 else 0.0


def compute_embedding_distance_gpu(model, seq_ref: str, seq_var: str) -> float:
    """Calcula distancia coseno entre embeddings de ref e variante."""
    import torch

    emb_ref = get_evo2_embeddings(model, seq_ref)
    emb_var = get_evo2_embeddings(model, seq_var)

    # Mean pooling sobre posicoes
    vec_ref = emb_ref.squeeze(0).mean(dim=0)
    vec_var = emb_var.squeeze(0).mean(dim=0)

    # Distancia coseno
    cos_sim = torch.nn.functional.cosine_similarity(
        vec_ref.unsqueeze(0), vec_var.unsqueeze(0)
    ).item()

    return 1.0 - cos_sim  # Distancia (0 = identico, 2 = oposto)


# =====================================================================
# MODO LIGHTWEIGHT - Scoring baseado em features de sequencia
# =====================================================================

def compute_kmer_frequencies(sequence: str, k: int = 4) -> dict:
    """Calcula frequencias de k-mers normalizadas."""
    seq = sequence.upper()
    kmers = Counter()
    total = len(seq) - k + 1
    if total <= 0:
        return {}
    for i in range(total):
        kmer = seq[i:i + k]
        if all(b in "ACGT" for b in kmer):
            kmers[kmer] += 1
    return {kmer: count / total for kmer, count in kmers.items()}


def compute_dinucleotide_bias(sequence: str) -> dict:
    """Calcula bias de dinucleotideos (CpG, etc.)."""
    seq = sequence.upper()
    observed = Counter()
    mono = Counter()

    for i in range(len(seq) - 1):
        di = seq[i:i + 2]
        if all(b in "ACGT" for b in di):
            observed[di] += 1
    for b in seq:
        if b in "ACGT":
            mono[b] += 1

    total = sum(mono.values())
    bias = {}
    for di, obs_count in observed.items():
        b1, b2 = di[0], di[1]
        expected = (mono[b1] / total) * (mono[b2] / total) * (total - 1)
        if expected > 0:
            bias[di] = obs_count / expected
    return bias


def compute_codon_adaptation_index(sequence: str) -> float:
    """
    Calcula indice simplificado de adaptacao de codons.
    Codons raros em bacterias indicam possivel perda de funcao.
    """
    # Tabela simplificada de codons frequentes em E. coli / K. pneumoniae
    # Valores relativos: 1.0 = mais frequente, <0.5 = raro
    freq_codons = {
        "GCG": 0.36, "GCC": 0.27, "GCA": 0.21, "GCT": 0.16,  # Ala
        "TGC": 0.56, "TGT": 0.44,  # Cys
        "GAC": 0.63, "GAT": 0.37,  # Asp
        "GAG": 0.31, "GAA": 0.69,  # Glu
        "TTC": 0.43, "TTT": 0.57,  # Phe
        "GGC": 0.40, "GGT": 0.34, "GGA": 0.11, "GGG": 0.15,  # Gly
        "ATG": 1.00,  # Met (start)
        "AAA": 0.76, "AAG": 0.24,  # Lys
        "CTG": 0.50, "TTA": 0.13, "TTG": 0.13, "CTC": 0.10, "CTT": 0.10, "CTA": 0.04,  # Leu
        "ATT": 0.51, "ATC": 0.42, "ATA": 0.07,  # Ile
        "TAA": 0.64, "TGA": 0.30, "TAG": 0.07,  # Stop
    }

    seq = sequence.upper()
    scores = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if len(codon) == 3 and all(b in "ACGT" for b in codon):
            score = freq_codons.get(codon, 0.5)  # Default 0.5 se nao mapeado
            scores.append(score)

    if not scores:
        return 0.5
    return sum(scores) / len(scores)


def compute_sequence_complexity(sequence: str) -> float:
    """Calcula complexidade linguistica (entropia de Shannon normalizada)."""
    seq = sequence.upper()
    if len(seq) == 0:
        return 0.0
    counts = Counter(seq)
    total = sum(counts.values())
    entropy = 0.0
    for count in counts.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log2(p)
    max_entropy = math.log2(min(4, len(set(seq))))
    return entropy / max_entropy if max_entropy > 0 else 0.0


def compute_functional_features(sequence: str) -> dict:
    """Computa todas as features funcionais de uma sequencia."""
    return {
        "gc_content": gc_content(sequence),
        "length": len(sequence),
        "complexity": compute_sequence_complexity(sequence),
        "codon_adaptation": compute_codon_adaptation_index(sequence),
        "kmer_4": compute_kmer_frequencies(sequence, 4),
        "dinuc_bias": compute_dinucleotide_bias(sequence),
    }


def compute_kmer_cosine_distance(kmer_a: dict, kmer_b: dict) -> float:
    """Calcula distancia coseno entre dois perfis de k-mers."""
    all_kmers = set(kmer_a.keys()) | set(kmer_b.keys())
    if not all_kmers:
        return 1.0

    dot = sum(kmer_a.get(k, 0) * kmer_b.get(k, 0) for k in all_kmers)
    norm_a = math.sqrt(sum(v ** 2 for v in kmer_a.values())) or 1e-10
    norm_b = math.sqrt(sum(v ** 2 for v in kmer_b.values())) or 1e-10

    cos_sim = dot / (norm_a * norm_b)
    return 1.0 - cos_sim


def compute_functional_distance_lightweight(features_ref: dict, features_var: dict) -> dict:
    """
    Calcula distancia funcional usando features de sequencia.
    Retorna score composto e componentes individuais.
    """
    # 1. Distancia de GC content
    gc_dist = abs(features_ref["gc_content"] - features_var["gc_content"])

    # 2. Distancia de complexidade
    complexity_dist = abs(features_ref["complexity"] - features_var["complexity"])

    # 3. Distancia de adaptacao de codons
    cai_dist = abs(features_ref["codon_adaptation"] - features_var["codon_adaptation"])

    # 4. Distancia de k-mers (proxy para embedding distance)
    kmer_dist = compute_kmer_cosine_distance(
        features_ref["kmer_4"], features_var["kmer_4"]
    )

    # 5. Distancia de dinucleotideos
    dinuc_a = features_ref["dinuc_bias"]
    dinuc_b = features_var["dinuc_bias"]
    all_di = set(dinuc_a.keys()) | set(dinuc_b.keys())
    dinuc_dist = 0.0
    if all_di:
        dinuc_dist = sum(abs(dinuc_a.get(d, 1.0) - dinuc_b.get(d, 1.0)) for d in all_di) / len(all_di)

    # Score composto (0 = identico, 1 = muito diferente)
    # Pesos calibrados para importancia relativa
    composite = (
        0.15 * min(gc_dist * 5, 1.0) +      # GC shift
        0.10 * min(complexity_dist * 5, 1.0) + # Complexity change
        0.25 * min(cai_dist * 3, 1.0) +       # Codon adaptation (funcional)
        0.35 * min(kmer_dist * 2, 1.0) +       # K-mer profile (proxy embedding)
        0.15 * min(dinuc_dist, 1.0)             # Dinucleotide bias
    )

    return {
        "functional_distance": round(composite, 4),
        "gc_shift": round(gc_dist, 4),
        "complexity_change": round(complexity_dist, 4),
        "codon_adaptation_shift": round(cai_dist, 4),
        "kmer_distance": round(kmer_dist, 4),
        "dinucleotide_shift": round(dinuc_dist, 4),
    }


def predict_functional_impact(distance: dict) -> dict:
    """
    Prediz impacto funcional baseado na distancia.
    Retorna classificacao e score de confianca.
    """
    fd = distance["functional_distance"]

    if fd < 0.02:
        impact = "conserved"
        label = "Funcao conservada"
        confidence = min(0.95, 1.0 - fd * 10)
        resistance_maintained = True
    elif fd < 0.08:
        impact = "likely_conserved"
        label = "Provavelmente conservada"
        confidence = 0.75 - (fd - 0.02) * 3
        resistance_maintained = True
    elif fd < 0.20:
        impact = "uncertain"
        label = "Impacto incerto"
        confidence = 0.50
        resistance_maintained = None  # Incerto
    elif fd < 0.40:
        impact = "likely_disrupted"
        label = "Provavelmente alterada"
        confidence = 0.60 + (fd - 0.20) * 1.5
        resistance_maintained = False
    else:
        impact = "disrupted"
        label = "Funcao provavelmente perdida"
        confidence = min(0.95, 0.70 + fd * 0.5)
        resistance_maintained = False

    return {
        "impact": impact,
        "label": label,
        "confidence": round(confidence, 3),
        "resistance_maintained": resistance_maintained,
        "functional_distance": fd,
    }


# =====================================================================
# PIPELINE DE SCORING
# =====================================================================

def fetch_sequence_ncbi(accession: str) -> str | None:
    """Baixa sequencia do NCBI."""
    import requests
    url = f"{NCBI_BASE_URL}/efetch.fcgi"
    params = {
        "db": "nucleotide", "id": accession,
        "rettype": "fasta", "retmode": "text", "email": NCBI_EMAIL,
    }
    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        content = resp.text.strip()
        if content.startswith(">"):
            lines = content.split("\n")
            return "".join(l.strip() for l in lines[1:] if not l.startswith(">")).upper()
    except Exception as e:
        print(f"    Erro fetch {accession}: {e}")
    return None


def load_cached_sequence(gene_name: str) -> str | None:
    """Tenta carregar sequencia do cache local."""
    from utils import parse_fasta
    fasta_path = os.path.join(SEQUENCES_DIR, f"{gene_name}.fasta")
    if os.path.exists(fasta_path):
        seqs = parse_fasta(fasta_path)
        if seqs:
            return list(seqs.values())[0]
    return None


def run_evo2_scoring(mode: str = "auto"):
    """
    Executa pipeline completo de scoring funcional.
    mode: 'gpu', 'lightweight', 'auto' (detecta automaticamente)
    """
    print("=" * 70)
    print("SCORE FUNCIONAL - Evo 2 Inspired Variant Analysis")
    print("SmartLab BacEnd | IA para Medicos")
    print("=" * 70)

    os.makedirs(EVO2_REPORTS_DIR, exist_ok=True)

    # Detectar modo
    use_gpu = False
    evo2_model = None
    if mode == "auto":
        use_gpu = check_evo2_available()
    elif mode == "gpu":
        use_gpu = check_evo2_available()
        if not use_gpu:
            print("[EVO2] GPU solicitada mas nao disponivel. Usando lightweight.")

    if use_gpu:
        evo2_model = load_evo2_model()
        print(f"[EVO2] Modo: GPU (Evo 2 {EVO2_MODEL_NAME})")
    else:
        print("[EVO2] Modo: Lightweight (features de sequencia)")

    # Carregar alvos e variantes
    targets = {}
    with open(os.path.join(BASE_DIR, "targets_brazil.csv")) as f:
        for row in csv.DictReader(f):
            targets[row["name"]] = dict(row)

    variants = []
    variants_csv = os.path.join(BASE_DIR, "targets_brazil_variants.csv")
    if os.path.exists(variants_csv):
        with open(variants_csv) as f:
            for row in csv.DictReader(f):
                variants.append(dict(row))

    print(f"\n[EVO2] Alvos de referencia: {len(targets)}")
    print(f"[EVO2] Variantes a analisar: {len(variants)}")

    # Mapear familias base (primeiras 12 entradas sao referencias)
    reference_genes = {name: data for name, data in list(targets.items())[:12]}

    # Carregar/baixar sequencias de referencia
    ref_sequences = {}
    ref_features = {}
    print(f"\n[EVO2] Carregando sequencias de referencia...")
    for gene_name, data in reference_genes.items():
        seq = load_cached_sequence(gene_name)
        if not seq:
            print(f"  Baixando {gene_name} ({data['gene_accession']})...")
            seq = fetch_sequence_ncbi(data["gene_accession"])
            time.sleep(0.4)
        if seq:
            ref_sequences[gene_name] = seq
            ref_features[gene_name] = compute_functional_features(seq)
            print(f"  {gene_name}: {len(seq)}bp, GC={gc_content(seq)*100:.1f}%")
        else:
            print(f"  {gene_name}: FALHA no fetch")

    # Scoring de variantes
    results = []
    print(f"\n[EVO2] Analisando variantes...")

    # Mapear variante -> gene de referencia
    family_map = {
        "mecA1": "mecA", "mecA2": "mecA",
        "blaKPC-2": "blaKPC", "blaKPC-3": "blaKPC", "blaKPC-4": "blaKPC",
        "blaKPC-5": "blaKPC", "blaKPC-11": "blaKPC", "blaKPC-30": "blaKPC",
        "blaKPC-31": "blaKPC",
        "blaNDM-1": "blaNDM", "blaNDM-2": "blaNDM", "blaNDM-5": "blaNDM",
        "blaNDM-7": "blaNDM",
        "vanA": "vanA",
        "mcr-1": "mcr-1", "mcr-1.1": "mcr-1", "mcr-5": "mcr-1",
        "blaOXA-48": "blaOXA-48", "blaOXA-181": "blaOXA-48", "blaOXA-232": "blaOXA-48",
        "blaVIM-1": "blaVIM", "blaVIM-2": "blaVIM", "blaVIM-4": "blaVIM",
        "blaIMP-1": "blaIMP", "blaIMP-6": "blaIMP",
        "blaGES-1": "blaGES", "blaGES-5": "blaGES",
        "blaCTX-M-2": "blaCTX-M-15", "blaCTX-M-8": "blaCTX-M-15",
        "blaCTX-M-9": "blaCTX-M-15", "blaCTX-M-14": "blaCTX-M-15",
        "blaCTX-M-27": "blaCTX-M-15",
        "qnrS1": "qnrS", "qnrS2": "qnrS",
        "armA": "armA",
    }

    for variant in variants:
        var_name = variant["name"]
        var_acc = variant["gene_accession"]
        var_family = variant.get("gene_family", "")

        # Encontrar referencia
        ref_gene = family_map.get(var_name, var_family)
        if ref_gene not in ref_sequences:
            print(f"  {var_name}: sem referencia ({ref_gene}), pulando")
            continue

        print(f"\n  [{var_name}] vs [{ref_gene}]")

        # Buscar sequencia da variante
        var_seq = fetch_sequence_ncbi(var_acc)
        time.sleep(0.4)

        if not var_seq:
            print(f"    FALHA no fetch")
            results.append({
                "variant": var_name, "reference": ref_gene, "accession": var_acc,
                "status": "fetch_failed",
            })
            continue

        print(f"    Sequencia: {len(var_seq)}bp")

        if use_gpu and evo2_model:
            # Modo GPU: embeddings reais do Evo 2
            try:
                emb_dist = compute_embedding_distance_gpu(
                    evo2_model, ref_sequences[ref_gene], var_seq
                )
                ll_ref = compute_log_likelihood_gpu(evo2_model, ref_sequences[ref_gene])
                ll_var = compute_log_likelihood_gpu(evo2_model, var_seq)

                distance = {
                    "functional_distance": round(emb_dist, 4),
                    "log_likelihood_ref": round(ll_ref, 4),
                    "log_likelihood_var": round(ll_var, 4),
                    "ll_ratio": round(ll_var - ll_ref, 4),
                    "mode": "gpu",
                }
            except Exception as e:
                print(f"    Erro GPU: {e}, fallback lightweight")
                var_features = compute_functional_features(var_seq)
                distance = compute_functional_distance_lightweight(
                    ref_features[ref_gene], var_features
                )
                distance["mode"] = "lightweight"
        else:
            # Modo lightweight
            var_features = compute_functional_features(var_seq)
            distance = compute_functional_distance_lightweight(
                ref_features[ref_gene], var_features
            )
            distance["mode"] = "lightweight"

        impact = predict_functional_impact(distance)

        result = {
            "variant": var_name,
            "reference": ref_gene,
            "accession": var_acc,
            "status": "scored",
            **distance,
            **impact,
        }
        results.append(result)

        icon = {"conserved": "=", "likely_conserved": "~", "uncertain": "?",
                "likely_disrupted": "!", "disrupted": "X"}
        print(f"    {icon.get(impact['impact'], '?')} {impact['label']} "
              f"(dist={distance['functional_distance']:.4f}, conf={impact['confidence']:.2f})")

    # Salvar resultados
    save_scoring_results(results)
    save_scoring_report(results, ref_sequences)
    save_scoring_json(results)

    # Resumo
    print(f"\n{'='*70}")
    print("SCORING FUNCIONAL CONCLUIDO")
    print(f"{'='*70}")

    scored = [r for r in results if r["status"] == "scored"]
    conserved = sum(1 for r in scored if r["impact"] in ("conserved", "likely_conserved"))
    uncertain = sum(1 for r in scored if r["impact"] == "uncertain")
    disrupted = sum(1 for r in scored if r["impact"] in ("likely_disrupted", "disrupted"))

    print(f"  Total analisadas: {len(scored)}")
    print(f"  Funcao conservada: {conserved}")
    print(f"  Impacto incerto: {uncertain}")
    print(f"  Funcao alterada: {disrupted}")
    print(f"  Modo: {'GPU (Evo 2)' if use_gpu else 'Lightweight'}")

    return results


def save_scoring_results(results: list):
    """Salva resultados em TSV."""
    filepath = os.path.join(EVO2_REPORTS_DIR, "functional_scores.tsv")
    scored = [r for r in results if r["status"] == "scored"]
    if not scored:
        return

    fieldnames = [
        "variant", "reference", "accession", "functional_distance",
        "impact", "label", "confidence", "resistance_maintained",
        "gc_shift", "kmer_distance", "codon_adaptation_shift", "mode",
    ]

    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for r in scored:
            writer.writerow(r)

    print(f"\n[EVO2] Scores salvos: {filepath}")


def save_scoring_report(results: list, ref_sequences: dict):
    """Gera relatorio textual."""
    filepath = os.path.join(EVO2_REPORTS_DIR, "functional_analysis_report.txt")

    with open(filepath, "w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("RELATORIO DE SCORING FUNCIONAL - Evo 2 Inspired\n")
        f.write(f"SmartLab BacEnd - {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write("=" * 70 + "\n\n")

        scored = [r for r in results if r["status"] == "scored"]
        failed = [r for r in results if r["status"] != "scored"]

        f.write(f"Variantes analisadas: {len(scored)}\n")
        f.write(f"Falhas no fetch: {len(failed)}\n\n")

        # Agrupar por referencia
        by_ref = {}
        for r in scored:
            ref = r["reference"]
            if ref not in by_ref:
                by_ref[ref] = []
            by_ref[ref].append(r)

        for ref, variants in sorted(by_ref.items()):
            f.write(f"\n{'─'*50}\n")
            ref_len = len(ref_sequences.get(ref, ""))
            f.write(f"[{ref}] Referencia: {ref_len}bp\n")
            f.write(f"{'Variante':<20} {'Dist.Func.':<12} {'Impacto':<25} {'Conf.':<8} {'Resist.'}\n")
            f.write("-" * 80 + "\n")

            for v in sorted(variants, key=lambda x: x["functional_distance"]):
                resist = "SIM" if v["resistance_maintained"] else ("?" if v["resistance_maintained"] is None else "NAO")
                f.write(f"{v['variant']:<20} {v['functional_distance']:<12.4f} "
                        f"{v['label']:<25} {v['confidence']:<8.3f} {resist}\n")

        # Resumo de priorizacao
        f.write(f"\n{'='*70}\n")
        f.write("PRIORIZACAO DINAMICA (baseada em impacto funcional)\n")
        f.write(f"{'='*70}\n\n")

        # Variantes com funcao alterada = podem escapar deteccao
        disrupted = [r for r in scored if r["impact"] in ("likely_disrupted", "disrupted")]
        if disrupted:
            f.write("ATENCAO - Variantes com funcao possivelmente alterada:\n")
            f.write("(podem ter perdido resistencia, reduzindo relevancia clinica)\n\n")
            for v in disrupted:
                f.write(f"  {v['variant']} (ref: {v['reference']}) - "
                        f"dist={v['functional_distance']:.4f} - {v['label']}\n")
        else:
            f.write("Todas as variantes analisadas mantem funcao conservada ou provavel.\n")
            f.write("O painel de deteccao cobre variantes clinicamente relevantes.\n")

    print(f"[EVO2] Relatorio salvo: {filepath}")


def save_scoring_json(results: list):
    """Salva resultados em JSON para consumo pelo frontend."""
    filepath = os.path.join(EVO2_REPORTS_DIR, "functional_scores.json")
    scored = [r for r in results if r["status"] == "scored"]

    output = {}
    for r in scored:
        output[r["variant"]] = {
            "reference": r["reference"],
            "functional_distance": r["functional_distance"],
            "impact": r["impact"],
            "label": r["label"],
            "confidence": r["confidence"],
            "resistance_maintained": r["resistance_maintained"],
            "mode": r.get("mode", "lightweight"),
            "scored_at": datetime.now().isoformat(),
        }

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False, indent=2)

    print(f"[EVO2] JSON salvo: {filepath}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Scoring funcional de variantes AMR")
    parser.add_argument("--mode", choices=["auto", "gpu", "lightweight"],
                        default="auto", help="Modo de scoring")
    args = parser.parse_args()
    run_evo2_scoring(mode=args.mode)
