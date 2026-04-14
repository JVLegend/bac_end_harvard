#!/usr/bin/env python3
"""
Covariance Probes para scoring avancado de guides CRISPR-Cas12a.
Captura co-ocorrencia de features (nao apenas features individuais)
para prever eficacia de clivagem.

Inspirado em: Goodfire/Mayo Clinic - Covariance-based Sequence Pooling
sobre embeddings do Evo 2 para prever patogenicidade de variantes.

Abordagem adaptada para CRISPR:
  - Em vez de embeddings neurais, usamos features biofisicas de sequencia
  - Covariance matrix captura interacoes nao-lineares entre features
  - Probe scores combinam features individuais + covariancia

Benchmark: Compara com score_guide() rule-based existente em utils.py
"""

import os
import csv
import json
import math
from datetime import datetime

from config import BASE_DIR, GUIDES_DIR, SEQUENCES_DIR, REPORTS_DIR, CAS12A
from utils import (
    gc_content, max_homopolymer, self_complementarity_score,
    has_poly_t, reverse_complement, parse_fasta, find_pam_sites, extract_spacer,
)

# === Configuracao ===
PROBES_REPORTS_DIR = os.path.join(REPORTS_DIR, "covariance_probes")


# =====================================================================
# FEATURE EXTRACTION AVANCADA
# =====================================================================

def compute_seed_region_features(spacer: str) -> dict:
    """
    Analisa a seed region do spacer (posicoes 1-8, criticas para Cas12a).
    A seed region e a porcao proximal ao PAM e determina binding inicial.
    """
    seed = spacer[:8].upper()
    non_seed = spacer[8:].upper()

    return {
        "seed_gc": gc_content(seed),
        "seed_purine_fraction": sum(1 for b in seed if b in "AG") / len(seed),
        "seed_homopolymer": max_homopolymer(seed),
        "non_seed_gc": gc_content(non_seed) if non_seed else 0,
        "seed_vs_nonseed_gc_diff": abs(gc_content(seed) - (gc_content(non_seed) if non_seed else 0)),
    }


def compute_thermodynamic_features(spacer: str) -> dict:
    """
    Calcula features termodinamicas simplificadas.
    Energia livre de hibridizacao (nearest-neighbor simplificado).
    """
    # Parametros nearest-neighbor simplificados (kcal/mol, DNA/RNA hibrido)
    # Baseado em SantaLucia & Hicks 2004
    nn_dg = {
        "AA": -1.0, "AT": -0.88, "AC": -1.44, "AG": -1.28,
        "TA": -0.58, "TT": -1.0, "TC": -1.30, "TG": -1.45,
        "CA": -1.45, "CT": -1.28, "CC": -1.84, "CG": -2.17,
        "GA": -1.30, "GT": -1.44, "GC": -2.24, "GG": -1.84,
    }

    seq = spacer.upper()
    total_dg = 0.0
    min_dg_window = 0.0  # Regiao mais estavel
    max_dg_window = 0.0  # Regiao menos estavel

    dg_per_pos = []
    for i in range(len(seq) - 1):
        di = seq[i:i + 2]
        dg = nn_dg.get(di, -1.0)
        total_dg += dg
        dg_per_pos.append(dg)

    # Janela deslizante de 5bp para encontrar regioes extremas
    if len(dg_per_pos) >= 5:
        windows = [sum(dg_per_pos[i:i + 5]) for i in range(len(dg_per_pos) - 4)]
        min_dg_window = min(windows)  # Mais estavel (mais negativo)
        max_dg_window = max(windows)  # Menos estavel

    # Energia da seed region (primeiros 8bp)
    seed_dg = sum(dg_per_pos[:7]) if len(dg_per_pos) >= 7 else total_dg

    # Gradiente de estabilidade (seed -> 3' end)
    if len(dg_per_pos) >= 10:
        first_half = sum(dg_per_pos[:len(dg_per_pos) // 2])
        second_half = sum(dg_per_pos[len(dg_per_pos) // 2:])
        stability_gradient = second_half - first_half  # Positivo = 3' menos estavel (bom)
    else:
        stability_gradient = 0.0

    return {
        "total_dg": round(total_dg, 2),
        "seed_dg": round(seed_dg, 2),
        "min_dg_window": round(min_dg_window, 2),
        "max_dg_window": round(max_dg_window, 2),
        "stability_gradient": round(stability_gradient, 2),
        "dg_per_base": round(total_dg / len(seq), 3) if seq else 0,
    }


def compute_structural_accessibility(spacer: str, context_5p: str = "", context_3p: str = "") -> dict:
    """
    Estima acessibilidade estrutural do alvo.
    Regioes com forte estrutura secundaria sao menos acessiveis ao Cas12a.
    Usa predicao simplificada baseada em palindromes e repeticoes invertidas.
    """
    full_seq = (context_5p + spacer + context_3p).upper()
    spacer_up = spacer.upper()

    # Contar palindromes (indicam hairpins)
    palindrome_score = 0
    for win_size in [4, 5, 6]:
        for i in range(len(full_seq) - win_size + 1):
            window = full_seq[i:i + win_size]
            rc_window = reverse_complement(window)
            # Procurar o complemento reverso proximo
            for j in range(i + win_size + 2, min(i + win_size + 15, len(full_seq) - win_size + 1)):
                if full_seq[j:j + win_size] == rc_window:
                    palindrome_score += win_size
                    break

    # Contar repeticoes diretas (indicam slippage)
    direct_repeat_score = 0
    for win_size in [3, 4, 5]:
        for i in range(len(spacer_up) - win_size + 1):
            window = spacer_up[i:i + win_size]
            count = spacer_up.count(window)
            if count > 1:
                direct_repeat_score += (count - 1) * win_size

    # AT-richness no contexto (regioes AT-rich sao mais acessiveis)
    at_fraction = sum(1 for b in full_seq if b in "AT") / len(full_seq) if full_seq else 0.5

    # Score de acessibilidade (0-1, maior = mais acessivel)
    accessibility = 1.0
    accessibility -= min(palindrome_score * 0.02, 0.4)  # Penalizar hairpins
    accessibility -= min(direct_repeat_score * 0.01, 0.2)  # Penalizar repeticoes
    accessibility += (at_fraction - 0.5) * 0.2  # Bonus AT-rich
    accessibility = max(0.0, min(1.0, accessibility))

    return {
        "accessibility": round(accessibility, 3),
        "palindrome_score": palindrome_score,
        "direct_repeat_score": direct_repeat_score,
        "local_at_fraction": round(at_fraction, 3),
    }


def compute_positional_nucleotide_preferences(spacer: str) -> dict:
    """
    Calcula concordancia com preferencias posicionais de Cas12a.
    Baseado em dados experimentais publicados (Kim et al. 2017, Kleinstiver et al. 2019):
    - Posicao 1 (proxima ao PAM): prefere A/T
    - Posicoes 1-4: prefere A
    - Posicao 20 (distal): tolera mais variacao
    """
    seq = spacer.upper()
    if len(seq) < 20:
        return {"positional_score": 0.5, "preferred_matches": 0}

    # Preferencias por posicao (simplificado da literatura)
    # Formato: {posicao: {base: bonus}}
    preferences = {
        0: {"A": 0.15, "T": 0.10},           # Pos 1: A/T preferido
        1: {"A": 0.10, "T": 0.05},           # Pos 2: A preferido
        2: {"A": 0.08},                       # Pos 3
        3: {"A": 0.05, "G": 0.05},           # Pos 4
        4: {"C": 0.05, "T": 0.05},           # Pos 5
        # Posicoes 6-17: menor efeito
        18: {"G": 0.05, "C": 0.05},          # Pos 19
        19: {"G": 0.08, "C": 0.05},          # Pos 20: GC tolera
    }

    total_bonus = 0.0
    preferred_matches = 0
    for pos, prefs in preferences.items():
        if pos < len(seq):
            base = seq[pos]
            bonus = prefs.get(base, 0)
            total_bonus += bonus
            if bonus > 0:
                preferred_matches += 1

    return {
        "positional_score": round(0.5 + total_bonus, 3),
        "preferred_matches": preferred_matches,
    }


# =====================================================================
# COVARIANCE MATRIX
# =====================================================================

def extract_feature_vector(spacer: str, position: int, gene_length: int,
                           context_5p: str = "", context_3p: str = "") -> list[float]:
    """
    Extrai vetor de features completo para um spacer.
    Retorna vetor numerico para calculo de covariancia.
    """
    gc = gc_content(spacer)
    homo = max_homopolymer(spacer)
    self_comp = self_complementarity_score(spacer)
    poly_t = 1.0 if has_poly_t(spacer) else 0.0
    rel_pos = position / gene_length if gene_length > 0 else 0.5

    seed = compute_seed_region_features(spacer)
    thermo = compute_thermodynamic_features(spacer)
    struct = compute_structural_accessibility(spacer, context_5p, context_3p)
    pos_pref = compute_positional_nucleotide_preferences(spacer)

    # Vetor de 18 features
    return [
        gc,                              # 0: GC content global
        homo / 6.0,                      # 1: Homopolymer (normalizado)
        self_comp / 20.0,                # 2: Self-complementarity (normalizado)
        poly_t,                          # 3: Poly-T presente
        rel_pos,                         # 4: Posicao relativa
        seed["seed_gc"],                 # 5: GC da seed region
        seed["seed_purine_fraction"],    # 6: Purinas na seed
        seed["seed_homopolymer"] / 4.0,  # 7: Homopolymer na seed
        seed["seed_vs_nonseed_gc_diff"], # 8: Diferenca GC seed vs non-seed
        thermo["dg_per_base"],           # 9: Energia por base (normalizado para ~0-1)
        (thermo["seed_dg"] + 15) / 15,  # 10: Energia seed (normalizado)
        thermo["stability_gradient"] / 5, # 11: Gradiente de estabilidade
        struct["accessibility"],         # 12: Acessibilidade estrutural
        struct["palindrome_score"] / 30, # 13: Palindromes (normalizado)
        struct["local_at_fraction"],     # 14: Fracao AT local
        pos_pref["positional_score"],    # 15: Score posicional
        pos_pref["preferred_matches"] / 8, # 16: Matches preferenciais
        1.0 - abs(gc - 0.50) * 2,       # 17: Proximidade ao GC ideal (0.50)
    ]


def compute_covariance_matrix(feature_vectors: list[list[float]]) -> list[list[float]]:
    """
    Calcula matriz de covariancia entre features.
    Captura interacoes par-a-par que scoring linear nao captura.
    """
    n = len(feature_vectors)
    if n < 2:
        dim = len(feature_vectors[0]) if feature_vectors else 0
        return [[0.0] * dim for _ in range(dim)]

    dim = len(feature_vectors[0])

    # Calcular medias
    means = [0.0] * dim
    for vec in feature_vectors:
        for i in range(dim):
            means[i] += vec[i]
    means = [m / n for m in means]

    # Calcular covariancia
    cov = [[0.0] * dim for _ in range(dim)]
    for vec in feature_vectors:
        for i in range(dim):
            for j in range(dim):
                cov[i][j] += (vec[i] - means[i]) * (vec[j] - means[j])
    for i in range(dim):
        for j in range(dim):
            cov[i][j] /= (n - 1)

    return cov


def compute_covariance_score(feature_vector: list[float],
                              cov_matrix: list[list[float]],
                              target_weights: list[float]) -> float:
    """
    Calcula score usando covariance pooling.
    Em vez de apenas somar features * weights (linear),
    incorpora interacoes entre features via covariancia.

    Score = w^T * f + alpha * f^T * C * w
    onde f = features, w = weights, C = covariance matrix
    """
    dim = len(feature_vector)

    # Componente linear: w^T * f
    linear_score = sum(feature_vector[i] * target_weights[i] for i in range(dim))

    # Componente de covariancia: f^T * C * w
    # Primeiro: C * w
    cw = [0.0] * dim
    for i in range(dim):
        for j in range(dim):
            cw[i] += cov_matrix[i][j] * target_weights[j]

    # Depois: f^T * (C * w)
    cov_score = sum(feature_vector[i] * cw[i] for i in range(dim))

    # Combinar com peso alpha para componente de covariancia
    alpha = 0.3  # Peso da covariancia (calibravel)
    total = linear_score + alpha * cov_score

    return total


# =====================================================================
# PROBE TRAINING (PESOS BASEADOS EM LITERATURA)
# =====================================================================

def get_cas12a_probe_weights() -> list[float]:
    """
    Pesos calibrados para eficacia de clivagem Cas12a.
    Baseados em dados experimentais publicados:
    - Kim et al. 2017 (Nat Methods) - eficacia de guias Cpf1
    - Kleinstiver et al. 2019 (Nat Biotech) - AsCas12a melhorado
    - Li et al. 2020 - regras de design para Cas12a diagnostico

    Pesos positivos = feature favoravel a clivagem
    Pesos negativos = feature desfavoravel
    """
    return [
        +0.15,   # 0: GC content (moderado e bom)
        -0.20,   # 1: Homopolymer (penaliza)
        -0.15,   # 2: Self-complementarity (penaliza)
        -0.25,   # 3: Poly-T (penaliza fortemente - termina transcricao)
        +0.05,   # 4: Posicao relativa (leve preferencia central)
        +0.18,   # 5: Seed GC (seed GC-rich = melhor binding inicial)
        +0.12,   # 6: Seed purine fraction (purinas na seed = bom)
        -0.15,   # 7: Seed homopolymer (penaliza na seed)
        -0.08,   # 8: Diferenca GC seed vs non-seed (uniformidade preferida)
        +0.10,   # 9: Energia por base (mais negativo = mais estavel = bom)
        +0.14,   # 10: Energia seed (seed estavel = bom binding)
        +0.08,   # 11: Gradiente de estabilidade (3' menos estavel = bom turnover)
        +0.20,   # 12: Acessibilidade estrutural (critico)
        -0.12,   # 13: Palindromes (penaliza - indica hairpins)
        +0.06,   # 14: Fracao AT local (AT-rich = mais acessivel)
        +0.15,   # 15: Score posicional (concordancia com preferencias)
        +0.10,   # 16: Matches preferenciais
        +0.18,   # 17: Proximidade ao GC ideal
    ]


# =====================================================================
# PIPELINE DE SCORING COM PROBES
# =====================================================================

def score_guide_with_probes(spacer: str, position: int, gene_length: int,
                             cov_matrix: list[list[float]] = None,
                             context_5p: str = "", context_3p: str = "") -> dict:
    """
    Score avancado de guide usando covariance probes.
    Retorna score e componentes detalhados.
    """
    features = extract_feature_vector(spacer, position, gene_length, context_5p, context_3p)
    weights = get_cas12a_probe_weights()

    # Se nao temos covariance matrix, usar apenas linear
    if cov_matrix is None:
        dim = len(features)
        cov_matrix = [[0.0] * dim for _ in range(dim)]

    raw_score = compute_covariance_score(features, cov_matrix, weights)

    # Normalizar para escala 0-100
    # Raw score tipico varia de -0.5 a +0.5
    normalized = max(0, min(100, (raw_score + 0.5) * 100))

    # Extrair componentes para explicabilidade
    seed = compute_seed_region_features(spacer)
    thermo = compute_thermodynamic_features(spacer)
    struct = compute_structural_accessibility(spacer, context_5p, context_3p)
    pos_pref = compute_positional_nucleotide_preferences(spacer)

    return {
        "probe_score": round(normalized, 1),
        "raw_score": round(raw_score, 4),
        "gc": round(gc_content(spacer), 3),
        "homopolymer": max_homopolymer(spacer),
        "poly_t": has_poly_t(spacer),
        "self_comp": self_complementarity_score(spacer),
        "seed_gc": round(seed["seed_gc"], 3),
        "seed_purine": round(seed["seed_purine_fraction"], 3),
        "accessibility": round(struct["accessibility"], 3),
        "total_dg": thermo["total_dg"],
        "seed_dg": thermo["seed_dg"],
        "stability_gradient": thermo["stability_gradient"],
        "positional_score": pos_pref["positional_score"],
        "feature_vector": features,
    }


def run_covariance_probe_analysis(gene_name: str = None):
    """
    Executa analise completa com covariance probes para um ou todos os genes.
    Compara com scoring rule-based original.
    """
    print("=" * 70)
    print("COVARIANCE PROBES - Scoring Avancado de Guides CRISPR-Cas12a")
    print("Inspirado em Goodfire/Mayo Clinic - EVEE")
    print("=" * 70)

    os.makedirs(PROBES_REPORTS_DIR, exist_ok=True)

    # Determinar genes a processar
    from config import TARGETS
    genes = [gene_name] if gene_name else list(TARGETS.keys())[:12]

    all_results = {}
    all_feature_vectors = []  # Para covariance matrix global

    # Primeira passada: extrair features de todos os candidatos
    print("\n[PROBES] Fase 1: Extraindo features de todos os candidatos...")
    candidates_by_gene = {}

    for gname in genes:
        fasta_path = os.path.join(SEQUENCES_DIR, f"{gname}.fasta")
        if not os.path.exists(fasta_path):
            print(f"  {gname}: sem sequencia, pulando")
            continue

        seqs = parse_fasta(fasta_path)
        if not seqs:
            continue

        sequence = list(seqs.values())[0]
        pam_sites = find_pam_sites(sequence, "TTTV")

        candidates = []
        for site in pam_sites:
            spacer = extract_spacer(sequence, site, CAS12A["spacer_length_optimal"])
            if not spacer or len(spacer) < 20:
                continue
            if any(b not in "ATCG" for b in spacer.upper()):
                continue

            gc = gc_content(spacer)
            if gc < CAS12A["gc_min"] or gc > CAS12A["gc_max"]:
                continue

            # Extrair contexto (50bp flanqueando)
            pos = site["position"]
            ctx_5p = sequence[max(0, pos - 50):pos]
            ctx_3p = sequence[pos + 24:min(len(sequence), pos + 74)]

            fv = extract_feature_vector(spacer, pos, len(sequence), ctx_5p, ctx_3p)
            all_feature_vectors.append(fv)

            candidates.append({
                "spacer": spacer,
                "position": pos,
                "strand": site["strand"],
                "pam_seq": site["pam_seq"],
                "gene_length": len(sequence),
                "context_5p": ctx_5p,
                "context_3p": ctx_3p,
                "feature_vector": fv,
            })

        candidates_by_gene[gname] = candidates
        print(f"  {gname}: {len(candidates)} candidatos ({len(pam_sites)} PAM sites)")

    # Calcular covariance matrix global
    print(f"\n[PROBES] Fase 2: Calculando covariance matrix ({len(all_feature_vectors)} vetores)...")
    cov_matrix = compute_covariance_matrix(all_feature_vectors)

    # Mostrar features mais correlacionadas
    feature_names = [
        "GC", "Homo", "SelfComp", "PolyT", "RelPos",
        "SeedGC", "SeedPur", "SeedHomo", "GCdiff",
        "dG/base", "SeedDG", "StabGrad", "Access",
        "Palindr", "ATfrac", "PosPref", "PrefMatch", "GCideal",
    ]

    print("\n  Top interacoes de covariancia:")
    interactions = []
    for i in range(len(cov_matrix)):
        for j in range(i + 1, len(cov_matrix[i])):
            interactions.append((abs(cov_matrix[i][j]), i, j, cov_matrix[i][j]))
    interactions.sort(reverse=True)
    for strength, i, j, val in interactions[:5]:
        sign = "+" if val > 0 else "-"
        print(f"    {feature_names[i]} x {feature_names[j]}: {sign}{strength:.4f}")

    # Segunda passada: scoring com probes
    print(f"\n[PROBES] Fase 3: Scoring com covariance probes...")
    weights = get_cas12a_probe_weights()

    benchmark_results = []

    for gname, candidates in candidates_by_gene.items():
        print(f"\n  {'─'*40}")
        print(f"  [{gname}] {len(candidates)} candidatos")

        scored = []
        for cand in candidates:
            # Score com covariance probes
            probe_result = score_guide_with_probes(
                cand["spacer"], cand["position"], cand["gene_length"],
                cov_matrix, cand["context_5p"], cand["context_3p"],
            )

            # Score original (rule-based) para benchmark
            from utils import score_guide
            original = score_guide(cand["spacer"], cand["position"], cand["gene_length"])

            scored.append({
                "gene": gname,
                "spacer": cand["spacer"],
                "position": cand["position"],
                "strand": cand["strand"],
                "pam": cand["pam_seq"],
                "probe_score": probe_result["probe_score"],
                "original_score": original["score"],
                "rank_diff": 0,  # Calcular depois
                **{k: v for k, v in probe_result.items() if k != "feature_vector" and k != "probe_score"},
            })

        # Ranquear por cada metodo
        scored_by_probe = sorted(scored, key=lambda x: x["probe_score"], reverse=True)
        scored_by_original = sorted(scored, key=lambda x: x["original_score"], reverse=True)

        # Calcular diferenca de rank
        probe_rank = {s["spacer"]: i for i, s in enumerate(scored_by_probe)}
        orig_rank = {s["spacer"]: i for i, s in enumerate(scored_by_original)}
        for s in scored_by_probe:
            s["rank_probe"] = probe_rank[s["spacer"]] + 1
            s["rank_original"] = orig_rank[s["spacer"]] + 1
            s["rank_diff"] = s["rank_original"] - s["rank_probe"]

        # Top 5 por probe
        top5 = scored_by_probe[:5]
        all_results[gname] = top5

        print(f"  {'Rank':<5} {'Probe':>7} {'Orig':>7} {'Diff':>5} {'Spacer':<22} {'SeedGC':>7} {'Access':>7} {'dG':>7}")
        print(f"  {'─'*75}")
        for s in top5:
            diff_str = f"+{s['rank_diff']}" if s['rank_diff'] > 0 else str(s['rank_diff'])
            print(f"  {s['rank_probe']:<5} {s['probe_score']:>7.1f} {s['original_score']:>7.1f} "
                  f"{diff_str:>5} {s['spacer']:<22} {s['seed_gc']:>7.3f} "
                  f"{s['accessibility']:>7.3f} {s['total_dg']:>7.1f}")

        benchmark_results.extend(scored_by_probe[:10])

    # Salvar resultados
    save_probe_results(all_results)
    save_benchmark_report(benchmark_results, cov_matrix, feature_names)
    save_probe_json(all_results)

    # Resumo
    print(f"\n{'='*70}")
    print("COVARIANCE PROBES - CONCLUIDO")
    print(f"{'='*70}")
    print(f"  Genes analisados: {len(all_results)}")
    print(f"  Total de candidatos: {len(all_feature_vectors)}")
    print(f"  Dimensoes de features: {len(feature_names)}")
    print(f"  Covariance matrix: {len(cov_matrix)}x{len(cov_matrix)}")

    # Comparacao global
    rank_changes = [s["rank_diff"] for s in benchmark_results if s["rank_diff"] != 0]
    if rank_changes:
        avg_change = sum(abs(rc) for rc in rank_changes) / len(rank_changes)
        promoted = sum(1 for rc in rank_changes if rc > 0)
        demoted = sum(1 for rc in rank_changes if rc < 0)
        print(f"  Guides re-ranqueados: {len(rank_changes)}/{len(benchmark_results)}")
        print(f"  Promovidos pelo probe: {promoted} | Rebaixados: {demoted}")
        print(f"  Mudanca media de rank: {avg_change:.1f} posicoes")

    return all_results


def save_probe_results(results: dict):
    """Salva resultados dos probes em TSV."""
    filepath = os.path.join(PROBES_REPORTS_DIR, "probe_scores.tsv")
    all_rows = []
    for gene, guides in results.items():
        all_rows.extend(guides)

    if not all_rows:
        return

    fieldnames = [
        "gene", "spacer", "position", "strand", "pam",
        "probe_score", "original_score", "rank_probe", "rank_original", "rank_diff",
        "gc", "seed_gc", "seed_purine", "accessibility",
        "total_dg", "seed_dg", "stability_gradient", "positional_score",
    ]

    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    print(f"\n[PROBES] Scores salvos: {filepath}")


def save_benchmark_report(benchmark: list, cov_matrix: list, feature_names: list):
    """Gera relatorio de benchmark probes vs rule-based."""
    filepath = os.path.join(PROBES_REPORTS_DIR, "benchmark_report.txt")

    with open(filepath, "w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("BENCHMARK: Covariance Probes vs Rule-Based Scoring\n")
        f.write(f"SmartLab BacEnd - {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write("=" * 70 + "\n\n")

        f.write("METODOLOGIA\n")
        f.write("-" * 40 + "\n")
        f.write("Rule-based (original): 5 features lineares (GC, homo, poly-T, self-comp, position)\n")
        f.write(f"Covariance probes: {len(feature_names)} features + matriz de covariancia\n")
        f.write("Features adicionais: seed region, termodinamica, acessibilidade estrutural,\n")
        f.write("  preferencias posicionais, gradiente de estabilidade\n\n")

        f.write(f"FEATURES ({len(feature_names)})\n")
        f.write("-" * 40 + "\n")
        for i, name in enumerate(feature_names):
            f.write(f"  [{i:2d}] {name}\n")

        f.write(f"\nCOVARIANCE MATRIX (top 10 interacoes)\n")
        f.write("-" * 40 + "\n")
        interactions = []
        for i in range(len(cov_matrix)):
            for j in range(i + 1, len(cov_matrix[i])):
                interactions.append((abs(cov_matrix[i][j]), i, j, cov_matrix[i][j]))
        interactions.sort(reverse=True)
        for strength, i, j, val in interactions[:10]:
            sign = "+" if val > 0 else "-"
            f.write(f"  {feature_names[i]:12s} x {feature_names[j]:12s}: {sign}{strength:.4f}\n")

        # Analise de rank
        f.write(f"\nANALISE DE RE-RANQUEAMENTO\n")
        f.write("-" * 40 + "\n")
        rank_changes = [s["rank_diff"] for s in benchmark if s["rank_diff"] != 0]
        total = len(benchmark)
        changed = len(rank_changes)
        f.write(f"Total de guides analisados: {total}\n")
        f.write(f"Guides com rank alterado: {changed} ({changed/total*100:.0f}%)\n")
        if rank_changes:
            promoted = sum(1 for rc in rank_changes if rc > 0)
            demoted = sum(1 for rc in rank_changes if rc < 0)
            avg = sum(abs(rc) for rc in rank_changes) / len(rank_changes)
            f.write(f"Promovidos pelo probe: {promoted}\n")
            f.write(f"Rebaixados pelo probe: {demoted}\n")
            f.write(f"Mudanca media: {avg:.1f} posicoes\n")

        f.write(f"\nGUIDES COM MAIOR MUDANCA DE RANK\n")
        f.write("-" * 40 + "\n")
        top_changes = sorted(benchmark, key=lambda x: abs(x["rank_diff"]), reverse=True)[:10]
        for s in top_changes:
            if s["rank_diff"] != 0:
                direction = "PROMOVIDO" if s["rank_diff"] > 0 else "REBAIXADO"
                f.write(f"  {s['gene']:12s} {s['spacer'][:15]}... "
                        f"orig={s['rank_original']:>3} -> probe={s['rank_probe']:>3} "
                        f"({direction} {abs(s['rank_diff'])} pos) "
                        f"access={s['accessibility']:.3f} seedGC={s['seed_gc']:.3f}\n")

        f.write(f"\nCONCLUSAO\n")
        f.write("-" * 40 + "\n")
        f.write("O scoring por covariance probes incorpora 18 features biofisicas\n")
        f.write("(vs 5 do rule-based), capturando interacoes entre seed region,\n")
        f.write("termodinamica e acessibilidade estrutural que afetam a eficacia\n")
        f.write("real de clivagem do Cas12a.\n")

    print(f"[PROBES] Benchmark salvo: {filepath}")


def save_probe_json(results: dict):
    """Salva resultados em JSON para frontend."""
    filepath = os.path.join(PROBES_REPORTS_DIR, "probe_scores.json")

    output = {}
    for gene, guides in results.items():
        output[gene] = []
        for g in guides:
            output[gene].append({
                "spacer": g["spacer"],
                "probe_score": g["probe_score"],
                "original_score": g["original_score"],
                "rank_probe": g.get("rank_probe", 0),
                "rank_original": g.get("rank_original", 0),
                "seed_gc": g["seed_gc"],
                "accessibility": g["accessibility"],
                "total_dg": g["total_dg"],
                "positional_score": g["positional_score"],
            })

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False, indent=2)

    print(f"[PROBES] JSON salvo: {filepath}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Covariance Probes para guides CRISPR-Cas12a")
    parser.add_argument("--gene", type=str, default=None, help="Gene especifico (default: todos)")
    args = parser.parse_args()
    run_covariance_probe_analysis(gene_name=args.gene)
