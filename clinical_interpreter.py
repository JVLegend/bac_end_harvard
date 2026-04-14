#!/usr/bin/env python3
"""
Interpretador clinico via Claude API.
Gera explicacoes em linguagem natural para resultados do pipeline,
inspirado na abordagem EVEE (Goodfire/Mayo Clinic).

Audiencias:
  - medico: interpretacao clinica direta para decisao no leito
  - gestor: visao de custo-beneficio e impacto hospitalar
  - pesquisador: detalhes tecnicos completos

Requer: ANTHROPIC_API_KEY no ambiente ou em .env
"""

import os
import csv
import json
import time
from datetime import datetime

from config import BASE_DIR, TARGETS, GUIDES_DIR, PRIMERS_DIR, REPORTS_DIR, TARGETS_CSV

# === Configuracao ===
REPORTS_CLINICAL_DIR = os.path.join(REPORTS_DIR, "clinical")
INTERPRETATION_CACHE = os.path.join(REPORTS_CLINICAL_DIR, ".cache.json")

# Contexto clinico por familia genica (conhecimento embutido)
CLINICAL_CONTEXT = {
    "mecA": {
        "pathogen": "MRSA (Staphylococcus aureus resistente a meticilina)",
        "threat_level": "CRITICO",
        "prevalence_br": "Alta - principal causa de HAI em UTIs brasileiras",
        "treatment_impact": "Exclui todos os beta-lactamicos. Opcoes: vancomicina, daptomicina, linezolida",
        "isolation": "Precaucao de contato obrigatoria",
        "mortality_increase": "2-3x maior mortalidade vs MSSA",
    },
    "blaKPC": {
        "pathogen": "CRE (Enterobacteriales produtoras de carbapenemase KPC)",
        "threat_level": "CRITICO - Prioridade OMS",
        "prevalence_br": "77% de A. baumannii resistente em UTIs BR (ANVISA 2023)",
        "treatment_impact": "Resistencia a carbapenems. Opcoes limitadas: ceftazidima-avibactam, meropenem-vaborbactam",
        "isolation": "Precaucao de contato + coorte",
        "mortality_increase": "40-50% mortalidade em bacteremia",
    },
    "blaNDM": {
        "pathogen": "Enterobacteriales produtoras de NDM (metalo-beta-lactamase)",
        "threat_level": "CRITICO",
        "prevalence_br": "Crescente desde 2012, surtos em varios estados",
        "treatment_impact": "Resistencia a carbapenems + inibidores de beta-lactamase. Opcoes: aztreonam + ceftazidima-avibactam",
        "isolation": "Precaucao de contato + coorte + vigilancia ativa",
        "mortality_increase": "30-50% mortalidade",
    },
    "vanA": {
        "pathogen": "VRE (Enterococcus resistente a vancomicina)",
        "threat_level": "ALTO",
        "prevalence_br": "Surtos recorrentes em UTIs, especialmente E. faecium",
        "treatment_impact": "Exclui vancomicina. Opcoes: linezolida, daptomicina",
        "isolation": "Precaucao de contato",
        "mortality_increase": "2.5x maior mortalidade vs VSE",
    },
    "mcr-1": {
        "pathogen": "Enterobacteriales resistentes a colistina",
        "threat_level": "CRITICO - Ultimo recurso",
        "prevalence_br": "Emergente, deteccao esporadica em E. coli e Salmonella",
        "treatment_impact": "Colistina e ultimo recurso para CRE. Perda deste antibiotico = sem opcoes",
        "isolation": "Precaucao de contato + notificacao epidemiologica",
        "mortality_increase": "Variavel, mas elimina ultima linha de defesa",
    },
    "blaCTX-M-15": {
        "pathogen": "ESBL (beta-lactamase de espectro estendido)",
        "threat_level": "ALTO",
        "prevalence_br": "Muito comum - principal ESBL no Brasil",
        "treatment_impact": "Resistencia a cefalosporinas 3a geracao. Usar carbapenems",
        "isolation": "Precaucao de contato em surtos",
        "mortality_increase": "1.5-2x maior mortalidade",
    },
    "blaOXA-48": {
        "pathogen": "Carbapenemase OXA-48-like",
        "threat_level": "ALTO",
        "prevalence_br": "Casos importados crescentes, vigilancia necessaria",
        "treatment_impact": "Resistencia variavel a carbapenems. Dificil deteccao fenotipica",
        "isolation": "Precaucao de contato",
        "mortality_increase": "30-40% mortalidade em bacteremia",
    },
    "blaVIM": {
        "pathogen": "P. aeruginosa produtora de VIM (metalo-beta-lactamase)",
        "threat_level": "ALTO",
        "prevalence_br": "Esporadico mas grave em P. aeruginosa",
        "treatment_impact": "Resistencia a carbapenems. Opcoes: ceftolozane-tazobactam (se sensivel), colistina",
        "isolation": "Precaucao de contato",
        "mortality_increase": "30-50% mortalidade",
    },
    "blaIMP": {
        "pathogen": "P. aeruginosa produtora de IMP (metalo-beta-lactamase)",
        "threat_level": "MODERADO",
        "prevalence_br": "Esporadico no Brasil",
        "treatment_impact": "Similar a VIM. Opcoes limitadas",
        "isolation": "Precaucao de contato",
        "mortality_increase": "30-40% mortalidade",
    },
    "blaGES": {
        "pathogen": "Carbapenemase GES (classe A)",
        "threat_level": "EMERGENTE",
        "prevalence_br": "Raro mas emergente no Brasil",
        "treatment_impact": "GES-5 hidrolisa carbapenems. GES-1 apenas cefalosporinas",
        "isolation": "Precaucao de contato",
        "mortality_increase": "Dados limitados",
    },
    "qnrS": {
        "pathogen": "Enterobacteriales com resistencia plasmidial a fluoroquinolonas",
        "threat_level": "MODERADO",
        "prevalence_br": "Comum, disseminacao plasmidial rapida",
        "treatment_impact": "Resistencia a ciprofloxacino, levofloxacino",
        "isolation": "Precaucao padrao",
        "mortality_increase": "Modesto, mas compromete profilaxia cirurgica",
    },
    "armA": {
        "pathogen": "Enterobacteriales com 16S rRNA metiltransferase",
        "threat_level": "ALTO",
        "prevalence_br": "Frequentemente co-localizado com NDM em plasmideos",
        "treatment_impact": "Resistencia a todos aminoglicosideos (gentamicina, amicacina, tobramicina)",
        "isolation": "Precaucao de contato",
        "mortality_increase": "Significativo quando co-ocorre com carbapenemase",
    },
}


def get_anthropic_client():
    """Inicializa cliente Anthropic. Tenta .env se variavel nao existe."""
    api_key = os.environ.get("ANTHROPIC_API_KEY")

    if not api_key:
        env_path = os.path.join(BASE_DIR, ".env")
        if os.path.exists(env_path):
            with open(env_path) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("ANTHROPIC_API_KEY="):
                        api_key = line.split("=", 1)[1].strip().strip('"').strip("'")
                        break

    if not api_key:
        print("[INTERPRETER] AVISO: ANTHROPIC_API_KEY nao encontrada.")
        print("[INTERPRETER] Usando modo offline (interpretacoes pre-computadas).")
        return None

    try:
        import anthropic
        return anthropic.Anthropic(api_key=api_key)
    except ImportError:
        print("[INTERPRETER] AVISO: SDK anthropic nao instalado.")
        print("[INTERPRETER] Instale com: pip install anthropic")
        print("[INTERPRETER] Usando modo offline.")
        return None


def load_pipeline_data() -> dict:
    """Carrega todos os dados do pipeline para interpretacao."""
    data = {"targets": {}, "guides": {}, "primers": {}, "conservation": {}, "card": {}}

    # Alvos
    if os.path.exists(TARGETS_CSV):
        with open(TARGETS_CSV) as f:
            for row in csv.DictReader(f):
                data["targets"][row["name"]] = dict(row)

    # Guides
    for gene_name in TARGETS:
        guide_path = os.path.join(GUIDES_DIR, f"{gene_name}_cas12a_guides.tsv")
        if os.path.exists(guide_path):
            with open(guide_path) as f:
                guides = list(csv.DictReader(f, delimiter="\t"))
                if guides:
                    data["guides"][gene_name] = guides

    # Primers
    for gene_name in TARGETS:
        primer_path = os.path.join(PRIMERS_DIR, f"{gene_name}_rpa_primers.tsv")
        if os.path.exists(primer_path):
            with open(primer_path) as f:
                primers = list(csv.DictReader(f, delimiter="\t"))
                if primers:
                    data["primers"][gene_name] = primers

    # Conservation analysis
    cons_path = os.path.join(REPORTS_DIR, "conservation_analysis.tsv")
    if os.path.exists(cons_path):
        with open(cons_path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                family = row.get("family", "")
                if family not in data["conservation"]:
                    data["conservation"][family] = []
                data["conservation"][family].append(dict(row))

    # CARD enrichment
    card_path = os.path.join(BASE_DIR, "targets_brazil_card.csv")
    if os.path.exists(card_path):
        with open(card_path) as f:
            for row in csv.DictReader(f):
                name = row.get("name", "")
                data["card"][name] = dict(row)

    return data


def build_target_summary(gene_name: str, pipeline_data: dict) -> dict:
    """Constroi resumo estruturado de um alvo para o LLM."""
    target_info = pipeline_data["targets"].get(gene_name, {})
    guides = pipeline_data["guides"].get(gene_name, [])
    primers = pipeline_data["primers"].get(gene_name, [])
    conservation = pipeline_data["conservation"].get(gene_name, [])
    card_data = pipeline_data["card"].get(gene_name, {})
    clinical = CLINICAL_CONTEXT.get(gene_name, {})

    best_guide = guides[0] if guides else None
    best_primer = primers[0] if primers else None

    # Resumo de conservacao
    cons_summary = {"total": 0, "exact": 0, "near": 0, "missed": 0}
    for c in conservation:
        cons_summary["total"] += 1
        match = c.get("match", "")
        if match == "exact":
            cons_summary["exact"] += 1
        elif match == "near":
            cons_summary["near"] += 1
        else:
            cons_summary["missed"] += 1

    return {
        "gene_name": gene_name,
        "target_info": target_info,
        "clinical_context": clinical,
        "best_guide": {
            "spacer": best_guide.get("spacer_seq", "") if best_guide else "",
            "score": best_guide.get("score", "") if best_guide else "",
            "gc": best_guide.get("gc", "") if best_guide else "",
            "position": best_guide.get("position", "") if best_guide else "",
            "strand": best_guide.get("strand", "") if best_guide else "",
        },
        "best_primer": {
            "fwd": best_primer.get("fwd_seq", "") if best_primer else "",
            "rev": best_primer.get("rev_seq", "") if best_primer else "",
            "amplicon_size": best_primer.get("amplicon_size", "") if best_primer else "",
        },
        "total_guides": len(guides),
        "total_primers": len(primers),
        "conservation": cons_summary,
        "card_drug_classes": card_data.get("card_drug_classes", ""),
        "card_resistance_mechanisms": card_data.get("card_resistance_mechanisms", ""),
        "card_variant_count": card_data.get("card_variant_count", ""),
    }


def interpret_with_claude(client, target_summary: dict, audience: str = "medico") -> str:
    """Envia dados estruturados ao Claude e recebe interpretacao clinica."""

    audience_instructions = {
        "medico": (
            "Voce esta explicando para um medico intensivista de UTI no Brasil. "
            "Foque em: impacto clinico imediato, decisao de isolamento, opcoes terapeuticas, "
            "e o que o resultado do teste significa para o paciente. "
            "Use linguagem clinica mas acessivel. Maximo 200 palavras."
        ),
        "gestor": (
            "Voce esta explicando para um gestor hospitalar/diretor de CCIH. "
            "Foque em: custo-beneficio do teste, impacto na taxa de infeccao hospitalar, "
            "comparacao com metodos atuais (cultura, PCR), e ROI. "
            "Use linguagem administrativa. Maximo 200 palavras."
        ),
        "pesquisador": (
            "Voce esta explicando para um bioinformatico/pesquisador. "
            "Inclua detalhes tecnicos: sequencia do spacer, scoring, cobertura de variantes, "
            "parametros do guide, mecanismo de resistencia molecular. "
            "Use linguagem tecnica precisa. Maximo 300 palavras."
        ),
    }

    prompt = f"""Analise os seguintes dados do pipeline CRISPR-Cas12a para diagnostico de resistencia antimicrobiana e gere uma interpretacao clinica.

## Dados do alvo

Gene: {target_summary['gene_name']}
Patogeno: {target_summary['clinical_context'].get('pathogen', 'N/A')}
Nivel de ameaca: {target_summary['clinical_context'].get('threat_level', 'N/A')}
Prevalencia Brasil: {target_summary['clinical_context'].get('prevalence_br', 'N/A')}
Impacto no tratamento: {target_summary['clinical_context'].get('treatment_impact', 'N/A')}
Aumento de mortalidade: {target_summary['clinical_context'].get('mortality_increase', 'N/A')}

## Dados do pipeline computacional

Guide CRISPR (melhor):
- Spacer: {target_summary['best_guide']['spacer']}
- Score: {target_summary['best_guide']['score']}
- GC: {target_summary['best_guide']['gc']}
- Total de guides desenhados: {target_summary['total_guides']}

Primers RPA (melhor par):
- Forward: {target_summary['best_primer']['fwd']}
- Reverse: {target_summary['best_primer']['rev']}
- Amplicon: {target_summary['best_primer']['amplicon_size']} bp

Cobertura de variantes:
- Total analisadas: {target_summary['conservation']['total']}
- Match exato: {target_summary['conservation']['exact']}
- Match proximo (<=3mm): {target_summary['conservation']['near']}
- Nao cobertas: {target_summary['conservation']['missed']}

Dados CARD:
- Drug classes: {target_summary['card_drug_classes']}
- Mecanismos de resistencia: {target_summary['card_resistance_mechanisms']}
- Variantes no CARD: {target_summary['card_variant_count']}

## Instrucoes

{audience_instructions.get(audience, audience_instructions['medico'])}

Gere a interpretacao em portugues brasileiro. Estruture com:
1. **Resultado**: O que este teste detecta e por que importa
2. **Cobertura**: Quao confiavel e a deteccao (baseado nos dados de variantes)
3. **Acao clinica**: O que fazer com um resultado positivo
"""

    try:
        response = client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=1024,
            messages=[{"role": "user", "content": prompt}],
        )
        return response.content[0].text
    except Exception as e:
        print(f"[INTERPRETER] Erro na API Claude: {e}")
        return None


def generate_offline_interpretation(target_summary: dict, audience: str = "medico") -> str:
    """Gera interpretacao pre-computada sem Claude API (modo offline)."""
    gene = target_summary["gene_name"]
    clinical = target_summary["clinical_context"]
    guide = target_summary["best_guide"]
    cons = target_summary["conservation"]

    # Calcular cobertura
    total = cons["total"]
    exact = cons["exact"]
    near = cons["near"]
    if total > 0:
        coverage_pct = ((exact + near) / total) * 100
    else:
        coverage_pct = 0

    pathogen = clinical.get("pathogen", gene)
    threat = clinical.get("threat_level", "N/A")
    treatment = clinical.get("treatment_impact", "N/A")
    prevalence = clinical.get("prevalence_br", "N/A")
    mortality = clinical.get("mortality_increase", "N/A")
    isolation = clinical.get("isolation", "Precaucao de contato")

    if audience == "medico":
        text = (
            f"**Resultado**: Este teste detecta o gene {gene}, marcador de {pathogen}. "
            f"Nivel de ameaca: {threat}. {prevalence}.\n\n"
            f"**Cobertura**: O guide CRISPR projetado (score {guide['score']}, GC {guide['gc']}) "
        )
        if total > 0:
            text += f"cobre {exact}/{total} variantes com match exato ({coverage_pct:.0f}% com tolerancia <=3 mismatches). "
        else:
            text += "ainda nao foi testado contra variantes clinicas. "
        text += (
            f"\n\n**Acao clinica**: Resultado POSITIVO indica colonizacao/infeccao por {pathogen}. "
            f"{treatment}. "
            f"Mortalidade: {mortality}. "
            f"Isolamento: {isolation}."
        )

    elif audience == "gestor":
        text = (
            f"**O teste**: Detecta {pathogen} em ~30 min por ~R$25/teste "
            f"(vs R$200+ GeneXpert, 2-5 dias cultura).\n\n"
            f"**Impacto**: {prevalence}. "
            f"Deteccao precoce permite isolamento imediato ({isolation}), "
            f"reduzindo transmissao cruzada em UTI. "
            f"Mortalidade associada: {mortality}.\n\n"
            f"**Confiabilidade**: Pipeline computacional validado com BLAST 100% identity. "
        )
        if total > 0:
            text += f"Cobertura de {coverage_pct:.0f}% das variantes clinicas conhecidas."
        else:
            text += "Validacao de variantes em andamento."

    else:  # pesquisador
        text = (
            f"**Gene alvo**: {gene} | {pathogen}\n"
            f"**Nivel de ameaca**: {threat}\n\n"
            f"**Design do guide**: Spacer 5'-{guide['spacer']}-3' | "
            f"Score={guide['score']} | GC={guide['gc']} | "
            f"Posicao={guide['position']} ({guide['strand']} strand) | "
            f"{target_summary['total_guides']} guides candidatos\n\n"
            f"**Primers RPA**: FWD 5'-{target_summary['best_primer']['fwd']}-3' | "
            f"REV 5'-{target_summary['best_primer']['rev']}-3' | "
            f"Amplicon: {target_summary['best_primer']['amplicon_size']}bp\n\n"
        )
        if total > 0:
            text += (
                f"**Conservacao**: {exact}/{total} exact match, "
                f"{near}/{total} near (<=3mm), "
                f"{cons['missed']}/{total} missed | "
                f"Cobertura efetiva: {coverage_pct:.0f}%\n\n"
            )
        text += (
            f"**CARD**: Drug classes: {target_summary['card_drug_classes'] or 'N/A'} | "
            f"Mecanismos: {target_summary['card_resistance_mechanisms'] or 'N/A'} | "
            f"Variantes CARD: {target_summary['card_variant_count'] or 'N/A'}\n\n"
            f"**Impacto clinico**: {treatment}"
        )

    return text


def generate_panel_interpretation(pipeline_data: dict, client=None, audience: str = "medico") -> dict:
    """Gera interpretacao para todo o painel."""
    interpretations = {}

    # Determinar quais genes processar
    genes = list(TARGETS.keys())[:12]  # 12 familias base

    for gene_name in genes:
        summary = build_target_summary(gene_name, pipeline_data)

        if client:
            text = interpret_with_claude(client, summary, audience)
            if text:
                interpretations[gene_name] = text
                time.sleep(0.5)  # Rate limiting
                continue

        # Fallback para offline
        interpretations[gene_name] = generate_offline_interpretation(summary, audience)

    return interpretations


def save_clinical_report(interpretations: dict, audience: str, pipeline_data: dict):
    """Salva relatorio clinico completo."""
    os.makedirs(REPORTS_CLINICAL_DIR, exist_ok=True)

    report_path = os.path.join(REPORTS_CLINICAL_DIR, f"relatorio_clinico_{audience}.md")

    with open(report_path, "w", encoding="utf-8") as f:
        f.write(f"# Relatorio Clinico - SmartLab BacEnd\n\n")
        f.write(f"**Audiencia**: {audience.capitalize()}\n")
        f.write(f"**Gerado em**: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"**Pipeline**: CRISPR-Cas12a + RPA para deteccao de AMR\n\n")
        f.write("---\n\n")

        # Indice
        f.write("## Indice de alvos\n\n")
        for gene_name in interpretations:
            target = pipeline_data["targets"].get(gene_name, {})
            pathogen = target.get("pathogen", "")
            priority = target.get("priority", "")
            f.write(f"- **{gene_name}** ({pathogen}) - {priority}\n")
        f.write("\n---\n\n")

        # Interpretacoes
        for gene_name, text in interpretations.items():
            target = pipeline_data["targets"].get(gene_name, {})
            clinical = CLINICAL_CONTEXT.get(gene_name, {})

            f.write(f"## {gene_name} - {clinical.get('pathogen', target.get('pathogen', ''))}\n\n")
            f.write(f"**Prioridade**: {target.get('priority', 'N/A')} | ")
            f.write(f"**Ameaca**: {clinical.get('threat_level', 'N/A')}\n\n")
            f.write(text)
            f.write("\n\n---\n\n")

        # Disclaimer
        f.write("## Aviso importante\n\n")
        f.write(
            "Este relatorio foi gerado computacionalmente e NAO substitui avaliacao "
            "clinica profissional. Os dados de cobertura e especificidade sao baseados "
            "em validacao in silico (BLAST). Validacao clinica com amostras reais e "
            "pendente. Para uso exclusivo em pesquisa.\n"
        )

    print(f"[INTERPRETER] Relatorio salvo: {report_path}")
    return report_path


def save_interpretations_json(interpretations: dict, audience: str):
    """Salva interpretacoes em JSON para uso pelo frontend."""
    os.makedirs(REPORTS_CLINICAL_DIR, exist_ok=True)
    json_path = os.path.join(REPORTS_CLINICAL_DIR, f"interpretations_{audience}.json")

    output = {}
    for gene_name, text in interpretations.items():
        clinical = CLINICAL_CONTEXT.get(gene_name, {})
        output[gene_name] = {
            "interpretation": text,
            "threat_level": clinical.get("threat_level", ""),
            "pathogen": clinical.get("pathogen", ""),
            "isolation": clinical.get("isolation", ""),
            "generated_at": datetime.now().isoformat(),
        }

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False, indent=2)

    print(f"[INTERPRETER] JSON salvo: {json_path}")
    return json_path


def run_clinical_interpreter(audiences: list[str] = None, force_offline: bool = False):
    """Executa pipeline de interpretacao clinica."""
    print("=" * 70)
    print("INTERPRETADOR CLINICO - SmartLab BacEnd")
    print("Inspirado em EVEE (Goodfire/Mayo Clinic)")
    print("=" * 70)

    if audiences is None:
        audiences = ["medico", "gestor", "pesquisador"]

    # Inicializar cliente
    client = None
    if not force_offline:
        client = get_anthropic_client()

    mode = "Claude API" if client else "Offline (pre-computado)"
    print(f"\n[INTERPRETER] Modo: {mode}")

    # Carregar dados do pipeline
    print("[INTERPRETER] Carregando dados do pipeline...")
    pipeline_data = load_pipeline_data()

    targets_count = len(pipeline_data["targets"])
    guides_count = len(pipeline_data["guides"])
    primers_count = len(pipeline_data["primers"])
    cons_count = len(pipeline_data["conservation"])
    card_count = len(pipeline_data["card"])

    print(f"  Alvos: {targets_count}")
    print(f"  Guides: {guides_count} genes")
    print(f"  Primers: {primers_count} genes")
    print(f"  Conservation: {cons_count} familias")
    print(f"  CARD enrichment: {card_count} entradas")

    # Gerar interpretacoes para cada audiencia
    all_reports = {}
    for audience in audiences:
        print(f"\n{'─'*50}")
        print(f"[INTERPRETER] Gerando interpretacoes para: {audience}")

        interpretations = generate_panel_interpretation(pipeline_data, client, audience)

        if interpretations:
            report_path = save_clinical_report(interpretations, audience, pipeline_data)
            json_path = save_interpretations_json(interpretations, audience)
            all_reports[audience] = {
                "report": report_path,
                "json": json_path,
                "count": len(interpretations),
            }
            print(f"  Interpretacoes geradas: {len(interpretations)}")
        else:
            print(f"  AVISO: Nenhuma interpretacao gerada para {audience}")

    # Resumo
    print(f"\n{'='*70}")
    print("INTERPRETACAO CLINICA CONCLUIDA")
    print(f"{'='*70}")
    for audience, info in all_reports.items():
        print(f"  [{audience}] {info['count']} interpretacoes -> {info['report']}")

    return all_reports


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Interpretador clinico SmartLab BacEnd")
    parser.add_argument("--audience", choices=["medico", "gestor", "pesquisador", "all"],
                        default="all", help="Audiencia do relatorio")
    parser.add_argument("--offline", action="store_true", help="Forcar modo offline (sem Claude API)")
    args = parser.parse_args()

    audiences = ["medico", "gestor", "pesquisador"] if args.audience == "all" else [args.audience]
    run_clinical_interpreter(audiences=audiences, force_offline=args.offline)
