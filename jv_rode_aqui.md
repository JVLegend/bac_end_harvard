# JV - Rode Aqui no Novo Computador

Guia completo para rodar o pipeline do zero neste computador.

---

## 1. Pre-requisitos

- Python 3.10+ (recomendado) ou 3.8+
- Internet (para NCBI e CARD)
- Opcional: ANTHROPIC_API_KEY para interpretacoes via Claude

---

## 2. Instalar dependencias

```bash
pip install -r requirements.txt
```

Instala: `biopython`, `requests`, `anthropic`

---

## 3. (Opcional) Configurar Claude API

Crie um arquivo `.env` na raiz do projeto:

```
ANTHROPIC_API_KEY=sk-ant-...sua-chave-aqui...
```

Se nao tiver a chave, tudo funciona no modo offline automaticamente.

---

## 4. Rodar o pipeline completo

### Opcao A - Tudo de uma vez (recomendado)

```bash
# Pipeline base (12 genes, sem BLAST para ser rapido)
python run_batch.py 12 --skip-blast

# Depois rodar os modulos novos:
python card_integration.py
python covariance_probes.py
python evo2_scoring.py
python clinical_interpreter.py --offline
```

### Opcao B - Passo a passo (mais controle)

```bash
# 1. Baixar sequencias dos 12 genes AMR do NCBI (~2-5 min)
python fetch_sequences.py

# 2. Projetar guides CRISPR-Cas12a (~30 seg)
python design_guides.py

# 3. Projetar primers RPA (~30 seg)
python design_primers.py

# 4. Validar especificidade via BLAST (~20-40 min, pode pular)
python specificity_check.py

# 5. Montar painel final + folha de oligos (~5 seg)
python multiplex_panel.py

# 6. Analisar cobertura de variantes (~5 min, requer internet)
python conservation_analysis.py
```

### Opcao C - Com BLAST completo (mais lento, mais completo)

```bash
python run_batch.py 12
```

---

## 5. Rodar os modulos novos (melhorias)

Rodar APOS o pipeline base acima ter terminado:

```bash
# Fase 1 - Integracao CARD (~5-10 min, baixa ~500MB)
python card_integration.py

# Fase 3 - Score funcional de variantes (~5 min, CPU puro)
python evo2_scoring.py

# Fase 5 - Covariance probes para guides (~1 min, CPU puro)
python covariance_probes.py

# Fase 2 - Interpretacoes clinicas
python clinical_interpreter.py --offline        # sem API key
python clinical_interpreter.py                  # com API key no .env
```

---

## 6. Ver os resultados

Todos os outputs ficam em:

```
reports/
├── panel_report.txt              # relatorio final do painel CRISPR
├── oligo_order.tsv               # folha de encomenda de oligos
├── specificity_report.tsv        # resultados BLAST
├── conservation_analysis.tsv     # cobertura de variantes
├── card_new_variants.csv         # novas variantes descobertas via CARD
├── card_integration_report.txt   # relatorio CARD completo
├── clinical/
│   ├── relatorio_clinico_medico.md
│   ├── relatorio_clinico_gestor.md
│   ├── relatorio_clinico_pesquisador.md
│   └── interpretations_*.json
├── evo2_scoring/
│   ├── functional_scores.tsv
│   ├── functional_scores.json
│   └── functional_analysis_report.txt
└── covariance_probes/
    ├── probe_scores.tsv
    ├── probe_scores.json
    └── benchmark_report.txt

targets_brazil_card.csv           # CSV principal enriquecido com CARD
```

---

## 7. Abrir o site/dashboard

Abra no navegador (sem servidor necessario):

```
public/index.html      # site vitrine
public/dashboard.html  # dashboard interativo com filtros
```

---

## 8. Resumo dos modulos do projeto

| Script | O que faz | Precisa de internet? |
|--------|-----------|----------------------|
| `fetch_sequences.py` | Baixa sequencias NCBI | Sim |
| `design_guides.py` | Projeta crRNA Cas12a | Nao |
| `design_primers.py` | Projeta primers RPA | Nao |
| `specificity_check.py` | Valida via BLAST | Sim |
| `multiplex_panel.py` | Monta painel final | Nao |
| `conservation_analysis.py` | Cobertura de variantes | Sim |
| `card_integration.py` | Integra CARD database | Sim (1a vez) |
| `evo2_scoring.py` | Score funcional variantes | Sim (fetch NCBI) |
| `covariance_probes.py` | Scoring avancado guides | Nao |
| `clinical_interpreter.py` | Interpretacao clinica IA | Opcional (Claude API) |

---

## 9. Problemas comuns

**NCBI lento ou com erro de rate limit:**
```bash
# Adicionar sua chave de email no config.py
# NCBI_EMAIL = "seu@email.com"
```

**CARD demora muito no download:**
```bash
# Usar cache forcado se ja baixou antes
python card_integration.py --force
```

**Erro de importacao:**
```bash
pip install --upgrade biopython requests anthropic
```

**Python 3.8 - erro de type hints (`str | None`):**
```bash
# Atualizar para Python 3.10+ recomendado
# Ou instalar: pip install typing_extensions
```

---

## 10. Estrutura do projeto

```
bac_end_harvard/
├── config.py                   # parametros centrais
├── utils.py                    # funcoes bioinformaticas
├── fetch_sequences.py          # etapa 1: NCBI
├── design_guides.py            # etapa 2: crRNA
├── design_primers.py           # etapa 4: RPA primers
├── specificity_check.py        # etapa 5: BLAST
├── multiplex_panel.py          # etapa 7: painel final
├── conservation_analysis.py    # analise de variantes
├── tracking.py                 # tracking de progresso
├── run_batch.py                # orquestrador batch
├── card_integration.py         # NOVO: integracao CARD
├── evo2_scoring.py             # NOVO: score funcional
├── covariance_probes.py        # NOVO: probes avancadas
├── clinical_interpreter.py     # NOVO: interpretacao IA
├── targets_brazil.csv          # 12 genes de referencia
├── targets_brazil_variants.csv # 47 variantes clinicas
├── requirements.txt            # pip install -r requirements.txt
├── public/
│   ├── index.html              # site vitrine
│   └── dashboard.html          # dashboard interativo
├── README.md
├── TODO.md
├── MELHORIAS.md
└── jv_rode_aqui.md             # este arquivo
```
