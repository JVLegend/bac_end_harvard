# MELHORIAS - Roadmap de Evolucao do BacEnd

Documento de melhorias planejadas, inspirado na pesquisa Goodfire/Mayo Clinic (EVEE - Evo Variant Effect Explorer, abril 2026) e no estado da arte em genomica computacional.

Ultima atualizacao: 2026-04-14

---

## Contexto: O que aprendemos com o EVEE (Goodfire + Mayo Clinic)

A Goodfire publicou o EVEE, que usa o modelo genomico fundacional **Evo 2** (Arc Institute) para prever patogenicidade de 4.2M variantes geneticas do ClinVar com **0.997 AUROC**. Principais inovacoes:

1. **Covariance probes** sobre embeddings do Evo 2 (captura co-ocorrencia de features)
2. **Perfis de disrupcao funcional** (splice sites, elementos regulatorios, dominios estruturais)
3. **Explicacoes em linguagem natural** via LLM (Claude Sonnet) a partir dos perfis
4. **Database aberta** com scores para todas as variantes ClinVar

**Como isso se aplica ao BacEnd**: Embora o EVEE foque em variantes humanas, o Evo 2 opera sobre DNA de todos os dominios da vida. A metodologia pode ser adaptada para prever impacto funcional de variantes de genes AMR bacterianos.

---

## Fase 1: Integracao CARD Database

**Prioridade**: ALTA | **Esforco**: Medio | **Impacto**: Alto

### Problema atual
- Variantes curadas manualmente em CSVs estaticos
- Novas variantes AMR nao sao descobertas automaticamente
- Falta metadados de mecanismo de resistencia

### Solucao
Integrar o **CARD** (Comprehensive Antibiotic Resistance Database) como fonte primaria:
- API/download automatico de variantes AMR
- Metadados ricos: mecanismo de resistencia (efflux, target modification, inactivation), organismo, prevalencia
- Anotacoes ARO (Antibiotic Resistance Ontology)
- Atualizacao periodica do catalogo de alvos

### Entregaveis
- `card_integration.py` - modulo de integracao
- `targets_brazil_variants.csv` enriquecido com dados CARD
- Auto-descoberta de novas variantes relevantes
- Secao no site mostrando fontes de dados

### Referencia
- CARD: https://card.mcmaster.ca/
- RGI (Resistance Gene Identifier): ferramenta de anotacao do CARD

---

## Fase 2: Explicacoes em Linguagem Natural (Claude API)

**Prioridade**: ALTA | **Esforco**: Baixo | **Impacto**: Alto

### Problema atual
- Reports sao tecnicos: TSVs com scores, sequencias, posicoes BLAST
- Medicos e gestores hospitalares nao conseguem interpretar diretamente
- Sem contexto clinico nos resultados

### Solucao
Usar a **Claude API** para gerar interpretacoes clinicas, inspirado na abordagem EVEE:

1. **Input**: Dados do pipeline (guide score, BLAST results, conservation analysis)
2. **Processamento**: Claude recebe contexto estruturado + prompt clinico
3. **Output**: Interpretacao em portugues para diferentes audiencias:
   - **Medico**: "O guide para mecA cobre 100% das variantes MRSA prevalentes no Brasil. Sensibilidade estimada: 10-100 copias/uL."
   - **Gestor**: "Este teste detecta as 3 bacterias resistentes mais comuns em UTIs do SUS por R$25/teste."
   - **Pesquisador**: "Spacer posicao 847-867 no gene mecA, 0 mismatches com variantes mecA1/mecA2. GC=55%, score=92.3."

### Entregaveis
- `clinical_interpreter.py` - modulo de interpretacao via Claude API
- Relatorio clinico automatico em portugues
- Secao no site com interpretacoes geradas por IA
- Perfis de disrupcao funcional (estilo EVEE simplificado)

### Referencia
- EVEE explanation pipeline: disruption profiles -> natural language via LLM
- Claude API: https://docs.anthropic.com/

---

## Fase 3: Score Funcional via Evo 2

**Prioridade**: ALTA | **Esforco**: Alto | **Impacto**: Muito alto

### Problema atual
- Conservation analysis compara por mismatch count (0-5) - puramente baseado em sequencia
- Nao distingue entre mutacao silenciosa e mutacao que altera funcao de resistencia
- Priorizacao (P1/P2/P3) e manual e estatica

### Solucao
Integrar o **Evo 2** (modelo fundacional genomico da Arc Institute, open source):

1. **Embeddings**: Gerar representacoes densas de cada variante AMR
2. **Distancia funcional**: Comparar variantes por distancia no espaco de embeddings (nao por mismatches)
3. **Score de impacto**: Prever se variante mantem/perde funcao de resistencia
4. **Priorizacao dinamica**: Re-ranquear alvos baseado em impacto funcional real

### Entregaveis
- `evo2_scoring.py` - modulo de scoring funcional
- Score de impacto funcional por variante
- Priorizacao automatica de alvos
- Comparacao: scoring atual (rule-based) vs Evo 2 (embedding-based)

### Referencia
- Evo 2: https://github.com/ArcInstitute/evo2
- Preprint Goodfire: biorxiv.org/content/10.64898/2026.04.10.717844v1

---

## Fase 4: Dashboard Interativo

**Prioridade**: MEDIA | **Esforco**: Medio | **Impacto**: Alto

### Problema atual
- Site e uma vitrine estatica
- Sem interatividade com os dados do pipeline
- Sem visualizacao da posicao dos guides/primers no gene

### Solucao
Dashboard interativo inspirado no EVEE Explorer:

1. **Mapa genico**: Visualizacao do gene com dominios funcionais, posicao do guide e primers
2. **Tabela interativa**: Filtros por familia, prioridade, mecanismo, score
3. **Perfil por variante**: Score de cobertura, mismatches, interpretacao IA
4. **Comparativo**: Antes/depois da integracao de cada melhoria

### Entregaveis
- Dashboard em HTML/JS (sem framework pesado, manter simplicidade)
- Visualizacao de genes com anotacoes
- Filtros e busca interativa
- Exportacao de relatorios

---

## Fase 5: Covariance Probes para Design de Guias

**Prioridade**: MEDIA | **Esforco**: Muito alto | **Impacto**: Diferencial competitivo

### Problema atual
- Guide scoring e rule-based (GC%, homopolimeros, self-complementarity)
- Nao considera contexto genomico alem da sequencia local
- Nao preve eficacia real de clivagem

### Solucao
Adaptar a abordagem de **covariance probes** da Goodfire:

1. **Covariance pooling**: Substituir mean pooling por covariance-based pooling nos embeddings
2. **Probes de eficacia**: Treinar sondas que preveem eficacia de clivagem Cas12a
3. **Contexto genomico**: Considerar estrutura secundaria, acessibilidade do alvo
4. **Benchmark**: Comparar com scoring atual usando dados experimentais publicados

### Entregaveis
- `covariance_probes.py` - implementacao de probes
- Novo scoring de guides baseado em embeddings
- Benchmark comparativo
- Publicacao dos resultados

### Referencia
- Covariance-based Sequence Pooling (Goodfire Research, abril 2026)
- Dados experimentais de eficacia Cas12a: PubMed

---

## Tabela Resumo

| Fase | Melhoria | Esforco | Impacto | Status |
|------|----------|---------|---------|--------|
| 1 | Integracao CARD database | Medio | Alto | **Concluida** |
| 2 | Explicacoes via Claude API | Baixo | Alto | **Concluida** |
| 3 | Score funcional Evo 2 | Alto | Muito alto | **Concluida** |
| 4 | Dashboard interativo | Medio | Alto | **Concluida** |
| 5 | Covariance probes | Muito alto | Diferencial | **Concluida** |

---

## Principios de implementacao

1. **Modular**: Cada fase e um modulo independente que se integra ao pipeline existente
2. **Incremental**: Cada fase entrega valor isoladamente
3. **Backward-compatible**: Pipeline existente continua funcionando sem as melhorias
4. **Testavel**: Cada modulo pode ser validado independentemente
5. **Open-source**: Todo codigo sera aberto no GitHub
