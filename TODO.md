# TODO - SmartLab BacEnd

Tarefas pendentes e em andamento do projeto.
Ultima atualizacao: 2026-04-14

---

## Concluido

- [x] Pipeline base: fetch, guides, primers, BLAST, panel (5 etapas)
- [x] 12 familias genicas de referencia processadas com BLAST 100%
- [x] 30 variantes clinicas processadas
- [x] Conservation analysis para cobertura de variantes
- [x] Sistema de tracking CSV
- [x] Batch mode (run_batch.py)
- [x] Site vitrine (public/index.html)
- [x] Documentacao: README.md, TODO.md, MELHORIAS.md

---

## Em andamento

### Fase 1 - Integracao CARD Database [CONCLUIDA]
- [x] Criar modulo `card_integration.py` para buscar dados do CARD (Comprehensive Antibiotic Resistance Database)
- [x] Mapear variantes CARD -> alvos existentes do pipeline (12 familias)
- [x] Auto-descoberta de novas variantes AMR relevantes para o Brasil
- [x] Enriquecer dados com metadados CARD (mecanismo de resistencia, drug class, ontologia ARO)
- [x] Gerar CSV enriquecido (`targets_brazil_card.csv`) e relatorio de descobertas
- [x] Atualizar site com secao CARD

### Fase 2 - Explicacoes em Linguagem Natural (Claude API) [CONCLUIDA]
- [x] Criar modulo `clinical_interpreter.py` com integracao Claude API
- [x] Gerar interpretacao clinica para cada alvo/guide/primer (3 audiencias: medico, gestor, pesquisador)
- [x] Contexto clinico embutido para 12 familias AMR (mortalidade, isolamento, tratamento, prevalencia BR)
- [x] Modo online (Claude API) + modo offline (pre-computado)
- [x] Relatorio clinico em .md + JSON para frontend
- [x] Atualizar site com secao "Interpretacao Clinica por IA"

### Fase 3 - Score Funcional via Evo 2 [CONCLUIDA]
- [x] Criar modulo `evo2_scoring.py` com modo GPU (Evo 2 real) e lightweight (features de sequencia)
- [x] 5 dimensoes de analise: k-mer profile, codon adaptation, GC shift, dinucleotideo, complexidade
- [x] Calcular distancia funcional composta entre variantes e referencia
- [x] Classificacao de impacto: conserved, likely_conserved, uncertain, likely_disrupted, disrupted
- [x] Priorizacao dinamica baseada em impacto funcional
- [x] Outputs: TSV, JSON, relatorio textual
- [x] Atualizar site com secao "Score Funcional de Variantes"

### Fase 4 - Dashboard Interativo [CONCLUIDA]
- [x] Criar `public/dashboard.html` - pagina interativa completa
- [x] Tabela sortavel com 8 colunas (gene, classe, prioridade, organismo, variantes, BLAST, impacto, acao)
- [x] Filtros por prioridade (P1/P2/P3), classe de resistencia e busca por texto
- [x] Painel de detalhes expandivel com 3 colunas: interpretacao clinica, dados pipeline, dados CARD
- [x] Mapa visual do gene mostrando posicao do guide e primers
- [x] Barra de impacto funcional por alvo
- [x] Exportacao CSV
- [x] Botao "Abrir Dashboard" no hero do site principal
- [x] Navegacao entre dashboard e site principal

### Fase 5 - Covariance Probes para Design de Guias [CONCLUIDA]
- [x] Criar modulo `covariance_probes.py` com 18 features biofisicas
- [x] Features novas: seed region (GC, purinas, homopolymer), termodinamica (dG nearest-neighbor, energia seed, gradiente estabilidade), acessibilidade estrutural (palindromes, repeticoes invertidas), preferencias posicionais (dados experimentais Cas12a)
- [x] Implementar covariance matrix sobre todos os candidatos
- [x] Score composto: w^T*f + alpha * f^T*C*w (linear + covariancia)
- [x] Benchmark automatico contra scoring rule-based original
- [x] Outputs: TSV, JSON, relatorio de benchmark detalhado
- [x] Atualizar site com secao "Covariance Probes" + pipeline atualizado para 7 etapas

---

## Backlog (futuro)

- [ ] Integracao com dados epidemiologicos em tempo real (ANVISA/BR-GLASS API)
- [ ] Suporte a novos patogenos alem de AMR (ex: virus respiratorios)
- [ ] App mobile para leitura de fluorescencia via camera
- [ ] Validacao clinica com amostras reais (parceria HC-FMUSP)
- [ ] Submissao regulatoria ANVISA (IVD, RDC 830/2023)
- [ ] Multi-idioma (ingles/espanhol) para adocao internacional
