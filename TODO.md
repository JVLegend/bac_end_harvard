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

### Fase 1 - Integracao CARD Database
- [ ] Criar modulo `card_integration.py` para buscar dados do CARD (Comprehensive Antibiotic Resistance Database)
- [ ] Mapear variantes CARD -> alvos existentes do pipeline
- [ ] Auto-descoberta de novas variantes AMR relevantes para o Brasil
- [ ] Enriquecer `targets_brazil_variants.csv` com metadados CARD (mecanismo de resistencia, organismo, etc.)
- [ ] Atualizar site com secao CARD

### Fase 2 - Explicacoes em Linguagem Natural (Claude API)
- [ ] Criar modulo `clinical_interpreter.py` com integracao Claude API
- [ ] Gerar interpretacao clinica para cada alvo/guide/primer
- [ ] Perfil de disruptacao funcional por variante (inspirado no EVEE/Goodfire)
- [ ] Relatorio clinico em portugues para medicos nao-bioinformaticos
- [ ] Atualizar site com secao de interpretacao IA

### Fase 3 - Score Funcional via Evo 2
- [ ] Integrar modelo Evo 2 (Arc Institute) para embeddings genomicos
- [ ] Calcular distancia funcional entre variantes (nao apenas sequencial)
- [ ] Score de patogenicidade/virulencia computacional por variante
- [ ] Priorizacao dinamica de alvos baseada em impacto funcional
- [ ] Atualizar site com scores funcionais

### Fase 4 - Dashboard Interativo
- [ ] Mapa de gene com dominios funcionais anotados
- [ ] Posicao do guide/primer no contexto funcional
- [ ] Score de cobertura por variante com indicador de confianca
- [ ] Visualizacao estilo EVEE para cada alvo
- [ ] Filtros por familia genica, prioridade, mecanismo

### Fase 5 - Covariance Probes para Design de Guias
- [ ] Implementar covariance-based sequence pooling (vs mean pooling)
- [ ] Treinar probes sobre embeddings Evo 2 para prever eficacia de clivagem
- [ ] Considerar estrutura secundaria do DNA alvo
- [ ] Benchmark contra scoring rule-based atual

---

## Backlog (futuro)

- [ ] Integracao com dados epidemiologicos em tempo real (ANVISA/BR-GLASS API)
- [ ] Suporte a novos patogenos alem de AMR (ex: virus respiratorios)
- [ ] App mobile para leitura de fluorescencia via camera
- [ ] Validacao clinica com amostras reais (parceria HC-FMUSP)
- [ ] Submissao regulatoria ANVISA (IVD, RDC 830/2023)
- [ ] Multi-idioma (ingles/espanhol) para adocao internacional
