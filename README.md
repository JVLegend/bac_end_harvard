# SmartLab BacEnd - CRISPR-Cas12a Diagnostic Pipeline for AMR Detection

Pipeline computacional de IA agentica para design automatizado de diagnosticos CRISPR-Cas12a, focado em resistencia antimicrobiana (AMR) em UTIs brasileiras.

**Projeto**: IA para Medicos / SmartLab BacEnd
**Instituicoes**: HC-FMUSP, Harvard Medical School (Broad Institute)
**Status**: Pesquisa ativa - validacao clinica pendente

---

## O que faz

Um profissional de saude coleta amostra fecal, aplica em dispositivo de papel, e obtem resultado fluorescente em **~30 minutos**. O dispositivo detecta genes de resistencia antimicrobiana sem necessidade de termociclador.

```
Amostra -> Lise (fibra de vidro) -> Extracao (fluxo lateral em papel)
-> Amplificacao RPA (37C, 20 min) -> Trans-clivagem Cas12a -> Leitura fluorescente
```

### Numeros do projeto

| Metrica | Valor |
|---------|-------|
| Alvos AMR processados | 42 (12 familias + 30 variantes) |
| BLAST identity (referencias) | 100% em todas as 12 familias |
| Custo por teste | ~R$25 (vs R$200+ GeneXpert) |
| Tempo de resultado | ~30 minutos |
| Bases de dados | NCBI RefSeq, AMRFinderPlus, BLAST nt |

---

## Painel de Alvos

| Familia genica | Classe de resistencia | Variantes cobertas | Relevancia clinica |
|----------------|-----------------------|---------------------|--------------------|
| **mecA** | Meticilina (beta-lactamicos) | mecA, mecA1, mecA2 | MRSA - principal HAI em UTIs |
| **blaKPC** | Carbapenems | KPC-2, KPC-3, KPC-4, KPC-5, KPC-11, KPC-30, KPC-31 | CRE - prioridade critica OMS |
| **blaNDM** | Carbapenems (MBL) | NDM-1, NDM-2, NDM-5, NDM-7 | Crescente no Brasil desde 2012 |
| **vanA** | Vancomicina | vanA (2 seq. ref.) | VRE - surtos recorrentes |
| **mcr** | Colistina | mcr-1, mcr-1.1, mcr-5 | Ultimo recurso terapeutico |
| **blaCTX-M** | Cefalosporinas (ESBL) | CTX-M-2, 8, 9, 14, 15, 27 | ESBL mais comum no Brasil |
| **blaOXA-48** | Carbapenems (classe D) | OXA-48, OXA-181, OXA-232 | Casos importados crescentes |
| **blaVIM** | Carbapenems (MBL) | VIM-1, VIM-2, VIM-4 | P. aeruginosa resistente |
| **blaIMP** | Carbapenems (MBL) | IMP-1, IMP-6 | Esporadico em P. aeruginosa |
| **blaGES** | Carbapenems (classe A) | GES-1, GES-5 | Emergente |
| **qnrS** | Fluoroquinolonas | qnrS1, qnrS2 | Resistencia por plasmideo |
| **armA** | Aminoglicosideos | armA (2 seq. ref.) | Co-localizado com NDM |

---

## Como usar

### Pre-requisitos

```bash
pip install -r requirements.txt
```

Dependencias: `biopython >= 1.83`, `requests >= 2.31.0`

### Execucao sequencial

```bash
# 1. Buscar sequencias DNA dos genes AMR no NCBI
python fetch_sequences.py

# 2. Projetar crRNA guides Cas12a (scan PAM TTTV + scoring)
python design_guides.py

# 3. Projetar primers RPA (30-35nt, amplicon 100-200bp)
python design_primers.py

# 4. Validar especificidade via BLAST (requer internet)
python specificity_check.py

# 5. Montar painel final + folha de encomenda de oligos
python multiplex_panel.py
```

### Execucao em batch

```bash
# Processar 4 genes de uma vez
python run_batch.py 4

# Pular BLAST (mais rapido, para iteracao)
python run_batch.py 4 --skip-blast
```

### Analise de conservacao de variantes

```bash
python conservation_analysis.py
```

---

## Arquitetura do Pipeline

```
config.py (parametros centrais: alvos, Cas12a, RPA)
    |
utils.py (funcoes compartilhadas: rev_comp, GC%, Tm, PAM scan)
    |
    v
fetch_sequences.py --> design_guides.py --> design_primers.py
                                                    |
                                                    v
                                        specificity_check.py --> multiplex_panel.py
                                                    |
                                        conservation_analysis.py
```

### Estrutura de diretorios

```
bac_end_harvard/
|-- config.py                    # Definicoes de alvos, parametros Cas12a/RPA
|-- utils.py                     # Funcoes bioinformaticas compartilhadas
|-- fetch_sequences.py           # Download de sequencias NCBI Entrez
|-- design_guides.py             # Design de crRNA (scan PAM + scoring)
|-- design_primers.py            # Design de primers RPA
|-- specificity_check.py         # Validacao BLAST
|-- multiplex_panel.py           # Montagem do painel + folha de oligos
|-- conservation_analysis.py     # Analise de cobertura de variantes
|-- tracking.py                  # Sistema de tracking do pipeline
|-- run_batch.py                 # Orquestrador batch
|-- targets_brazil.csv           # 12 alvos de referencia
|-- targets_brazil_variants.csv  # 47 variantes clinicas
|-- requirements.txt             # Dependencias Python
|-- public/index.html            # Frontend/vitrine do projeto
|-- TODO.md                      # Tarefas pendentes
|-- MELHORIAS.md                 # Roadmap de melhorias futuras
|-- sequences/                   # Sequencias FASTA baixadas
|-- guides/                      # crRNA guides projetados (.tsv)
|-- primers/                     # Primers RPA projetados (.tsv)
|-- reports/                     # Relatorios finais e folhas de oligos
```

---

## Tecnologia CRISPR-Cas12a

- **Enzima**: LbCas12a / AsCas12a
- **PAM**: TTTV (T-T-T-[A/C/G]) na fita nao-alvo, 5' do spacer
- **Spacer**: 20-24 nt (otimo: 20 nt)
- **GC ideal**: 40-60%
- **Amplificacao**: RPA isotermica a 37C, 20 min, amplicon 100-200bp
- **Reporter**: ssDNA-FQ (5'-6-FAM / 3'-BHQ-1), trans-clivagem
- **Leitura**: Fluorescencia UV ou camera de smartphone

---

## Fontes de Dados

- **NCBI Entrez API** - Sequencias RefSeq de genes AMR
- **NCBI BLAST** - Validacao contra 124M genomas (banco nt)
- **AMRFinderPlus** - Catalogo de referencia (11.251 entradas)
- **ANVISA/BR-GLASS** - Dados epidemiologicos brasileiros para priorizacao

---

## Time

| Nome | Papel | Afiliacao |
|------|-------|-----------|
| **Joao Victor Pacheco Dias** | CTO & AI Architect | PhD candidate HC-FMUSP, CTO WingsAI, ITU/WHO AI for Health |
| **Dr. Gustavo Sakuno** | Clinical & Scientific Lead | Postdoc Harvard Medical School, Broad Institute, PhD USP |
| **Raul Primo** | Software Engineer | Engenheiro de Software, Pipeline development |

---

## Roadmap

Veja [MELHORIAS.md](MELHORIAS.md) para o roadmap completo de evolucao do projeto, incluindo integracoes com modelos genomicos fundacionais (Evo 2), CARD database, e explicacoes por IA.

Veja [TODO.md](TODO.md) para tarefas pendentes e em andamento.

---

## Licenca e Uso

Projeto de pesquisa - uso nao-diagnostico. Validacao clinica pendente.
Estrategia regulatoria: Fase 1 (pesquisa) -> Fase 2 (submissao IVD ANVISA, RDC 830/2023).

**Contato**: [IA para Medicos](https://www.iaparamedicos.com.br/)
