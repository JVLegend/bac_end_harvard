/* SmartSepsis — i18n (PT-BR / EN) */
(function () {
  const T = {
    pt: {
      /* ── NAV COMPARTILHADO ── */
      'nav.tab.smartsepsis': 'SmartSepsis',
      'nav.tab.crowdfunding': 'Crowdfunding Médico',
      'nav.tab.sequence': 'Sequenciar em Casa',
      'nav.link.clinical': 'Dor Clínica',
      'nav.link.pipeline': 'Pipeline',
      'nav.link.team': 'Time',
      'nav.link.pangenome': 'Pangenoma',
      'nav.link.dashboard': 'Dashboard',
      'nav.link.paper': 'Papers',

      /* ── INDEX.HTML ── */
      'index.title': 'SmartSepsis-Oph — AI-First Molecular Diagnostics for Ophthalmic Infections',
      'index.hero.badge': 'SmartSepsis-Oph — Programa AI-First de Diagnóstico Molecular Ocular',
      'index.hero.h1': 'Endoftalmite bacteriana pode cegar um olho em 48 horas.<br>Cultura demora 48 a 72.<br><span class="hl-green">Estamos projetando o teste molecular que cabe em uma amostra de nanolitros.</span>',
      'index.hero.quote': '&ldquo;O laboratório molecular que cabe na sala cirúrgica oftalmológica.&rdquo;',
      'index.hero.desc': 'SmartSepsis-Oph é um programa de design AI-first para diagnóstico molecular point-of-care em infecções oftalmológicas — a amostra mais difícil da medicina molecular: volumes sub-microlitro, baixa carga bacteriana, decisão em minutos. Estamos atualmente em <strong>fase de design computacional (Fase 0)</strong>, construindo bibliotecas de guias-RNA, modelos de especificidade e arquitetura de ensaio em tira de papel para o primeiro alvo declarado (<em>S. aureus</em> / mecA). <strong>Buscamos abertamente parceiros de wet-lab, colaboradores clínicos e financiadores</strong> para avançar à validação experimental.',
      'index.stat1.label': 'Guias candidatos desenhados in silico',
      'index.stat2.label': 'Famílias gênicas modeladas',
      'index.stat3.label': 'Genomas comparados (BLAST)',
      'index.stat4.label': 'Estágio atual',
      'index.stat5.label': 'Primeiro alvo declarado',
      'index.cta.support': 'Seja nosso parceiro',
      'index.cta.dashboard': 'Ver biblioteca de design',

      'index.dor.title': 'A dor clínica que motiva o projeto',
      'index.dor.sub': 'Cenários reais de oftalmologia onde a demora microbiológica se traduz em perda visual.',
      'index.dor.s1': '<strong>Cenário 1 — Endoftalmite pós-facoemulsificação:</strong> POI 3 dias. Hipopio, hialite, queda BCVA. Você faz tap vítreo + injeta vanco+ceftaz empírico. Cultura volta em 72h: <em>S. aureus</em> mecA-positivo. Vanco cobria, mas se fosse <em>P. aeruginosa</em> com blaVIM você teria perdido o olho. Cada hora de espera = mais inflamação, mais perda visual permanente.',
      'index.dor.s2': '<strong>Cenário 2 — Queratite por lente de contato:</strong> Paciente jovem, úlcera centro-periférica, dor 8/10. Suspeita <em>Pseudomonas</em> vs <em>Acanthamoeba</em>. Raspado pra cultura, fluoroquinolona tópica empirica. 48h depois, sem melhora. Era <em>Pseudomonas</em> com qnrS — resistente. Tempo perdido = perfuração corneana iminente.',
      'index.dor.s3': '<strong>Cenário 3 — CCIH oftalmológica:</strong> Surto de endoftalmite pós-injeção intraocular numa clínica de retina. 4 casos em 2 semanas. Vigilância ativa de profilaxia exige cultura semanal de superfície ocular — resultado em 5 dias. Antes disso, você já agendou 30 novas injeções. Risco invisivel.',
      'index.dor.solution': '<strong>Objetivo de design do SmartSepsis-Oph:</strong> um teste point-of-care em tira de papel que devolva informação molecular sobre genes de resistência conhecidos em fração do tempo da cultura, a partir de amostra ocular sub-microlitro. Em fase de design computacional. <em>Métricas de performance só serão reportadas após validação experimental com parceiro wet-lab.</em>',
      'index.dor.why': '<strong>Por que microbiologia ocular brasileira?</strong> O perfil epidemiológico ocular brasileiro difere de painéis comerciais europeus: <em>S. epidermidis</em> e <em>S. aureus</em> mecA dominantes em endoftalmite pós-catarata, <em>Pseudomonas</em> em queratite por lente de contato, ESBL emergente em úlcera neonatal, <em>Klebsiella</em> com KPC em endoftalmite endógena. Projeto liderado por João Victor Pacheco Dias (IA para Médicos / doutorando HC-FMUSP) com validação clínica do Dr. Gustavo Sakuno (postdoc Mass Eye and Ear / Harvard).',

      'index.docu.title': 'Docusérie (em planejamento): &ldquo;O Laboratório que Cabe na Sala Cirúrgica&rdquo;',
      'index.docu.sub': '<strong>Em planejamento — ainda não gravada.</strong> Pretendemos documentar publicamente cada fase do programa, do problema clínico ao design computacional e (uma vez firmada parceria) à validação experimental. Estrutura prevista de 5 episódios:',
      'index.docu.ep1.title': 'EP 01 — O problema clínico ocular',
      'index.docu.ep1.desc': 'Por que endoftalmite, queratite e infecções pós-injeção intravítrea continuam sendo tratadas empiricamente. Por que cultura não resolve. Por que o volume de amostra ocular é o gargalo molecular.',
      'index.docu.ep2.title': 'EP 02 — Por dentro do design AI-first',
      'index.docu.ep2.desc': 'Como agentes de IA desenham e priorizam guias CRISPR para um alvo declarado (S. aureus / mecA). Demonstração do pipeline computacional rodando.',
      'index.docu.ep3.title': 'EP 03 — A bancada (a confirmar)',
      'index.docu.ep3.desc': 'Episódio dependente de parceria de wet-lab confirmada. Documentaremos a primeira validação experimental do design — sem encenação, sem fingir resultado que ainda não temos.',
      'index.docu.ep4.title': 'EP 04 — Da bancada à tira',
      'index.docu.ep4.desc': 'Da reação molecular validada ao protótipo de leitura em papel. Cronograma e desafios reais de prototipagem.',
      'index.docu.ep5.title': 'EP 05 — Caminho regulatório e clínico',
      'index.docu.ep5.desc': 'O que separa um protótipo validado de um diagnóstico aprovado. RDC 830/2023, IRB, estudos prospectivos. O caminho honesto.',
      'index.docu.founders': 'Cronograma de gravação será divulgado após confirmação de parceria experimental.',
      'index.docu.why': '<strong>Por que documentar?</strong> Porque ciência AI-first em estágio inicial costuma se vender com efeitos visuais sem substância. Queremos fazer o oposto: documentar o que existe (design computacional), o que falta (validação experimental) e o caminho honesto entre os dois. Distribuição prevista: YouTube + IA para Médicos.',

      'index.diff.title': 'Por que esta abordagem',
      'index.diff.sub': 'CRISPR-Dx em papel existe desde 2017 (SHERLOCK/DETECTR). Nossa diferenciação não está na bioquímica — está em atacar o pior caso de amostra (volumes oculares sub-microlitro) com pipeline AI-first de design de guias, modelagem de especificidade in silico contra o microbioma ocular, e iteração rápida sobre alvos epidemiologicamente relevantes para o contexto brasileiro.',

      'index.dash.title': 'Biblioteca Computacional de Design (Fase 0)',
      'index.dash.sub': 'Bibliotecas de guias e primers desenhados in silico para 12 famílias gênicas de resistência relevantes em infecção ocular. Cada alvo passou por design CRISPR-Cas12a, primers RPA isotérmicos e comparação BLAST contra repositórios públicos. <strong>Resultados computacionais, ainda não validados experimentalmente.</strong>',

      'index.tech.title': 'Como o Design Funciona',
      'index.tech.sub': 'Pipeline computacional de 7 etapas, do fetch de sequências de referência ao scoring funcional e interpretação por IA. Todo o output é in silico — o passo experimental é o que estamos buscando colaboração para executar.',

      'index.valid.title': 'O que afirmamos — e o que não',
      'index.valid.sub': 'Honestidade científica sobre o estágio atual do programa.',

      'index.team.title': 'Equipe e Conselheiros',
      'index.team.sub': 'Time enxuto em fase de design. Buscando ativamente conselheiros seniores e parceiros de bancada.',

      'index.footer.tagline': 'Programa AI-first de diagnóstico molecular para infecções oftalmológicas',
      'index.footer.desc': 'Fase de design computacional. Pipeline agêntico de desenho de guias CRISPR-Cas12a + modelagem de especificidade + análise estrutural. Não é dispositivo médico. Não é para uso clínico.',
      'index.footer.location': 'São Paulo, Brasil · Programa de pesquisa por IA para Médicos · Aberto a colaborações',

      /* ── SMARTWETLAB.HTML ── */
      'wetlab.title': 'SmartSepsis — Crowdfunding entre Médicos | IA para Médicos',
      'wetlab.hero.eyebrow': 'Crowdfunding Médico · Fase Semente',
      'wetlab.hero.h1': 'Médicos que <em>entendem</em> o problema financiam a solução.',
      'wetlab.hero.sub': 'SmartSepsis-Oph: diagnóstico point-of-care de infecções oculares com CRISPR-Cas12a em 30 min. Sem termociclador. Desenhado por oftalmologistas, para oftalmologistas.',

      'wetlab.pain.title': 'Você prescreve antibiótico<br>às cegas. E o paciente paga.',
      'wetlab.pain.sub': 'Cultura microbiológica leva 48-72h. Em endoftalmite, cada hora importa.',

      'wetlab.why.title': 'Por que médicos financiam ciência',
      'wetlab.tiers.title': 'Estrutura de Apoio',
      'wetlab.tiers.sub': 'Modelos transparentes de participação — do simbólico ao estratégico.',
      'wetlab.mission.title': 'Nossa Missão',
      'wetlab.team.title': 'Time',
      'wetlab.faq.title': 'FAQ',

      /* ── SEQUENCIAR-EM-CASA.HTML ── */
      'seq.title': 'Sequenciar o Próprio Genoma em Casa — Versão Brasileira | SmartLab',
      'seq.hero.kicker': 'Guia Prático · Versão Brasileira',
      'seq.hero.h1': 'Quero sequenciar meu <em>genoma inteiro</em> na cozinha de casa.',
      'seq.hero.lede': 'Sim, é possível. Sim, no Brasil. Sim, com imposto. Este é o passo-a-passo honesto — equipamento, reagentes, preços em real, armadilhas de importação, e como conectar tudo isso à visão do SmartLab. Adaptado livremente do <a href="https://iwantosequencemygenomeathome.com/" target="_blank" rel="noopener">post original de Seth Showes</a>.',
      'seq.meta.difficulty': '<strong>Dificuldade:</strong> Médio-alto (requer paciência de pipetagem)',
      'seq.meta.time': '<strong>Tempo total:</strong> 3 dias úteis',
      'seq.meta.cost': '<strong>Custo:</strong> ~R$ 15k por run + R$ 40k de setup',
      'seq.s1.title': 'Por que sequenciar o próprio genoma?',
      'seq.s1.callout.label': 'Aviso importante',
      'seq.s2.title': 'Como a mágica acontece: MinION em 90 segundos',
      'seq.s3.title': 'Equipamento: o que comprar, onde, e por quanto',
      'seq.s3.callout.label': 'Dica brasileira',
      'seq.s4.title': 'Reagentes: kits, onde achar, e como driblar o "pacote de 24"',
      'seq.s4.callout.label': 'Realidade da importação',
      'seq.s5.title': 'O workflow de 3 dias',
      'seq.s5.callout.label': 'Dica: adaptive sampling é seu amigo',
      'seq.s6.title': 'O custo total em reais',
      'seq.s7.title': 'O que pode dar errado (e provavelmente vai)',
      'seq.s7.callout.label': 'Não é consultoria médica',
      'seq.s8.title': 'Por que esse post está no site do SmartLab?',
      'seq.s9.title': 'Quero começar. Por onde?',
      'seq.s9.callout.label': 'Precisa de ajuda?',
      'seq.table.item': 'Item',
      'seq.table.where': 'Onde comprar',
      'seq.table.price': 'Preço (R$)',
      'seq.table.reagent': 'Reagente',
      'seq.table.size': 'Tamanho / Fornecedor',

      /* ── DASHBOARD.HTML ── */
      'dash.title': 'SmartLab BacEnd — Dashboard Interativo',
      'dash.nav.back': 'Voltar ao site',
      'dash.subtitle': 'Análise computacional de alvos AMR',
      'dash.header.title': 'Biblioteca Computacional de Design — Explorador Interativo',
      'dash.header.desc': 'Explore os 42 alvos AMR desenhados in silico pelo pipeline: filtre por prioridade, família, impacto funcional predito. Clique em um gene para ver guides candidatos, primers RPA e scoring. Resultados computacionais — validação experimental pendente.',
      'dash.esm.label': '🧬 Inferencia ESM-2 (Fase 7) + Phenotype Probe (Fase 9)',
      'dash.esm.waiting': 'aguardando dados...',
      'dash.esm.loading': 'Carregando metricas...',
      'dash.slider.label': '🎚️ Multi-Objective Design Slider (Fase 8)',
      'dash.slider.efficacy': 'Eficácia (Cov)',
      'dash.slider.specificity': 'Especificidade',
      'dash.slider.coverage': 'Cobertura',
      'dash.slider.rpa': 'RPA (Tm/GC)',
      'dash.slider.cost': 'Custo (Oligo)',
      'dash.btn.reset': 'Reset',
      'dash.btn.save': 'Salvar',
      'dash.filter.priority': 'Prioridade',
      'dash.filter.class': 'Classe de resistencia',
      'dash.filter.search': 'Buscar gene',
      'dash.filter.all': 'Todas',
      'dash.filter.p1': 'P1 - Critica',
      'dash.filter.p2': 'P2 - Alta',
      'dash.filter.p3': 'P3 - Moderada',
      'dash.stat.targets': 'Alvos',
      'dash.stat.visible': 'Visivel',
      'dash.stat.blast': 'BLAST 100%',
      'dash.th.gene': 'Gene <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.class': 'Classe <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.priority': 'Prioridade <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.organism': 'Organismo <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.variants': 'Variantes <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.dnaimp': 'Impacto DNA (Evo2) <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.protimp': 'Impacto Prot (ESM2) <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.action': 'Acao',

      /* ── PANGENOME.HTML ── */
      'pangen.title': 'SmartSepsis — Pangenoma AMR',
      'pangen.nav.back': '← SmartSepsis',
      'pangen.hero.title': 'Análise de Pangenoma AMR',
      'pangen.header.title': 'Pangenoma de patógenos AMR brasileiros',
      'pangen.header.sub': '21 genomas (K. pneumoniae + E. coli) RefSeq · panaroo 1.6.1 · clean-mode sensitive',
      'pangen.kpi.genomes': 'Genomas',
      'pangen.kpi.genes': 'Genes (todos)',
      'pangen.kpi.shell': 'Shell genes',
      'pangen.kpi.cloud': 'Cloud genes',
      'pangen.s1.title': 'Distribuição de genes por número de cepas',
      'pangen.s1.desc': 'Histograma: quantos genes estão presentes em N das 21 cepas. Picos em valores baixos = cloud (raros); picos altos = conservados.',
      'pangen.s2.title': 'Top 30 genes mais conservados',
      'pangen.s2.desc': 'Genes presentes na maior fração de cepas. Genes "group_*" são clusters sem nome funcional anotado (prodigal não anotou função).',
      'pangen.s3.title': 'Resumo metodológico',
      'pangen.table.rank': '#',
      'pangen.table.gene': 'Gene',
      'pangen.table.present': 'Presente em',
      'pangen.table.frac': 'Fração',

      /* ── STRUCTURE.HTML ── */
      'struct.title': 'SmartSepsis — Estrutura 3D das Variantes AMR',
      'struct.hero.title': 'Estruturas 3D',
      'struct.sidebar.variants': 'Variantes AMR',
      'struct.loading': 'Carregando...',
      'struct.select.hint': 'Selecione uma variante',
      'struct.sidebar.details': 'Detalhes da Variante',
      'struct.details.hint': 'Selecione uma variante para ver score ESM-2, phenotype probe e estatísticas estruturais.',
      'struct.btn.surface': 'Superficie',

      /* ── PAPER.HTML ── */
      'paper.title': 'SmartSepsis — Propostas de Paper (Oftalmologia / Oculomics)',
    },

    en: {
      /* ── SHARED NAV ── */
      'nav.tab.smartsepsis': 'SmartSepsis',
      'nav.tab.crowdfunding': 'Medical Crowdfunding',
      'nav.tab.sequence': 'Sequence at Home',
      'nav.link.clinical': 'Clinical Problem',
      'nav.link.pipeline': 'Pipeline',
      'nav.link.team': 'Team',
      'nav.link.pangenome': 'Pangenome',
      'nav.link.dashboard': 'Dashboard',
      'nav.link.paper': 'Papers',

      /* ── INDEX.HTML ── */
      'index.title': 'SmartSepsis-Oph — AI-First Molecular Diagnostics for Ophthalmic Infections',
      'index.hero.badge': 'SmartSepsis-Oph — AI-First Molecular Diagnostics for Ophthalmic Infections',
      'index.hero.h1': 'Bacterial endophthalmitis can blind an eye in 48 hours.<br>Cultures take 48 to 72.<br><span class="hl-green">We are designing the molecular test that fits in a nanoliter sample.</span>',
      'index.hero.quote': '&ldquo;The molecular lab that fits inside the operating room.&rdquo;',
      'index.hero.desc': 'SmartSepsis-Oph is an AI-first design program for CRISPR-based point-of-care diagnostics in ophthalmic infections — the hardest sample in molecular medicine: sub-microliter volumes, low bacterial load, minutes to act. We are currently in <strong>computational design phase (Phase 0)</strong>, building guide-RNA libraries, specificity models, and a paper-strip assay architecture for our first declared target (<em>S. aureus</em> / mecA). <strong>We are openly seeking wet-lab partners, clinical collaborators, and funders</strong> to advance to experimental validation.',
      'index.stat1.label': 'Computationally designed guides',
      'index.stat2.label': 'Gene families modeled',
      'index.stat3.label': 'Genomes BLAST-compared',
      'index.stat4.label': 'Current stage',
      'index.stat5.label': 'First declared target',
      'index.cta.support': 'Partner with us',
      'index.cta.dashboard': 'View design library',

      'index.dor.title': 'The clinical pain that motivates the project',
      'index.dor.sub': 'Real ophthalmology scenarios where microbiological delay translates into vision loss.',
      'index.dor.s1': '<strong>Scenario 1 — Post-phacoemulsification endophthalmitis:</strong> POD 3. Hypopyon, vitritis, BCVA drop. You perform vitreous tap + inject empirical vanco+ceftaz. Culture returns at 72h: <em>S. aureus</em> mecA-positive. Vanco covered it, but if it were <em>P. aeruginosa</em> with blaVIM, the patient would have lost the eye. Every hour of waiting = more inflammation, more permanent visual loss.',
      'index.dor.s2': '<strong>Scenario 2 — Contact lens keratitis:</strong> Young patient, central-peripheral ulcer, pain 8/10. Suspicion of <em>Pseudomonas</em> vs <em>Acanthamoeba</em>. Scraping for culture, empirical topical fluoroquinolone. 48h later, no improvement. It was <em>Pseudomonas</em> with qnrS — resistant. Time lost = imminent corneal perforation.',
      'index.dor.s3': '<strong>Scenario 3 — Ophthalmic infection control:</strong> Outbreak of endophthalmitis post-intravitreal injection at a retina clinic. 4 cases in 2 weeks. Active prophylaxis surveillance requires weekly ocular surface culture — results in 5 days. By then, 30 new injections have already been scheduled. Invisible risk.',
      'index.dor.solution': '<strong>SmartSepsis-Oph design objective:</strong> a paper-strip point-of-care test that returns molecular information about known resistance genes in a fraction of culture time, from sub-microliter ocular samples. In computational design phase. <em>Performance metrics will only be reported after experimental validation with a wet-lab partner.</em>',
      'index.dor.why': '<strong>Why Brazilian ocular microbiology?</strong> The Brazilian ocular epidemiologic profile differs from European commercial panels: <em>S. epidermidis</em> and <em>S. aureus</em> mecA dominant in post-cataract endophthalmitis, <em>Pseudomonas</em> in contact-lens keratitis, emerging ESBL in neonatal ulcers, <em>Klebsiella</em> KPC in endogenous endophthalmitis. Project led by João Victor Pacheco Dias (AI for Physicians / HC-FMUSP doctoral candidate) with clinical validation from Dr. Gustavo Sakuno (postdoc, Mass Eye and Ear / Harvard).',

      'index.docu.title': 'Documentary Series (Planned): &ldquo;The Lab that Fits in the OR&rdquo;',
      'index.docu.sub': '<strong>In planning — not yet filmed.</strong> We intend to publicly document every phase of the program: the clinical problem, the computational design, and (once a partnership is in place) experimental validation. Planned 5-episode structure:',
      'index.docu.ep1.title': 'EP 01 — The ocular clinical problem',
      'index.docu.ep1.desc': 'Why endophthalmitis, keratitis, and post-injection infections are still treated empirically. Why culture does not solve it. Why ocular sample volume is the molecular bottleneck.',
      'index.docu.ep2.title': 'EP 02 — Inside the AI-first design',
      'index.docu.ep2.desc': 'How AI agents design and prioritize CRISPR guides for a declared target (S. aureus / mecA). Computational pipeline running on screen.',
      'index.docu.ep3.title': 'EP 03 — The bench (TBC)',
      'index.docu.ep3.desc': 'Episode contingent on confirmed wet-lab partnership. We will document the first experimental validation honestly — no staging, no faking results we do not yet have.',
      'index.docu.ep4.title': 'EP 04 — From bench to strip',
      'index.docu.ep4.desc': 'From a validated molecular reaction to a paper-strip readout prototype. Real timeline and prototyping challenges.',
      'index.docu.ep5.title': 'EP 05 — Regulatory and clinical pathway',
      'index.docu.ep5.desc': 'What separates a validated prototype from an approved diagnostic. ANVISA RDC 830/2023, IRB, prospective studies. The honest road.',
      'index.docu.founders': 'Filming schedule will be announced after experimental partnership is confirmed.',
      'index.docu.why': '<strong>Why document?</strong> Because early-stage AI-first science is too often sold with visuals that outrun substance. We want the opposite: document what exists (computational design), what is missing (experimental validation), and the honest path between them. Planned distribution: YouTube + AI for Physicians.',

      'index.diff.title': 'Why this approach',
      'index.diff.sub': 'Paper-based CRISPR-Dx has existed since 2017 (SHERLOCK/DETECTR). Our differentiation is not the biochemistry — it is attacking the worst-case sample (sub-microliter ocular volumes) with an AI-first guide-design pipeline, in silico specificity modeling against the ocular microbiome, and fast iteration on targets epidemiologically relevant to the Brazilian context.',

      'index.dash.title': 'Computational Design Library (Phase 0)',
      'index.dash.sub': 'Guide and primer libraries computationally designed for 12 resistance-gene families relevant to ocular infection. Each target underwent CRISPR-Cas12a guide design, isothermal RPA primer design, and BLAST comparison against public repositories. <strong>Computational results — not yet experimentally validated.</strong>',

      'index.tech.title': 'How the Design Works',
      'index.tech.sub': '7-step computational pipeline, from fetching reference sequences to functional scoring and AI interpretation. All output is in silico — the experimental step is what we are seeking collaboration to perform.',

      'index.valid.title': 'What we claim — and what we do not',
      'index.valid.sub': 'Scientific honesty about the program\'s current stage.',

      'index.team.title': 'Team and Advisors',
      'index.team.sub': 'Lean team in design phase. Actively seeking senior scientific advisors and wet-lab partners.',

      'index.footer.tagline': 'AI-first molecular diagnostics program for ophthalmic infections',
      'index.footer.desc': 'Computational design phase. Agentic pipeline for CRISPR-Cas12a guide design + specificity modeling + structural analysis. Not a medical device. Not for clinical use.',
      'index.footer.location': 'São Paulo, Brazil · A research program by AI for Physicians · Open to collaboration',

      /* ── SMARTWETLAB.HTML ── */
      'wetlab.title': 'SmartSepsis — Medical Crowdfunding | AI for Physicians',
      'wetlab.hero.eyebrow': 'Medical Crowdfunding · Seed Stage',
      'wetlab.hero.h1': 'Physicians who <em>understand</em> the problem fund the solution.',
      'wetlab.hero.sub': 'SmartSepsis-Oph: point-of-care diagnosis of ocular infections with CRISPR-Cas12a in 30 min. No thermocycler. Designed by ophthalmologists, for ophthalmologists.',

      'wetlab.pain.title': 'You prescribe antibiotics<br>blindly. And the patient pays.',
      'wetlab.pain.sub': 'Microbiological culture takes 48-72h. In endophthalmitis, every hour matters.',

      'wetlab.why.title': 'Why physicians fund science',
      'wetlab.tiers.title': 'Support Structure',
      'wetlab.tiers.sub': 'Transparent participation models — from symbolic to strategic.',
      'wetlab.mission.title': 'Our Mission',
      'wetlab.team.title': 'Team',
      'wetlab.faq.title': 'FAQ',

      /* ── SEQUENCIAR-EM-CASA.HTML ── */
      'seq.title': 'Sequencing Your Own Genome at Home — Brazilian Guide | SmartLab',
      'seq.hero.kicker': 'Brazilian Practical Guide',
      'seq.hero.h1': 'I want to sequence my entire <em>genome</em> at home.',
      'seq.hero.lede': 'Yes, it\'s possible. Yes, in Brazil. Yes, with import taxes. This is the honest step-by-step — equipment, reagents, prices in BRL, import pitfalls, and how to connect it all to the SmartLab vision. Freely adapted from the <a href="https://iwantosequencemygenomeathome.com/" target="_blank" rel="noopener">original post by Seth Showes</a>.',
      'seq.meta.difficulty': '<strong>Difficulty:</strong> Medium-high (requires pipetting patience)',
      'seq.meta.time': '<strong>Total time:</strong> 3 business days',
      'seq.meta.cost': '<strong>Cost:</strong> ~R$ 15k per run + R$ 40k setup',
      'seq.s1.title': 'Why sequence your own genome?',
      'seq.s1.callout.label': 'Important notice',
      'seq.s2.title': 'How the magic works: MinION in 90 seconds',
      'seq.s3.title': 'Equipment: what to buy, where, and for how much',
      'seq.s3.callout.label': 'Brazilian tip',
      'seq.s4.title': 'Reagents: kits, where to find them, and how to work around the "pack of 24"',
      'seq.s4.callout.label': 'Import reality',
      'seq.s5.title': 'The 3-day workflow',
      'seq.s5.callout.label': 'Tip: adaptive sampling is your friend',
      'seq.s6.title': 'Total cost in reais',
      'seq.s7.title': 'What can go wrong (and probably will)',
      'seq.s7.callout.label': 'Not medical advice',
      'seq.s8.title': 'Why is this post on the SmartLab site?',
      'seq.s9.title': 'I want to start. Where do I begin?',
      'seq.s9.callout.label': 'Need help?',
      'seq.table.item': 'Item',
      'seq.table.where': 'Where to buy',
      'seq.table.price': 'Price (R$)',
      'seq.table.reagent': 'Reagent',
      'seq.table.size': 'Size / Supplier',

      /* ── DASHBOARD.HTML ── */
      'dash.title': 'SmartLab BacEnd — Interactive Dashboard',
      'dash.nav.back': 'Back to site',
      'dash.subtitle': 'Computational analysis of AMR targets',
      'dash.header.title': 'Computational Design Library — Interactive Explorer',
      'dash.header.desc': 'Explore the 42 in silico–designed AMR targets: filter by priority, family, predicted functional impact. Click a gene to see candidate guides, RPA primers, and scoring. Computational results — experimental validation pending.',
      'dash.esm.label': '🧬 ESM-2 Inference (Phase 7) + Phenotype Probe (Phase 9)',
      'dash.esm.waiting': 'awaiting data...',
      'dash.esm.loading': 'Loading metrics...',
      'dash.slider.label': '🎚️ Multi-Objective Design Slider (Phase 8)',
      'dash.slider.efficacy': 'Efficacy (Cov)',
      'dash.slider.specificity': 'Specificity',
      'dash.slider.coverage': 'Coverage',
      'dash.slider.rpa': 'RPA (Tm/GC)',
      'dash.slider.cost': 'Cost (Oligo)',
      'dash.btn.reset': 'Reset',
      'dash.btn.save': 'Save',
      'dash.filter.priority': 'Priority',
      'dash.filter.class': 'Resistance class',
      'dash.filter.search': 'Search gene',
      'dash.filter.all': 'All',
      'dash.filter.p1': 'P1 - Critical',
      'dash.filter.p2': 'P2 - High',
      'dash.filter.p3': 'P3 - Moderate',
      'dash.stat.targets': 'Targets',
      'dash.stat.visible': 'Visible',
      'dash.stat.blast': 'BLAST 100%',
      'dash.th.gene': 'Gene <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.class': 'Class <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.priority': 'Priority <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.organism': 'Organism <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.variants': 'Variants <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.dnaimp': 'DNA Impact (Evo2) <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.protimp': 'Prot Impact (ESM2) <span class=\"sort-icon\">&#x25B4;&#x25BE;</span>',
      'dash.th.action': 'Action',

      /* ── PANGENOME.HTML ── */
      'pangen.title': 'SmartSepsis — AMR Pangenome',
      'pangen.nav.back': '← SmartSepsis',
      'pangen.hero.title': 'AMR Pangenome Analysis',
      'pangen.header.title': 'AMR pangenome of Brazilian pathogens',
      'pangen.header.sub': '21 genomes (K. pneumoniae + E. coli) RefSeq · panaroo 1.6.1 · clean-mode sensitive',
      'pangen.kpi.genomes': 'Genomes',
      'pangen.kpi.genes': 'Genes (all)',
      'pangen.kpi.shell': 'Shell genes',
      'pangen.kpi.cloud': 'Cloud genes',
      'pangen.s1.title': 'Gene distribution by number of strains',
      'pangen.s1.desc': 'Histogram: how many genes are present in N of the 21 strains. Peaks at low values = cloud (rare); peaks at high values = conserved.',
      'pangen.s2.title': 'Top 30 most conserved genes',
      'pangen.s2.desc': 'Genes present in the largest fraction of strains. "group_*" entries are clusters without annotated function (prodigal did not annotate).',
      'pangen.s3.title': 'Methodological summary',
      'pangen.table.rank': '#',
      'pangen.table.gene': 'Gene',
      'pangen.table.present': 'Present in',
      'pangen.table.frac': 'Fraction',

      /* ── STRUCTURE.HTML ── */
      'struct.title': 'SmartSepsis — 3D Structure of AMR Variants',
      'struct.hero.title': '3D Structures',
      'struct.sidebar.variants': 'AMR Variants',
      'struct.loading': 'Loading...',
      'struct.select.hint': 'Select a variant',
      'struct.sidebar.details': 'Variant Details',
      'struct.details.hint': 'Select a variant to view ESM-2 score, phenotype probe, and structural statistics.',
      'struct.btn.surface': 'Surface',

      /* ── PAPER.HTML ── */
      'paper.title': 'SmartSepsis — Paper Proposals (Ophthalmology / Oculomics)',
    },
  };

  /* ─── Engine ─── */
  function getLang() {
    return localStorage.getItem('ss_lang') || 'pt';
  }

  function setLang(lang) {
    localStorage.setItem('ss_lang', lang);
    applyLang(lang);
    updateButtons(lang);
  }

  function applyLang(lang) {
    const dict = T[lang] || T['pt'];

    /* document.title */
    document.querySelectorAll('[data-i18n-title]').forEach(el => {
      const key = el.dataset.i18nTitle;
      if (dict[key]) document.title = dict[key];
    });

    /* textContent */
    document.querySelectorAll('[data-i18n]').forEach(el => {
      const key = el.dataset.i18n;
      if (dict[key] !== undefined) el.textContent = dict[key];
    });

    /* innerHTML */
    document.querySelectorAll('[data-i18n-html]').forEach(el => {
      const key = el.dataset.i18nHtml;
      if (dict[key] !== undefined) el.innerHTML = dict[key];
    });

    /* html lang attr */
    document.documentElement.lang = lang === 'en' ? 'en' : 'pt-BR';
  }

  function updateButtons(lang) {
    document.querySelectorAll('.lang-btn').forEach(btn => {
      btn.classList.toggle('active', btn.dataset.lang === lang);
    });
  }

  function buildSwitcher() {
    const switcher = document.createElement('div');
    switcher.className = 'lang-switcher';
    switcher.innerHTML =
      '<button class="lang-btn" data-lang="pt" title="Português (BR)" onclick="SS_Lang.set(\'pt\')">🇧🇷</button>' +
      '<button class="lang-btn" data-lang="en" title="English (US)"  onclick="SS_Lang.set(\'en\')">🇺🇸</button>';
    return switcher;
  }

  function init() {
    /* Inject CSS */
    const style = document.createElement('style');
    style.textContent = `
      .lang-switcher{display:flex;align-items:center;gap:2px;margin-left:6px;}
      .lang-btn{background:none;border:none;cursor:pointer;font-size:1.25rem;line-height:1;
        padding:2px 4px;border-radius:3px;opacity:.45;transition:opacity .2s,transform .15s;
        display:flex;align-items:center;}
      .lang-btn:hover{opacity:.85;}
      .lang-btn.active{opacity:1;transform:scale(1.12);}
    `;
    document.head.appendChild(style);

    /* Inject switcher into every nav */
    document.querySelectorAll('nav').forEach(nav => {
      nav.appendChild(buildSwitcher());
    });

    const lang = getLang();
    applyLang(lang);
    updateButtons(lang);
  }

  /* Public API */
  window.SS_Lang = { set: setLang, get: getLang };

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
