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
      'index.title': 'SmartSepsis — O laboratório que cabe no PS | IA para Médicos',
      'index.hero.badge': 'SmartSepsis-Oph — Microbiologia Inteligente em Oftalmologia',
      'index.hero.h1': 'Endoftalmite pós-cirurgia. Queratite com hidrogel.<br>Você prescreve vanco+ceftaz empírico e espera 48h.<br><span class="hl-green">E se você soubesse em 30 minutos?</span>',
      'index.hero.quote': '&ldquo;O laboratório molecular que cabe no centro cirúrgico oftalmológico.&rdquo;',
      'index.hero.desc': 'Diagnóstico point-of-care de patógenos oculares com CRISPR-Cas12a em papel. Endoftalmite, queratite, conjuntivite refratária, profilaxia perioperatória. Em vez de esperar a cultura de 48-72h após tap vítreo ou raspado de córnea, você sabe em 30 minutos se o patógeno carrega gene de resistência — mecA (MRSA pós-LASIK), blaKPC/blaNDM (Klebsiella endoftalmite), vanA, mcr, blaCTX-M (ESBL) e mais 6 famílias relevantes a infecções oculares. Sem termociclador. Custo: ~R$25 por teste. Desenhado com Dr. Sakuno (Mass Eye and Ear / Harvard) e equipe HC-FMUSP.',
      'index.stat1.label': 'Resultado à beira do leito',
      'index.stat2.label': 'Custo por teste (vs R$200+ PCR)',
      'index.stat3.label': 'Alvos AMR clínicos cobertos',
      'index.stat4.label': 'Antibiograma vs SmartSepsis',
      'index.stat5.label': 'Famílias gênicas validadas',
      'index.cta.support': 'Apoiar como médico fundador',
      'index.cta.dashboard': 'Ver dashboard técnico',

      'index.dor.title': 'A dor clínica que nos move',
      'index.dor.sub': 'Não é um problema de bancada. É uma decisão que você toma plantão após plantão.',
      'index.dor.s1': '<strong>Cenário 1 — Endoftalmite pós-facoemulsificação:</strong> POI 3 dias. Hipopio, hialite, queda BCVA. Você faz tap vítreo + injeta vanco+ceftaz empírico. Cultura volta em 72h: <em>S. aureus</em> mecA-positivo. Vanco cobria, mas se fosse <em>P. aeruginosa</em> com blaVIM você teria perdido o olho. Cada hora de espera = mais inflamação, mais perda visual permanente.',
      'index.dor.s2': '<strong>Cenário 2 — Queratite por lente de contato:</strong> Paciente jovem, úlcera centro-periférica, dor 8/10. Suspeita <em>Pseudomonas</em> vs <em>Acanthamoeba</em>. Raspado pra cultura, fluoroquinolona tópica empirica. 48h depois, sem melhora. Era <em>Pseudomonas</em> com qnrS — resistente. Tempo perdido = perfuração corneana iminente.',
      'index.dor.s3': '<strong>Cenário 3 — CCIH oftalmológica:</strong> Surto de endoftalmite pós-injeção intraocular numa clínica de retina. 4 casos em 2 semanas. Vigilância ativa de profilaxia exige cultura semanal de superfície ocular — resultado em 5 dias. Antes disso, você já agendou 30 novas injeções. Risco invisivel.',
      'index.dor.solution': '<strong>SmartSepsis-Oph muda isso:</strong> Em 30 minutos, no centro cirúrgico ou no consultório, sem termociclador, você sabe se o patógeno ocular carrega gene de resistência. Não substitui cultura — mas <strong>antecipa a decisão em 47h30min</strong>. Antibiótico intraocular mais preciso. Profilaxia perioperatória racional. Vigilância ativa em centros de retina.',
      'index.dor.why': '<strong>Por que oftalmologistas brasileiros + Mass Eye and Ear?</strong> Painel desenhado para a microbiologia ocular do contexto brasileiro: <em>S. epidermidis</em> e <em>S. aureus</em> mecA dominantes em endoftalmite pós-catarata, <em>Pseudomonas</em> em queratite por LC, ESBL emergente em úlcera neonatal, <em>Klebsiella</em> com KPC em endoftalmite endógena. Liderado pelo Dr. Sakuno (Mass Eye and Ear / Harvard, postdoc em <em>oculomica e biomarcadores</em>) e equipe HC-FMUSP.',

      'index.docu.title': 'Docusérie: &ldquo;O Laboratório que Cabe no PS&rdquo;',
      'index.docu.sub': 'Em paralelo à validação clínica, estamos filmando uma docusérie médica em 5 episódios documentando a construção do SmartSepsis — do problema clínico real ao MVP no leito.',
      'index.docu.ep1.title': 'A noite no PS oftalmológico',
      'index.docu.ep1.desc': 'Plantão real em emergência ocular. Tap vítreo, raspado de córnea, decisão de injeção intraocular sem dado microbiológico. A dor concreta de prescrever sem saber o patógeno.',
      'index.docu.ep2.title': 'Bancada do HC-FMUSP',
      'index.docu.ep2.desc': 'Karine no wet lab. CRISPR-Cas12a clivando DNA em fluorescência visualizável. A bioquímica que vira diagnóstico.',
      'index.docu.ep3.title': 'Pipeline de IA',
      'index.docu.ep3.desc': '5 agentes de IA desenhando guides CRISPR em 5 minutos. O que antes levava semanas de bioinformática.',
      'index.docu.ep4.title': 'MVP no centro cirúrgico',
      'index.docu.ep4.desc': 'Primeiro teste em tap vítreo real. Dr. Sakuno conduz. 30 min cronometrados, durante a injeção intraocular. Resultado fluorescente antes do paciente sair da maca.',
      'index.docu.ep5.title': 'ROI clínico',
      'index.docu.ep5.desc': 'Custo de um olho perdido por endoftalmite mal-tratada vs custo do teste. BCVA preservada. Re-injeções evitadas. O número que importa em oftalmologia.',
      'index.docu.founders': 'Apoiadores do crowdfunding recebem acesso antecipado aos episódios + behind-the-scenes da bancada e do PS.',
      'index.docu.why': '<strong>Por que docusérie?</strong> Porque ciência que fica em paper não muda prática clínica. Queremos que infectologista, intensivista e CCIH de qualquer hospital brasileiro <em>vejam</em> como o SmartSepsis foi construído — e cobrem a tecnologia de seus diretores técnicos. Distribuição: YouTube + IA para Médicos + recortes Instagram/LinkedIn.',

      'index.diff.title': 'Por que SmartSepsis é diferente',
      'index.diff.sub': 'CRISPR em papel existe desde 2017 (SHERLOCK/DETECTR). Nosso diferencial não é a biologia — é combinar IA clínica + epidemiologia brasileira + custo de PS.',

      'index.dash.title': 'Cobertura do Pipeline',
      'index.dash.sub': '42 alvos AMR processados computacionalmente: 12 famílias gênicas de referência + 30 variantes clínicas prevalentes em UTIs brasileiras. Cada alvo passou por design de guide CRISPR-Cas12a, primers RPA e validação BLAST contra 124 milhões de genomas.',

      'index.tech.title': 'Como Funciona',
      'index.tech.sub': 'Pipeline computacional de 7 etapas, incluindo scoring avançado e interpretação por IA.',

      'index.valid.title': 'Validação e Posicionamento',
      'index.valid.sub': 'Transparência sobre o que já validamos e o que ainda não.',

      'index.team.title': 'Time',
      'index.team.sub': 'Expertise combinada em IA, medicina e engenharia.',

      'index.footer.tagline': 'O laboratório molecular que cabe no centro cirúrgico oftalmológico',
      'index.footer.desc': 'Diagnóstico point-of-care de patógenos oculares com CRISPR-Cas12a + foundation models de proteína (ESM-2 + ProtT5) + pipeline agêntico',
      'index.footer.location': 'São Paulo, Brasil | HC-FMUSP × Mass Eye and Ear (Harvard) | 43 variantes AMR · 12 famílias gênicas oculares',

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
      'dash.header.title': 'Dashboard Interativo — Painel AMR',
      'dash.header.desc': 'Explore os 42 alvos do pipeline: filtre por prioridade, família, impacto funcional. Clique em um gene para ver detalhes.',
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
      'index.title': 'SmartSepsis — The Lab that Fits in the OR | AI for Physicians',
      'index.hero.badge': 'SmartSepsis-Oph — Intelligent Microbiology in Ophthalmology',
      'index.hero.h1': 'Post-surgery endophthalmitis. Hydrogel keratitis.<br>You prescribe vanco+ceftaz empirically and wait 48h.<br><span class="hl-green">What if you knew in 30 minutes?</span>',
      'index.hero.quote': '&ldquo;The molecular lab that fits in the ophthalmic surgical center.&rdquo;',
      'index.hero.desc': 'Point-of-care diagnosis of ocular pathogens with paper-based CRISPR-Cas12a. Endophthalmitis, keratitis, refractory conjunctivitis, perioperative prophylaxis. Instead of waiting 48-72h for culture after vitreous tap or corneal scraping, you know in 30 minutes whether the pathogen carries a resistance gene — mecA (MRSA post-LASIK), blaKPC/blaNDM (Klebsiella endophthalmitis), vanA, mcr, blaCTX-M (ESBL) and 6 more families relevant to ocular infections. No thermocycler. Cost: ~R$25 per test. Designed with Dr. Sakuno (Mass Eye and Ear / Harvard) and the HC-FMUSP team.',
      'index.stat1.label': 'Bedside result',
      'index.stat2.label': 'Cost per test (vs R$200+ PCR)',
      'index.stat3.label': 'Clinical AMR targets covered',
      'index.stat4.label': 'Antibiogram vs SmartSepsis',
      'index.stat5.label': 'Validated gene families',
      'index.cta.support': 'Support as a founding physician',
      'index.cta.dashboard': 'View technical dashboard',

      'index.dor.title': 'The clinical pain that drives us',
      'index.dor.sub': 'This is not a bench problem. It is a decision you make shift after shift.',
      'index.dor.s1': '<strong>Scenario 1 — Post-phacoemulsification endophthalmitis:</strong> POD 3. Hypopyon, vitritis, BCVA drop. You perform vitreous tap + inject empirical vanco+ceftaz. Culture returns at 72h: <em>S. aureus</em> mecA-positive. Vanco covered it, but if it were <em>P. aeruginosa</em> with blaVIM, the patient would have lost the eye. Every hour of waiting = more inflammation, more permanent visual loss.',
      'index.dor.s2': '<strong>Scenario 2 — Contact lens keratitis:</strong> Young patient, central-peripheral ulcer, pain 8/10. Suspicion of <em>Pseudomonas</em> vs <em>Acanthamoeba</em>. Scraping for culture, empirical topical fluoroquinolone. 48h later, no improvement. It was <em>Pseudomonas</em> with qnrS — resistant. Time lost = imminent corneal perforation.',
      'index.dor.s3': '<strong>Scenario 3 — Ophthalmic infection control:</strong> Outbreak of endophthalmitis post-intravitreal injection at a retina clinic. 4 cases in 2 weeks. Active prophylaxis surveillance requires weekly ocular surface culture — results in 5 days. By then, 30 new injections have already been scheduled. Invisible risk.',
      'index.dor.solution': '<strong>SmartSepsis-Oph changes this:</strong> In 30 minutes, in the OR or clinic, without a thermocycler, you know whether the ocular pathogen carries a resistance gene. It does not replace culture — but <strong>anticipates the decision by 47h30min</strong>. More precise intravitreal antibiotic. Rational perioperative prophylaxis. Active surveillance at retina centers.',
      'index.dor.why': '<strong>Why Brazilian ophthalmologists + Mass Eye and Ear?</strong> Panel designed for Brazilian ocular microbiology: <em>S. epidermidis</em> and <em>S. aureus</em> mecA dominant in post-cataract endophthalmitis, <em>Pseudomonas</em> in CL keratitis, emerging ESBL in neonatal ulcers, <em>Klebsiella</em> KPC in endogenous endophthalmitis. Led by Dr. Sakuno (Mass Eye and Ear / Harvard, postdoc in <em>oculomics and biomarkers</em>) and the HC-FMUSP team.',

      'index.docu.title': 'Documentary Series: &ldquo;The Lab that Fits in the ER&rdquo;',
      'index.docu.sub': 'In parallel with clinical validation, we are filming a 5-episode medical documentary series chronicling the construction of SmartSepsis — from the real clinical problem to the bedside MVP.',
      'index.docu.ep1.title': 'A night in the ophthalmic ER',
      'index.docu.ep1.desc': 'Real shift in ocular emergency. Vitreous tap, corneal scraping, intravitreal injection decision without microbiological data. The concrete pain of prescribing without knowing the pathogen.',
      'index.docu.ep2.title': 'HC-FMUSP wet lab bench',
      'index.docu.ep2.desc': 'Karine in the wet lab. CRISPR-Cas12a cleaving DNA with visible fluorescence. The biochemistry that becomes a diagnostic.',
      'index.docu.ep3.title': 'AI pipeline',
      'index.docu.ep3.desc': '5 AI agents designing CRISPR guides in 5 minutes. What used to take weeks of bioinformatics.',
      'index.docu.ep4.title': 'MVP in the surgical center',
      'index.docu.ep4.desc': 'First test on real vitreous tap. Dr. Sakuno leads. 30 minutes timed during intravitreal injection. Fluorescent result before the patient leaves the table.',
      'index.docu.ep5.title': 'Clinical ROI',
      'index.docu.ep5.desc': 'Cost of an eye lost to poorly treated endophthalmitis vs cost of the test. BCVA preserved. Re-injections avoided. The number that matters in ophthalmology.',
      'index.docu.founders': 'Crowdfunding supporters receive early access to episodes + behind-the-scenes from the bench and ER.',
      'index.docu.why': '<strong>Why a documentary series?</strong> Because science that stays in papers does not change clinical practice. We want infectious disease specialists, intensivists, and infection control teams at any Brazilian hospital to <em>see</em> how SmartSepsis was built — and demand the technology from their technical directors. Distribution: YouTube + AI for Physicians + Instagram/LinkedIn clips.',

      'index.diff.title': 'Why SmartSepsis is different',
      'index.diff.sub': 'Paper-based CRISPR has existed since 2017 (SHERLOCK/DETECTR). Our differentiator is not the biology — it is combining clinical AI + Brazilian epidemiology + ER-level cost.',

      'index.dash.title': 'Pipeline Coverage',
      'index.dash.sub': '42 AMR targets computationally processed: 12 reference gene families + 30 clinical variants prevalent in Brazilian ICUs. Each target underwent CRISPR-Cas12a guide design, RPA primers, and BLAST validation against 124 million genomes.',

      'index.tech.title': 'How It Works',
      'index.tech.sub': '7-step computational pipeline including advanced scoring and AI interpretation.',

      'index.valid.title': 'Validation and Positioning',
      'index.valid.sub': 'Transparency about what we have and have not yet validated.',

      'index.team.title': 'Team',
      'index.team.sub': 'Combined expertise in AI, medicine, and engineering.',

      'index.footer.tagline': 'The molecular lab that fits in the ophthalmic surgical center',
      'index.footer.desc': 'Point-of-care diagnosis of ocular pathogens with CRISPR-Cas12a + protein foundation models (ESM-2 + ProtT5) + agentic pipeline',
      'index.footer.location': 'São Paulo, Brazil | HC-FMUSP × Mass Eye and Ear (Harvard) | 43 AMR variants · 12 ocular gene families',

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
      'dash.header.title': 'Interactive Dashboard — AMR Panel',
      'dash.header.desc': 'Explore the 42 pipeline targets: filter by priority, family, functional impact. Click a gene to see details.',
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
