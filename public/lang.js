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
      'seq.hero.kicker': 'Guia Prático Brasileiro',
      'seq.hero.h1': 'Sequenciar o próprio genoma <em>em casa</em>.',
      'seq.intro.title': 'O que é possível hoje',
      'seq.cost.title': 'Custos Reais (2025)',
      'seq.kit.title': 'Kits disponíveis no Brasil',
      'seq.analysis.title': 'O que fazer com os dados',
      'seq.legal.title': 'Questões Legais e Éticas',
      'seq.nanopore.title': 'Sequenciamento Nanopore DIY',

      /* ── DASHBOARD.HTML ── */
      'dash.title': 'SmartLab BacEnd — Dashboard Interativo',
      'dash.nav.back': 'Voltar ao site',
      'dash.subtitle': 'Análise computacional de alvos AMR',

      /* ── PANGENOME.HTML ── */
      'pangen.title': 'SmartSepsis — Pangenoma AMR',
      'pangen.nav.back': '← SmartSepsis',
      'pangen.hero.title': 'Análise de Pangenoma AMR',

      /* ── STRUCTURE.HTML ── */
      'struct.title': 'SmartSepsis — Estrutura 3D das Variantes AMR',
      'struct.hero.title': 'Estrutura 3D das Variantes AMR',

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
      'seq.hero.h1': 'Sequencing your own genome <em>at home</em>.',
      'seq.intro.title': 'What is possible today',
      'seq.cost.title': 'Real Costs (2025)',
      'seq.kit.title': 'Available kits in Brazil',
      'seq.analysis.title': 'What to do with the data',
      'seq.legal.title': 'Legal and Ethical Questions',
      'seq.nanopore.title': 'DIY Nanopore Sequencing',

      /* ── DASHBOARD.HTML ── */
      'dash.title': 'SmartLab BacEnd — Interactive Dashboard',
      'dash.nav.back': 'Back to site',
      'dash.subtitle': 'Computational analysis of AMR targets',

      /* ── PANGENOME.HTML ── */
      'pangen.title': 'SmartSepsis — AMR Pangenome',
      'pangen.nav.back': '← SmartSepsis',
      'pangen.hero.title': 'AMR Pangenome Analysis',

      /* ── STRUCTURE.HTML ── */
      'struct.title': 'SmartSepsis — 3D Structure of AMR Variants',
      'struct.hero.title': '3D Structure of AMR Variants',

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
