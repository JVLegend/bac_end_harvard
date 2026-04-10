"""
Funções utilitárias para o pipeline CRISPR-Cas12a.
"""

COMPLEMENT = str.maketrans("ATCGatcg", "TAGCtagc")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    seq = seq.upper()
    if len(seq) == 0:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return gc / len(seq)


def max_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    seq = seq.upper()
    max_run = 1
    current_run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1
    return max_run


def tm_basic(seq: str) -> float:
    """Cálculo básico de Tm (Wallace rule para <14nt, ajustado para maiores)."""
    seq = seq.upper()
    n = len(seq)
    gc = sum(1 for b in seq if b in "GC")
    at = sum(1 for b in seq if b in "AT")
    if n < 14:
        return 2 * at + 4 * gc
    # Fórmula salt-adjusted simplificada (50mM Na+)
    return 64.9 + 41.0 * (gc - 16.4) / n


def self_complementarity_score(seq: str) -> int:
    """Score simples de auto-complementaridade (3' end)."""
    seq = seq.upper()
    rc = reverse_complement(seq)
    score = 0
    # Check 3' end (últimos 6 nt) contra toda a sequência
    tail = seq[-6:]
    for i in range(len(rc) - 5):
        matches = sum(1 for a, b in zip(tail, rc[i : i + 6]) if a == b)
        if matches >= 4:
            score += matches
    return score


def has_poly_t(seq: str, min_run: int = 4) -> bool:
    """Verifica se há poly-T (termina transcrição in vitro)."""
    return "T" * min_run in seq.upper()


def parse_fasta(filepath: str) -> dict[str, str]:
    """Parse arquivo FASTA, retorna dict {header: sequence}."""
    sequences = {}
    current_header = None
    current_seq = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            elif line:
                current_seq.append(line)
    if current_header is not None:
        sequences[current_header] = "".join(current_seq)
    return sequences


def write_fasta(filepath: str, header: str, sequence: str, line_width: int = 70):
    """Escreve sequência em formato FASTA."""
    with open(filepath, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i : i + line_width] + "\n")


def find_pam_sites(sequence: str, pam_pattern: str = "TTTV") -> list[dict]:
    """
    Encontra sites PAM TTTV em ambas as fitas.
    Para Cas12a, PAM fica upstream (5') do protospacer na fita não-alvo.

    Retorna lista de dicts com: position, strand, pam_seq, spacer_start
    """
    seq = sequence.upper()
    rc = reverse_complement(seq)
    v_bases = {"A", "C", "G"}  # V = não T
    sites = []

    # Fita sense: procurar TTT[V] seguido do spacer
    for i in range(len(seq) - 24):
        if seq[i : i + 3] == "TTT" and seq[i + 3] in v_bases:
            sites.append(
                {
                    "position": i,
                    "strand": "+",
                    "pam_seq": seq[i : i + 4],
                    "spacer_start": i + 4,
                }
            )

    # Fita antisense: procurar na reverse complement
    for i in range(len(rc) - 24):
        if rc[i : i + 3] == "TTT" and rc[i + 3] in v_bases:
            # Converter posição para coordenada na fita sense
            orig_pos = len(seq) - i - 4
            sites.append(
                {
                    "position": orig_pos,
                    "strand": "-",
                    "pam_seq": rc[i : i + 4],
                    "spacer_start": i + 4,
                }
            )

    return sites


def extract_spacer(sequence: str, site: dict, length: int = 20) -> str:
    """Extrai sequência do spacer a partir de um site PAM."""
    seq = sequence.upper()
    if site["strand"] == "+":
        start = site["spacer_start"]
        return seq[start : start + length]
    else:
        rc = reverse_complement(seq)
        start = site["spacer_start"]
        return rc[start : start + length]


def score_guide(spacer: str, position: int, gene_length: int) -> dict:
    """
    Pontua um guide RNA candidato.
    Score mais alto = melhor guide.
    """
    gc = gc_content(spacer)
    homo = max_homopolymer(spacer)
    self_comp = self_complementarity_score(spacer)
    poly_t = has_poly_t(spacer)

    score = 100.0

    # GC content scoring
    if 0.40 <= gc <= 0.60:
        score += 20  # ideal
    elif 0.30 <= gc <= 0.70:
        score += 5  # aceitável
    else:
        score -= 30  # ruim

    # Penalizar homopolímeros
    if homo >= 5:
        score -= 40
    elif homo >= 4:
        score -= 15

    # Penalizar poly-T (termina transcrição)
    if poly_t:
        score -= 25

    # Penalizar auto-complementaridade
    score -= self_comp * 2

    # Preferir posições centrais no gene
    rel_pos = position / gene_length
    if 0.2 <= rel_pos <= 0.8:
        score += 10
    elif 0.1 <= rel_pos <= 0.9:
        score += 5

    return {
        "score": round(score, 1),
        "gc": round(gc, 3),
        "homopolymer": homo,
        "self_comp": self_comp,
        "poly_t": poly_t,
        "rel_position": round(rel_pos, 3),
    }
