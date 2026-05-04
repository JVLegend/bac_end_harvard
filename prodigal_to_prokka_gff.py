#!/usr/bin/env python3
"""
Converte GFF do prodigal em GFF prokka-compativel para panaroo.

Uso: python prodigal_to_prokka_gff.py prodigal.gff genome.fna prefix output.gff

Diferencas:
- IDs: prodigal usa "ID=1_1" -> prokka usa "ID=PREFIX_00001"
- Adiciona inference, locus_tag, product (placeholder), gene (vazio)
- Source virá como "Prodigal:002006" em vez de "prodigal-2.6.3"
- Anexa FASTA da sequencia ao final com separador ##FASTA
"""

import sys
from pathlib import Path


def main():
    if len(sys.argv) != 5:
        print("Uso: prodigal_to_prokka_gff.py prodigal.gff genome.fna prefix output.gff")
        sys.exit(1)

    prodigal_gff, genome_fa, prefix, out_gff = sys.argv[1:]
    # Sanitiza prefix: panaroo nao gosta de pontos/hifens em locus_tag prefix
    sanitized = "".join(c if c.isalnum() else "" for c in prefix)[:8].upper() or "GENOME"

    contig_seqs = {}
    cur = None
    with open(genome_fa) as f:
        for line in f:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                contig_seqs[cur] = []
            elif cur and line.strip():
                contig_seqs[cur].append(line.strip())

    out_lines = ["##gff-version 3"]
    for contig in contig_seqs:
        seq = "".join(contig_seqs[contig])
        out_lines.append(f"##sequence-region {contig} 1 {len(seq)}")

    counter = 0
    with open(prodigal_gff) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            counter += 1
            locus_tag = f"{sanitized}_{counter:05d}"
            attrs = (
                f"ID={locus_tag};"
                f"locus_tag={locus_tag};"
                f"inference=ab initio prediction:Prodigal:002006;"
                f"product=hypothetical protein"
            )
            out_lines.append(
                "\t".join(
                    [
                        cols[0],   # contig
                        "Prodigal:002006",
                        "CDS",
                        cols[3],   # start
                        cols[4],   # end
                        cols[5],   # score
                        cols[6],   # strand
                        cols[7],   # phase
                        attrs,
                    ]
                )
            )

    out_lines.append("##FASTA")
    for contig, parts in contig_seqs.items():
        seq = "".join(parts)
        out_lines.append(f">{contig}")
        for i in range(0, len(seq), 60):
            out_lines.append(seq[i : i + 60])

    Path(out_gff).write_text("\n".join(out_lines) + "\n")
    print(f"  -> {out_gff} ({counter} CDS)")


if __name__ == "__main__":
    main()
