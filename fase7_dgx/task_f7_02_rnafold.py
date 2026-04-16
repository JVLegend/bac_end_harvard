"""Fase 7 · Task 2 — Secondary structure analysis dos guides CRISPR via RNAfold.

Por que: guides com muita estrutura secundaria (hairpins) podem nao hibridizar
eficientemente com o alvo. RNAfold calcula a energia minima (MFE) e estrutura
dot-bracket.

Estrategia: roda RNAfold (ViennaRNA, via pip `ViennaRNA` OU fallback
pela biblioteca `RNA` / subprocess). Se nao disponivel, usa aproximacao por
contagem de palindromos internos (como proxy).

Output:
- rnafold_guides.csv [gene, guide_id, spacer_seq, mfe, structure, palindromes, classificacao]
- rnafold_summary.json
"""
import os
import subprocess
import time
import pathlib
import pandas as pd
from _fase7_utils import save_csv, save_json, emit_result, EXP_DIR

GUIDES_DIR = pathlib.Path("/home/oftalmousp/jv-teste/harvard_bacend/guides")


def try_viennarna(seq: str):
    """Tenta RNA.fold() da biblioteca ViennaRNA, se instalada."""
    try:
        import RNA
        structure, mfe = RNA.fold(seq)
        return structure, float(mfe)
    except ImportError:
        pass
    # Fallback: chama binario externo RNAfold se existir
    try:
        res = subprocess.run(["RNAfold", "--noPS"], input=seq.encode(),
                             capture_output=True, timeout=30)
        lines = res.stdout.decode().split("\n")
        if len(lines) >= 2:
            # linha 1: sequencia, linha 2: "structure (mfe)"
            s_line = lines[1].strip()
            # parse "((((.....)))) (-5.30)"
            if " " in s_line:
                structure, mfe_part = s_line.rsplit(" ", 1)
                mfe = float(mfe_part.strip("()"))
                return structure, mfe
    except Exception:
        pass
    return None, None


def count_palindromes(seq: str, min_len: int = 4) -> int:
    """Fallback: conta palindromos (self-complementaridade) como proxy de hairpins."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A", "N": "N"}
    n = 0
    for k in range(min_len, min(12, len(seq)//2) + 1):
        for i in range(len(seq) - k + 1):
            sub = seq[i:i+k]
            rc = "".join(comp.get(b, "N") for b in reversed(sub))
            # procura rc em pos >= i+k (nao overlap)
            if rc in seq[i+k:]:
                n += 1
    return n


def classify(mfe, palindromes) -> str:
    if mfe is None:
        if palindromes == 0:
            return "likely_linear"
        elif palindromes <= 2:
            return "weak_structure"
        else:
            return "strong_structure"
    # com MFE real: kcal/mol. -10+ = estrutura forte
    if mfe >= -3:
        return "linear_ok"
    elif mfe >= -7:
        return "weak_hairpin"
    elif mfe >= -12:
        return "moderate_hairpin"
    else:
        return "strong_hairpin_risk"


def parse_guide_tsv(path: pathlib.Path):
    """Pipeline salva guides em .tsv com colunas: rank, gene, spacer_seq, pam_seq, ..."""
    try:
        df = pd.read_csv(path, sep="\t")
        return df
    except Exception as e:
        print(f"skip {path}: {e}")
        return None


def main():
    t0 = time.time()
    rows = []
    files = list(GUIDES_DIR.glob("*.tsv"))
    print(f"arquivos guides: {len(files)}")
    for f in files:
        gene = f.stem.replace("_guides", "").replace("_design", "")
        df = parse_guide_tsv(f)
        if df is None or df.empty:
            continue
        spacer_col = next((c for c in df.columns if "spacer" in c.lower() or "seq" in c.lower()), None)
        if not spacer_col:
            continue
        for idx, row in df.iterrows():
            spacer = str(row[spacer_col]).upper()
            if not spacer or any(b not in "ACGTU" for b in spacer):
                continue
            struct, mfe = try_viennarna(spacer)
            palins = count_palindromes(spacer)
            cls = classify(mfe, palins)
            rows.append({
                "gene": gene,
                "guide_idx": int(idx),
                "spacer": spacer,
                "length": len(spacer),
                "mfe_kcal_mol": round(mfe, 2) if mfe is not None else None,
                "structure": struct if struct else "",
                "palindromes": palins,
                "classification": cls,
            })

    if not rows:
        print("nenhum guide valido encontrado")
        return

    df = pd.DataFrame(rows)
    save_csv(df, "rnafold_guides.csv")

    summary = {
        "n_guides": len(df),
        "with_mfe": int(df["mfe_kcal_mol"].notna().sum()),
        "class_counts": df["classification"].value_counts().to_dict(),
        "mfe_stats": {
            "mean": float(df["mfe_kcal_mol"].mean()) if df["mfe_kcal_mol"].notna().any() else None,
            "min": float(df["mfe_kcal_mol"].min()) if df["mfe_kcal_mol"].notna().any() else None,
            "max": float(df["mfe_kcal_mol"].max()) if df["mfe_kcal_mol"].notna().any() else None,
        },
        "palindromes_stats": {
            "mean": round(float(df["palindromes"].mean()), 2),
            "max": int(df["palindromes"].max()),
        },
    }
    save_json(summary, "rnafold_summary.json")
    emit_result("n_guides", len(df))
    emit_result("mean_palindromes", round(float(df["palindromes"].mean()), 2))
    emit_result("strong_hairpin_count",
                int((df["classification"] == "strong_hairpin_risk").sum()))
    emit_result("seconds", f"{time.time()-t0:.1f}")


if __name__ == "__main__":
    main()
