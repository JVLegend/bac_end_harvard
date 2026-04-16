"""Fase 7 BacEnd — utilidades compartilhadas para análises GPU-intensivas na DGX."""
import os
import json
import time
import pathlib
from typing import Iterable

EXP_DIR = os.environ.get(
    "FASE7_DIR",
    "/home/oftalmousp/jv-teste/harvard_bacend/fase7_results",
)
os.makedirs(EXP_DIR, exist_ok=True)

PROJECT_ROOT = pathlib.Path("/home/oftalmousp/jv-teste/harvard_bacend")

# Mapa das 12 familias AMR (referencias NCBI do pipeline existente)
AMR_FAMILIES = [
    "mecA", "blaKPC", "blaNDM", "vanA", "mcr-1", "blaCTX-M",
    "blaOXA-48", "blaVIM", "blaIMP", "blaGES", "qnrS", "armA",
]

DNA_ALPHA = {"A", "C", "G", "T"}


def save_json(obj, filename: str, pretty: bool = True) -> str:
    path = os.path.join(EXP_DIR, filename)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2 if pretty else None, default=str)
    print(f"SAVED: {path}")
    return path


def save_text(content: str, filename: str) -> str:
    path = os.path.join(EXP_DIR, filename)
    with open(path, "w") as f:
        f.write(content)
    print(f"SAVED: {path}")
    return path


def save_csv(df, filename: str) -> str:
    path = os.path.join(EXP_DIR, filename)
    df.to_csv(path, index=False)
    print(f"SAVED: {path}")
    return path


def emit_result(key: str, value) -> None:
    print(f"RESULT: {key}={value}")


def read_fasta(path) -> dict:
    """Le um FASTA simples → {header: sequence}."""
    result = {}
    current_header = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    result[current_header] = "".join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            result[current_header] = "".join(current_seq)
    return result


def list_family_seqs() -> dict:
    """Retorna {familia: [caminhos_fasta]} varrendo PROJECT_ROOT/sequences."""
    seqs_dir = PROJECT_ROOT / "sequences"
    out = {f: [] for f in AMR_FAMILIES}
    if not seqs_dir.exists():
        return out
    for f in seqs_dir.iterdir():
        if not f.is_file():
            continue
        name_lower = f.name.lower()
        for fam in AMR_FAMILIES:
            if fam.lower().replace("-", "") in name_lower.replace("-", "").replace("_", ""):
                out[fam].append(str(f))
                break
    return out


def timer(fn):
    def wrap(*a, **k):
        t0 = time.time()
        result = fn(*a, **k)
        print(f"  [{fn.__name__}] {time.time()-t0:.1f}s")
        return result
    return wrap
