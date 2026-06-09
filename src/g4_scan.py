"""Поиск потенциальных G-квадруплексов в геноме.

Паттерн из задания: (G{3,5}[ATGC]{1,7}){3,}G{3,5}
Ищется на обоих стрендах:
  * стренд (+): G-богатый паттерн;
  * стренд (-): C-богатый паттерн на прямой цепи — это эквивалент
    G-квадруплекса на обратной цепи, координаты сразу в системе (+).

Выход: BED6 (results/g4.bed), 0-based полуинтервал.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# Регистр приводим к верхнему до поиска ("Внимание к регистру" в задании).
G4_PLUS = re.compile(r"(?:G{3,5}[ATGC]{1,7}){3,}G{3,5}")
G4_MINUS = re.compile(r"(?:C{3,5}[ATGC]{1,7}){3,}C{3,5}")


def iter_fasta(path: Path):
    """Простой потоковый парсер FASTA -> (имя, последовательность в верхнем регистре)."""
    name = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).upper()


def scan(fasta: Path, out_bed: Path) -> dict[str, int]:
    counts = {"+": 0, "-": 0}
    with open(out_bed, "w") as out:
        for name, seq in iter_fasta(fasta):
            for strand, pattern in (("+", G4_PLUS), ("-", G4_MINUS)):
                for m in pattern.finditer(seq):
                    start, end = m.start(), m.end()
                    counts[strand] += 1
                    idx = counts[strand]
                    # score = длина структуры (для BED колонки 5)
                    out.write(
                        f"{name}\t{start}\t{end}\tG4_{strand}_{idx}\t{end - start}\t{strand}\n"
                    )
            print(f"  {name}: готово", file=sys.stderr)
    return counts


def main() -> None:
    ap = argparse.ArgumentParser(description="Поиск G-квадруплексов (regex, оба стренда)")
    ap.add_argument("--genome", default="data/genome.fna", type=Path)
    ap.add_argument("--out", default="results/g4.bed", type=Path)
    args = ap.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    counts = scan(args.genome, args.out)
    total = counts["+"] + counts["-"]
    print(f"G-квадруплексы: всего {total}  (+ {counts['+']}, - {counts['-']})")
    print(f"Записано в {args.out}")


if __name__ == "__main__":
    main()
