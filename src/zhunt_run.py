"""Поиск Z-ДНК алгоритмом Z-Hunt (пакет zhuntrs — Rust-порт оригинального zhunt).

Для каждой последовательности считаются локально-оптимальные окна Z-ДНК.
Окна с ZH-score > порога (по заданию — 400) сливаются в непрерывные участки;
за участком сохраняется максимальный score.

Считается параллельно по последовательностям (multiprocessing).
Выход: BED6 (results/zhunt.bed), 0-based полуинтервал, strand = '.'.
"""
from __future__ import annotations

import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

import zhuntrs

# Параметры окна Z-ДНК в динуклеотидах.
# maxdn=8 (окна 6..16 п.н.) — баланс скорости и чувствительности:
# ловит те же максимальные score, что и maxdn=12, но на порядок быстрее.
MINDN = 3
MAXDN = 8
THRESHOLD = 400.0  # по заданию: z-score > 400


def iter_fasta(path: Path):
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


def merge_windows(windows: list[tuple[int, int, float]]) -> list[tuple[int, int, float]]:
    """Слить перекрывающиеся/смежные окна, сохранив максимальный score."""
    if not windows:
        return []
    windows.sort()
    merged: list[list] = [list(windows[0])]
    for start, end, score in windows[1:]:
        last = merged[-1]
        if start <= last[1]:
            last[1] = max(last[1], end)
            last[2] = max(last[2], score)
        else:
            merged.append([start, end, score])
    return [(s, e, sc) for s, e, sc in merged]


def process_seq(item: tuple[str, str]) -> tuple[str, list[tuple[int, int, float]]]:
    """Обработать одну последовательность: предсказание Z-ДНК -> слитые участки."""
    name, seq = item
    start, end, score, _seq, _conf = zhuntrs.predict(
        seq.encode(), MINDN, MAXDN, THRESHOLD, False
    )
    windows = [(start[i], end[i], score[i]) for i in range(len(start))]
    regions = merge_windows(windows)
    print(f"  {name}: {len(regions)} участков Z-ДНК", file=sys.stderr)
    return name, regions


def main() -> None:
    ap = argparse.ArgumentParser(description="Поиск Z-ДНК (zhunt, z-score > 400)")
    ap.add_argument("--genome", default="data/genome.fna", type=Path)
    ap.add_argument("--out", default="results/zhunt.bed", type=Path)
    ap.add_argument(
        "--only",
        default=None,
        help="ограничить список последовательностей (через запятую)",
    )
    ap.add_argument("--procs", default=8, type=int, help="число процессов")
    args = ap.parse_args()

    only = set(args.only.split(",")) if args.only else None
    args.out.parent.mkdir(parents=True, exist_ok=True)

    seqs = [
        (name, seq)
        for name, seq in iter_fasta(args.genome)
        if not only or name in only
    ]
    print(f"Последовательностей к обработке: {len(seqs)}", file=sys.stderr)

    with Pool(args.procs) as pool:
        results = pool.map(process_seq, seqs)

    total = 0
    with open(args.out, "w") as out:
        for name, regions in results:
            for i, (s, e, sc) in enumerate(regions, 1):
                out.write(f"{name}\t{s}\t{e}\tZhunt_{name}_{i}\t{sc:.1f}\t.\n")
            total += len(regions)

    print(f"Z-ДНК (zhunt, score>{THRESHOLD:.0f}, maxdn={MAXDN}): {total} участков")
    print(f"Записано в {args.out}")


if __name__ == "__main__":
    main()
