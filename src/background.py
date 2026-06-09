"""Сравнение распределения структур с фоном (рандомизация).

Для каждого набора структур координаты случайно перемешиваются по геному
(bedtools shuffle, сохраняются число и длины интервалов). На каждой реплике
пересчитывается приоритетное отнесение к участкам (как в Таблице 1).
Среднее по репликам — ожидаемая доля при случайном размещении.

Выход: results/background.csv c наблюдаемой и ожидаемой долей и log2(obs/exp).
"""
from __future__ import annotations

import argparse
import math
import subprocess
from pathlib import Path

import pandas as pd

FEATURES = [
    ("exons", "Экзоны"),
    ("introns", "Интроны"),
    ("promoters", "Промоторы (1000 п.н. до TSS)"),
    ("downstream", "Downstream (200 п.н.)"),
    ("intergenic", "Межгенные участки"),
]
STRUCTURES = [("g4", "Квадруплексы"), ("zhunt", "Zhunt"), ("zdnabert", "ZDNABERT")]


def sh(cmd: str) -> str:
    return subprocess.run(
        cmd, shell=True, check=True, capture_output=True, text=True
    ).stdout


def count_lines(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    return int(sh(f"wc -l < {path}").strip())


def assign_fractions(struct_bed: Path, feat_sorted: dict[str, Path],
                     tmp: Path, tag: str) -> dict[str, float]:
    """Приоритетное отнесение структур к участкам -> доли по участкам."""
    total = count_lines(struct_bed)
    fractions: dict[str, float] = {}
    remaining = struct_bed
    for fkey, _ in FEATURES:
        inside = tmp / f"{tag}_{fkey}_in.bed"
        rest = tmp / f"{tag}_{fkey}_rest.bed"
        sh(f"bedtools intersect -u -a {remaining} -b {feat_sorted[fkey]} > {inside}")
        sh(f"bedtools intersect -v -a {remaining} -b {feat_sorted[fkey]} > {rest}")
        fractions[fkey] = count_lines(inside) / total if total else 0.0
        remaining = rest
    return fractions


def main() -> None:
    ap = argparse.ArgumentParser(description="Сравнение с фоном (рандомизация)")
    ap.add_argument("--results", default="results", type=Path)
    ap.add_argument("--sizes", default="data/genome.chrom.sizes", type=Path)
    ap.add_argument("--replicates", default=100, type=int)
    args = ap.parse_args()
    res = args.results
    tmp = res / "_tmp_bg"
    tmp.mkdir(parents=True, exist_ok=True)

    feat_sorted: dict[str, Path] = {}
    for fkey, _ in FEATURES:
        s = tmp / f"feat_{fkey}.bed"
        sh(f"sort -k1,1 -k2,2n {res/'features'/(fkey+'.bed')} > {s}")
        feat_sorted[fkey] = s

    rows: list[dict] = []
    for skey, sname in STRUCTURES:
        sbed = res / f"{skey}.bed"
        if not sbed.exists() or sbed.stat().st_size == 0:
            print(f"[пропуск] {sbed} отсутствует")
            continue
        struct_sorted = tmp / f"struct_{skey}.bed"
        sh(f"sort -k1,1 -k2,2n {sbed} > {struct_sorted}")

        observed = assign_fractions(struct_sorted, feat_sorted, tmp, f"obs_{skey}")

        # Реплики случайного размещения
        exp_acc: dict[str, list[float]] = {f[0]: [] for f in FEATURES}
        for i in range(args.replicates):
            shuf = tmp / f"shuf_{skey}_{i}.bed"
            sh(
                f"bedtools shuffle -i {struct_sorted} -g {args.sizes} "
                f"-seed {i + 1} -noOverlapping 2>/dev/null | "
                f"sort -k1,1 -k2,2n > {shuf}"
            )
            frac = assign_fractions(shuf, feat_sorted, tmp, f"sh_{skey}")
            for fkey in exp_acc:
                exp_acc[fkey].append(frac[fkey])
            shuf.unlink(missing_ok=True)

        for fkey, fname in FEATURES:
            exp_vals = exp_acc[fkey]
            exp_mean = sum(exp_vals) / len(exp_vals)
            n = len(exp_vals)
            sd = math.sqrt(sum((x - exp_mean) ** 2 for x in exp_vals) / n) if n else 0
            obs = observed[fkey]
            log2fc = math.log2(obs / exp_mean) if obs > 0 and exp_mean > 0 else float("nan")
            z = (obs - exp_mean) / sd if sd > 0 else float("nan")
            rows.append(
                {
                    "Структура": sname,
                    "Участок": fname,
                    "Доля_набл": round(obs, 4),
                    "Доля_ожид": round(exp_mean, 4),
                    "SD_фон": round(sd, 4),
                    "log2(набл/ожид)": round(log2fc, 3),
                    "z": round(z, 2),
                }
            )
        print(f"{sname}: фон посчитан ({args.replicates} реплик)")

    df = pd.DataFrame(rows)
    df.to_csv(res / "background.csv", index=False)
    print(df.to_string(index=False))
    print(f"\nЗаписано в {res/'background.csv'}")


if __name__ == "__main__":
    main()
