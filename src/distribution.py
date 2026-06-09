"""Распределение вторичных структур ДНК по участкам генома.

Строит две таблицы (по заданию, п. 7):
  * Таблица 1 — число и доля СТРУКТУР, попавших в каждый тип участка.
    Каждая структура относится ровно к одному участку по приоритету:
    Экзон > Интрон > Промотор > Downstream > Межгенный.
  * Таблица 2 — число и доля УЧАСТКОВ, в которых есть хотя бы одна структура.

Использует bedtools intersect.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd

# Участки в порядке приоритета для отнесения структуры (Таблица 1).
FEATURES = [
    ("exons", "Экзоны"),
    ("introns", "Интроны"),
    ("promoters", "Промоторы (1000 п.н. до TSS)"),
    ("downstream", "Downstream (200 п.н.)"),
    ("intergenic", "Межгенные участки"),
]
STRUCTURES = [
    ("g4", "Квадруплексы"),
    ("zhunt", "Zhunt"),
    ("zdnabert", "ZDNABERT"),
]


def sh(cmd: str) -> str:
    return subprocess.run(
        cmd, shell=True, check=True, capture_output=True, text=True
    ).stdout


def sort_bed(src: Path, dst: Path) -> None:
    sh(f"sort -k1,1 -k2,2n {src} > {dst}")


def count_lines(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    return int(sh(f"wc -l < {path}").strip())


def main() -> None:
    ap = argparse.ArgumentParser(description="Таблицы распределения структур")
    ap.add_argument("--results", default="results", type=Path)
    args = ap.parse_args()
    res = args.results
    feat_dir = res / "features"
    tmp = res / "_tmp"
    tmp.mkdir(parents=True, exist_ok=True)

    # Отсортированные BED участков
    feat_sorted: dict[str, Path] = {}
    for key, _ in FEATURES:
        s = tmp / f"feat_{key}.bed"
        sort_bed(feat_dir / f"{key}.bed", s)
        feat_sorted[key] = s

    table1_rows: list[dict] = []
    table2_rows: list[dict] = []

    for skey, sname in STRUCTURES:
        sbed = res / f"{skey}.bed"
        if not sbed.exists() or sbed.stat().st_size == 0:
            print(f"[пропуск] {sbed} отсутствует")
            continue
        struct_sorted = tmp / f"struct_{skey}.bed"
        sort_bed(sbed, struct_sorted)
        total_struct = count_lines(struct_sorted)

        # --- Таблица 1: приоритетное отнесение каждой структуры ---
        remaining = struct_sorted
        for fkey, fname in FEATURES:
            inside = tmp / f"{skey}_in_{fkey}.bed"
            sh(
                f"bedtools intersect -u -a {remaining} -b {feat_sorted[fkey]} "
                f"> {inside}"
            )
            rest = tmp / f"{skey}_not_{fkey}.bed"
            sh(
                f"bedtools intersect -v -a {remaining} -b {feat_sorted[fkey]} "
                f"> {rest}"
            )
            n = count_lines(inside)
            table1_rows.append(
                {
                    "Участок": fname,
                    "Структура": sname,
                    "Число": n,
                    "Доля": round(n / total_struct, 4) if total_struct else 0.0,
                }
            )
            remaining = rest

        # --- Таблица 2: участки с >=1 структурой ---
        for fkey, fname in FEATURES:
            total_feat = count_lines(feat_sorted[fkey])
            with_struct = tmp / f"{fkey}_has_{skey}.bed"
            sh(
                f"bedtools intersect -u -a {feat_sorted[fkey]} -b {struct_sorted} "
                f"> {with_struct}"
            )
            n = count_lines(with_struct)
            table2_rows.append(
                {
                    "Участок": fname,
                    "Структура": sname,
                    "Число участков": n,
                    "Всего участков": total_feat,
                    "Доля участков": round(n / total_feat, 4) if total_feat else 0.0,
                }
            )

    # Сводим в широкие таблицы (структуры — колонки)
    t1 = pd.DataFrame(table1_rows)
    t2 = pd.DataFrame(table2_rows)
    if not t1.empty:
        wide1 = t1.pivot(index="Участок", columns="Структура",
                         values=["Число", "Доля"])
        wide1 = wide1.reindex([f[1] for f in FEATURES])
        wide1.to_csv(res / "distribution_table1.csv")
        print("=== Таблица 1: число и доля структур по участкам ===")
        print(wide1.to_string())
    if not t2.empty:
        wide2 = t2.pivot(index="Участок", columns="Структура",
                         values=["Число участков", "Доля участков"])
        wide2 = wide2.reindex([f[1] for f in FEATURES])
        wide2.to_csv(res / "distribution_table2.csv")
        print("\n=== Таблица 2: доля участков со структурой ===")
        print(wide2.to_string())

    print(f"\nЗаписано: {res/'distribution_table1.csv'}, {res/'distribution_table2.csv'}")


if __name__ == "__main__":
    main()
