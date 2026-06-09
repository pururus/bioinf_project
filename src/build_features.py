"""Извлечение участков генома из GFF-аннотации.

Строит пять BED-файлов в {outdir}:
  * exons.bed       — экзоны
  * introns.bed     — интроны (тело гена минус экзоны)
  * promoters.bed   — 1000 п.н. вверх от TSS (с учётом стренда)
  * downstream.bed  — 200 п.н. вниз от конца гена (с учётом стренда)
  * intergenic.bed  — всё остальное (геном минус ген ∪ промотор ∪ downstream)

GFF читаем построчно (pyranges.read_gff3 неприемлемо медленный на больших
GFF). Интервальные операции (merge/subtract) делаются через bedtools.
Пишет также data/{organism}/chrom.sizes (нужен для bedtools).
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

PROMOTER_UP = 1000
DOWNSTREAM_LEN = 200


def chrom_sizes(fai: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with open(fai) as fh:
        for line in fh:
            f = line.split("\t")
            sizes[f[0]] = int(f[1])
    return sizes


def parse_features(gff: Path) -> tuple[list, list]:
    """Прочитать GFF -> (genes [(seqid,start,end,strand)], exons [(seqid,start,end)])."""
    genes: list[tuple[str, int, int, str]] = []
    exons: list[tuple[str, int, int]] = []
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            ftype = f[2]
            if ftype == "gene":
                # BED 0-based: GFF 1-based -> start-1, end остаётся
                genes.append((f[0], int(f[3]) - 1, int(f[4]), f[6]))
            elif ftype == "exon":
                exons.append((f[0], int(f[3]) - 1, int(f[4])))
    return genes, exons


def sort_uniq_bed(rows: list, path: Path) -> None:
    rows = sorted(set((r[0], r[1], r[2]) for r in rows))
    with open(path, "w") as fh:
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")


def bedtools_merge(src: Path, dst: Path) -> int:
    """Слить перекрывающиеся интервалы."""
    cmd = f"sort -k1,1 -k2,2n {src} | bedtools merge -i - > {dst}"
    subprocess.run(cmd, shell=True, check=True)
    return int(
        subprocess.run(["wc", "-l", str(dst)], capture_output=True, text=True).stdout.split()[0]
    )


def bedtools_subtract(a: Path, b: Path, dst: Path) -> int:
    cmd = f"bedtools subtract -a {a} -b {b} | sort -k1,1 -k2,2n > {dst}"
    subprocess.run(cmd, shell=True, check=True)
    return int(
        subprocess.run(["wc", "-l", str(dst)], capture_output=True, text=True).stdout.split()[0]
    )


def main() -> None:
    ap = argparse.ArgumentParser(description="BED-файлы участков генома из GFF")
    ap.add_argument("--gff", default="data/annotation.gff", type=Path)
    ap.add_argument("--fai", default="data/genome.fna.fai", type=Path)
    ap.add_argument("--outdir", default="results/features", type=Path)
    ap.add_argument("--sizes", default="data/genome.chrom.sizes", type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    sizes = chrom_sizes(args.fai)
    with open(args.sizes, "w") as fh:
        for name, length in sizes.items():
            fh.write(f"{name}\t{length}\n")

    print("Парсинг GFF...", flush=True)
    genes_raw, exons_raw = parse_features(args.gff)
    print(f"  gene: {len(genes_raw)}, exon: {len(exons_raw)}", flush=True)

    # Промоторы (1000 п.н. вверх от TSS, стренд) и downstream (200 п.н. за концом)
    promoters_raw: list[tuple[str, int, int]] = []
    downstream_raw: list[tuple[str, int, int]] = []
    for seqid, start, end, strand in genes_raw:
        chrlen = sizes.get(seqid, end)
        if strand == "+":
            p_s, p_e = max(0, start - PROMOTER_UP), start
            d_s, d_e = end, min(chrlen, end + DOWNSTREAM_LEN)
        else:
            p_s, p_e = end, min(chrlen, end + PROMOTER_UP)
            d_s, d_e = max(0, start - DOWNSTREAM_LEN), start
        if p_s < p_e:
            promoters_raw.append((seqid, p_s, p_e))
        if d_s < d_e:
            downstream_raw.append((seqid, d_s, d_e))

    tmp = Path(tempfile.mkdtemp(prefix="build_features_"))
    try:
        # Сырые BED-ы
        sort_uniq_bed([(g[0], g[1], g[2]) for g in genes_raw], tmp / "genes.raw.bed")
        sort_uniq_bed(exons_raw, tmp / "exons.raw.bed")
        sort_uniq_bed(promoters_raw, tmp / "promoters.raw.bed")
        sort_uniq_bed(downstream_raw, tmp / "downstream.raw.bed")

        # Финальные слитые
        n_exons = bedtools_merge(tmp / "exons.raw.bed", args.outdir / "exons.bed")
        n_genes_merged = bedtools_merge(tmp / "genes.raw.bed", tmp / "genes.merged.bed")
        n_introns = bedtools_subtract(
            tmp / "genes.merged.bed", args.outdir / "exons.bed", args.outdir / "introns.bed"
        )
        n_proms = bedtools_merge(tmp / "promoters.raw.bed", args.outdir / "promoters.bed")
        n_down = bedtools_merge(tmp / "downstream.raw.bed", args.outdir / "downstream.bed")

        # Межгенник: геном - (гены ∪ промоторы ∪ downstream)
        with open(tmp / "occupied.raw.bed", "w") as fh:
            for path in (
                tmp / "genes.merged.bed",
                args.outdir / "promoters.bed",
                args.outdir / "downstream.bed",
            ):
                fh.write(path.read_text())
        bedtools_merge(tmp / "occupied.raw.bed", tmp / "occupied.bed")

        # Полный геном как BED
        with open(tmp / "whole.bed", "w") as fh:
            for name, length in sizes.items():
                fh.write(f"{name}\t0\t{length}\n")
        sort_uniq = subprocess.run(
            f"sort -k1,1 -k2,2n {tmp / 'whole.bed'} > {tmp / 'whole.sorted.bed'}",
            shell=True, check=True,
        )
        n_inter = bedtools_subtract(
            tmp / "whole.sorted.bed", tmp / "occupied.bed", args.outdir / "intergenic.bed"
        )

        print(f"exons       : {n_exons}")
        print(f"introns     : {n_introns}")
        print(f"promoters   : {n_proms}")
        print(f"downstream  : {n_down}")
        print(f"intergenic  : {n_inter}")
        print(f"genes_merged: {n_genes_merged} (для справки)")
    finally:
        subprocess.run(["rm", "-rf", str(tmp)], check=False)


if __name__ == "__main__":
    main()
