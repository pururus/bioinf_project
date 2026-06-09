"""Поиск эпигенетических генов: hmmscan протеома по локальной базе Pfam.

Шаги:
  1. hmmscan: protein.faa против data/pfam/epigenetics.hmm -> domtblout;
  2. фильтр попаданий по E-value полной последовательности;
  3. привязка белок -> ген (по GFF: CDS.protein_id -> mRNA -> gene);
  4. таблица results/epigenetic_genes.csv.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from epigenetic_families import ALL_FAMILIES

EVALUE_MAX = 1e-5


def parse_gff_protein_to_gene(gff: Path) -> dict[str, dict]:
    """protein_id -> {gene_id, gene_name, seqid, start, end, strand}.

    Поддерживает два стиля GFF:
    * NCBI/RefSeq — у CDS есть атрибут protein_id (XP_*);
    * funannotate — protein_id отсутствует, идентификатор белка совпадает
      с ID мРНК (FUN_*-T1).
    """
    mrna_to_gene: dict[str, str] = {}
    gene_info: dict[str, dict] = {}
    cds_protein: dict[str, str] = {}  # mRNA_id -> protein_id (только NCBI)
    mrnas_seen: set[str] = set()

    attr_re = re.compile(r"(\w+)=([^;]+)")
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            attrs = dict(attr_re.findall(f[8]))
            if ftype == "gene":
                gid = attrs.get("ID", "")
                gene_info[gid] = {
                    "gene_id": gid,
                    "gene_name": attrs.get("gene", attrs.get("Name", gid)),
                    "seqid": f[0],
                    "start": int(f[3]),
                    "end": int(f[4]),
                    "strand": f[6],
                }
            elif ftype in ("mRNA", "transcript", "lnc_RNA"):
                mid = attrs.get("ID", "")
                parent = attrs.get("Parent", "")
                if mid and parent:
                    mrna_to_gene[mid] = parent
                    mrnas_seen.add(mid)
            elif ftype == "CDS":
                pid = attrs.get("protein_id")
                parent = attrs.get("Parent", "")
                if pid and parent:
                    cds_protein.setdefault(parent, pid)

    # Fallback: если у CDS нет protein_id (funannotate),
    # считаем ID мРНК идентификатором белка.
    for mid in mrnas_seen:
        cds_protein.setdefault(mid, mid)

    protein_to_gene: dict[str, dict] = {}
    for mid, pid in cds_protein.items():
        gid = mrna_to_gene.get(mid)
        if gid and gid in gene_info:
            protein_to_gene[pid] = gene_info[gid]
    return protein_to_gene


def run_hmmscan(hmm_db: Path, proteome: Path, domtbl: Path, cpu: int) -> None:
    domtbl.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "hmmscan",
        "--cpu", str(cpu),
        "--domtblout", str(domtbl),
        "-E", str(EVALUE_MAX),
        str(hmm_db),
        str(proteome),
    ]
    print("Запуск:", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def parse_domtbl(domtbl: Path) -> pd.DataFrame:
    """Прочитать domtblout hmmscan, оставить лучшее попадание на (белок, семейство)."""
    rows = []
    with open(domtbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 22:
                continue
            family_acc = f[1].split(".")[0]  # 'PF00125.27' -> 'PF00125'
            rows.append(
                {
                    "family_pfam": family_acc,
                    "family_name": ALL_FAMILIES.get(family_acc, f[0]),
                    "protein_id": f[3],
                    "evalue": float(f[6]),
                    "score": float(f[7]),
                }
            )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df[df["evalue"] < EVALUE_MAX]
    df = df.sort_values("evalue").drop_duplicates(["protein_id", "family_pfam"])
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description="Поиск эпигенетических генов (hmmscan)")
    ap.add_argument("--hmm", default="data/pfam/epigenetics.hmm", type=Path)
    ap.add_argument("--proteome", default="data/proteome.faa", type=Path)
    ap.add_argument("--gff", default="data/annotation.gff", type=Path)
    ap.add_argument("--domtbl", default="results/hmmscan.domtbl", type=Path)
    ap.add_argument("--out", default="results/epigenetic_genes.csv", type=Path)
    ap.add_argument("--cpu", default=4, type=int)
    ap.add_argument("--skip-hmmscan", action="store_true", help="использовать готовый domtbl")
    args = ap.parse_args()

    if not args.skip_hmmscan:
        run_hmmscan(args.hmm, args.proteome, args.domtbl, args.cpu)

    hits = parse_domtbl(args.domtbl)
    if hits.empty:
        print("Попаданий не найдено", file=sys.stderr)
        return

    p2g = parse_gff_protein_to_gene(args.gff)
    for col in ("gene_id", "gene_name", "seqid", "start", "end", "strand"):
        hits[col] = hits["protein_id"].map(
            lambda p, c=col: p2g.get(p, {}).get(c, "")
        )

    hits = hits.sort_values(["family_pfam", "evalue"])
    cols = [
        "family_pfam", "family_name", "gene_id", "gene_name",
        "protein_id", "seqid", "start", "end", "strand", "evalue", "score",
    ]
    hits[cols].to_csv(args.out, index=False)

    n_fam = hits["family_pfam"].nunique()
    print(f"Найдено попаданий: {len(hits)}; семейств с генами: {n_fam}")
    summary = hits.groupby(["family_pfam", "family_name"]).size()
    for (acc, name), n in summary.items():
        print(f"  {acc}  {name}: {n} ген(ов)")
    print(f"Таблица записана в {args.out}")


if __name__ == "__main__":
    main()
