"""Сборка локальной HMM-базы Pfam только из нужных семейств.

HMM-профили скачиваются по одному с InterPro, объединяются в один файл и
индексируются hmmpress. Так не нужно качать полную Pfam-A (~1.5 ГБ).
"""
from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).parent))
from epigenetic_families import ALL_FAMILIES

INTERPRO_HMM = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{acc}?annotation=hmm"
OUT_DIR = Path("data/pfam")
DB_PATH = OUT_DIR / "epigenetics.hmm"


def fetch_hmm(acc: str) -> bytes:
    """Скачать один HMM-профиль Pfam (InterPro отдаёт gzip)."""
    resp = requests.get(INTERPRO_HMM.format(acc=acc), timeout=60)
    resp.raise_for_status()
    data = resp.content
    if data[:2] == b"\x1f\x8b":  # gzip-сигнатура
        data = gzip.decompress(data)
    if not data.startswith(b"HMMER"):
        raise ValueError(f"{acc}: ответ не похож на HMM-профиль")
    return data


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ok = 0
    with open(DB_PATH, "wb") as db:
        for acc in ALL_FAMILIES:
            try:
                hmm = fetch_hmm(acc)
            except (requests.HTTPError, ValueError) as exc:
                print(f"  ПРОПУСК {acc}: {exc}", file=sys.stderr)
                continue
            db.write(hmm)
            if not hmm.endswith(b"\n"):
                db.write(b"\n")
            ok += 1
            print(f"  {acc}: {ALL_FAMILIES[acc]}")
    print(f"HMM-база собрана: {DB_PATH} ({ok}/{len(ALL_FAMILIES)} семейств)")

    subprocess.run(["hmmpress", "-f", str(DB_PATH)], check=True)
    print("hmmpress: индексы построены")


if __name__ == "__main__":
    main()
