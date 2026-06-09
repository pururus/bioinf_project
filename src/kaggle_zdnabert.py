"""Оркестрация запуска Z-DNABERT на Kaggle через API (без ручного копирования).

Подкоманды:
  dataset  — упаковать геном и залить/обновить как приватный датасет Kaggle;
  push     — собрать kernel-metadata.json и запустить ядро (script, GPU);
  status   — показать статус ядра;
  fetch    — скачать результат (zdnabert.bed) в results/.

Требуется настроенный Kaggle API (~/.kaggle/kaggle.json либо переменные
окружения KAGGLE_USERNAME / KAGGLE_KEY).
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

NOTEBOOK = Path("notebooks")


def paths_for(organism: str) -> dict:
    """Стандартные пути и slug-и для конкретного организма."""
    return {
        "dataset_slug": f"acropora-{organism}-genome",
        "kernel_slug": f"zdnabert-acropora-{organism}",
        "genome": Path(f"data/{organism}/genome.fna"),
        "build": Path(f"data/{organism}/kaggle_dataset"),
        "out_dir": Path(f"results/{organism}"),
    }


def kaggle_username() -> str:
    """Имя пользователя Kaggle: env, kaggle.json или kaggle config view."""
    import os

    if os.environ.get("KAGGLE_USERNAME"):
        return os.environ["KAGGLE_USERNAME"]
    for cand in (Path.home() / ".kaggle/kaggle.json",
                 Path(os.environ.get("KAGGLE_CONFIG_DIR", "")) / "kaggle.json"):
        if cand.is_file():
            return json.loads(cand.read_text())["username"]
    out = subprocess.run(
        ["kaggle", "config", "view"], capture_output=True, text=True, check=True
    ).stdout
    for line in out.splitlines():
        if "username:" in line:
            return line.split(":", 1)[1].strip()
    raise SystemExit("Не удалось определить имя пользователя Kaggle")


def run(cmd: list[str]) -> None:
    print("$", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


def cmd_dataset(args: argparse.Namespace) -> None:
    user = kaggle_username()
    p = paths_for(args.organism)
    p["build"].mkdir(parents=True, exist_ok=True)
    gz = p["build"] / "genome.fna.gz"
    if not gz.exists() or args.force:
        print("Упаковка генома в gzip...", file=sys.stderr)
        with open(p["genome"], "rb") as src, gzip.open(gz, "wb") as dst:
            shutil.copyfileobj(src, dst, length=1 << 20)
    meta = {
        "title": f"Acropora {args.organism} genome",
        "id": f"{user}/{p['dataset_slug']}",
        "licenses": [{"name": "CC0-1.0"}],
    }
    (p["build"] / "dataset-metadata.json").write_text(json.dumps(meta, indent=2))

    exists = subprocess.run(
        ["kaggle", "datasets", "status", f"{user}/{p['dataset_slug']}"],
        capture_output=True, text=True,
    ).returncode == 0
    if exists:
        run(["kaggle", "datasets", "version", "-p", str(p["build"]),
             "-m", "update genome", "--dir-mode", "zip"])
    else:
        run(["kaggle", "datasets", "create", "-p", str(p["build"]), "--dir-mode", "zip"])
    print(f"Датасет: https://www.kaggle.com/datasets/{user}/{p['dataset_slug']}")


WEIGHTS_DATASET = "zdnabert-hg-kouzine-weights"


def cmd_push(args: argparse.Namespace) -> None:
    user = kaggle_username()
    p = paths_for(args.organism)
    meta = {
        "id": f"{user}/{p['kernel_slug']}",
        "title": f"ZDNABERT Acropora {args.organism}",
        "code_file": "zdnabert_kaggle.py",
        "language": "python",
        "kernel_type": "script",
        "is_private": True,
        "enable_gpu": True,
        "enable_internet": True,
        "dataset_sources": [
            f"{user}/{p['dataset_slug']}",
            f"{user}/{WEIGHTS_DATASET}",
        ],
        "competition_sources": [],
        "kernel_sources": [],
    }
    (NOTEBOOK / "kernel-metadata.json").write_text(json.dumps(meta, indent=2))
    run(["kaggle", "kernels", "push", "-p", str(NOTEBOOK)])
    print(f"Ядро: https://www.kaggle.com/code/{user}/{p['kernel_slug']}")


def cmd_status(args: argparse.Namespace) -> None:
    user = kaggle_username()
    p = paths_for(args.organism)
    run(["kaggle", "kernels", "status", f"{user}/{p['kernel_slug']}"])


def cmd_fetch(args: argparse.Namespace) -> None:
    user = kaggle_username()
    p = paths_for(args.organism)
    p["out_dir"].mkdir(parents=True, exist_ok=True)
    run(["kaggle", "kernels", "output", f"{user}/{p['kernel_slug']}",
         "-p", str(p["out_dir"])])
    print(f"Результат скачан в {p['out_dir']}/ (ожидается zdnabert.bed)")


def main() -> None:
    ap = argparse.ArgumentParser(description="Запуск Z-DNABERT на Kaggle")
    ap.add_argument(
        "--organism",
        required=True,
        help="имя организма-подкаталога (например cervicornis, palmata)",
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("dataset", help="залить геном датасетом")
    p.add_argument("--force", action="store_true", help="пересобрать gzip")
    p.set_defaults(func=cmd_dataset)
    sub.add_parser("push", help="запустить ядро").set_defaults(func=cmd_push)
    sub.add_parser("status", help="статус ядра").set_defaults(func=cmd_status)
    sub.add_parser("fetch", help="скачать результат").set_defaults(func=cmd_fetch)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
