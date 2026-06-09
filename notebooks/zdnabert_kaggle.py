from __future__ import annotations

import glob
import gzip
import os
import subprocess
import sys

_SENTINEL = "__torch_reinstalled__"


def _setup_torch():
    if os.environ.get(_SENTINEL):
        return
    try:
        import torch as _t
        if _t.cuda.is_available():
            _t.zeros(1, device="cuda").sum().item()
            print(
                f"GPU: {_t.cuda.get_device_name(0)} "
                f"sm_{''.join(map(str, _t.cuda.get_device_capability(0)))} OK",
                flush=True,
            )
            return
    except Exception as exc:
        print(f"CUDA self-test FAILED: {exc}", flush=True)

    print("Переустановка torch 2.7.1+cu126...", flush=True)
    subprocess.run(
        [sys.executable, "-m", "pip", "install", "-q",
         "torch==2.7.1", "torchvision==0.22.1",
         "--index-url", "https://download.pytorch.org/whl/cu126"],
        check=True,
    )
    os.environ[_SENTINEL] = "1"
    os.execv(sys.executable, [sys.executable, *sys.argv])


_setup_torch()

import numpy as np
import torch
from scipy import ndimage
from transformers import BertForTokenClassification, BertTokenizer

print(f"После настройки: torch={torch.__version__}, cuda={torch.cuda.is_available()}", flush=True)

# --- параметры (как в оригинальном ноутбуке) ---
THRESHOLD = 0.5          # порог уверенности модели
MIN_LEN = 10             # минимальная длина участка Z-ДНК (в токенах)
K = 6                    # размер k-mer (DNABERT 6-mer)
WIN = 512                # длина окна в токенах
PAD = 16                 # перекрытие окон
BATCH = 64               # размер батча для GPU


def find_input(patterns: list[str]) -> str:
    for pat in patterns:
        hits = sorted(glob.glob(pat, recursive=True))
        if hits:
            return hits[0]
    print("Дамп /kaggle/input:")
    for root, dirs, files in os.walk("/kaggle/input"):
        for f in files:
            print(" ", os.path.join(root, f))
    raise FileNotFoundError(f"не найдено: {patterns}")


def iter_fasta(path: str):
    """Простой потоковый парсер FASTA: yields (name, sequence_upper)."""
    opener = gzip.open if path.endswith(".gz") else open
    name = None
    chunks: list[str] = []
    with opener(path, "rt") as fh:
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


def seq2kmer(seq: str, k: int = K) -> list[str]:
    return [seq[x:x + k] for x in range(len(seq) + 1 - k)]


def split_windows(kmers: list[str], length: int = WIN, pad: int = PAD):
    """Окна по `length` токенов с перекрытием `pad`."""
    for st in range(0, len(kmers), length - pad):
        yield kmers[st:min(st + length, len(kmers))]


def stitch(pieces: list[np.ndarray], pad: int = PAD) -> np.ndarray:
    """Склейка предсказаний из окон в один массив. O(n), не O(n²).

    Промежуточный буфер аллоцируется как sum(len(p)) — это безопасная верхняя
    оценка; финальный размер sum-pad*(M-1) только если все окна не короче pad.
    Возвращаем срез до фактической позиции `off`.
    """
    if not pieces:
        return np.empty(0, dtype=np.float32)
    max_size = sum(len(p) for p in pieces)
    res = np.empty(max_size, dtype=pieces[0].dtype)
    off = 0
    for i, piece in enumerate(pieces):
        if i > 0:
            off = max(0, off - pad)
        res[off:off + len(piece)] = piece
        off += len(piece)
    return res[:off]


def predict_sequence(kmers, tokenizer, model, device) -> np.ndarray:
    """Вероятность Z-ДНК для каждого токена (6-mer) последовательности.

    Токенизируем напрямую через vocab-словарь (БезBertTokenizer.encode) —
    в 50-100 раз быстрее, иначе на 31-Mb хромосоме упирается в CPU.
    """
    vocab = tokenizer.get_vocab()
    unk_id = vocab.get(tokenizer.unk_token, 0)
    all_ids = np.fromiter(
        (vocab.get(k, unk_id) for k in kmers),
        dtype=np.int64,
        count=len(kmers),
    )

    pad_id = tokenizer.pad_token_id
    preds: list[np.ndarray] = []
    window_starts = list(range(0, len(all_ids), WIN - PAD))

    for bs in range(0, len(window_starts), BATCH):
        batch_starts = window_starts[bs:bs + BATCH]
        lengths = [min(WIN, len(all_ids) - st) for st in batch_starts]
        maxlen = max(lengths)
        arr = np.full((len(batch_starts), maxlen), pad_id, dtype=np.int64)
        mask = np.zeros_like(arr)
        for i, (st, ln) in enumerate(zip(batch_starts, lengths)):
            arr[i, :ln] = all_ids[st : st + ln]
            mask[i, :ln] = 1
        with torch.no_grad():
            out = model(
                torch.from_numpy(arr).to(device),
                attention_mask=torch.from_numpy(mask).to(device),
            )[-1]
            prob = torch.softmax(out, dim=-1)[:, :, 1].cpu().numpy()
        for i, ln in enumerate(lengths):
            preds.append(prob[i, :ln])
    return stitch(preds)


def regions_to_bed(name: str, prob: np.ndarray, idx_start: int):
    """Связные участки prob>порога -> строки BED. O(n) через find_objects."""
    labeled, n = ndimage.label(prob > THRESHOLD)
    slices = ndimage.find_objects(labeled)
    lines = []
    for sl_tuple in slices:
        if sl_tuple is None:
            continue
        sl = sl_tuple[0]
        length = sl.stop - sl.start
        if length <= MIN_LEN:
            continue
        start = int(sl.start)
        end = int(sl.stop - 1) + K
        score = float(prob[sl].max())
        idx_start += 1
        lines.append(
            f"{name}\t{start}\t{end}\tZDNABERT_{name}_{idx_start}\t{score:.3f}\t."
        )
    return lines, idx_start


def main() -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}", flush=True)

    model_dir = find_input([
        "/kaggle/input/**/pytorch_model.bin",
    ])
    model_dir = os.path.dirname(model_dir)
    print(f"веса модели: {model_dir}", flush=True)
    print(f"  файлы: {sorted(os.listdir(model_dir))}", flush=True)

    genome = find_input([
        "/kaggle/input/**/*.fna.gz",
        "/kaggle/input/**/*.fna",
        "/kaggle/input/**/*.fa.gz",
        "/kaggle/input/**/*.fa",
    ])
    print(f"геном: {genome}", flush=True)

    tokenizer = BertTokenizer.from_pretrained(model_dir)
    model = BertForTokenClassification.from_pretrained(model_dir).to(device)
    model.eval()

    out_path = "/kaggle/working/zdnabert.bed"
    counter = 0
    import time
    with open(out_path, "w") as out:
        for name, seq in iter_fasta(genome):
            if len(seq) < K:
                continue
            t0 = time.time()
            print(f"  {name}: длина {len(seq)} нт, обработка...", flush=True)
            kmers = seq2kmer(seq)
            prob = predict_sequence(kmers, tokenizer, model, device)
            lines, counter = regions_to_bed(name, prob, counter)
            if lines:
                out.write("\n".join(lines) + "\n")
            out.flush()
            print(
                f"  {name}: {len(lines)} участков Z-ДНК ({time.time() - t0:.0f}s)",
                flush=True,
            )

    print(f"Готово. Всего участков: {counter}. Файл: {out_path}", flush=True)


if __name__ == "__main__":
    main()
