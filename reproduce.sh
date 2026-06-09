#!/bin/bash
# Воспроизведение результатов: Эпигенетика и вторичные структуры ДНК у *Acropora palmata*

set -e  # Exit on error

echo "=== 0. Зависимости ==="

apt-get install -y bedtools hmmer samtools seqkit
wget -c https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
chmod +x datasets
sudo mv datasets /usr/local/bin/

uv sync

echo "=== 1. Скачать геном + аннотацию ==="
~/.local/bin/datasets download genome accession GCF_964030605.1 \
    --include genome,gff3,protein --filename data/palmata.zip
unzip -o data/palmata.zip -d data/palmata_raw
mv data/palmata_raw/ncbi_dataset/data/GCF_964030605.1/*genomic.fna data/genome.fna
mv data/palmata_raw/ncbi_dataset/data/GCF_964030605.1/genomic.gff data/annotation.gff
mv data/palmata_raw/ncbi_dataset/data/GCF_964030605.1/protein.faa data/proteome.faa
samtools faidx data/genome.fna

echo "=== 2. Локальная HMM-база (16 семейств) ==="
uv run python src/build_pfam_db.py

echo "=== 3. Эпигенетические гены (HMMER) ==="
uv run python src/hmmer_epigenetics.py \
    --proteome data/proteome.faa \
    --gff data/annotation.gff \
    --domtbl results/hmmscan.domtbl \
    --out results/epigenetic_genes.csv

echo "=== 4. G-квадруплексы (regex) и Z-ДНК (zhunt) ==="
uv run python src/g4_scan.py --genome data/genome.fna --out results/g4.bed
uv run python src/zhunt_run.py --genome data/genome.fna --out results/zhunt.bed --procs 8

echo "=== 5. Z-DNABERT на Kaggle (GPU ~2-3 ч) ==="
uv run python src/kaggle_zdnabert.py --organism palmata dataset
uv run python src/kaggle_zdnabert.py --organism palmata push
uv run python src/kaggle_zdnabert.py --organism palmata fetch

echo "=== 6. Разметка участков + таблицы + фон ==="
uv run python src/build_features.py --gff data/annotation.gff \
    --fai data/genome.fna.fai \
    --outdir results/features --sizes data/chrom.sizes
uv run python src/distribution.py --results results
uv run python src/background.py --results results \
    --sizes data/chrom.sizes --replicates 50

echo "=== Готово! ==="
