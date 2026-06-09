# data/

## Геном *Acropora palmata*

GitHub не позволяет хранить файлы >100 МБ, поэтому сам FASTA-файл генома не закоммичен. Скачать его можно командой ниже:

```bash
~/.local/bin/datasets download genome accession GCF_964030605.1 \
    --include genome,gff3,protein --filename data/palmata.zip
unzip -o data/palmata.zip -d data/palmata_raw
mv data/palmata_raw/ncbi_dataset/data/GCF_964030605.1/*genomic.fna data/genome.fna
mv data/palmata_raw/ncbi_dataset/data/GCF_964030605.1/genomic.gff data/annotation.gff
mv data/palmata_raw/ncbi_dataset/data/GCF_964030605.1/protein.faa data/proteome.faa
samtools faidx data/genome.fna
```

## Что уже лежит в репозитории

| Файл | Что это | Размер |
|---|---|---|
| `annotation.gff.gz` | GFF от NCBI RefSeq (gzip) | 11 МБ |
| `proteome.faa.gz` | белки от NCBI RefSeq (gzip) | 7.5 МБ |
| `chrom.sizes` | длины контигов (`samtools faidx ...`) | <10 КБ |
| `pfam/` | локальная HMM-база (16 семейств, собрана `src/build_pfam_db.py`) | 2 МБ |

Перед запуском пайплайна распаковать gzip:

```bash
gunzip data/annotation.gff.gz data/proteome.faa.gz
```

## NCBI metadata

- **Accession:** GCF_964030605.1 (RefSeq)
- **Уровень сборки:** chromosome (244 скаффолда, самый длинный 31.8 Мб)
- **Размер:** 334 Мб; N50 = 22.6 Мб
- **Аннотация:** 36 946 генов / 34 126 белков
- **GC-содержание:** 39 %
- **Ссылка:** https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_964030605.1/
