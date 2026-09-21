# ESDEG

ESDEG tests motif enrichment of differential-expression genes or a gene set against a GC-matched promoter background. Motif conversion and scanning use [Mimosa](https://github.com/ubercomrade/mimosa); bundled HOCOMOCO and JASPAR database models are supported.

## Install

```bash
uv sync --group dev
```

## Enrichment

```bash
esdeg enrichment promoters.fa results.csv enrichment.tsv \
  --format deg --n-permutations 10000 --seed 0 --match-ratio 5
```

Input DEG tables require `id`, `log2FoldChange`, and `padj`. CSV and TSV are accepted. For a one-ID-per-line set, use `--format set`.

## Annotation

```bash
esdeg annotation enrichment.tsv deg.csv counts.csv genes.gtf annotated.tsv \
  --filter --min-auc 0.5
```

Annotation uses ROC AUC enrichment and keeps missing expression values as `NaN`.

## Development

```bash
uv run pytest
uv run ruff check .
```
