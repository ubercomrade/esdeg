# ESDEG

ESDEG is a Python command-line tool for testing transcription-factor motif enrichment in promoters of differentially expressed genes (DEGs) or an arbitrary gene set. It scans promoters with [Mimosa](https://github.com/ubercomrade/mimosa), compares foreground promoters with a background matched on GC content and CpG depletion, and reports ROC AUC enrichment with analytical Mann–Whitney p-values and FDR correction. Bundled HOCOMOCO v12 and JASPAR 2024 models are supported.

## Background

Transcription factors regulate gene expression by binding transcription-factor binding sites (TFBSs) in regulatory DNA. TFBSs are variable and are commonly represented by position weight matrices (PWMs). If a TF contributes to an expression response, its motif can be over-represented among promoters of response genes.

Promoter GC content and CpG depletion affect motif scores, especially for GC-rich or CpG-containing motifs. ESDEG therefore compares a foreground with a promoter background matched on both, rather than with an arbitrary set of genes.

## Installation

ESDEG requires Python 3.12 or 3.13 and [uv](https://docs.astral.sh/uv/).

```bash
git clone https://github.com/ubercomrade/esdeg.git
cd esdeg
uv sync --locked
```

For development tools as well:

```bash
uv sync --locked --group dev
```

Run `uv run esdeg --help` to see the complete command-line reference.

## CLI quickstart

Test HOCOMOCO motifs for enrichment in a DEG table:

```bash
uv run esdeg enrichment \
  example/promoters.p400m100.fa \
  example/E-MTAB-6598.degs.csv \
  enrichment.tsv \
  --motifs hocomoco \
  --seed 0 --match-ratio 5
```

Test a one-ID-per-line gene set instead:

```bash
uv run esdeg enrichment \
  example/promoters.p400m100.fa \
  example/HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
  set-enrichment.tsv \
  --format set --motifs hocomoco
```

For JASPAR, select its database and a taxonomic group:

```bash
uv run esdeg enrichment promoters.fa results.csv enrichment.tsv \
  --motifs jaspar --taxon plants
```

## Inputs

### Promoters

Provide a non-empty FASTA file with one unique promoter per gene. ESDEG uses the text before the first colon in each FASTA header as the gene ID, so `>ENSG000001:coordinates` becomes `ENSG000001`. IDs must match the input table or gene set. Every sequence must be at least as long as the longest selected motif; use promoters defined by the same genomic rule for meaningful comparisons.

### Differential-expression table

The `deg` format accepts CSV or TSV with these required columns:

| Column | Meaning |
| --- | --- |
| `id` | Gene ID matching the promoter FASTA. |
| `log2FoldChange` | Estimated expression log2 fold change. |
| `padj` | Adjusted differential-expression p-value. |

By default, foreground genes satisfy `padj <= 0.05` and `abs(log2FoldChange) >= 1`. The `--regulated` option selects `all`, `up`, or `down` foreground genes. The initial background pool contains non-DEGs with `padj > 0.05` and `abs(log2FoldChange) <= log2(5/4)`.

### Gene set

The `set` format is a plain text file with one gene ID per line. Blank lines and duplicate IDs are ignored. The background pool is all promoter IDs not present in the set.

## Enrichment

```text
esdeg enrichment PROMOTERS INPUT OUTPUT [options]
```

`PROMOTERS` is the FASTA file, `INPUT` is a DEG table or gene-set file, and `OUTPUT` is the enrichment TSV.

| Option | Default | Purpose |
| --- | --- | --- |
| `--format {deg,set}` | `deg` | Choose a DEG table or a one-ID-per-line gene set. |
| `--motifs {jaspar,hocomoco}` | `hocomoco` | Use bundled HOCOMOCO v12 or load JASPAR 2024 models. |
| `--model PATH` | — | Load a custom Mimosa model; repeat for multiple models. Cannot be combined with `--motifs`. |
| `--taxon` | `vertebrates` | JASPAR taxonomic group: `plants`, `vertebrates`, `insects`, `urochordates`, `nematodes`, or `fungi`. |
| `--regulated {all,up,down}` | `all` | Select the DEG direction. |
| `--pvalue` | `0.05` | Foreground DEG `padj` threshold. |
| `--log2fc-deg` | `1.0` | Absolute foreground log2 fold-change threshold. |
| `--log2fc-back` | `log2(5/4)` | Absolute background log2 fold-change ceiling. |
| `--match-ratio` | `5` | Maximum background-to-foreground size ratio after matching. |
| `--seed` | `0` | Reproducible seed for background matching. |
| `--background-output PATH` | — | Write the IDs of the selected background promoters to a TSV. |
| `--visualization PATH` | — | Write an interactive HTML scatter plot. |
| `--report PATH` | — | Write an interactive HTML results table. |
| `--xlsx PATH` | — | Write an XLSX table with available motif logos. |

## Method and statistics

For each motif, ESDEG scans every promoter on the best strand and retains the maximum PWM score across positions. A high score means that a promoter contains a stronger match to that motif.

The background used for testing has at most `match-ratio × N_foreground` promoters. ESDEG describes every promoter by its GC fraction and by its log2 CpG observed/expected ratio, `log2((CG + 1) / (C·G / L + 1))`. Each foreground promoter in turn takes the nearest unused candidate from the background pool by Mahalanobis distance on these two features, and this repeats until the target size is reached, so no promoter is used twice. Matching is performed once for the run. The log reports the mean of each feature in the foreground, the pool, and the selected background, together with the standardized mean difference (SMD); a warning is issued when `|SMD| > 0.1`.

For each motif, the foreground and selected background score distributions produce the effect size `auc_roc`, the ROC AUC. A value of `0.5` indicates no score separation; values above `0.5` indicate enrichment of high motif scores in the foreground. It equals the Mann–Whitney statistic `U / (N_foreground · N_background)`.

Significance is the one-sided (enrichment only) p-value of the Mann–Whitney test, from the normal approximation of `U` with tie and continuity corrections. The same method is used for every group size, so p-values for very small groups are approximate. The null hypothesis is that foreground and background scores come from the same distribution, so p-values are conditional on the selected background and are valid only to the extent that the background matches the foreground in promoter composition.

ESDEG applies the Benjamini–Hochberg procedure across all tested motifs. The enrichment TSV contains every motif; it does not automatically remove non-significant results. The optional plot displays motifs with `p_value_roc_adj < 0.05`.

## Results

The enrichment TSV is sorted by decreasing `auc_roc` and contains:

| Column | Meaning |
| --- | --- |
| `motif_id` | Motif model identifier. |
| `tf_name`, `tf_class`, `tf_family` | Motif annotation. |
| `auc_roc` | ROC AUC effect size. |
| `p_value_roc` | One-sided Mann–Whitney p-value. |
| `p_value_roc_adj` | Benjamini–Hochberg adjusted p-value. |
| `jaspar_cluster` | JASPAR motif cluster, included only for bundled JASPAR runs. |

## Expression annotation

Add TF expression information to an enrichment table:

```bash
uv run esdeg annotation enrichment.tsv deg.csv counts.csv genes.gtf annotated.tsv \
  --filter --min-auc 0.5
```

The DEG table needs `id`, `log2FoldChange`, and `padj`; the counts table needs `id` and `counts`; the GTF supplies `gene_id` to `gene_name` mapping. Missing expression values remain `NaN`.

With `--filter`, a row must satisfy `p_value_roc_adj < --me-padj` (default `0.05`) and `auc_roc >= --min-auc` (default `0.5`). Its TF must also be differentially expressed according to `--de-padj` and `--log-fc`, or have counts at least `--ncounts`. Use `--best` to retain the best motif per TF after annotation.

## Development

```bash
uv run pytest -q
uv run ruff check .
```

## References

Oshchepkov D, Chadaeva I, Kozhemyakina R, et al. Transcription Factors as Important Regulators of Changes in Behavior through Domestication of Gray Rats: Quantitative Data from RNA Sequencing. *International Journal of Molecular Sciences*. 2022;23(20):12269. [https://doi.org/10.3390/ijms232012269](https://doi.org/10.3390/ijms232012269)

## License

MIT. See [LICENSE](LICENSE).
