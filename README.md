# ESDEG

ESDEG is a Python command-line tool for testing transcription-factor motif enrichment in promoters of differentially expressed genes (DEGs) or an arbitrary gene set. It scans promoters with [Mimosa](https://github.com/ubercomrade/mimosa), compares foreground promoters with a background matched on GC content and CpG depletion, and reports auPRC enrichment with permutation p-values calibrated to the matched design and FDR correction. Bundled HOCOMOCO v12 and JASPAR 2024 models are supported.

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
  --seed 0
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
| `--match-ratio` | `5` | Maximum background-to-foreground size ratio. It is lowered one step at a time (down to 1) while the standardized mean difference of GC or CpG exceeds 0.1. |
| `--seed` | `0` | Reproducible seed for background matching, tie-breaking and permutations. |
| `--n-permutations` | `100000` | Size of the shared permutation sample whose moments define the null; p-values come from a Cornish–Fisher fit to it, so they are not limited by `1 / n`. |
| `--background-output PATH` | — | Write the IDs of the selected background promoters to a TSV. |
| `--visualization PATH` | — | Write an interactive HTML scatter plot. |
| `--report PATH` | — | Write an interactive HTML results table. |
| `--xlsx PATH` | — | Write an XLSX table with available motif logos. |

## Method and statistics

For each motif, ESDEG scans every promoter on the best strand and scores each promoter by the maximum motif score over all positions.

The background used for testing has at most `match-ratio × N_foreground` promoters. ESDEG describes every promoter by its GC fraction and by its log2 CpG observed/expected ratio, `log2((CG + 1) / (C·G / L + 1))`. Each foreground promoter in turn takes the nearest unused candidate from the background pool by Mahalanobis distance on these two features, and this repeats until the target size is reached, so no promoter is used twice. Matching is performed once for the run. The log reports the mean of each feature in the foreground, the pool, and the selected background, together with the standardized mean difference (SMD); if `|SMD| > 0.1` the ratio is lowered and matching repeats, and a warning is issued if the ratio of 1 is still unbalanced. A larger background ranked the motif of the factor behind a gene set higher in the GTRD benchmark: the mean reciprocal rank of the own-TF motif was 0.28 (ratio 1), 0.34 (ratio 3) and 0.36 (ratio 5) with the `auprc` p-value.

For each motif, the pooled foreground and selected background promoters are ranked by the summarized score. The enrichment statistic is `auprc`, the average precision of that ranking: the mean, over foreground promoters, of the fraction of foreground promoters among all promoters ranked at or above them. Its chance level is `N_foreground / (N_foreground + N_background)`, which is `1/6` at the default ratio of 5; larger values indicate enrichment of high motif scores in the foreground, and the value depends on the ratio, so compare `auprc` only between runs with the same background size.

Significance is a one-sided (enrichment only) permutation p-value, computed for `auprc`. Without tied scores, the foreground labels of a motif read in descending score order are a uniformly random arrangement under the null hypothesis, and that null depends only on the group sizes. One sample of `--n-permutations` such arrangements therefore serves every motif. The upper tail of that sample is not counted but modelled: its mean, standard deviation, skewness and excess kurtosis define a Cornish–Fisher expansion that converts a statistic into a normal quantile, so p-values below `1 / n` are available and the sample can stay small. The background was matched, not drawn at random, so labels are exchangeable only inside a matched group: a foreground promoter and the background promoters assigned to it. For each motif the mean and standard deviation of `auprc` under 2000 permutations that choose one foreground promoter per group give its calibrated null; the observed value is mapped from that scale to the scale of the large sample before the tail count. On random foreground sets of 25–1000 genes the fraction of `p < 0.05` is 4.5–5.2 % (2.8–3.8 % without the calibration). Against 4·10⁷-permutation samples the Cornish–Fisher tail is within about 0.5 orders of magnitude at `p = 1e-6` for group sizes 25 to 1000 and ratios 1 and 3, and against exact matched-group permutations the calibrated p-values agree within about 20 % down to `p ~ 3e-3` and are conservative by at most a factor of about 3 near `1e-5`. These checks cover foreground groups of at least 25 promoters; for smaller groups the statistics are strongly discrete and p-values are approximate. Ties in scores are broken by a random order fixed by `--seed`; a motif with all scores equal gets `p = 1`. The null hypothesis is that scores are unrelated to the labels within matched groups, so p-values are conditional on the selected background and are valid only to the extent that the background matches the foreground in promoter composition.

ESDEG applies the Benjamini–Hochberg procedure across all tested motifs. The enrichment TSV contains every motif; it does not automatically remove non-significant results. The optional plot displays motifs with `p_value_auprc_adj < 0.05`.

## Results

The enrichment TSV is sorted by decreasing `auprc` and contains:

| Column | Meaning |
| --- | --- |
| `motif_id` | Motif model identifier. |
| `tf_name`, `tf_class`, `tf_family` | Motif annotation. |
| `auprc` | auPRC: area under the precision–recall curve of foreground against background, computed as average precision. |
| `p_value_auprc` | One-sided calibrated permutation p-value for `auprc`. |
| `p_value_auprc_adj` | Benjamini–Hochberg adjusted p-value. |
| `jaspar_cluster` | JASPAR motif cluster, included only for bundled JASPAR runs. |

## Expression annotation

Add TF expression information to an enrichment table:

```bash
uv run esdeg annotation enrichment.tsv deg.csv counts.csv genes.gtf annotated.tsv \
  --filter
```

The DEG table needs `id`, `log2FoldChange`, and `padj`; the counts table needs `id` and `counts`; the GTF supplies `gene_id` to `gene_name` mapping. Missing expression values remain `NaN`.

With `--filter`, a row must satisfy `p_value_auprc_adj < --me-padj` (default `0.05`). Its TF must also be differentially expressed according to `--de-padj` and `--log-fc`, or have counts at least `--ncounts`. Use `--best` to retain the best motif per TF after annotation.

## Development

```bash
uv run pytest -q
uv run ruff check .
```

## References

Oshchepkov D, Chadaeva I, Kozhemyakina R, et al. Transcription Factors as Important Regulators of Changes in Behavior through Domestication of Gray Rats: Quantitative Data from RNA Sequencing. *International Journal of Molecular Sciences*. 2022;23(20):12269. [https://doi.org/10.3390/ijms232012269](https://doi.org/10.3390/ijms232012269)

## License

MIT. See [LICENSE](LICENSE).
