#!/bin/bash
# JASPAR 2026 plant models on up- and down-regulated genes of an Arabidopsis DEG table.
# The inputs are not bundled: bash example/arf/example_run.sh <promoters.fa> <deg_table.csv>
set -euo pipefail
cd "$(dirname "$0")"

promoters=${1:?promoter FASTA path}
table=${2:?DEG table path}

time uv run esdeg enrichment \
$promoters \
$table \
./auxin_meristem_simonini_05h_DEGs.up.tsv \
-f deg \
-m jaspar \
-t plants \
--xlsx ./auxin_meristem_simonini_05h_DEGs.up.report.xlsx \
--regulated up


time uv run esdeg enrichment \
$promoters \
$table \
./auxin_meristem_simonini_05h_DEGs.down.tsv \
-f deg \
-m jaspar \
-t plants \
--xlsx ./auxin_meristem_simonini_05h_DEGs.down.report.xlsx \
--regulated down
