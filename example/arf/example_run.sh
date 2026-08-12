#!/bin/bash

promoters=/home/anton/Documents/genomes/promoters/at/promoters.p1500m0.fa
table=/home/anton/Documents/PhD/ESDEG/ARF/auxin_meristem_simonini_05h_DEGs.csv
#time uv run esdeg enrichment \
#$promoters \
#/home/anton/Documents/PhD/ESDEG/ARF/auxin_apex_mogzova_55h_DEGs.csv \
#./auxin_apex_mogzova_55h.all.tsv \
#-f deg \
#-m jaspar \
#-t plants \
#--xlsx ./auxin_apex_mogzova_55h.all.report.xlsx \
#--regulated all \
#--nproc 6


time uv run esdeg enrichment \
$promoters \
$table \
./auxin_meristem_simonini_05h_DEGs.up.tsv \
-f deg \
-m jaspar \
-t plants \
--xlsx ./auxin_meristem_simonini_05h_DEGs.up.report.xlsx \
--regulated up \
--nproc 6



time uv run esdeg enrichment \
$promoters \
$table \
./auxin_meristem_simonini_05h_DEGs.down.tsv \
-f deg \
-m jaspar \
-t plants \
--xlsx ./auxin_meristem_simonini_05h_DEGs.down.report.xlsx \
--regulated down \
--nproc 6


