#!/bin/bash

promoters=./promoters.p400m100.fa
gtf=./Homo_sapiens.GRCh38.112.gtf.gz
if [[ ! -e $promoters ]]; then
    unzip ./promoters.p400m100.zip
fi

# if [[ ! -e $gtf ]]; then
#     wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
# fi


time uv run esdeg enrichment \
$promoters \
./E-MTAB-6598.degs.csv \
./stat_irf.montecarlo.enrichment.tsv \
-f deg \
-m hocomoco \
--visualization ./stat_irf.montecarlo.all.html \
--report ./stat_irf.montecarlo.all.report.html \
--xlsx ./stat_irf.montecarlo.all.report.xlsx \
--regulated all \
--nproc 8


time uv run esdeg enrichment \
$promoters \
./E-MTAB-6598.degs.csv \
./stat_irf.montecarlo.up.tsv \
-f deg \
-m hocomoco \
--visualization ./stat_irf.montecarlo.up.report.html \
--report ./stat_irf.montecarlo.up.report.html \
--xlsx ./stat_irf.montecarlo.up.report.xlsx \
--regulated up \
--nproc 8


time uv run esdeg enrichment \
$promoters \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./nfkb1.montecarlo.set.tsv \
-f set \
-m hocomoco \
--visualization ./nfkb1.montecarlo.set.report.html \
--report ./nfkb1.montecarlo.set.report.html \
--xlsx ./nfkb1.montecarlo.set.report.xlsx \
--nproc 8

# time ESDEG annotation \
# ./stat_irf.montecarlo.enrichment.tsv \
# ./E-MTAB-6598.degs.csv \
# ./E-MTAB-6598.fpkm.csv \
# $gtf \
# ./stat_irf.montecarlo.enrichment.filtered.tsv \
# --xlsx ./stat_irf.montecarlo.enrichment.filtered.xlsx \
# --filter
