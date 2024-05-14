#!/bin/bash

promoters=./promoters.p400m100.fa
gtf=./Homo_sapiens.GRCh38.112.gtf.gz
if [[ ! -e $promoters ]]; then
    unzip ./promoters.p400m100.zip
fi

if [[ ! -e $gtf ]]; then
    wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
fi


# comment if you want use JASPAR
time ESDEG hocomoco \
$promoters \
./matrices_db \
--nproc 6

# uncomment if you want use JASPAR
#time ESDEG jaspar \
#vertebrates \
#$promoters \
#./matrices_db \
#--nproc 4


time ESDEG deg \
./E-MTAB-6598.degs.csv \
./matrices_db \
./stat_irf.montecarlo.enrichment.tsv \
--visualization ./stat_irf.montecarlo.enrichment.html \
--report ./stat_irf.montecarlo.enrichment.report.html \
--xlsx ./stat_irf.montecarlo.enrichment.report.xlsx \
--parameter enrichment \
--regulated all

time ESDEG deg \
./E-MTAB-6598.degs.csv \
./matrices_db \
./stat_irf.montecarlo.fraction.tsv \
--visualization ./stat_irf.montecarlo.fraction.html \
--parameter fraction \
--regulated up


time ESDEG set \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./matrices_db \
./nfkb1.montecarlo.enrichment.tsv \
--visualization ./nfkb1.montecarlo.enrichment.html \
--parameter enrichment


time ESDEG set \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./matrices_db \
./nfkb1.montecarlo.fraction.tsv \
--visualization ./nfkb1.montecarlo.fraction.html \
--parameter fraction

time ESDEG annotation \
./stat_irf.montecarlo.enrichment.tsv \
./E-MTAB-6598.degs.csv \
./E-MTAB-6598.fpkm.csv \
$gtf \
./stat_irf.montecarlo.enrichment.filtered.tsv \
--xlsx ./stat_irf.montecarlo.enrichment.filtered.xlsx \
--filter
