#!/bin/bash

xz -d ./hs.ensembl.promoters.fa.xz
promoters=./hs.ensembl.promoters.fa

time ESDEG.py preparation \
vertebrates \
$promoters \
./matrices_db
 
time ESDEG.py deg \
./E-GEOD-48230-query-results.csv \
./matrices_db \
./ovol1.montecarlo.enrichment.tsv \
--parameter enrichment \
--regulated down 

time ESDEG.py deg \
./E-GEOD-48230-query-results.csv \
./matrices_db \
./ovol1.montecarlo.fraction.tsv \
--parameter fraction \
--regulated up


time ESDEG.py set \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./matrices_db \
./nfkb1.montecarlo.enrichment.tsv \
--parameter enrichment


time ESDEG.py set \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./matrices_db \
./nfkb1.montecarlo.fraction.tsv \
--parameter fraction
