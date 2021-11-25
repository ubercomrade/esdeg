#!/bin/bash

enREST.py \
./E-GEOD-48230-query-results.csv \
./OVOL1_HUMAN.H11MO.0.C.pcm \
hg38 \
./ovol1.montecarlo.enrichment.tsv \
--method montecarlo \
--parameter enrichment \
--format hocomoco


enREST.py \
./E-GEOD-48230-query-results.csv \
./OVOL1_HUMAN.H11MO.0.C.pcm \
hg38 \
./ovol1.montecarlo.fraction.tsv \
--method montecarlo \
--parameter fraction \
--format hocomoco


enREST.py \
./E-GEOD-48230-query-results.csv \
./OVOL1_HUMAN.H11MO.0.C.pcm \
hg38 \
./ovol1.binom.enrichment.tsv \
--method binom \
--parameter enrichment \
--format hocomoco


enREST.py \
./E-GEOD-48230-query-results.csv \
./OVOL1_HUMAN.H11MO.0.C.pcm \
hg38 \
./ovol1.hypergeom.fraction.tsv \
--method hypergeom \
--parameter fraction \
--format hocomoco
