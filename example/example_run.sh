#!/bin/bash

enREST.py deg \
./E-GEOD-48230-query-results.csv \
./OVOL1_HUMAN.H11MO.0.C.pcm \
hg38 \
./ovol1.montecarlo.enrichment.tsv \
--parameter enrichment \
--format hocomoco


enREST.py deg \
./E-GEOD-48230-query-results.csv \
./OVOL1_HUMAN.H11MO.0.C.pcm \
hg38 \
./ovol1.montecarlo.fraction.tsv \
--parameter fraction \
--format hocomoco