#!/bin/bash

enREST.py deg \
./E-GEOD-48230-query-results.csv \
./OVOL1_HUMAN.H11MO.0.C.pcm \
hg38 \
./ovol1.montecarlo.enrichment \
--parameter enrichment \
--format hocomoco


enREST.py deg \
./E-GEOD-48230-query-results.csv \
./OVOL1_HUMAN.H11MO.0.C.pcm \
hg38 \
./ovol1.montecarlo.fraction \
--parameter fraction \
--format hocomoco


enREST.py set \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./NFKB1_HUMAN.H11MO.1.B.pcm \
hg38 \
./nfkb1.montecarlo.enrichment \
--parameter enrichment \
--format hocomoco


enREST.py set \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./NFKB1_HUMAN.H11MO.1.B.pcm \
hg38 \
./nfkb1.montecarlo.fraction \
--parameter fraction \
--format hocomoco