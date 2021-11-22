#!/bin/bash

enREST.py matrix \
./E-GEOD-48230-query-results.tsv \
./OVOL1_HUMAN.H11MO.0.C.pcm \
hg38 \
./ovol1.enrichment.tsv \
-m enrichment

enREST.py matrix \
./E-GEOD-48230-query-results.tsv \
./OVOL1_HUMAN.H11MO.0.C.pcm \
hg38 \
./ovol1.fraction.tsv \
-m fraction

enREST.py hocomoco \
./E-GEOD-48230-query-results.tsv \
./HOCOMOCOv11_core_pcms_HUMAN_mono_example.txt \
hg38 \
./hocomoco.enrichment.tsv \
-m enrichment

enREST.py hocomoco \
./E-GEOD-48230-query-results.tsv \
./HOCOMOCOv11_core_pcms_HUMAN_mono_example.txt \
hg38 \
./hocomoco.fraction.tsv \
-m fraction
