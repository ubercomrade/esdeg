#!/bin/bash
time ESDEG.py pwm_preparation \
./several_matrices_for_deg_testing.txt \
hg38 \
./matrices_for_deg_testing \
--format hocomoco
 
time ESDEG.py pwm_preparation \
./several_matrices_for_set_testing.txt \
hg38 \
./matrices_for_set_testing \
--format hocomoco

time ESDEG.py bamm_preparation \
./several_bamms_for_set_testing \
hg38 \
./bamms_for_set_testing

time ESDEG.py bamm_preparation \
./several_bamms_for_deg_testing \
hg38 \
./bamms_for_deg_testing

#PWM
time ESDEG.py deg \
./E-GEOD-48230-query-results.csv \
./matrices_for_deg_testing \
hg38 \
./ovol1.montecarlo.enrichment.tsv \
--parameter enrichment \
--regulated down 

time ESDEG.py deg \
./E-GEOD-48230-query-results.csv \
./matrices_for_deg_testing \
hg38 \
./ovol1.montecarlo.fraction.tsv \
--parameter fraction \
--regulated up


time ESDEG.py set \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./matrices_for_set_testing \
hg38 \
./nfkb1.montecarlo.enrichment.tsv \
--parameter enrichment


time ESDEG.py set \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./matrices_for_set_testing \
hg38 \
./nfkb1.montecarlo.fraction.tsv \
--parameter fraction

#BaMM
time ESDEG.py deg \
./E-GEOD-48230-query-results.csv \
./bamms_for_deg_testing \
hg38 \
./ovol1.bamms.montecarlo.enrichment.tsv \
--parameter enrichment \
--regulated down 

time ESDEG.py deg \
./E-GEOD-48230-query-results.csv \
./bamms_for_deg_testing \
hg38 \
./ovol1.bamms.montecarlo.fraction.tsv \
--parameter fraction \
--regulated up


time ESDEG.py set \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./bamms_for_set_testing \
hg38 \
./nfkb1.bamms.montecarlo.enrichment.tsv \
--parameter enrichment


time ESDEG.py set \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./bamms_for_set_testing \
hg38 \
./nfkb1.bamms.montecarlo.fraction.tsv \
--parameter fraction
 
