#!/bin/bash

tar -xjvf fasta.tar.bz

time enREST.py deg \
./E-GEOD-48230-query-results.csv \
./several_matrices_for_deg_testing.txt \
hg38 \
./ovol1.montecarlo.enrichment \
--parameter enrichment \
--format hocomoco 

time enREST.py deg \
./E-GEOD-48230-query-results.csv \
./several_matrices_for_deg_testing.txt \
hg38 \
./ovol1.montecarlo.fraction \
--parameter fraction \
--format hocomoco


time enREST.py set \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./several_matrices_for_set_testing.txt \
hg38 \
./nfkb1.montecarlo.enrichment \
--parameter enrichment \
--format hocomoco 


time enREST.py set \
./HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt \
./several_matrices_for_set_testing.txt \
hg38 \
./nfkb1.montecarlo.fraction \
--parameter fraction \
--format hocomoco 


time enREST.py fasta \
./foreground.seq \
./background.seq \
./tair10.meme \
tair10 \
./tair10.montecarlo.enrichment \
--parameter enrichment \
--format meme 


time enREST.py fasta \
./foreground.seq \
./background.seq \
./tair10.meme \
tair10 \
./tair10.montecarlo.fraction \
--parameter fraction \
--format meme 
