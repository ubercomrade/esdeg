# ESDEG

## Introduction
We developed approach that estimates the enrichment of motifs respecting transcription factor binding sites (TFBS) in promoters of differentially expressed genes (DEGs) derived from RNA-seq experiment. The application of our tool requires designation of species (options hg38, mm10 and tair10 respect human, mouse and Arabidopsis), an input list of gene ID of GEGs is required (Ensembl ID for human/mouse, TAIR ID for Arabidopsis)

## Python requirements

  * numpy
  * scipy
  * pythran
  * pandas
  * statsmodels

  You can use pip to install dependences `pip install numpy, scipy, statsmodels, pythran, pandas`,
  Or you can use conda `conda install numpy, scipy, statsmodels, pythran, pandas`

## Installation

```  
git clone https://github.com/ubercomrade/enrest.git  
cd enrest/  
pip install -e .  
```

## Usage
The command `ESDEG.py -h` return:

```
usage: ESDEG.py [-h] {pwm_preparation,bamm_preparation,deg,set} ...

positional arguments:
  {pwm_preparation,bamm_preparation,deg,set}
                        Available commands:
    pwm_preparation     Run data base preparation (PWM models)
    bamm_preparation    Run data base preparation (BaMM models, experimental!)
    deg                 Run test on DEGs
    set                 Run test on SET of genes

options:
  -h, --help            show this help message and exit

```

Command `enREST.py deg` is used for analisys DEGs (input of DEseq2 is used).

Command `enREST.py set` is used for analisys SET of genes (input of gene list is used).

Command `enREST.py pwm_preparation` is used for preparing motif (PWM) data base for next analisys (input of motifs in HOCOMOCO/MEME format is used).

Command `enREST.py bamm_preparation` is used for preparing motif (BaMM) data base for next analisys.

## Preparation (PWM)

```
usage: ESDEG.py pwm_preparation [-h] [-f FORMAT] matrices N output

positional arguments:
  matrices              Path to matrices in HOCOMOCO (PCM) or in MEME (PFM) format
  N                     promoters of organism (hg38, mm10, tair10, rnor6)
  output                Name of directory for output files

options:
  -h, --help            show this help message and exit
  -f FORMAT, --format FORMAT
                        Format of file with matrices (meme or hocomoco), default= meme
```
#### Required arguments description

**First positional argument** `matrices`:

It's PATH to file with matrices in HOCOMOCO [3] (PCMs) format or MEME [4] (PFMs) format. _It's should be noted, that file can also contain only one model_

Example (HOCOMOCO):

```
>AHR_HUMAN.H11MO.0.B
40.51343240527031 18.259112547756697  56.41253757072521 38.77363485291994
10.877470982533044  11.870876719950774  34.66312982331297 96.54723985087516
21.7165707818416  43.883079837598544  20.706746561638717  67.6523201955933
2.5465132509466635  1.3171620263517245  145.8637051322628 4.231336967110781
0.0 150.35847450464382  1.4927836298652875  2.1074592421627525
3.441039751299748 0.7902972158110341  149.37613720253387  0.3512432070271259
0.0 3.441039751299748 0.7024864140542533  149.81519121131782
0.0 0.0 153.95871737667187  0.0
43.07922333291745 66.87558226865211 16.159862546986584  27.844049228115868
>AIRE_HUMAN.H11MO.0.C
16.428571428571484  10.795918367346953  5.1632653061224625  8.918367346938803
7.51020408163268  7.51020408163268  5.632653061224489 20.65306122448985
5.632653061224489 5.632653061224489 8.448979591836776 21.591836734693906
1.877551020408166 0.0 36.612244897959265  2.8163265306122534
0.0 0.0 33.79591836734698 7.51020408163268
13.14285714285717 1.877551020408166 0.938775510204083 25.346938775510242
15.959183673469415  3.755102040816336 0.938775510204083 20.65306122448985
15.02040816326536 5.632653061224489 5.632653061224489 15.02040816326536
14.081632653061265  2.8163265306122534  5.632653061224489 18.775510204081698
21.591836734693945  3.755102040816336 5.632653061224489 10.326530612244925
15.959183673469415  1.877551020408166 2.8163265306122534  20.65306122448985
9.38775510204083  5.632653061224489 2.8163265306122534  23.469387755102094
0.0 0.938775510204083 36.612244897959265  3.755102040816336
0.0 0.0 40.36734693877556 0.938775510204083
9.38775510204083  3.755102040816336 10.326530612244925  17.836734693877606
2.8163265306122534  7.51020408163268  7.51020408163268  23.469387755102094
18.77551020408166 0.938775510204083 6.571428571428585 15.02040816326536
16.663265306122497  11.030612244898006  4.45918367346938  9.153061224489816
> ... etc.
```

Example (MEME):

`````
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF MA0004.1 Arnt
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
 0.200000  0.800000  0.000000  0.000000
 0.950000  0.000000  0.050000  0.000000
 0.000000  1.000000  0.000000  0.000000
 0.000000  0.000000  1.000000  0.000000
 0.000000  0.000000  0.000000  1.000000
 0.000000  0.000000  1.000000  0.000000
URL http://jaspar2018.genereg.net/matrix/MA0004.1

MOTIF MA0006.1 Ahr::Arnt
letter-probability matrix: alength= 4 w= 6 nsites= 24 E= 0
 0.125000  0.333333  0.083333  0.458333
 0.000000  0.000000  0.958333  0.041667
 0.000000  0.958333  0.000000  0.041667
 0.000000  0.000000  0.958333  0.041667
 0.000000  0.000000  0.000000  1.000000
 0.000000  0.000000  1.000000  0.000000
URL http://jaspar2018.genereg.net/matrix/MA0006.1
`````

**Second positional argument**  `N`:

Options for `N` are _hg38_, _mm10_, _tair10_ or _rnor6_. Depend on organism usage in research

**Third positional argument** `output`:

Directory to write prepared database of motifs.

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-f; --format` :

Options for `-f/--format ` are  _meme_ and _hocomoco_. You should choose value of parameter based on your input format data. The default value is _meme_.

## Preparation (BaMM, experimental!)

```
usage: ESDEG.py bamm_preparation [-h] [-o ORDER] bamms N output

positional arguments:
  bamms                 Path to directory with list of subdirectories contained BaMM models. BaMM files (TAG_motif_1.ihbcp, TAG.hbcp) must have the same TAG as the subdirectory name where files are placed
  N                     organism (hg38, mm10, tair10, rnor6)
  output                Name of directory to write output files

options:
  -h, --help            show this help message and exit
  -o ORDER, --order ORDER
                        Order of BaMMs. Order have to be common for all models. Default is 2.
```
#### Required arguments description

**First positional argument** `bamms`:

It's PATH to DIRECTORY with list of subdirectories contained BaMM models [5]. Each subdirectory contains BaMM files (TAG_motif_1.ihbcp, TAG.hbcp). They must have the same TAG as the subdirectory name where files are placed. Prepared BaMM models can be downloaded from https://bammmotif.soedinglab.org (http://wwwuser.gwdg.de/~compbiol/bamm)

Example:
```
DIRECTORY
    │
    ├── ERR1
    │    ├── ERR1.hbcp
    │    └── ERR1_motif_1.ihbcp
    ├── FOXA1
    │    ├── FOXA1.hbcp
    │    └── FOXA1_motif_1.ihbcp
    ├── FOXA2
    │    ├── FOXA2.hbcp
    │    └── FOXA2_motif_1.ihbcp
    ├── GATA-1
    │    ├── GATA-1.hbcp
    │    └── GATA-1_motif_1.ihbcp
    └── GATA-3
         ├── GATA-3.hbcp
         └── GATA-3_motif_1.ihbcp
```

**Second positional argument**  `N`:

Options for `N` are _hg38_, _mm10_, _tair10_ or _rnor6_. Depend on organism usage in research

**Third positional argument** `output`:

Directory to write prepared database of motifs.

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-o; --order`:

Order of BaMMs. Order have to be common for all models.


## DEG case

``````
usage: ESDEG.py deg [-h] [-p PARAMETER] [-r N] [-P PVALUE] [-l LOG2FC_DEG] [-L LOG2FC_BACK] [-c CONTENT] deg matrices N output

positional arguments:
  deg                   TSV file with DEG with ..., The NAME column must contain ensemble gene IDS
  matrices              Path to prepared data base of matrices
  N                     Organism (hg38, mm10, tair10)
  output                Path to write table with results

options:
  -h, --help            show this help message and exit
  -p PARAMETER, --parameter PARAMETER
                        Parameter estimated in test (enrichment or fraction), default= enrichment
  -r N, --regulated N   The parameter is used to choose up/down/all DEGs, default= all
  -P PVALUE, --pvalue PVALUE
                        The pvalue is used as threshold to choose DEGs, default= 0.05
  -l LOG2FC_DEG, --log2fc_deg LOG2FC_DEG
                        The absolute value of log2FoldChange used as threshold to choose DEGs promoters (DEGs >= thr OR DEGs <= -thr), default= 1
  -L LOG2FC_BACK, --log2fc_back LOG2FC_BACK
                        The absolute value of log2FoldChange used as threshold to choose background promoters (-thr <= BACK <= thr), default= log2(5/4)
  -c CONTENT, --content CONTENT
                        The maximal GC content difference between promoters of foreground and background in Monte Carlo algorithm. Range of possible threshold [0.01 .. 1.0]. If threshold is equal to 1.0 then GC content is not taken into
                        account. In this case (thr = 1.0) algorithm works faster. Default= 0.3.
``````

#### Required arguments description

**First positional argument** `deg` :

It's PATH to your Tab-separated file with full list of DEGs with columns: id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj. Same table can be generated by DESeq2 [1] (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or IRIS [2] (https://bmbls.bmi.osumc.edu/IRIS/). The first column (with name id) must contain ensemble ids for genes.

**Second positional argument** `matrices`:

It's PATH to directory with motifs database prepared by `ESDEG.py preparation`

**Third positional argument**  `N`:

Options for `N` are _hg38_ or _mm10_. Depend on organism usage in research. It has to be the same as for motifs database.

**Fourth positional argument** `output`:

Path to write result table.


#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-p; --parameter` :

The value of `-p; --parameter ` could be  _enrichment_ or _fraction_. If you choose _enrichment_, statistics are calculated based on number of TFBS in DEGs promoters. In this case number of TFBS in each promoters may play role and influences the result. If you choose _fraction_, statistics are calculated based on number of DEGs with TFBS.  In this case number of TFBS in promoters doesn't matter. The default value is _enrichment_.

**Third optional argument** `-r; --regulated` :

The argument `-r/--regulated` are used to choose type of DEGs in analisys. It could be `down` -> promoters of down regulated genes will be used in analisys; `up` -> promoters of up regulated genes will be used in analisys; `all` -> promoters of  up and down regulated genes will be used in analisys. The default value is _all_.

**Fourth optional argument** `-p; --pvalue` :

The argument  `-p; --pvalue ` is pvalue cutoff for DEGs choosing (DEGs <= pvalue). The default value is _0.05_.

**Fifth optional argument** `-l; --log2fc_deg` :

The argument  `-l; --log2fc_deg ` is Log2FoldChange cutoff for DEGs choosing (DEGs >= Log2FoldChange OR DEGs <= -Log2FoldChange). The default value is _1_.

**Sixth optional argument** `-l; --log2fc_back` :

The argument  `-l; --log2fc_back ` is Log2FoldChange cutoff for background choosing (-Log2FoldChange <= BACKGROUND <= Log2FoldChange). The default value is _log2(5/4)_.

**Seventh optional argument** `-c; --content` :

The argument `-c; --content` is used to set threshold of GC content for generating background.


## SET case

````
sage: ESDEG.py set [-h] [-p PARAMETER] [-f FORMAT] [-c CONTENT] set matrices N output

positional arguments:
  set                   File with list of genes. Genes must be in Ensemble format (ensemble gene IDS)
  matrices              Path to prepared data base of matrices
  N                     Organism (hg38, mm10, tair10)
  output                Path to write table with results

options:
  -h, --help            show this help message and exit
  -p PARAMETER, --parameter PARAMETER
                        Parameter estimated in test (enrichment or fraction), default= enrichment
  -c CONTENT, --content CONTENT
                        The maximal GC content difference between promoters of foreground and background in Monte Carlo algorithm. Range of possible threshold [0.01 .. 1.0]. If threshold is equal to 1.0 then GC content is not taken into
                        account. In this case (thr = 1.0) algorithm works faster. Default= 0.3.
````

#### Required arguments description

**First positional argument** `set` :

It's PATH to your SET of genes. Example:

```
ENSG00000072415
ENSG00000094796
ENSG00000143217
ENSG00000038427
ENSG00000186480
...
```

**Second positional argument** `matrices`:

It's PATH to directory with motifs database prepared by `ESDEG.py preparation`

**Third positional argument**  `N`:

Options for `N` are _hg38_ or _mm10_. Depend on organism usage in research. It has to be the same as for motifs database.

**Fourth positional argument** `output`:

Path to write result table.

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-p; --parameter` :

Options for `-p; --parameter ` are  _enrichment_ and _fraction_. If you choose _enrichment_, statistics are calculated based on number of TFBS in DEGs promoters. In this case number of TFBS in each promoters may play role and influences the result. If you choose _fraction_, statistics are calculated based on number of DEGs with TFBS.  In this case number of TFBS in promoters doesn't matter. The default value is _enrichment_.

**Third optional argument** `-c; --content` :

The argument `-c; --content` is used to set threshold of GC content for generating background.\


## Example run

Bash script with examples and data are located in `./example/example_run.sh` . You should run this script in `./example` directory.

```
ESDEG.py pwm_preparation \
./several_matrices_for_deg_testing.txt \
hg38 \
./matrices_for_deg_testing \
--format hocomoco
```

```
ESDEG.py deg \
./E-GEOD-48230-query-results.csv \
./matrices_for_deg_testing \
hg38 \
./ovol1.montecarlo.enrichment.tsv \
--parameter enrichment \
--regulated down
```

## Output file format
ESDEG generates output file in tsv format.
Here is example of file:
```
motif  log(or) distance  pval  adj.pval genes
OVOL1_HUMAN.H11MO.0.C 0.366649526240003 3.78217128582708  0.000278219178978296  0.000834657536934888  ENSG00000000971;ENSG00000008311;ENSG00000009694...
AHR_HUMAN.H11MO.0.B -0.468957662848219  -2.86315741842022 0.00311620230691013 0.00623240461382027 ENSG00000009694;ENSG00000019549;ENSG00000021300...
AIRE_HUMAN.H11MO.0.C  0.166577225854898 1.42032679890637  0.134777945940099 0.161733535128118 ENSG00000008311;ENSG00000009694;ENSG00000019549...
ALX1_HUMAN.H11MO.0.B  0.533935023187338 98.3543133675692  4.53256741484101E-99  2.71954044890461E-98  ENSG00000000971;ENSG00000008311;ENSG00000009694...
PBX1_HUMAN.H11MO.0.A  -0.0796380533521957 -1.1967039140339  0.274421836462296 0.274421836462296 ENSG00000000971;ENSG00000008311;ENSG00000009694...
PBX2_HUMAN.H11MO.0.C  0.243876166666887 1.77057549835486  0.0482711130305087  0.0724066695457631  ENSG00000000971;ENSG00000008311;ENSG00000009694...
```
Where:
**motif** - name of motif from data base. Names depends on preparation step.

**log(or)** - it's -log2 transforamation of  odds ratio (OR). It could be defined as $OR = N_f / N_b$, where $N_f$ - is the number of foreground promoters with predicted sites and $N_b$ - is a mean value of number of background promoters with predicted sites estimated by Monte-Carlo approach (fraction approach). Also It could be defined as $OR = N_f / N_b$, where $N_f$ - is a number of predicted sites in foreground and $N_b$ - is the mean value of number of predicted sites in background estimated by Monte-Carlo approach (enrichment approach)

**distance** - it's the euclidean distance of each factor calculated by using p-value and OR [6]. 

**pval** - it's combined p-value culculated by Hartung method [7]. 

**adj.pval** - adjasted p-value by Benjamini-Hochberg FDR correction. `statsmodels` is used to apply this correction. 

**genes** - list of genes with predicted sites (threshold(ERR) = 0.0005). 

## Reference

1. Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, *15*(12), 550. https://doi.org/10.1186/s13059-014-0550-8
2. Monier, B., McDermaid, A., Wang, C., Zhao, J., Miller, A., Fennell, A., & Ma, Q. (2019). IRIS-EDA: An integrated RNA-Seq interpretation system for gene expression data analysis. *PLOS Computational Biology*, *15*(2), e1006792. https://doi.org/10.1371/journal.pcbi.1006792
3. Kulakovskiy, I. v., Vorontsov, I. E., Yevshin, I. S., Sharipov, R. N., Fedorova, A. D., Rumynskiy, E. I., Medvedeva, Y. A., Magana-Mora, A., Bajic, V. B., Papatsenko, D. A., Kolpakov, F. A., & Makeev, V. J. (2018). HOCOMOCO: Towards a complete collection of transcription factor binding models for human and mouse via large-scale ChIP-Seq analysis. *Nucleic Acids Research*, *46*(D1), D252–D259. https://doi.org/10.1093/nar/gkx1106
4. Machanick, P., & Bailey, T. L. (2011). MEME-ChIP: motif analysis of large DNA datasets. *Bioinformatics*, *27*(12), 1696–1697. https://doi.org/10.1093/bioinformatics/btr189
5. Ge, W., Meier, M., Roth, C., & Söding, J. (2021). Bayesian Markov models improve the prediction of binding motifs beyond first order. NAR genomics and bioinformatics, 3(2), lqab026. https://doi.org/10.1093/nargab/lqab026
6. Puente-Santamaria, L., Wasserman, W. W., & Del Peso, L. (2019). TFEA.ChIP: a tool kit for transcription factor binding site enrichment analysis capitalizing on ChIP-seq datasets. Bioinformatics (Oxford, England), 35(24), 5339–5340. https://doi.org/10.1093/bioinformatics/btz573
7. Hartung, J. (1999). A note on combining dependent tests of significance. Biometrical Journal, 41(7), 849-855.

## License

Copyright (c) 2021 Anton Tsukanov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
