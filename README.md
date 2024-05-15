# ESDEG

## Introduction

We developed tool ESDEG that estimates the enrichment of motifs respecting transcription factor binding sites (TFBS) in promoters of differentially expressed genes (DEGs) derived from RNA-seq experiment. We applied this tool to perform analysis on promoters of the genes upregulated and downregulated in the brain tissue samples of tame rats [^a].
To identify enriched motifs, we first performed motif recognition in the given set of promoter (upstream regions of the same length), utilizing each from a library of nucleotide frequency matrices taken from the [JASPAR](https://jaspar.elixir.no/). Next, a table of recognition thresholds was compiled as described in works[^b][^c], where expected recognition rates (_ERRs_) for each threshold were calculated as probabilities of site prediction in the whole-genome dataset of promoters. Given that the enrichment of sites depends on a chosen threshold, subsequent calculations were performed for 20 recognition thresholds (equidistant on a logarithmic scale of _ERRs_) in the range of 5E−5 to 1E−3. Next, all promoters were categorized into either a foreground or background group. The foreground group contained the promoters of the significantly DEGs: either only upregulated or only downregulated or both upregulated and downregulated genes. The parameters (_adjusted p-value_ and _log2(fold change)_) were used to determine whether a gene was differentially expressed. DEGs are defined as follows: adjusted _p-value_ < 0.05 and |_log2(fold change)_| > 1. The background group included the same number of promoters as the foreground group, but they were randomly selected from none-DEGs according to given parameters (_adjusted p-value_, _log2(fold change)_ and _G/C content_) and have the same _G/C content_ as foreground. non-DEGs are defined as follows: _adjusted p-value_ > 0.05 and |_log2(fold change)_| <  _log2(5/4)_. Next, we calculated the frequency of a motif in each of the gene promoter groups. To estimate the significance of the motif enrichment, the Monte Carlo approach was applied. For a given motif and a recognition threshold, site frequency was calculated for the foreground group (_AVFOR_), whereas a background group was generated many times to estimate the site frequency distribution. Next, the average (_AVBACK_) and standard deviation (_SDBACK_) for the frequencies in the background group set were calculated. Then, a _Z-score_ was calculated according to the formula _(AVFOR − AVBACK)/SDBACK_. This _Z-score_ allowed us to compute the significance of enrichment (_p-value_) by fitting a site’s frequency distribution to a normal distribution. To confirm enrichment for a given motif, multiple _p-values_ were combined, referring all recognition thresholds into a single unified _p-value_ according to the Bonferroni method [^5]. Thus, a unified _p-value_ was calculated for each motif. Finally, an adjusted p-value (_padj_) was calculated from the set of unified p-values for all motifs with the Benjamini–Hochberg correction for multiple comparisons; in this way, our analysis involved multiple simultaneous statistical tests for various motifs.

* You can find more details in article _Oshchepkov, Dmitry, Irina Chadaeva, Rimma Kozhemyakina, Svetlana Shikhevich, Ekaterina Sharypova, Ludmila Savinkova, Natalya V. Klimova, Anton Tsukanov, Victor G. Levitsky, and Arcady L. Markel. 2022. "Transcription Factors as Important Regulators of Changes in Behavior through Domestication of Gray Rats: Quantitative Data from RNA Sequencing" International Journal of Molecular Sciences 23, no. 20: 12269. https://doi.org/10.3390/ijms232012269_

## Requirements

ESDEG is a command-line toolkit and Python package for estimating the enrichment of motifs in promoters. ESDEG runs on [Python](https://www.python.org/) 3.7 and later. Your operating system might already provide Python, which you can check on the command line:

```
python --version
```

### Python dependencies

If you haven't already satisfied these dependencies on your system, install
these Python packages via ``pip``:

* numpy
* scipy
* pandas
* biopython
* pyjaspar
* xlsxwriter
* panel

```
pip install numpy scipy plotly pandas biopython pyjaspar xlsxwriter panel
```

## Installation

```
git clone https://github.com/ubercomrade/esdeg.git  
cd esdeg  
pip install -e .  
```

## Usage

The command `ESDEG -h` return:

```
usage: ESDEG [-h] {jaspar,hocomoco,deg,set,annotation} ...

positional arguments:
  {jaspar,hocomoco,deg,set,annotation}
                        Available commands:
    jaspar              Run preparation DB based on JASPAR 2024 motifs (motifs
                        are selected based on taxon)
    hocomoco            Run preparation DB based on HOCOMOCO 2024 motifs (all
                        core motifs are used [M. musculus and H. sapiens])
    deg                 Run test on DEGs
    set                 Run test on SET of genes
    annotation          Add expression level and different expression data for
                        each TF/motif to filter ESDEG results

options:
  -h, --help            show this help message and exit
```

ESDEG includes five commands: `jaspar`, `hocomoco`, `deg`, `set` and `annotaion`.

The first step is to prepare the database using the commands `jaspar` or `hocomoco` depending on motifs database. This step requires (a) the library of motifs, which are extracted from the [JASPAR](https://jaspar.elixir.no) or [HOCOMOCO](https://hocomoco12.autosome.org/), and (b) promoters of the same length in FASTA format [example](https://github.com/ubercomrade/esdeg/blob/main/example/hs.ensembl.promoters.fa.zip), which are set by the user. The results of the command `preparation` represent a database including the list of motifs and their annotations in `.npy` format. This database is used further as an input for `deg` and `set` commands.

The next step is depending on type of input data.

If you have the results of RNA-seq analysis of different gene expression ([example](https://github.com/ubercomrade/esdeg/blob/main/example/E-MTAB-6598.degs.csv)) you should use command `deg`.

If you have only set of genes without any additional information as log2(fold change) and p-value adjusted ([example](https://github.com/ubercomrade/esdeg/blob/main/example/HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt)) you should use command `set`.

In the final step, ESDEG results can be complemented with gene expression information to improve filtering and refine reliable TFs. Expression level values (FPKM) as well as expression level changes (log2(Fold Change)) and its significance (padj) are given for all TFs that are included in the list of motifs. For this step you have to have results of RNA-seq analysis of different gene expression ([example](https://github.com/ubercomrade/esdeg/blob/main/example/E-MTAB-6598.degs.csv)) and normalized counts (FPKM or other) ([example](https://github.com/ubercomrade/esdeg/blob/main/example/E-MTAB-6598.fpkm.csv)). The `annotation` command is used for this step.

## Preparation (jaspar)

This step requires (a) the library of motifs, which are extracted from the [JASPAR](https://jaspar.elixir.no), and (b) promoters of the same length in FASTA format [example](https://github.com/ubercomrade/esdeg/blob/main/example/hs.ensembl.promoters.fa.zip), which are set by the user. The results of the command `preparation` represent a database including the list of motifs and their annotations in `.npy` format. This database is used further as an input for `deg` and `set` commands.

_Timing 40–120 min_

```
usage: ESDEG preparation [-h]
                            {plants,vertebrates,insects,urochordates,nematodes,fungi}
                            N output

positional arguments:
  {plants,vertebrates,insects,urochordates,nematodes,fungi}
                        Prepare database for respective JASPAR CORE taxonomic
                        group of motifs. Possible options are plants,
                        vertebrates, insects, urochordates, nematodes, fungi.
                        For more details see https://jaspar.uio.no/ and
                        https://pyjaspar.readthedocs.io/en/latest/index.html
  promoters             Path to promoters in fasta format. All promoters have
                        to be with same length. After the symbol ">" unique gene ID has
                        to be written (>ENSG00000160072::1:1469886-1472284 or
                        >ENSG00000160072)
  output                Name of directory to write output files

options:
  -h, --help            show this help message and exit
  -p NPROC, --nproc NPROC
                        Number of processes to split the work between.
                        Default= 4
```

#### Required arguments description

**First positional argument** `matrices`:

It's name of taxon that is avaliable in JASPAR database[^1]. Possible options for taxon are _plants_, _vertebrates_, _insects_, _urochordates_, _nematodes_, _fungi_.
Only motifs from CORE COLLECTION and with length >= 8 are used for analysis.
For more details see https://jaspar.elixir.no and https://pyjaspar.readthedocs.io/en/latest/index.html

**Second positional argument** `promoters`:

Argument `promoters` is the path to file with promoters in FASTA format. Promoters should have same length. After ">" unique gene ID have to be written (>ENSG00000160072::1:1469886-1472284 or >ENSG00000160072)

[Example](https://github.com/ubercomrade/esdeg/blob/main/example/hs.ensembl.promoters.fa.zip):

```
>ENSG00000160072::1:1469886-1472284
ACATTCCACCATTGTGATTTGTTTCTGCCCCACCCTAG...
>ENSG00000142611::1:3067261-3069659
TCGATAGACCCTCGAAAGGACGGCAGGGAATGGGGCTG...
>ENSG00000157911::1:2412102-2414502
GGACCTGCCCTGAGCTGGGGACGGGAAGGGCTTGGGCG...
>ENSG00000142655::1:10472968-10475366
ATCGGAGAGGAAAGGGGAAATGCGACTCGGCACTGTTG...
```

This type of FASTA file can be generated by using `bedtools getfasta` (https://bedtools.readthedocs.io/en/latest/) [^2] with flag `-name+`

**Third positional argument** `output`:

Directory to write prepared database of motifs.

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-p; --nproc` :

The number of processes (threads) that will be used to prepare the database of motifs. By default, the number of processes is four. Increasing the number of processes allows you to speed up calculations; as the number of processes increases, the amount of RAM used increases.

## Preparation (hocomoco)

This step requires (a) the library of motifs, which are extracted from the [HOCOMOCO](https://hocomoco12.autosome.org/), and (b) promoters of the same length in FASTA format [example](https://github.com/ubercomrade/esdeg/blob/main/example/hs.ensembl.promoters.fa.zip), which are set by the user. The results of the command `preparation` represent a database including the list of motifs and their annotations in `.npy` format. This database is used further as an input for `deg` and `set` commands.

_Timing 40–120 min_

```
usage: ESDEG preparation [-h]
                            N output

positional arguments:
  promoters             Path to promoters in fasta format. All promoters have
                        to be with same length. After the symbol ">" unique gene ID has
                        to be written (>ENSG00000160072::1:1469886-1472284 or
                        >ENSG00000160072)
  output                Name of directory to write output files

options:
  -h, --help            show this help message and exit
  -p NPROC, --nproc NPROC
                        Number of processes to split the work between.
                        Default= 4
```

#### Required arguments description

**First positional argument** `promoters`:

Argument `promoters` is the path to file with promoters in FASTA format. Promoters should have same length. After ">" unique gene ID have to be written (>ENSG00000160072::1:1469886-1472284 or >ENSG00000160072)

[Example](https://github.com/ubercomrade/esdeg/blob/main/example/hs.ensembl.promoters.fa.zip):

```
>ENSG00000160072::1:1469886-1472284
ACATTCCACCATTGTGATTTGTTTCTGCCCCACCCTAG...
>ENSG00000142611::1:3067261-3069659
TCGATAGACCCTCGAAAGGACGGCAGGGAATGGGGCTG...
>ENSG00000157911::1:2412102-2414502
GGACCTGCCCTGAGCTGGGGACGGGAAGGGCTTGGGCG...
>ENSG00000142655::1:10472968-10475366
ATCGGAGAGGAAAGGGGAAATGCGACTCGGCACTGTTG...
```

This type of FASTA file can be generated by using `bedtools getfasta` (https://bedtools.readthedocs.io/en/latest/) [^2] with flag `-name+`

**Second positional argument** `output`:

Directory to write prepared database of motifs.

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-p; --nproc` :

The number of processes (threads) that will be used to prepare the database of motifs. By default, the number of processes is four. Increasing the number of processes allows you to speed up calculations; as the number of processes increases, the amount of RAM used increases.

## DEG case

This command is used when you have comma-separated table including results of differentially expressed genes ([example](https://github.com/ubercomrade/esdeg/blob/main/example/E-GEOD-48230-query-results.csv)). This command allows you to detect which motifs are enriched in differentially expressed genes (foreground) compared to a background. Differentially expressed genes are determined by the parameters: `--pvalue` and `--log2fc_deg`. Background is determined by the parameters: `--pvalue` and `--log2fc_back`. See below.

_Timing from 10 seconds to 5 minutes_

```
usage: ESDEG deg [-h] [-v VISUALIZATION] [-x XLSX] [-R REPORT] [-p PARAMETER] [-r N] [-P PVALUE] [-l LOG2FC_DEG] [-L LOG2FC_BACK]
                    [-c CONTENT]
                    deg matrices output

positional arguments:
  deg                   Input file in CSV format with results of RNA-seq analysis. File must contain next columns: id,
                        log2FoldChange, padj
  matrices              Directory with prepared database contained .npy files (e.g. /path/to/database)
  output                Output file in TSV format (e.g. /path/to/output/file.tsv)

options:
  -h, --help            show this help message and exit
  -v VISUALIZATION, --visualization VISUALIZATION
                        Path to write interactive picture in HTML format (path/to/pic.html). if '-v' is given, then ESDEG creates
                        picutre. By default it isn't used
  -x XLSX, --xlsx XLSX  Path to write table with results in XLSX format (path/to/table.xlsx). XLSX table contains logo of motifs. if
                        '-x' is given, then ESDEG creates XSLX table. By default it isn't used
  -R REPORT, --report REPORT
                        Path to write interactive table with results in HTML format (path/to/report.html). HTML report contains logo
                        of motifs. if '-r' is given, then ESDEG creates HTML report. By default it isn't used
  -p PARAMETER, --parameter PARAMETER
                        Parameter estimated in test (enrichment or fraction), default= enrichment
  -r N, --regulated N   The parameter is used to choose up/down/all DEGs, default= all
  -P PVALUE, --pvalue PVALUE
                        The pvalue is used as threshold to choose DEGs, default= 0.05
  -l LOG2FC_DEG, --log2fc_deg LOG2FC_DEG
                        The absolute value of log2FoldChange used as threshold (L2FC_THR) to choose DEGs promoters (actual L2FC >=
                        L2FC_THR OR actual L2FC <= -L2FC_THR), default= 1.0
  -L LOG2FC_BACK, --log2fc_back LOG2FC_BACK
                        The absolute value of log2FoldChange used as threshold (L2FC_BACK_THR) to choose background promoters
                        (-L2FC_BACK_THR <= actual L2FC <= L2FC_BACK_THR), default= log2(5/4)=0.321928...
  -c CONTENT, --content CONTENT
                        The maximal GC content difference between promoters of foreground and background in Monte Carlo algorithm.
                        Range of possible threshold [0.01 .. 1.0]. If threshold is equal to 1.0 then GC content is not taken into
                        account. In this case (thr = 1.0) algorithm works faster. Default= 0.3.
```

#### Required arguments description

**First positional argument** `deg` :

It is PATH to your comma-separated file (.csv) with full list of DEGs with required columns: id, log2FoldChange, padj ([example](https://github.com/ubercomrade/esdeg/blob/main/example/E-MTAB-6598.degs.csv)). Similar table can be generated by DESeq2[^3] (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or IRIS[^4] (https://bmbls.bmi.osumc.edu/IRIS/).

[Example](https://github.com/ubercomrade/esdeg/blob/main/example/E-MTAB-6598.degs.csv):

```
id,log2FoldChange,padj
ENSG00000000003,-0.2,0.472804281
ENSG00000000419,0.2,0.46039097
ENSG00000000457,-0.1,0.841129918
ENSG00000000460,0.2,0.534296601
ENSG00000000971,-1.4,1.06E-05
ENSG00000001036,0.1,0.81968747
ENSG00000001167,-0.2,0.444068007
ENSG00000001460,-0.1,0.808782152
```

![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) **IMPORTANT! You have to use the same type of gene IDs (ENSEMBL ID, NCBI Gene ID, HUGO Gene ID ...) as that used in the preparation step for promoters**

**Second positional argument** `matrices`:

It is PATH to directory with prepared database by `ESDEG.py preparation`

**Third positional argument** `output`:

Output file in tab-separated format (.tsv) e.g. /path/to/output/file.tsv.

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-v; --visualization` :

It's path to write HTML report contained dotplot (log2(OR) vs -log10(padj)) (e.g. /path/to/file/with/graphics.html). Plotly is used for visualization.

**Third optional argument** `-x; --xlsx` :

It's path to write results in XLSX format. For this format logo is ploted for each motif in last column.

**Fourth optional argument** `-R; --report` :

It's path to write results in HTML format. For this format logo is ploted for each motif in last column.

**Fifth optional argument** `-p; --parameter` :

The value of `-p; --parameter ` can be  _enrichment_ or _fraction_. If you choose _enrichment_ option, statistics is calculated based on the number of predicted sites in promoters of DEGs. In this case the number of predicted sites in each promoter influence the result. If you choose _fraction_ option, statistics is calculated based on the number of DEGs with predicted sites. In this case the number of predicted sites in each promoter is not so important. The default value is _enrichment_.

**Sixth optional argument** `-r; --regulated` :

The argument `-r/--regulated` are used to choose DEGs: `down` (promoters of down regulated genes), `up` (promoters of up regulated genes) and `all` (promoters of up and down regulated genes). The default value is _all_.

**Seventh optional argument** `-p; --pvalue` :

The argument  `-p; --pvalue ` is p-value threshold (_P_THR_) for DEGs selection (actual p-value of gene <= _P_THR_). The default value is _0.05_.

**Eighth optional argument** `-l; --log2fc_deg` :

The argument  `-l; --log2fc_deg ` is Log2FoldChange threshold (_L2FC_THR_) for DEGs selection (actual value of Log2FoldChange of gene >= _L2FC_THR_ OR actual Log2FoldChange of gene <= -_L2FC_THR_). The default value is _1_.

**Ninth optional argument** `-l; --log2fc_back` :

The argument  `-l; --log2fc_back ` is Log2FoldChange threshold for background (_L2FC_BACK_THR_) selection (-_L2FC_BACK_THR_ <= actual value of Log2FoldChange of gene <= _L2FC_BACK_THR_). The default value is _log2(5/4) = 0.376287495_.

**Tenth optional argument** `-c; --content` :

The argument `-c; --content` is used to set threshold of GC content for generating background.

## SET case

This command is used when you have only set of genes from arbitrary source without any information related of expression of genes ([example](https://github.com/ubercomrade/esdeg/blob/main/example/E-GEOD-48230-query-results.csv)). This command allows you to detect which motifs are enriched in a given set of genes (foreground) compared to a background. In this case background includes all genes from database except given as a foreground.

_Timing from 10 seconds to 5 minutes_

```
usage: ESDEG set [-h] [-v VISUALIZATION] [-x XLSX] [-R REPORT] [-p PARAMETER] [-c CONTENT] set matrices output

positional arguments:
  set                   File with list of genes.
  matrices              Path to prepared data base of matrices
  output                Path to write table with results

options:
  -h, --help            show this help message and exit
  -v VISUALIZATION, --visualization VISUALIZATION
                        Path to write interactive picture in HTML format (path/to/pic.html). if '-v' is given, then ESDEG creates
                        picutre. By default it isn't used
  -x XLSX, --xlsx XLSX  Path to write table with results in XLSX format (path/to/table.xlsx). XLSX table contains logo of motifs. if
                        '-x' is given, then ESDEG creates XSLX table. By default it isn't used
  -R REPORT, --report REPORT
                        Path to write interactive table with results in HTML format (path/to/report.html). HTML report contains logo
                        of motifs. if '-r' is given, then ESDEG creates HTML report. By default it isn't used
  -p PARAMETER, --parameter PARAMETER
                        Parameter estimated in test (enrichment or fraction), default= enrichment
  -c CONTENT, --content CONTENT
                        The maximal GC content difference between promoters of foreground and background in Monte Carlo algorithm.
                        Range of possible threshold [0.01 .. 1.0]. If threshold is equal to 1.0 then GC content is not taken into
                        account. In this case (thr = 1.0) algorithm works faster. Default= 0.3.
```

#### Required arguments description

**First positional argument** `set` :

It's PATH to file with your SET of genes. [Example](https://github.com/ubercomrade/esdeg/blob/main/example/HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt):

```
ENSG00000072415
ENSG00000094796
ENSG00000143217
ENSG00000038427
ENSG00000186480
...
```

![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) **IMPORTANT! You have to use the same type of gene IDs (ENSEMBL ID, NCBI Gene ID, HUGO Gene ID ...) as that used in the preparation step for promoters**

**Second positional argument** `matrices`:

It is PATH to directory with prepared database by `ESDEG.py preparation`

**Third positional argument** `output`:

Output file in tab-separated format (.tsv) e.g. /path/to/output/file.tsv.

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-v; --visualization` :

HTML output report contained graphical files (e.g. /path/to/file/with/graphics.html). Plotly is used for visualization.

**Third optional argument** `-x; --xlsx` :

It's path to write results in XLSX format. For this format logo is ploted for each motif in last column.

**Fourth optional argument** `-R; --report` :

It's path to write results in HTML format. For this format logo is ploted for each motif in last column.

**Fifth optional argument** `-p; --parameter` :

The value of `-p; --parameter ` can be  _enrichment_ or _fraction_. If you choose _enrichment_ option, statistics is calculated based on the number of predicted sites in promoters of DEGs. In this case the number of predicted sites in each promoter influence the result. If you choose _fraction_ option, statistics is calculated based on the number of DEGs with predicted sites. In this case the number of predicted sites in each promoter is not so important. The default value is _enrichment_.

**Sixth optional argument** `-c; --content` :

The argument `-c; --content` is used to set threshold of G/C content for generating background.\

## annotation command

ESDEG results can be complemented with gene expression information to improve filtering and refine reliable TFs. Expression level values (FPKM) as well as expression level changes (log2(Fold Change)) and its significance (padj) are given for all TFs that are included in the list of motifs. For this step you have to have results of RNA-seq analysis of different gene expression ([example](https://github.com/ubercomrade/esdeg/blob/main/example/E-MTAB-6598.degs.csv)) and normalized counts (FPKM or other) ([example](https://github.com/ubercomrade/esdeg/blob/main/example/E-MTAB-6598.fpkm.csv)).


_Timing from 1 minutes to 5 minutes_

```
usage: ESDEG annotation [-h] [-x XLSX] [-b] [-f] [-m ME_PADJ] [-l LOG_ODDS]
                        [-d DE_PADJ] [-L LOG_FC] [-n NCOUNTS]
                        esdeg deg counts gtf output

positional arguments:
  esdeg                 Path to esdeg result in TSV format
  deg                   Input file in CSV format with results of RNA-seq
                        analysis. File must contain next columns: id,
                        log2FoldChange, padj
  counts                Input file in CSV format with normalized counts of
                        read (FPKM or other). File must contain next columns:
                        id, counts
  gtf                   Input file in GTF format with genome features
                        annotation (ENSEMBL is prefered). It`s used to convert
                        gene names and IDs
  output                Path to write table with annotated ESDEG results.
                        Table includes information related to expression level
                        and different expression level of each TF/motif

options:
  -h, --help            show this help message and exit
  -x XLSX, --xlsx XLSX  Path to write table with results in XLSX format
                        (path/to/table.xlsx). XLSX table contains logo of
                        motifs. if '-x' is given, then ESDEG creates XSLX
                        table. By default it isn't used
  -b, --best            If this argument is used, then for each TF only the
                        motif with best enrichment will be left (Some times TF
                        have several motifs in DB)
  -f, --filter          If this argument is used, then filtration will be
                        applied to ESDEG table based on expression level of
                        TFs.
  -m ME_PADJ, --me_padj ME_PADJ
                        p-value threshold for motif enrichment (adj.pval
                        column in table. adj.pval will be renamed to me_padj).
                        Default value = 0.05. Applied only when flag --filter
                        is used.
  -l LOG_ODDS, --log_odds LOG_ODDS
                        log(odds ratio) threshold for motif enrichment
                        (log2(or) column in table). Default value = 1.0.
                        Applied only when flag --filter is used.
  -d DE_PADJ, --de_padj DE_PADJ
                        p-value threshold for DEGs. It`s used to found TF
                        among DEGs. Default value = 0.05. Applied only when
                        flag --filter is used.
  -L LOG_FC, --log_fc LOG_FC
                        log(fold change) threshold for DEGs. It`s used to
                        found TF among DEGs. Default value = 1.0. Applied only
                        when flag --filter is used.
  -n NCOUNTS, --ncounts NCOUNTS
                        Counts threshold (expression level threshold) for
                        genes. It`s used to remove TFs that are not expressed
                        (With the exception of TFs, which are DEGs). Default
                        value = 5.0. Applied only when flag --filter is used.
```

#### Required arguments description

**First positional argument** `esdeg` :

It's PATH to esdeg result in TSV format

**Second positional argument** `deg` :

It is PATH to your comma-separated file (.csv) with full list of DEGs with required columns: id, log2FoldChange, padj ([example](https://github.com/ubercomrade/esdeg/blob/main/example/E-MTAB-6598.degs.csv)). Similar table can be generated by DESeq2[^3] (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or IRIS[^4] (https://bmbls.bmi.osumc.edu/IRIS/).

[Example](https://github.com/ubercomrade/esdeg/blob/main/example/E-MTAB-6598.degs.csv):

```
id,log2FoldChange,padj
ENSG00000000003,-0.2,0.472804281
ENSG00000000419,0.2,0.46039097
ENSG00000000457,-0.1,0.841129918
ENSG00000000460,0.2,0.534296601
ENSG00000000971,-1.4,1.06E-05
ENSG00000001036,0.1,0.81968747
ENSG00000001167,-0.2,0.444068007
ENSG00000001460,-0.1,0.808782152
```

**Third positional argument** `counts`:

It is PATH to CSV file that contained normalized counts (FPKM or other) for each gene with required columns: id, counts.

[Example](https://github.com/ubercomrade/esdeg/blob/main/example/E-MTAB-6598.fpkm.csv):

```
id,counts
ENSG00000000003,0.0206177149643403
ENSG00000000005,0.0567211708870934
ENSG00000000419,7.76915277767639
ENSG00000000457,1.79968444256887
ENSG00000000460,6.17131000049197
ENSG00000000938,0.0310136079961313
ENSG00000000971,0.524588698013393
ENSG00000001036,59.3947651344486
ENSG00000001084,2.42592773904779
```

**Fourth positional argument** `gtf`:

It is PATH to GTF format with genome features annotation (ENSEMBL is prefered). It is used to convert gene names and IDs

**Fifth positional argument** `output`:

Output file in tab-separated format (.tsv) e.g. /path/to/output/file.tsv.

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-x; --xlsx` :

It's path to write results in XLSX format. For this format logo is ploted for each motif in last column.

**Third optional argument** `-b; --best` :

 If this argument is used, then filtration will be applied to ESDEG table based on enrichment of motifs and expression level of TFs. For filtration next parameters are used:

**Fourth optional argument** `-f; --filter` :

 If this argument is used, then filtration will be applied to ESDEG table based on enrichment of motifs and expression level of TFs. For filtration next parameters are used:

1. `--me_padj` - p-value threshold for motif enrichment
2. `--log_odds` - log(odds ratio) threshold for motif enrichment
3. `--de_padj` - p-value threshold for DEGs
4. `--log_fc` - log(fold change) threshold for DEGs
5. `--ncounts` - counts threshold (expression level threshold) for
   genes.

**Fifth optional argument** `-m; --me_padj` :

The p-value threshold for motif enrichment (adj.pval column in table. adj.pval will be renamed to me_padj). Default value = 0.05. Applied only when flag `--filter` is used.

**Sixth optional argument** `-l; --log_odds` :

The log(odds ratio) threshold for motif enrichment (log2(or) column in table). Default value = 1.0. Applied only when flag `--filter` is used.

**Seventh optional argument** `-d; --de_padj` :

The p-value threshold for DEGs. It is used to found TF among DEGs. Default value = 0.05. Applied only when flag `--filter` is used.

**Eighth optional argument** `-L; --log_fc` :

The log(fold change) threshold for DEGs. It is used to found TF among DEGs. Default value = 1.0. Applied only when flag `--filter` is used.

**Ninth optional argument** `-L; --log_fc` :

The counts threshold (expression level threshold) for genes. It is used to remove TFs that are not expressed (With the exception of TFs, which are DEGs). Default value = 5.0. Applied only when flag `--filter` is used.

## Example run

Example bash script (for Linux), example batch script (for Windows) and data are located in `./example`. You should run this scripts in `./example` directory.

```
ESDEG jaspar \
vertebrates \
promoters.p400m100.fa \
./matrices_db \
--nproc 4
```

```
ESDEG deg \
./E-MTAB-6598.degs.csv \
./matrices_db \
./stat_irf.montecarlo.enrichment.tsv \
--visualization ./stat_irf.montecarlo.enrichment.html \
--report ./stat_irf.montecarlo.enrichment.report.html \
--xlsx ./stat_irf.montecarlo.enrichment.report.xlsx \
--parameter enrichment \
--regulated all
```

## Output file format

### Example output file of `ESDEG deg` and `ESDEG set` commands

| motif_id | tf_name | tf_class           | tf_family                | jaspar_cluster | log2(or)          | log10(pval)       | log10(adj.pval)   | adj.pval             | genes                                               |
| -------- | ------- | ------------------ | ------------------------ | -------------- | ----------------- | ----------------- | ----------------- | -------------------- | --------------------------------------------------- |
| MA0680.2 | Pax7    | Paired box factors | Paired plus homeo domain | cluster_41     | 0.354260746038436 | -56.5585489729083 | -53.6399944423581 | 2.29089696905926E-54 | ENSG00000000971;ENSG00000008311;ENSG00000009694;... |

Where:
**motif_id** - jaspar id from data base.

**tf_name** - name of Transcription factor

**tf_class** - class of DBD's transcrition factor

**tf_family** - family of DBD's transcrition factor

**jaspar_cluster** - сluster to which the motif belongs (https://jaspar.genereg.net/matrix-clusters/) [writen only for JASPAR database]

**log(or)** - it's -log2 transformation of  odds ratio (OR). For the option `fraction` it can be defined as follow $OR = N_f / N_b$, where $N_f$ - is the number of foreground promoters with predicted sites and $N_b$ - is a mean value of number of background promoters with predicted sites estimated by Monte-Carlo approach. Alternatively for the option `enrichment` it can be defined as follow $OR = N_f / N_b$, where $N_f$ - is a number of predicted sites in foreground and $N_b$ - is the mean value of number of predicted sites in background estimated by Monte-Carlo approach. The best log(or) value is represented in table. It's chosen by minimal p-value.

**log10(pval)** - it's log10(combined p-value). Combined p-value is calculated by Hartung method[^5].

**log10(adj.pval)** - it's log10(adjusted p-value) obtained by using Benjamini-Hochberg FDR correction.

**adj.pval** - it's adjusted p-value obtained by using Benjamini-Hochberg FDR correction.

**genes** - list of genes with predicted sites for the best case (min p-value).

**logo** - for XLSX and HTML format files contain `logo` column instead of `genes`.

### Example output file of`ESDEG annotation` command

| motif_id            | tf_name | tf_class                   | tf_family | log2(or)         | me_padj               | counts           | lfc                 | de_padj           | genes                               |
| ------------------- | ------- | -------------------------- | --------- | ---------------- | --------------------- | ---------------- | ------------------- | ----------------- | ----------------------------------- |
| IRF3.H12CORE.0.PS.A | IRF3    | Tryptophan cluster factors | IRF       | 3.44067709622015 | 2.44862817869797E-122 | 17.7506790122365 | -0.0822566108747206 | 0.999587214534912 | ENSG00000002549;ENSG00000010030;... |

Where:
**motif_id** - jaspar id from data base.

**tf_name** - name of Transcription factor

**tf_class** - class of DBD's transcrition factor

**tf_family** - family of DBD's transcrition factor

**jaspar_cluster** - сluster to which the motif belongs (https://jaspar.genereg.net/matrix-clusters/) [writen only for JASPAR database]

**log(or)** - it's -log2 transformation of  odds ratio (OR). For the option `fraction` it can be defined as follow $OR = N_f / N_b$, where $N_f$ - is the number of foreground promoters with predicted sites and $N_b$ - is a mean value of number of background promoters with predicted sites estimated by Monte-Carlo approach. Alternatively for the option `enrichment` it can be defined as follow $OR = N_f / N_b$, where $N_f$ - is a number of predicted sites in foreground and $N_b$ - is the mean value of number of predicted sites in background estimated by Monte-Carlo approach. The best log(or) value is represented in table. It's chosen by minimal p-value.

**me_padj** - it's adjusted p-value (significance for motif enrichment) obtained by using Benjamini-Hochberg FDR correction.

**counts** - it's normalized number of reads (expression level) for each TF.

**lfc** - it's log2(Fold Change) for each TF.

**de_padj** - it's adjusted p-value (significance for different expression) for each TF.

**genes** - list of genes with predicted sites for the best case (min p-value). This column is not available for XLSX and HTML files.

**logo** - for XLSX and HTML format files contain `logo` column instead of `genes`.

## Results visualization

![image](example.png)

## Citing

Oshchepkov, D., Chadaeva, I., Kozhemyakina, R., Shikhevich, S., Sharypova, E., Savinkova, L., Klimova, N. V., Tsukanov, A., Levitsky, V. G., & Markel, A. L. (2022). Transcription Factors as Important Regulators of Changes in Behavior through Domestication of Gray Rats: Quantitative Data from RNA Sequencing. International journal of molecular sciences, 23(20), 12269. https://doi.org/10.3390/ijms232012269

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

[^a]: Oshchepkov, Dmitry, Irina Chadaeva, Rimma Kozhemyakina, Svetlana Shikhevich, Ekaterina Sharypova, Ludmila Savinkova, Natalya V. Klimova, Anton Tsukanov, Victor G. Levitsky, and Arcady L. Markel. 2022. "Transcription Factors as Important Regulators of Changes in Behavior through Domestication of Gray Rats: Quantitative Data from RNA Sequencing" International Journal of Molecular Sciences 23, no. 20: 12269. https://doi.org/10.3390/ijms232012269
[^b]: Levitsky, V., Zemlyanskaya, E., Oshchepkov, D., Podkolodnaya, O., Ignatieva, E., Grosse, I., Mironova, V., & Merkulova, T. (2019). A single ChIP-seq dataset is sufficient for comprehensive analysis of motifs co-occurrence with MCOT package. Nucleic acids research, 47(21), e139. https://doi.org/10.1093/nar/gkz800
[^c]: Tsukanov, A. V., Mironova, V. V., & Levitsky, V. G. (2022). Motif models proposing independent and interdependent impacts of nucleotides are related to high and low affinity transcription factor binding sites in Arabidopsis. Frontiers in plant science, 13, 938545. https://doi.org/10.3389/fpls.2022.938545
[^1]: Castro-Mondragon, J. A., Riudavets-Puig, R., Rauluseviciute, I., Lemma, R. B., Turchi, L., Blanc-Mathieu, R., Lucas, J., Boddie, P., Khan, A., Manosalva Pérez, N., Fornes, O., Leung, T. Y., Aguirre, A., Hammal, F., Schmelter, D., Baranasic, D., Ballester, B., Sandelin, A., Lenhard, B., Vandepoele, K., … Mathelier, A. (2022). JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles. Nucleic acids research, 50(D1), D165–D173. https://doi.org/10.1093/nar/gkab1113
[^2]: Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics (Oxford, England), 26(6), 841–842. https://doi.org/10.1093/bioinformatics/btq033
[^3]: Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, *15*(12), 550. https://doi.org/10.1186/s13059-014-0550-8
[^4]: Monier, B., McDermaid, A., Wang, C., Zhao, J., Miller, A., Fennell, A., & Ma, Q. (2019). IRIS-EDA: An integrated RNA-Seq interpretation system for gene expression data analysis. *PLOS Computational Biology*, *15*(2), e1006792. https://doi.org/10.1371/journal.pcbi.1006792
[^5]: Cinar, O., & Viechtbauer, W. (2022). The poolr package for combining independent and dependent p values. Journal of Statistical Software, 101, 1-42. https://doi.org/10.18637/jss.v101.i01
