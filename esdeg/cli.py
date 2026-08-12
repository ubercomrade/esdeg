import os
import sys
import argparse
from esdeg.writers import write_table, write_table_ann, \
write_xlsx, write_xlsx_ann, write_report, create_picture
from esdeg.annotation import annotation
from esdeg.functions import esdeg as run


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser_name', help='Available commands:')

    enrichment_parser = subparsers.add_parser(name='enrichment', help='Run Monte-Carlo for motif enrichment evaluation on data (1. main step)')
    annotation_parser = subparsers.add_parser(name='annotation', help='Adding expression data of TF (2. optional step)')

    enrichment_parser.add_argument('promoters', action='store', metavar='promoters', help='Path to promoters in fasta format. All promoters have to be with same length. After the symbol ">" unique gene ID has to be written (>ENSG00000160072::1:1469886-1472284 or >ENSG00000160072)')
    enrichment_parser.add_argument('input', action='store', help='Path to input file for chosen data type.')
    enrichment_parser.add_argument('output', action='store', help='Output file in TSV format (e.g. /path/to/output/file.tsv).')
    enrichment_parser.add_argument('-f', '--format', action='store', choices=['deg', 'set'], default='deg', help='Data format that you want to analyse: for the `deg` option, input file is in CSV format with the results of RNA-seq analysis. \
                                                                                                                The file must contain the following columns: id, log2FoldChange, and padj. \
                                                                                                                For the `set` option, the input file is a TXT file where only the list of genes is written. Default = `deg`')
    enrichment_parser.add_argument('-m', '--motifs', action='store', choices=['jaspar', 'hocomoco'],  default='hocomoco', help='Choose motif DB: jaspar or hocomoco. Default = `hocomoco`')
    enrichment_parser.add_argument('-t', '--taxon', action='store', choices=['plants', 'vertebrates', 'insects', 'urochordates', 'nematodes', 'fungi'], default='vertebrates',
        help='Prepare database for respective JASPAR CORE taxonomic group of motifs. Possible options are plants, vertebrates, insects, urochordates, nematodes, fungi. \
        For more detailes see https://jaspar.uio.no/ and https://pyjaspar.readthedocs.io/en/latest/index.html')
    enrichment_parser.add_argument('-r', '--regulated', action='store', choices=['all', 'up', 'down'], default='all', metavar='N',
                        help='The parameter is used to choose up/down/all DEGs, default= all (used only for `deg` data type).')
    enrichment_parser.add_argument('-P', '--pvalue', action='store', type=float, default=0.05,
                        help='The pvalue is used as threshold to choose DEGs, default= 0.05 (used only for `deg` data type).')
    enrichment_parser.add_argument('-l', '--log2fc_deg', action='store', type=float, default=1.,
                        help='The absolute value of log2FoldChange used as threshold (L2FC_THR) to choose DEGs promoters (actual L2FC >= L2FC_THR OR actual L2FC <= -L2FC_THR), default= 1.0 (used only for `deg` data type).')
    enrichment_parser.add_argument('-L', '--log2fc_back', action='store', type=float, default=0.32192809488736235,
                        help='The absolute value of log2FoldChange used as threshold (L2FC_BACK_THR) to choose background promoters (-L2FC_BACK_THR <= actual L2FC <= L2FC_BACK_THR), default= log2(5/4)=0.321928... (used only for `deg` data type).')
    enrichment_parser.add_argument('-c', '--content', action='store', type=float, default=0.3,
                        help='The maximal GC content difference between promoters of foreground and background in Monte Carlo algorithm. \
                        Range of possible threshold [0.01 .. 1.0]. If threshold is equal to 1.0 then GC content is not taken into account. \
                        In this case (thr = 1.0) algorithm works faster. Default= 0.3. (used only for `deg` data type).')
    enrichment_parser.add_argument('-v', '--visualization', action='store', type=str, default='None',
                            help="Path to write interactive picture in HTML format (path/to/pic.html). if '-v' is given, then ESDEG creates picutre. By default it isn't used.")
    enrichment_parser.add_argument('-x', '--xlsx', action='store', type=str, default='None',
                            help="Path to write table with results in XLSX format (path/to/table.xlsx). XLSX table contains logo of motifs. if '-x' is given, then ESDEG creates XSLX table. By default it isn't used.")
    enrichment_parser.add_argument('-R', '--report', action='store', type=str, default='None',
                            help="Path to write interactive table with results in HTML format (path/to/report.html). HTML report contains logo of motifs. if '-r' is given, then ESDEG creates HTML report. By default it isn't used.")
    enrichment_parser.add_argument('-p', '--nproc', action='store', type=int, default=4, help='Number of processes to split the work between. Default= 4')


    annotation_parser.add_argument('esdeg', action='store', help='Path to esdeg result in TSV format')
    annotation_parser.add_argument('deg', action='store', help='Input file in CSV format with results of RNA-seq analysis. File must contain next columns: id, log2FoldChange, padj')
    annotation_parser.add_argument('counts', action='store', help='Input file in CSV format with normalized counts of read (FPKM or other). File must contain next columns: id, counts')
    annotation_parser.add_argument('gtf', action='store', help='Input file in GTF format with genome features annotation (ENSEMBL is prefered). It`s used to convert gene names and IDs')
    annotation_parser.add_argument('output', action='store', help='Path to write table with annotated ESDEG results. Table includes information related to expression level and different expression level of each TF/motif')
    annotation_parser.add_argument('-x', '--xlsx', action='store', type=str, default='None',
                            help="Path to write table with results in XLSX format (path/to/table.xlsx). XLSX table contains logo of motifs. if '-x' is given, then ESDEG creates XSLX table. By default it isn't used")
    annotation_parser.add_argument('-b', '--best', action='store_true', dest='best',
                        required=False, help='If this argument is used, then for each TF only the motif with best enrichment will be left (Some times TF have several motifs in DB)')
    annotation_parser.add_argument('-f', '--filter', action='store_true', dest='filter',
                        required=False, help='If this argument is used, then filtration will be applied to ESDEG table based on expression level of TFs.')
    annotation_parser.add_argument('-m', '--me_padj', action='store', type=float, default=0.05,
                        help='p-value threshold for motif enrichment (adj.pval column in table. adj.pval will be renamed to me_padj). Default value = 0.05. Applied only when flag --filter is used.')
    annotation_parser.add_argument('-l', '--log_odds', action='store', type=float, default=1.,
                        help='log(odds ratio) threshold for motif enrichment (log2(or) column in table). Default value = 1.0. Applied only when flag --filter is used.')
    annotation_parser.add_argument('-d', '--de_padj', action='store', type=float, default=0.05,
                        help='p-value threshold for DEGs. It`s used to found TF among DEGs. Default value = 0.05. Applied only when flag --filter is used.')
    annotation_parser.add_argument('-L', '--log_fc', action='store', type=float, default=1.,
                        help='log(fold change) threshold for DEGs. It`s used to found TF among DEGs. Default value = 1.0. Applied only when flag --filter is used.')
    annotation_parser.add_argument('-n', '--ncounts', action='store', type=float, default=5.,
                    help='Counts threshold (expression level threshold) for genes.  It`s used to remove TFs that are not expressed (With the exception of TFs, which are DEGs). Default value = 5.0. Applied only when flag --filter is used.')


    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main_cli():
    args = parse_args()

    if args.subparser_name == 'enrichment':
        motif_db = args.motifs
        path_to_promoters = args.promoters
        taxon = args.taxon
        type_of_data = args.format
        path_to_data = args.input
        path_to_output = args.output
        path_to_vis = args.visualization
        path_to_xlsx = args.xlsx
        path_to_report = args.report
        condition = args.regulated
        padj_thr= args.pvalue
        log2fc_thr_deg = args.log2fc_deg
        log2fc_thr_background = args.log2fc_back
        gc_threshold = args.content
        nproc = args.nproc

        if motif_db == 'hocomoco':
            taxon = 'human'

        df = run(motif_db,
                  taxon,
                  gc_threshold,
                  path_to_promoters,
                  path_to_data,
                  nproc,
                  type_of_data=type_of_data,
                  condition=condition,
                  log2fc_thr_deg=log2fc_thr_deg,
                  log2fc_thr_background=log2fc_thr_background,
                  padj_thr=padj_thr)

        if path_to_vis != 'None':
            create_picture(df, path_to_vis)
        if path_to_report != 'None':
            write_report(df, taxon, path_to_report)
        if path_to_xlsx != 'None':
            write_xlsx(df, taxon, path_to_xlsx)
        write_table(df, path_to_output)

    elif args.subparser_name == 'annotation':

        path_to_output = args.output
        deg_path = args.deg
        counts_path = args.counts
        esdeg_path = args.esdeg
        gtf_path = args.gtf

        path_to_xlsx = args.xlsx
        best_flag = args.best
        filter_flag = args.filter
        me_padj_thr = args.me_padj
        de_padj_thr = args.de_padj
        lor_thr = args.log_odds
        lfc_thr = args.log_fc
        counts_filter = args.ncounts

        df = annotation(deg_path, counts_path, esdeg_path, gtf_path, filter_flag, best_flag,
               me_padj_thr, de_padj_thr, lor_thr, lfc_thr, counts_filter)

        if path_to_xlsx != 'None':
            write_xlsx_ann(df, path_to_xlsx)
        write_table_ann(df, path_to_output)

    pass


if __name__ == '__main__':
    main_cli()
