import os
import sys
import argparse
import pkg_resources
from esdeg.set import set_case
from esdeg.deg import deg_case
from esdeg.preparation_hocomoco import prepare_motif_db as hocomoco_db
from esdeg.preparation_jaspar import prepare_motif_db as jaspar_db
from esdeg.writer import write_table, write_table_ann, \
write_xlsx, write_xlsx_ann, write_report, create_picture
from esdeg.annotation import annotation


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser_name', help='Available commands:')
    jaspar_preparation_parser = subparsers.add_parser('jaspar', help='Run preparation DB based on JASPAR 2024 motifs (motifs are selected based on taxon)')
    hocomoco_preparation_parser = subparsers.add_parser('hocomoco', help='Run preparation DB based on HOCOMOCO 2024 motifs (all core motifs are used [M. musculus and H. sapiens])')
    deg_parser = subparsers.add_parser('deg', help='Run test on DEGs')
    set_parser = subparsers.add_parser('set', help='Run test on SET of genes')
    annotation_parser = subparsers.add_parser('annotation', help='Add expression level and different expression data for each TF/motif to filter ESDEG results')

    jaspar_preparation_parser.add_argument('taxon', action='store', choices=['plants', 'vertebrates', 'insects', 'urochordates', 'nematodes', 'fungi'],
        help='Prepare database for respective JASPAR CORE taxonomic group of motifs. Possible options are plants, vertebrates, insects, urochordates, nematodes, fungi. \
        For more detailes see https://jaspar.uio.no/ and https://pyjaspar.readthedocs.io/en/latest/index.html')
    jaspar_preparation_parser.add_argument('promoters', action='store', metavar='promoters',
         help='Path to promoters in fasta format. All promoters have to be with same length. After the symbol ">" unique gene ID has to be written (>ENSG00000160072::1:1469886-1472284 or >ENSG00000160072)')
    jaspar_preparation_parser.add_argument('output', action='store', help='Name of directory to write output files')
    jaspar_preparation_parser.add_argument('-p', '--nproc', action='store', type=int, default=4, help='Number of processes to split the work between. Default= 4')

    hocomoco_preparation_parser.add_argument('promoters', action='store', metavar='promoters',
         help='Path to promoters in fasta format. All promoters have to be with same length. After the symbol ">" unique gene ID has to be written (>ENSG00000160072::1:1469886-1472284 or >ENSG00000160072)')
    hocomoco_preparation_parser.add_argument('output', action='store', help='Name of directory to write output files')
    hocomoco_preparation_parser.add_argument('-p', '--nproc', action='store', type=int, default=4, help='Number of processes to split the work between. Default= 4')

    deg_parser.add_argument('deg', action='store', help='Input file in CSV format with results of RNA-seq analysis. File must contain next columns: id, log2FoldChange, padj')
    deg_parser.add_argument('matrices', action='store', help='Directory with prepared database contained .npy files (e.g. /path/to/database)')
    deg_parser.add_argument('output', action='store', help='Output file in TSV format (e.g. /path/to/output/file.tsv)')
    deg_parser.add_argument('-v', '--visualization', action='store', type=str, default='None',
                            help="Path to write interactive picture in HTML format (path/to/pic.html). if '-v' is given, then ESDEG creates picutre. By default it isn't used")
    deg_parser.add_argument('-x', '--xlsx', action='store', type=str, default='None',
                            help="Path to write table with results in XLSX format (path/to/table.xlsx). XLSX table contains logo of motifs. if '-x' is given, then ESDEG creates XSLX table. By default it isn't used")
    deg_parser.add_argument('-R', '--report', action='store', type=str, default='None',
                            help="Path to write interactive table with results in HTML format (path/to/report.html). HTML report contains logo of motifs. if '-r' is given, then ESDEG creates HTML report. By default it isn't used")
    deg_parser.add_argument('-p', '--parameter', action='store', choices=['enrichment', 'fraction'],
                        metavar='PARAMETER', type=str, default='enrichment',
                        help='Parameter estimated in test (enrichment or fraction), default= enrichment')
    deg_parser.add_argument('-r', '--regulated', action='store', choices=['all', 'up', 'down'], default='all', metavar='N',
                        help='The parameter is used to choose up/down/all DEGs, default= all')
    deg_parser.add_argument('-P', '--pvalue', action='store', type=float, default=0.05,
                        help='The pvalue is used as threshold to choose DEGs, default= 0.05')
    deg_parser.add_argument('-l', '--log2fc_deg', action='store', type=float, default=1.,
                        help='The absolute value of log2FoldChange used as threshold (L2FC_THR) to choose DEGs promoters (actual L2FC >= L2FC_THR OR actual L2FC <= -L2FC_THR), default= 1.0')
    deg_parser.add_argument('-L', '--log2fc_back', action='store', type=float, default=0.32192809488736235,
                        help='The absolute value of log2FoldChange used as threshold (L2FC_BACK_THR) to choose background promoters (-L2FC_BACK_THR <= actual L2FC <= L2FC_BACK_THR), default= log2(5/4)=0.321928...')
    deg_parser.add_argument('-c', '--content', action='store', type=float, default=0.3,
                        help='The maximal GC content difference between promoters of foreground and background in Monte Carlo algorithm. \
                        Range of possible threshold [0.01 .. 1.0]. If threshold is equal to 1.0 then GC content is not taken into account. \
                        In this case (thr = 1.0) algorithm works faster. Default= 0.3.')

    set_parser.add_argument('set', action='store', help='File with list of genes.')
    set_parser.add_argument('matrices', action='store', help='Path to prepared data base of matrices')
    set_parser.add_argument('output', action='store', help='Path to write table with results')
    set_parser.add_argument('-v', '--visualization', action='store', type=str, default='None',
                            help="Path to write interactive picture in HTML format (path/to/pic.html). if '-v' is given, then ESDEG creates picutre. By default it isn't used")
    set_parser.add_argument('-x', '--xlsx', action='store', type=str, default='None',
                            help="Path to write table with results in XLSX format (path/to/table.xlsx). XLSX table contains logo of motifs. if '-x' is given, then ESDEG creates XSLX table. By default it isn't used")
    set_parser.add_argument('-R', '--report', action='store', type=str, default='None',
                            help="Path to write interactive table with results in HTML format (path/to/report.html). HTML report contains logo of motifs. if '-r' is given, then ESDEG creates HTML report. By default it isn't used")
    set_parser.add_argument('-p', '--parameter', action='store', choices=['enrichment', 'fraction'],
                        metavar='PARAMETER', type=str, default='enrichment',
                        help='Parameter estimated in test (enrichment or fraction), default= enrichment')
    set_parser.add_argument('-c', '--content', action='store', type=float, default=0.3,
                        help='The maximal GC content difference between promoters of foreground and background in Monte Carlo algorithm. \
                        Range of possible threshold [0.01 .. 1.0]. If threshold is equal to 1.0 then GC content is not taken into account. \
                        In this case (thr = 1.0) algorithm works faster. Default= 0.3.')


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


def main():
    args = parse_args()
    if args.subparser_name == 'deg':
        path_to_deg = args.deg
        path_to_db = args.matrices
        path_to_output = args.output
        path_to_vis = args.visualization
        path_to_xlsx = args.xlsx
        path_to_report = args.report
        parameter = args.parameter
        condition = args.regulated
        padj_thr= args.pvalue
        log2fc_thr_deg = args.log2fc_deg
        log2fc_thr_background = args.log2fc_back
        gc_threshold = args.content

        df, taxon = deg_case(path_to_deg,
                 path_to_db,
                 parameter=parameter,
                 padj_thr=padj_thr,
                 log2fc_thr_deg=log2fc_thr_deg,
                 log2fc_thr_background=log2fc_thr_background,
                 gc_threshold=gc_threshold,
                 condition=condition)
        if path_to_vis != 'None':
            create_picture(df, path_to_vis)
        if path_to_report != 'None':
            write_report(df, taxon, path_to_report)
        if path_to_xlsx != 'None':
            write_xlsx(df, taxon, path_to_xlsx)
        write_table(df, path_to_output)

    elif args.subparser_name == 'set':
        path_to_set = args.set
        path_to_db = args.matrices
        path_to_output = args.output
        path_to_vis = args.visualization
        path_to_xlsx = args.xlsx
        path_to_report = args.report
        parameter = args.parameter
        gc_threshold = args.content

        df, taxon = set_case(path_to_set,
                 path_to_db,
                 parameter=parameter,
                 gc_threshold=gc_threshold)
        if path_to_vis != 'None':
            create_picture(df, path_to_vis)
        if path_to_report != 'None':
            write_report(df, taxon, path_to_report)
        if path_to_xlsx != 'None':
            write_xlsx(df, taxon, path_to_xlsx)
        write_table(df, path_to_output)

    elif args.subparser_name == 'jaspar':
        taxon = args.taxon
        output_dir = args.output
        path_to_promoters = args.promoters
        nproc = args.nproc
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        jaspar_db(output_dir,
            path_to_promoters,
            taxon,
            nproc)


    elif args.subparser_name == 'hocomoco':
        output_dir = args.output
        path_to_promoters = args.promoters
        nproc = args.nproc
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        hocomoco_db(output_dir,
            path_to_promoters,
            nproc)

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
    main()
