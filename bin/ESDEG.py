import os
import sys
import argparse
from esdeg.set import set_case
from esdeg.deg import deg_case
from esdeg.preparation import prepare_motif_db
from esdeg.writer import write_table, create_picture


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser_name', help='Available commands:')
    preparation_parser = subparsers.add_parser('preparation', help='Run data base preparation')
    deg_parser = subparsers.add_parser('deg', help='Run test on DEGs')
    set_parser = subparsers.add_parser('set', help='Run test on SET of genes')

    preparation_parser.add_argument('taxon', action='store', choices=['plants', 'vertebrates', 'insects', 'urochordates', 'nematodes', 'fungi'],
        help='Prepare database for respective JASPAR CORE taxonomic group of motifs. Possible options are plants, vertebrates, insects, urochordates, nematodes, fungi. \
        For more detailes see https://jaspar.uio.no/ and https://pyjaspar.readthedocs.io/en/latest/index.html')
    preparation_parser.add_argument('organism', action='store', choices=['mm', 'hs', 'dm', 'dr', 'rn', 'ce', 'tair10', 'rnor6', 'rnor6_ucsc'], metavar='N',
         help='Promoters of organism (hs - H. sapiens, mm - M. musculus, dm - D. melanogaster, dr - D. rerio, rn - R. norvegicus, ce - C. elegans, tair10)')
    preparation_parser.add_argument('output', action='store', help='Name of directory to write output files')

    deg_parser.add_argument('deg', action='store', help='TSV file with DEG with ..., The NAME column must contain ensemble gene IDS')
    deg_parser.add_argument('matrices', action='store', help='Path to prepared data base of matrices')
    deg_parser.add_argument('organism', action='store', choices=['mm', 'hs', 'dm', 'dr', 'rn', 'ce', 'tair10', 'rnor6', 'rnor6_ucsc'], metavar='N',
         help='Organism (hs - H. sapiens, mm - M. musculus, dm - D. melanogaster, dr - D. rerio, rn - R. norvegicus, ce - C. elegans, tair10)')
    deg_parser.add_argument('output', action='store', help='Path to write table with results')
    deg_parser.add_argument('-v', '--visualization', action='store', type=str, default='None',
                            help="Path to write interactive picture in HTML format (path/to/pic.html). if '--v' is given, then ESDEG creates picutre. By default it isn't used")
    deg_parser.add_argument('-p', '--parameter', action='store', choices=['enrichment', 'fraction'],
                        metavar='PARAMETER', type=str, default='enrichment',
                        help='Parameter estimated in test (enrichment or fraction), default= enrichment')
    deg_parser.add_argument('-r', '--regulated', action='store', choices=['all', 'up', 'down'], default='all', metavar='N',
                        help='The parameter is used to choose up/down/all DEGs, default= all')
    deg_parser.add_argument('-P', '--pvalue', action='store', type=float, default=0.05,
                        help='The pvalue is used as threshold to choose DEGs, default= 0.05')
    deg_parser.add_argument('-l', '--log2fc_deg', action='store', type=float, default=1.,
                        help='The absolute value of log2FoldChange used as threshold to choose DEGs promoters (DEGs >= thr OR DEGs <= -thr), default= 1')
    deg_parser.add_argument('-L', '--log2fc_back', action='store', type=float, default=0.32192809488736235,
                        help='The absolute value of log2FoldChange used as threshold to choose background promoters (-thr <= BACK <= thr), default= log2(5/4)')
    deg_parser.add_argument('-c', '--content', action='store', type=float, default=0.3,
                        help='The maximal GC content difference between promoters of foreground and background in Monte Carlo algorithm. \
                        Range of possible threshold [0.01 .. 1.0]. If threshold is equal to 1.0 then GC content is not taken into account. \
                        In this case (thr = 1.0) algorithm works faster. Default= 0.3.')

    set_parser.add_argument('set', action='store', help='File with list of genes. Genes must be in Ensemble format (ensemble gene IDS)')
    set_parser.add_argument('matrices', action='store', help='Path to prepared data base of matrices')
    set_parser.add_argument('organism', action='store', choices=['mm', 'hs', 'dm', 'dr', 'rn', 'ce', 'tair10', 'rnor6', 'rnor6_ucsc'], metavar='N',
         help='Organism (hs - H. sapiens, mm - M. musculus, dm - D. melanogaster, dr - D. rerio, rn - R. norvegicus, ce - C. elegans, tair10)')
    set_parser.add_argument('output', action='store', help='Path to write table with results')
    set_parser.add_argument('-v', '--visualization', action='store', type=str, default='None',
                            help="Path to write interactive picture in HTML format (path/to/pic.html). if '--v' is given, then ESDEG creates picutre. By default it isn't used")
    set_parser.add_argument('-p', '--parameter', action='store', choices=['enrichment', 'fraction'],
                        metavar='PARAMETER', type=str, default='enrichment',
                        help='Parameter estimated in test (enrichment or fraction), default= enrichment')
    set_parser.add_argument('-c', '--content', action='store', type=float, default=0.3,
                        help='The maximal GC content difference between promoters of foreground and background in Monte Carlo algorithm. \
                        Range of possible threshold [0.01 .. 1.0]. If threshold is equal to 1.0 then GC content is not taken into account. \
                        In this case (thr = 1.0) algorithm works faster. Default= 0.3.')

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
        organism = args.organism
        parameter = args.parameter
        condition = args.regulated
        padj_thr= args.pvalue
        log2fc_thr_deg = args.log2fc_deg
        log2fc_thr_background = args.log2fc_back
        gc_threshold = args.content

        df = deg_case(path_to_deg,
                 path_to_db,
                 organism,
                 parameter=parameter,
                 padj_thr=padj_thr,
                 log2fc_thr_deg=log2fc_thr_deg,
                 log2fc_thr_background=log2fc_thr_background,
                 gc_threshold=gc_threshold,
                 condition=condition)
        if path_to_vis != 'None':
            create_picture(df, path_to_vis)
        write_table(df, path_to_output)

    elif args.subparser_name == 'set':
        path_to_set = args.set
        path_to_db = args.matrices
        path_to_output = args.output
        path_to_vis = args.visualization
        organism = args.organism
        parameter = args.parameter
        gc_threshold = args.content

        df = set_case(path_to_set,
                 path_to_db,
                 organism,
                 parameter=parameter,
                 gc_threshold=gc_threshold)
        if path_to_vis != 'None':
            create_picture(df, path_to_vis)
        write_table(df, path_to_output)

    elif args.subparser_name == 'preparation':
        taxon = args.taxon
        output_dir = args.output
        organism = args.organism

        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        this_dir, this_filename = os.path.split(__file__)
        if organism == 'mm':
            path_to_promoters = pkg_resources.resource_filename('esdeg', 'data/mm.ensembl.promoters.fa.xz')
        elif organism == 'hs':
            path_to_promoters = pkg_resources.resource_filename('esdeg', 'data/hs.ensembl.promoters.fa.xz')
        elif organism == 'dm':
            path_to_promoters = pkg_resources.resource_filename('esdeg', 'data/dm.ensembl.promoters.fa.xz')
        elif organism == 'dr':
            path_to_promoters = pkg_resources.resource_filename('esdeg', 'data/dr.ensembl.promoters.fa.xz')
        elif organism == 'rn':
            path_to_promoters = pkg_resources.resource_filename('esdeg', 'data/rn.ensembl.promoters.fa.xz')
        elif organism == 'ce':
            path_to_promoters = pkg_resources.resource_filename('esdeg', 'data/ce.ensembl.promoters.fa.xz')
        elif organism == 'tair10':
            path_to_promoters = pkg_resources.resource_filename('esdeg', 'data/tair10.ensembl.promoters.fa.xz')
        elif organism == 'rnor6_ucsc':
            path_to_promoters = pkg_resources.resource_filename('esdeg', 'data/rnor6.ucsc.promoters.fa.xz')

        prepare_motif_db(output_dir,
            path_to_promoters,
            organism,
            taxon)
    pass


if __name__ == '__main__':
    main()
