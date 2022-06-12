import os
import sys
import argparse
from enrest.functions import *
from enrest.set import set_case
from enrest.deg import deg_case
from enrest.fasta import fasta_case

def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser_name', help='Available commands:')
    deg_parser = subparsers.add_parser('deg', help='Run test on DEGs')
    set_parser = subparsers.add_parser('set', help='Run test on SET of genes')
    fasta_parser = subparsers.add_parser('fasta', help='Run test on user FASTA files')

    deg_parser.add_argument('deg', action='store', help='TSV file with DEG with ..., The NAME column must contain ensemble gene IDS')
    deg_parser.add_argument('matrices', action='store', help='Path to matrices in HOCOMOCO (PCM) or in MEME (PFM) format')
    deg_parser.add_argument('promoters', action='store', choices=['mm10', 'hg38', 'tair10', 'rnor6'], metavar='N',
         help='promoters of organism (hg38, mm10, tair10)')
    deg_parser.add_argument('output', action='store', help='Name of directory for output files')
    deg_parser.add_argument('-p', '--parameter', action='store', choices=['enrichment', 'fraction'],
                        metavar='PARAMETER', type=str, default='enrichment', 
                        help='Parameter estimated in test (enrichment or fraction), default= enrichment')
    deg_parser.add_argument('-f', '--format', action='store', choices=['meme', 'hocomoco'],
                        metavar='FORMAT', type=str, default='meme', 
                        help='Format of file with matrices (meme or hocomoco), default= meme')
    deg_parser.add_argument('-P', '--pvalue', action='store', type=float, default=0.05, 
                        help='The pvalue is used as threshold to choose DEGs, default= 0.05')
    deg_parser.add_argument('-l', '--log2fc_deg', action='store', type=float, default=1., 
                        help='The absolute value of log2FoldChange used as threshold to choose DEGs promoters (DEGs >= thr OR DEGs <= -thr), default= 1')
    deg_parser.add_argument('-L', '--log2fc_back', action='store', type=float, default=0.32192809488736235, 
                        help='The absolute value of log2FoldChange used as threshold to choose background promoters (-thr <= BACK <= thr), default= log2(5/4)')
    deg_parser.add_argument('-c', '--content', action='store', type=float, default=0.1, 
                        help='The maximal GC content difference between promoters of interest and generated sample in Monte Carlo algorithm. default= 0.1')     

    set_parser.add_argument('set', action='store', help='File with list of genes. Genes must be in Ensemble format (ensemble gene IDS)')
    set_parser.add_argument('matrices', action='store', help='Path to matrices in HOCOMOCO (PCM) or in MEME (PFM) format')
    set_parser.add_argument('promoters', action='store', choices=['mm10', 'hg38', 'tair10', 'rnor6'], metavar='N',
         help='promoters of organism (hg38, mm10, tair10)')
    set_parser.add_argument('output', action='store', help='Name of directory for output files')
    set_parser.add_argument('-p', '--parameter', action='store', choices=['enrichment', 'fraction'],
                        metavar='PARAMETER', type=str, default='enrichment', 
                        help='Parameter estimated in test (enrichment or fraction), default= enrichment')
    set_parser.add_argument('-f', '--format', action='store', choices=['meme', 'hocomoco'],
                        metavar='FORMAT', type=str, default='meme', 
                        help='Format of file with matrices (meme or hocomoco), default= meme')
    set_parser.add_argument('-c', '--content', action='store', type=float, default=0.1, 
                        help='The maximal GC content difference between promoters of interest and generated sample in Monte Carlo algorithm. default= 0.1')     


    fasta_parser.add_argument('foreground', action='store', help='Fasta file with sequences are used as foreground')
    fasta_parser.add_argument('background', action='store', help='Fasta file with sequences are used as background')
    fasta_parser.add_argument('matrices', action='store', help='Path to matrices in HOCOMOCO (PCM) or in MEME (PFM) format')
    fasta_parser.add_argument('promoters', action='store', choices=['mm10', 'hg38', 'tair10', 'rnor6'], metavar='N',
         help='promoters of organism (hg38, mm10, tair10)')
    fasta_parser.add_argument('output', action='store', help='Name of directory for output files')
    fasta_parser.add_argument('-p', '--parameter', action='store', choices=['enrichment', 'fraction'],
                        metavar='PARAMETER', type=str, default='enrichment', 
                        help='Parameter estimated in test (enrichment or fraction), default= enrichment')
    fasta_parser.add_argument('-f', '--format', action='store', choices=['meme', 'hocomoco'],
                        metavar='FORMAT', type=str, default='meme', 
                        help='Format of file with matrices (meme or hocomoco), default= meme')
    fasta_parser.add_argument('-c', '--content', action='store', type=float, default=0.1, 
                        help='The maximal GC content difference between foreground and generated sample in Monte Carlo algorithm. default= 0.1')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()    
    if args.subparser_name == 'deg':
        path_to_deg = args.deg
        path_to_db = args.matrices
        output_dir = args.output
        promoters = args.promoters
        parameter = args.parameter
        file_format = args.format
        padj_thr= args.pvalue
        log2fc_thr_deg = args.log2fc_deg
        log2fc_thr_background = args.log2fc_back
        gc_threshold = args.content
        
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        this_dir, this_filename = os.path.split(__file__)
        if promoters == 'mm10':
            path_to_promoters = os.path.join(this_dir, "../data", "mm10.ensembl.promoters.fa.xz")
        elif promoters == 'hg38':
            path_to_promoters = os.path.join(this_dir, "../data", "hg38.ensembl.promoters.fa.xz")
        elif promoters == 'tair10':
            path_to_promoters = os.path.join(this_dir, "../data", "tair10.ensembl.promoters.fa.xz")
        elif promoters == 'rnor6':
            path_to_promoters = os.path.join(this_dir, "../data", "rnor6.ensembl.promoters.fa.xz")
            
        deg_case(path_to_deg, 
                 path_to_db,
                 output_dir, 
                 path_to_promoters, 
                 file_format=file_format, 
                 parameter=parameter, 
                 padj_thr=padj_thr,
                 log2fc_thr_deg=log2fc_thr_deg, 
                 log2fc_thr_background=log2fc_thr_background,
                 gc_threshold=gc_threshold)
        
    elif args.subparser_name == 'set':
        path_to_set = args.set
        path_to_db = args.matrices
        output_dir = args.output
        promoters = args.promoters
        parameter = args.parameter
        file_format = args.format
        gc_threshold = args.content

        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
    
        this_dir, this_filename = os.path.split(__file__)
        if promoters == 'mm10':
            path_to_promoters = os.path.join(this_dir, "../data", "mm10.ensembl.promoters.fa.xz")
        elif promoters == 'hg38':
            path_to_promoters = os.path.join(this_dir, "../data", "hg38.ensembl.promoters.fa.xz")
        elif promoters == 'tair10':
            path_to_promoters = os.path.join(this_dir, "../data", "tair10.ensembl.promoters.fa.xz")
        elif promoters == 'rnor6':
            path_to_promoters = os.path.join(this_dir, "../data", "rnor6.ensembl.promoters.fa.xz")
        set_case(path_to_set, 
                 path_to_db, 
                 output_dir, 
                 path_to_promoters, 
                 file_format=file_format, 
                 parameter=parameter,
                 gc_threshold=gc_threshold)

    elif args.subparser_name == 'fasta':
        path_to_foreground = args.foreground
        path_to_background = args.background
        path_to_db = args.matrices
        output_dir = args.output
        promoters = args.promoters
        parameter = args.parameter
        file_format = args.format
        gc_threshold = args.content

        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
    
        this_dir, this_filename = os.path.split(__file__)
        if promoters == 'mm10':
            path_to_promoters = os.path.join(this_dir, "../data", "mm10.ensembl.promoters.fa.xz")
        elif promoters == 'hg38':
            path_to_promoters = os.path.join(this_dir, "../data", "hg38.ensembl.promoters.fa.xz")
        elif promoters == 'tair10':
            path_to_promoters = os.path.join(this_dir, "../data", "tair10.ensembl.promoters.fa.xz")
        elif promoters == 'rnor6':
            path_to_promoters = os.path.join(this_dir, "../data", "rnor6.ensembl.promoters.fa.xz")
        fasta_case(path_to_foreground,
                 path_to_background,
                 path_to_db, 
                 output_dir, 
                 path_to_promoters, 
                 file_format=file_format, 
                 parameter=parameter,
                 gc_threshold=gc_threshold)      
    pass

    
if __name__ == '__main__':
    main()
