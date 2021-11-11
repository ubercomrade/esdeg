import os
import sys
import argparse
from scipy import stats
from enrest.functions import *
from enrest.parsers import matrix_parser, hocomoco_parser, fasta_parser
from enrest.scaner import scaner


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser_name', help='Available commands:')
    matrix = subparsers.add_parser('matrix', help='Run method for only one given matrix')
    hocomoco = subparsers.add_parser('hocomoco', help='Run method for all matrixes from HOCOMOCO')
    
    matrix.add_argument('deg', action='store', help='TSV file with DEG with ..., The NAME column must contain ensemble gene IDS')
    matrix.add_argument('matrix', action='store', help='path to matrix (PCM) in HOCOMOCO format')
    matrix.add_argument('promoters', action='store', choices=['mm10', 'hg38', 'rnor6'], metavar='N',
         help='promoters of organism (hg38, mm10)')
    matrix.add_argument('output', action='store', help='path to write table with results')
    matrix.add_argument('-t', '--tag', action='store', type=str, dest='tag',
                        required=False, default='matrix', help='Used as name of matrix for first column in output, def=matrix')
    matrix.add_argument('-m', '--method', action='store', choices=['enrichment', 'fraction'],
                        metavar='N', type=str, default='enrichment', 
                        help='Realisation of montecarlo approach (enrichment or fraction), default= enrichment')
    
    hocomoco.add_argument('deg', action='store', help='TSV file with DEG with ..., The NAME column must contain ensemble gene IDS')
    hocomoco.add_argument('db', action='store', help='Path HOCOMOCO DB')
    hocomoco.add_argument('promoters', action='store', choices=['mm10', 'hg38'], metavar='N',
         help='promoters of organism (hg38, mm10)')
    hocomoco.add_argument('output', action='store', help='path to write table with results')  
    hocomoco.add_argument('-m', '--method', action='store', choices=['enrichment', 'fraction'],
                        metavar='N', type=str, default='enrichment', 
                        help='Realisation of montecarlo approach (enrichment or fraction), default= enrichment')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def run_montecarlo(deg_scores, other_scores, threshold_table, method):
    results_montecarlo = []
    for index in range(0, len(threshold_table)):
        threshold, fpr = threshold_table[index]
        if method == "enrichment":
            z_score = montecarlo_enrichment(deg_scores, other_scores, threshold)
        elif method == "fraction":
            z_score = montecarlo_fraction(deg_scores, other_scores, threshold)
        pval = stats.norm.sf(abs(z_score))
        results_montecarlo.append(pval)
    return results_montecarlo


#def run_fisher(deg_scores, other_scores, threshold_table):
#    results_fisher = []
#    for index in range(0, len(threshold_table)):
#        threshold, fpr = threshold_table[index]
#        pval = fisher_test(deg_scores, other_scores, threshold)
#        results_fisher.append(pval)
#    return results_fisher


# def create_plot(results, wdir, tag, condition, method):
#     fpr = [i[0] for i in results]
#     pvals = [i[1] for i in results]
#     with plt.style.context('ggplot'):
#         plt.plot(fpr, pvals, color='b')
#         plt.hlines(0.05, xmin=0, xmax=max(fpr), color='r')
#         plt.xlim(xmin=0, xmax=max(fpr))
#         plt.ylim(ymin=0, ymax=1)
#         plt.xlabel('FPR')
#         plt.ylabel('P-value')
#         plt.savefig(f"{wdir}/{tag}_{condition}_{method}.tif",
#                format="tif", dpi=400, bbox_inches='tight', pad_inches = 0.1)
#         plt.clf()
#         plt.cla()
#     pass
    
    
def matrix_case(args):
    path_to_deg = args.deg
    path_to_matrix = args.matrix
    output_path = args.output
    promoters = args.promoters
    method = args.method
    tag = args.tag
    
    this_dir, this_filename = os.path.split(__file__)
    if promoters == 'mm10':
        path_to_promoters = os.path.join(this_dir, "../data", "mm10.ensembl.promoters.fa")
    elif promoters == 'hg38':
        path_to_promoters = os.path.join(this_dir, "../data", "hg38.ensembl.promoters.fa")
    elif promoters == 'rnor6':
        path_to_promoters = os.path.join(this_dir, "../data", "rnor6.ensembl.promoters.fa")
        
    print('-'*30)
    print('Read DEG table')
    deg_table = pd.read_csv(path_to_deg, sep='\t', comment='#')
    deg_table = deg_table[deg_table['padj'] <= 1]
    print('-'*30)
    print('Read promoters')
    promoters = fasta_parser(path_to_promoters)
    print('-'*30)
    print('Read PWM')
    matrix, matrix_length, middle_score = matrix_parser(path_to_matrix)
    print('-'*30)
    print('Scan promotrers')
    all_results = scaner(promoters, matrix)
    best_results = get_best_scores(all_results)
    print('-'*30)
    print('Calculate threshold table')
    all_scores_flatten = get_all_flatten_scores(all_results)
    threshold_table = get_threshold(all_scores_flatten)
    threshold_table = np.array(threshold_table)
    fprs_table = threshold_table[:,1]
    fprs_choosen = [calculate_fpr(i) for i in range(9)]
    indexes = np.searchsorted(fprs_table, fprs_choosen)
    threshold_table = threshold_table[indexes]
    print('-'*30)
    print('Run tests:')
    container = []
    container.append([tag])
    for index, condition in enumerate(['ALL', 'UP', 'DOWN'], 1):
        print(f'{index}. {condition} - condition')
        deg_ids = get_deg_gene_ids(deg_table, condition)
        other_ids = get_other_gene_ids(deg_table)
        if method == "enrichment":
            deg_scores, other_scores = split_scores_by_gene_ids(all_results, deg_ids, other_ids)
        elif method == "fraction":
            deg_scores, other_scores = split_scores_by_gene_ids(best_results, deg_ids, other_ids)
        results = run_montecarlo(deg_scores, other_scores, threshold_table, method)
        container[-1] += results
    print('-'*30)
    head = '\t'.join(['FPR->'] + list(map(str, fprs_choosen))*3) + '\n' 
    head += '\t'.join(['ID'] + ['ALL']*9 + ['UP']*9 + ['DOWN']*9) + '\n'
    write_table(head, container, output_path)
    print('All done. Exit')
    pass


def hocomoco_case(args):
    path_to_deg = args.deg
    path_to_db = args.db
    output_path = args.output
    promoters = args.promoters
    method = args.method
    
    this_dir, this_filename = os.path.split(__file__)
    if promoters == 'mm10':
        path_to_promoters = os.path.join(this_dir, "../data", "mm10.ensembl.promoters.fa")
    elif promoters == 'hg38':
        path_to_promoters = os.path.join(this_dir, "../data", "hg38.ensembl.promoters.fa")
    elif promoters == 'rnor6':
        path_to_promoters = os.path.join(this_dir, "../data", "rnor6.ensembl.promoters.fa")
        
    print('-'*30)
    print('Read DEG table')
    deg_table = pd.read_csv(path_to_deg, sep='\t', comment='#')
    deg_table = deg_table[deg_table['padj'] <= 1]
    print('-'*30)
    print('Read promoters')
    promoters = fasta_parser(path_to_promoters)
    print('-'*30)
    print('Read HOCOMOCO')
    matrices = hocomoco_parser(path_to_db)
    number_of_matrices = len(matrices)
    print(f'Number of matrices = {number_of_matrices}')
    print('-'*30)
    container = []
    for i in range(number_of_matrices):
        name, matrix, matrix_length, middle_score = matrices[i]
        container.append([name])
        print(f'{i+1}. {name}:')
        print('Scan promotrers')
        all_results = scaner(promoters, matrix)
        best_results = get_best_scores(all_results)
        print('Calculate threshold table')
        all_scores_flatten = get_all_flatten_scores(all_results)
        threshold_table = get_threshold(all_scores_flatten)
        threshold_table = np.array(threshold_table)
        fprs_table = threshold_table[:,1]
        fprs_choosen = [calculate_fpr(i) for i in range(9)]
        indexes = np.searchsorted(fprs_table, fprs_choosen)
        threshold_table = threshold_table[indexes]
        print('Run tests:')
        for index, condition in enumerate(['ALL', 'UP', 'DOWN'], 1):
            print(f'{index}. {condition} - condition')
            deg_ids = get_deg_gene_ids(deg_table, condition)
            other_ids = get_other_gene_ids(deg_table)
            if method == "enrichment":
                deg_scores, other_scores = split_scores_by_gene_ids(all_results, deg_ids, other_ids)
            elif method == "fraction":
                deg_scores, other_scores = split_scores_by_gene_ids(best_results, deg_ids, other_ids)
            results = run_montecarlo(deg_scores, other_scores, threshold_table, method)
            container[-1] += results
        print('-'*30)  
    head = '\t'.join(['FPR->'] + list(map(str, fprs_choosen))*3) + '\n'   
    head += '\t'.join(['ID'] + ['ALL']*9 + ['UP']*9 + ['DOWN']*9) + '\n'
    write_table(head, container, output_path)
    print('All done. Exit')
    pass


def main():
    args = parse_args()    
    if args.subparser_name == 'matrix':
        matrix_case(args)
    elif args.subparser_name == 'hocomoco':
        hocomoco_case(args)   
    pass

    
if __name__ == '__main__':
    main()