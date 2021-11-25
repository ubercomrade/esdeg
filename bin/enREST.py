import os
import sys
import argparse
from scipy import stats
from enrest.functions import *
from enrest.parsers import matrices_parser, fasta_parser
from enrest.scaner import scaner


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('deg', action='store', help='TSV file with DEG with ..., The NAME column must contain ensemble gene IDS')
    parser.add_argument('matrices', action='store', help='path to matrices in HOCOMOCO (PCM) or in MEME (PFM) format')
    parser.add_argument('promoters', action='store', choices=['mm10', 'hg38', 'tair10', 'rnor6'], metavar='N',
         help='promoters of organism (hg38, mm10, tair10)')
    parser.add_argument('output', action='store', help='path to write table with results')
    parser.add_argument('-m', '--method', action='store', choices=['montecarlo', 'binom', 'hypergeom'],
                        metavar='METHOD', type=str, default='enrichment', 
                        help='Realisation of montecarlo approach (enrichment or fraction), default= enrichment')
    parser.add_argument('-p', '--parametr', action='store', choices=['enrichment', 'fraction'],
                        metavar='METHOD', type=str, default='enrichment', 
                        help='Parametr to test (enrichment or fraction), default= enrichment')
    parser.add_argument('-f', '--format', action='store', choices=['meme', 'hocomoco'],
                        metavar='FORMAT', type=str, default='meme', 
                        help='Format of file with matrices (meme or hocomoco), default= meme')
    parser.add_argument('-P', '--pvalue', action='store', type=float, default=0.05, 
                        help='The pvalue is used as threshold to choose DEGs, default= 0.05')
    parser.add_argument('-l', '--log2fc_deg', action='store', type=float, default=1., 
                        help='The absolute value of log2FoldChange used as threshold to choose DEGs promoters (DEGs >= thr OR DEGs <= -thr), default= 1')
    parser.add_argument('-L', '--log2fc_back', action='store', type=float, default=0.32192809488736235, 
                        help='The absolute value of log2FoldChange used as threshold to choose background promoters (-thr <= BACK <= thr), default= log2(5/4)')    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


#OLD
def run_test(deg_scores, other_scores, threshold_table, method, parametr):
    results = []
    for index in range(0, len(threshold_table)):
        threshold, fpr = threshold_table[index]
        if method == 'montecarlo':
            if parametr == "enrichment":
                z_score = montecarlo_enrichment(deg_scores, other_scores, threshold)
            elif parametr == "fraction":
                z_score = montecarlo_fraction(deg_scores, other_scores, threshold)
            pval = stats.norm.sf(abs(z_score))
        else:
            if parametr == "enrichment":
                pval = stat_tests_enrichment(deg_scores, other_scores, threshold, method=method)
            elif parametr == "fraction":
                pval = stat_tests_fraction(deg_scores, other_scores, threshold, method=method)            
        results.append(pval)
    cpval = hartung(np.array(results))
    results.append(cpval)
    return results


#ANOTHER APPROACH
# def run_montecarlo(deg_scores, other_scores, threshold_table, method):
#     results_montecarlo = []
#     for line1, line2 in pairwise(threshold_table):
#         threshold_min, fpr_min = line1
#         threshold_max, fpr_max = line2
#         if method == "enrichment":
#             z_score = montecarlo_enrichment(deg_scores, other_scores, threshold_min, threshold_max)
#         elif method == "fraction":
#             z_score = montecarlo_fraction(deg_scores, other_scores, threshold_min, threshold_max)
#         pval = stats.norm.sf(abs(z_score))
#         results_montecarlo.append(pval)
#     st, cpval = stats.combine_pvalues(results_montecarlo)
#     results_montecarlo.append(cpval)
#     return results_montecarlo
    
    
def main():
    args = parse_args()
    path_to_deg = args.deg
    path_to_db = args.db
    output_path = args.output
    promoters = args.promoters
    parametr = args.parametr
    method = args.method
    file_format = args.format
    padj_thr= args.pvalue
    log2fc_thr_deg = args.log2fc_deg
    log2fc_thr_background = args.log2fc_back
    
    this_dir, this_filename = os.path.split(__file__)
    if promoters == 'mm10':
        path_to_promoters = os.path.join(this_dir, "../data", "mm10.ensembl.promoters.fa")
    elif promoters == 'hg38':
        path_to_promoters = os.path.join(this_dir, "../data", "hg38.ensembl.promoters.fa")
    elif promoters == 'rnor6':
        path_to_promoters = os.path.join(this_dir, "../data", "rnor6.ensembl.promoters.fa")
        
    print('-'*30)
    print('Read DEG table')
    deg_table = pd.read_csv(path_to_deg, sep=',', comment='#')
    deg_table = deg_table[deg_table['padj'] <= 1]
    print('-'*30)
    print('Read promoters')
    promoters = fasta_parser(path_to_promoters)
    print('-'*30)
    print('Read matrices')
    matrices = matrices_parser(path_to_db, f=file_format)
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
        fprs_choosen = np.array([0.0005, 0.00015, 0.00005]) # LOW, MIDDLE, HIGH
        indexes = np.searchsorted(fprs_table, fprs_choosen)
        threshold_table = threshold_table[indexes]
        print('Run tests:')
        for index, condition in enumerate(['ALL', 'UP', 'DOWN'], 1):
            print(f'{index}. {condition} - condition')
            deg_ids = get_deg_gene_ids(deg_table, condition, padj_thr=padj_thr, log2fc_thr=log2fc_thr_deg)
            other_ids = get_other_gene_ids(deg_table, padj_thr=padj_thr, log2fc_thr=log2fc_thr_background)
            if parametr == "enrichment":
                deg_scores, other_scores = split_scores_by_gene_ids(all_results, deg_ids, other_ids)
            elif parametr == "fraction":
                deg_scores, other_scores = split_scores_by_gene_ids(best_results, deg_ids, other_ids)
            results = run_test(deg_scores, other_scores, threshold_table, method, parametr)
            container[-1] += results
        print('-'*30)  
    head = '\t'.join(['FPR->'] + ['LOW', 'MIDDLE', 'HIGHT', 'COMMON']*3) + '\n' # new
    head += '\t'.join(['ID'] + ['ALL']*4 + ['UP']*4 + ['DOWN']*4) + '\n' #new
    write_table(head, container, output_path)
    print('All done. Exit')
    pass

    
if __name__ == '__main__':
    main()