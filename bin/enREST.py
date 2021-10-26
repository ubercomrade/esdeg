import os
import sys
import argparse
import matplotlib.pyplot as plt
from scipy import stats
from enrest.functions import *
from enrest.parsers import matrix_parser, fasta_parser
from enrest.scaner import scaner


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('deg', action='store', help='TSV file with DEG with ..., The NAME column must contain ensemble gene IDS')
    parser.add_argument('matrix', action='store', help='path to matrix (PCM) in HOCOMOCO format')
    parser.add_argument('promoters', action='store', choices=['mm10', 'hg38', 'tair10', 'b73'], metavar='N',
         help='promoters of organism (hg38, mm10, tair10, b73)')
    parser.add_argument('output', action='store', help='output dir')
    parser.add_argument('-t', '--tag', action='store', type=str, dest='tag',
                        required=False, default='results', help='prefix of output files -> $tag_$otherpart')
    parser.add_argument('-m', '--method', action='store', choices=['montecarlo', 'fisher'],
                        metavar='N', type=str, default='montecarlo', 
                        help='Method to calculate statistics of enrichment (montecarlo or fisher), default= montecarlo')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def run_montecarlo(deg_scores, other_scores, threshold):
    results_montecarlo = []
    for index in range(0, len(threshold_table), 3):
        threshold, fpr = threshold_table[index]
        z_score = montecarlo(deg_scores, other_scores, threshold)
        results_montecarlo.append((threshold, fpr, stats.norm.sf(abs(z_score))))
    return results_montecarlo


def run_fisher(deg_scores, other_scores, threshold):
    results_fisher = []
    for index in range(0, len(threshold_table), 3):
        threshold, fpr = threshold_table[index]
        pval = fisher_test(deg_scores, other_scores, threshold)
        results_fisher.append((threshold, fpr, pval))
    return results_fisher


def create_plot(results, wdir, tag, condition, method):
    fpr = [i[1] for i in results]
    pvals = [i[2] for i in results]
    with plt.style.context('ggplot'):
        plt.plot(fpr, pvals, color='b')
        plt.hlines(0.05, xmin=0, xmax=max(fpr), color='r')
        plt.xlim(xmin=0, xmax=max(fpr))
        plt.ylim(ymin=0, ymax=1)
        plt.xlabel('FPR')
        plt.ylabel('P-value')
        plt.savefig(f"{wdir}/{tag}_{condition}_{method}.tif",
               format="tif", dpi=400, bbox_inches='tight', pad_inches = 0.1)
        plt.clf()
        plt.cla()
    pass
    

def main():
    args = parse_args()
    path_to_deg = args.deg
    path_to_matrix = args.matrix
    wdir = args.output
    promoters = args.promoters
    tag = args.tag
    method = args.method
    
    if not os.path.isdir(wdir):
        os.mkdir(wdir)
    
    this_dir, this_filename = os.path.split(__file__)
    if promoters == 'mm10':
        path_to_promoters = os.path.join(this_dir, "data", "mm10.ensembl.promoters.fa")
    elif promoters == 'hg38':
        path_to_promoters = os.path.join(this_dir, "data", "hg38.ensembl.promoters.fa")
    # elif promoters == 'tair10':
    #     path_to_promoters = os.path.join(this_dir, "data", "tair10.fasta")
    # elif promoters == 'b73':
    #     path_to_promoters = os.path.join(this_dir, "data", "b73_v5.fasta")
        
    print('-'*30)
    print('Read DEG table')
    deg_table = pd.read_csv(path_to_deg, sep='\t', comment='#')
    deg_table = deg_table[deg_table['padj'] <= 1]
    print('-'*30)
    print('Read PWM')
    matrix, matrix_length, middle_score = matrix_parser(path_to_matrix)
    print('-'*30)
    print('Read promoters')
    promoters = fasta_parser(promoters_fasta_path)
    print('-'*30)
    print('Promoters scanning')
    scan_results = scan_promoters(promoters, matrix)
    best_results = get_best_scores(scan_results)
    print('-'*30)
    print('Calculate threshold table')
    all_scores_flatten = get_all_flatten_scores(scan_results)
    threshold_table = get_threshold(all_scores_flatten)
    threshold_table = threshold_table[::-1]  
    print('-'*30)
    
    print('Run tests:')
    for index, condition in enumerate(['ALL', 'UP', 'DOWN'], 1):
        print(f'{index}. {condition} - condition')
        deg_ids = get_deg_gene_ids(deg_table, cond)
        other_ids = get_other_gene_ids(deg_table)
        if method == 'montecarlo':
            deg_scores, other_scores = split_scores_by_gene_ids(scan_results,
                                                                    deg_ids,
                                                                    other_ids)
            results = run_montecarlo(deg_scores, other_scores, threshold)
            create_plot(results, wdir, tag, condition, method)
        elif method == 'fisher':
            deg_scores, other_scores = split_scores_by_gene_ids(best_results,
                                                                    deg_ids,
                                                                    other_ids)
            results = run_fisher(deg_scores, other_scores, threshold)
            create_plot(results, wdir, tag, condition, method)
        else:
            sys.stderr.write("I don't know that method, I exit, See You!\n")
            sys.exit(0)
    print('-'*30)
    print('All done. Exit')
    pass

    
if __name__ == '__main__':
    main()