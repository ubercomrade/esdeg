import os
import sys
import argparse
from scipy import stats
from enrest.functions import *
from enrest.parsers import matrices_parser, fasta_parser, read_set_of_genes
from enrest.scaner import scaner


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser_name', help='Available commands:')
    deg_parser = subparsers.add_parser('deg', help='Run test on DEGs')
    set_parser = subparsers.add_parser('set', help='Run test on SET of genes')

    deg_parser.add_argument('deg', action='store', help='TSV file with DEG with ..., The NAME column must contain ensemble gene IDS')
    deg_parser.add_argument('matrices', action='store', help='Path to matrices in HOCOMOCO (PCM) or in MEME (PFM) format')
    deg_parser.add_argument('promoters', action='store', choices=['mm10', 'hg38', 'tair10', 'rnor6'], metavar='N',
         help='promoters of organism (hg38, mm10, tair10)')
    deg_parser.add_argument('output', action='store', help='Name of directory for output files')
    deg_parser.add_argument('-p', '--parameter', action='store', choices=['enrichment', 'fraction'],
                        metavar='METHOD', type=str, default='enrichment', 
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

    set_parser.add_argument('set', action='store', help='File with list of genes. Genes must be in Ensemble format (ensemble gene IDS)')
    set_parser.add_argument('matrices', action='store', help='Path to matrices in HOCOMOCO (PCM) or in MEME (PFM) format')
    set_parser.add_argument('promoters', action='store', choices=['mm10', 'hg38', 'tair10', 'rnor6'], metavar='N',
         help='promoters of organism (hg38, mm10, tair10)')
    set_parser.add_argument('output', action='store', help='Name of directory for output files')
    set_parser.add_argument('-p', '--parameter', action='store', choices=['enrichment', 'fraction'],
                        metavar='METHOD', type=str, default='enrichment', 
                        help='Parameter estimated in test (enrichment or fraction), default= enrichment')
    set_parser.add_argument('-f', '--format', action='store', choices=['meme', 'hocomoco'],
                        metavar='FORMAT', type=str, default='meme', 
                        help='Format of file with matrices (meme or hocomoco), default= meme')
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())



def run_test(genes, set_scores, other_scores, threshold_table, parameter):
    results = {('LOW', 'log2Fold'): 0, ('LOW', 'pval'): 0,
               ('MIDDLE', 'log2Fold'): 0, ('MIDDLE', 'pval'): 0,
               ('HIGH', 'log2Fold'): 0, ('HIGH', 'pval'): 0,
               ('COMBINE', 'pval'): 0}
    pvals = []
    genes = np.array(genes)
    for index, level in zip(range(0, len(threshold_table)), ['LOW', 'MIDDLE', 'HIGH']):
        threshold, fpr = threshold_table[index]
        if parameter == "enrichment":
            z_score, fold = montecarlo_enrichment(set_scores, other_scores, threshold)
        elif parameter == "fraction":
            z_score, fold = montecarlo_fraction(set_scores, other_scores, threshold)
        pval = stats.norm.sf(z_score)          
        pvals.append(pval)
        results[(level, 'pval')] = pval
        results[(level, 'log2Fold')] = fold
        # OFF
        # genes_with_bs = genes[np.sum(np.greater_equal(set_scores, threshold), axis=1) > 0]
        # genes_without_bs = genes[np.sum(np.greater_equal(set_scores, threshold), axis=1) == 0]
        # results[(level, 'GeneIdsWithBS')] = '; '.join(genes_with_bs)
        # results[(level, 'GeneIdsWithoutBS')] = '; '.join(genes_without_bs)
    cpval = hartung(np.array(pvals))
    results[('COMBINE', 'pval')] = cpval
    return results


# BINOM and GEOM
# def run_test(deg_scores, other_scores, threshold_table, method, parameter):
#     results = []
#     for index in range(0, len(threshold_table)):
#         threshold, fpr = threshold_table[index]
#         if method == 'montecarlo':
#             if parameter == "enrichment":
#                 z_score = montecarlo_enrichment(deg_scores, other_scores, threshold)
#             elif parameter == "fraction":
#                 z_score = montecarlo_fraction(deg_scores, other_scores, threshold)
#             pval = stats.norm.sf(abs(z_score))
#         else:
#             if parameter == "enrichment":
#                 pval = stat_tests_enrichment(deg_scores, other_scores, threshold, method=method)
#             elif parameter == "fraction":
#                 pval = stat_tests_fraction(deg_scores, other_scores, threshold, method=method)            
#         results.append(pval)
#     cpval = hartung(np.array(results))
#     results.append(cpval)
#     return results


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


def set_case(path_to_set, path_to_db, output_dir, path_to_promoters, 
             file_format="meme", parameter="enrichment"):    
    # WORK IN PROGRESS? PICs
    # logos_dir = f"{output_dir}/logos/"
    # distributions_dir = f"{output_dir}/distributions/"
    # if not os.path.isdir(logos_dir):
    #     os.mkdir(logos_dir)
    # if not os.path.isdir(distributions_dir):
    #     os.mkdir(distributions_dir)
    print('-'*30)
    print('Read SET of genes')
    set_ids = read_set_of_genes(path_to_set)
    print('-'*30)
    print('Read promoters')
    promoters = fasta_parser(path_to_promoters)
    all_ids = set([i[0] for i in promoters])
    print('-'*30)
    print('Read matrices')
    matrices = matrices_parser(path_to_db, f=file_format)
    number_of_matrices = len(matrices)
    print(f'Number of matrices = {number_of_matrices}')
    print('-'*30)
    container = []
    for i in range(number_of_matrices):
        name, pwm, pfm, matrix_length, middle_score = matrices[i]
        line = {('', 'ID'): name}
        print(f'{i+1}. {name}:')
        print('Scan promotrers')
        all_results = scaner(promoters, pwm)
        best_results = get_best_scores(all_results)
        print('Calculate threshold table')
        all_scores_flatten = get_all_flatten_scores(all_results)
        threshold_table = get_threshold(all_scores_flatten)
        threshold_table = np.array(threshold_table)
        fprs_table = threshold_table[:,1]
        fprs_choosen = np.array([0.0005, 0.00015, 0.00005]) # LOW, MIDDLE, HIGH
        indexes = np.searchsorted(fprs_table, fprs_choosen)
        threshold_table = threshold_table[indexes]
        print('Run test')
        other_ids = get_other_gene_ids_for_set_case(set_ids, all_ids)
        if parameter == "enrichment":
            set_scores, other_scores, genes = split_scores_by_gene_ids(all_results, set_ids, other_ids)
        elif parameter == "fraction":
            set_scores, other_scores, genes = split_scores_by_gene_ids(best_results, set_ids, other_ids)
        results = run_test(genes, set_scores, other_scores, threshold_table, parameter)
        # WORK IN PROGRESS
        # write_distrubution = f"{distributions_dir}/{name}.jpg"
        # write_logo = f"{logos_dir}/{name}.jpg"
        # plot_bs_distribution(set_scores, threshold_table, write_distrubution, length=2000, window=20)
        # plot_logo(pfm, write_logo)
        # line.update({('COMBINE', 'Logo'): write_logo,
        #              ('COMBINE', 'BsDistrubution'): write_distrubution})
        line.update(results)
        container.append(line)
        print('-'*30)  
    df = pd.DataFrame(container, columns=container[0].keys())
    output_path = f"{output_dir}/table.all.tsv"
    df.to_csv(output_path, sep='\t', index=False)
    print('All done. Exit')
    return df, set_scores, threshold_table

    
def deg_case(path_to_deg, path_to_db, output_dir, path_to_promoters, 
             file_format='meme', parameter='enrichment', padj_thr=0.05,
             log2fc_thr_deg=1, log2fc_thr_background=np.log2(5/4)):    
    # WORK IN PROGRESS? PICs
    # logos_dir = f"{output_dir}/logos/"
    # distributions_dir = f"{output_dir}/distributions/"
    # if not os.path.isdir(logos_dir):
    #     os.mkdir(logos_dir)
    # if not os.path.isdir(distributions_dir):
    #     os.mkdir(distributions_dir)
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
    container = {'ALL': [],
                 'UP': [],
                 'DOWN': []}
    for i in range(number_of_matrices):
        name, pwm, pfm, matrix_length, middle_score = matrices[i]
        print(f'{i+1}. {name}:')
        print('Scan promotrers')
        all_results = scaner(promoters, pwm)
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
            line = {('', 'ID'): name}
            print(f' {index}. {condition}')
            deg_ids = get_deg_gene_ids(deg_table, condition, padj_thr=padj_thr, log2fc_thr=log2fc_thr_deg)
            other_ids = get_other_gene_ids_for_deg_case(deg_table, padj_thr=padj_thr, log2fc_thr=log2fc_thr_background)
            if parameter == "enrichment":
                deg_scores, other_scores, genes = split_scores_by_gene_ids(all_results, deg_ids, other_ids)
            elif parameter == "fraction":
                deg_scores, other_scores, genes = split_scores_by_gene_ids(best_results, deg_ids, other_ids)
            results = run_test(genes, deg_scores, other_scores, threshold_table, parameter)
            line.update(results)
            # WORK IN PROGRESS? (LOGO and DISTR)
            # write_distrubution = f"{distributions_dir}/{name}.{condition}.jpg"
            # write_logo = f"{logos_dir}/{name}.{condition}.jpg"
            # plot_bs_distribution(deg_scores, threshold_table, write_distrubution, length=2000, window=20)
            # plot_logo(pfm, write_logo)
            # line.update({('COMBINE', 'Logo'): write_logo})
            #             ('COMBINE', 'BsDistrubution'): write_distrubution})
            container[condition].append(line)
    print('-'*30)
    for index, condition in enumerate(['ALL', 'UP', 'DOWN'], 1):
        df = pd.DataFrame(container[condition], columns=container[condition][0].keys())
        condition = condition.lower()
        output_path = f"{output_dir}/table.{condition}.tsv"
        df.to_csv(output_path, sep='\t', index=False)
    print('All done. Exit')
    pass


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
        
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        this_dir, this_filename = os.path.split(__file__)
        if promoters == 'mm10':
            path_to_promoters = os.path.join(this_dir, "../data", "mm10.ensembl.promoters.fa")
        elif promoters == 'hg38':
            path_to_promoters = os.path.join(this_dir, "../data", "hg38.ensembl.promoters.fa")
        elif promoters == 'rnor6':
            path_to_promoters = os.path.join(this_dir, "../data", "rnor6.ensembl.promoters.fa")
            
        deg_case(path_to_deg, 
                 path_to_db,
                 output_dir, 
                 path_to_promoters, 
                 file_format=file_format, 
                 parameter=parameter, 
                 padj_thr=padj_thr,
                 log2fc_thr_deg=log2fc_thr_deg, 
                 log2fc_thr_background=log2fc_thr_background)
        
    elif args.subparser_name == 'set':
        path_to_set = args.set
        path_to_db = args.matrices
        output_dir = args.output
        promoters = args.promoters
        parameter = args.parameter
        file_format = args.format

        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
    
        this_dir, this_filename = os.path.split(__file__)
        if promoters == 'mm10':
            path_to_promoters = os.path.join(this_dir, "../data", "mm10.ensembl.promoters.fa")
        elif promoters == 'hg38':
            path_to_promoters = os.path.join(this_dir, "../data", "hg38.ensembl.promoters.fa")
        elif promoters == 'rnor6':
            path_to_promoters = os.path.join(this_dir, "../data", "rnor6.ensembl.promoters.fa")
        set_case(path_to_set, 
                 path_to_db, 
                 output_dir, 
                 path_to_promoters, 
                 file_format=file_format, 
                 parameter=parameter)   
    pass

    
if __name__ == '__main__':
    main()