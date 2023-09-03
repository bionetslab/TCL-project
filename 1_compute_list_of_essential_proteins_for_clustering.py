from melc import Clustering_FPremoval_CORRECTED_hierarchicalClustering
import argparse


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('adata_pickle_path', type=str, help='Please specify path to anndata pickle data!')
    # parser.add_argument('--user_selected_cluster_level', type=str, default='complete', help='Please specify cluster level to get essential proteins list from (optional).\nValues accepted: {any_integer_value, "complete"}.\nIf value is integer positive, that represents the cluster level. If value is integer  negative, that represents the cluster number in reverse order (for example: -1 represents (n-1)-th cluster). If value is "complete", that represents the last (n-th) cluster.\nFor all other values (i. negative integers; ii. fractional values; iii. any strings other than "complete"; iv. integer value is out of range of [1, n] where n is the total number of clusters): the default is n.\nPlease note: If your input is "complete", just type it in normally without inverted commas/ quotation marks -- the program is capable of reading strings directly without any quotation symbols.')
    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    # Clustering_FPremoval_CORRECTED_hierarchicalClustering.run(args.adata_pickle_path, args.user_selected_cluster_level)
    Clustering_FPremoval_CORRECTED_hierarchicalClustering.run(args.adata_pickle_path)
    