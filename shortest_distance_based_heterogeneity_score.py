from melc import Shortest_Distance_Based_Heterogeneity_Score
import argparse


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('adata_pickle_path', type=str, help='Please specify path to anndata pickle data!')
    parser.add_argument('--negative_class_name', type=str, default='0', help='Please specify negative_class_name (optional).\nValues accepted: any string value.\nIf no value is specified, the default is 0.')
    parser.add_argument('--positive_class_name', type=str, default='1', help='Please specify positive_class_name (optional).\nValues accepted: any string value.\nIf no value is specified, the default is 1.')
    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    Shortest_Distance_Based_Heterogeneity_Score.run(args.adata_pickle_path, args.negative_class_name, args.positive_class_name)
    # Clustering_FPremoval_CORRECTED_hierarchicalClustering.run(args.adata_pickle_path)