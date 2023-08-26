from melc import Clustering_FPremoval_CORRECTED_hierarchicalClustering
import argparse


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('adata_pickle_path', type=str, help='Please specify path to anndata pickle data!')
    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    Clustering_FPremoval_CORRECTED_hierarchicalClustering.run(args.adata_pickle_path)
    