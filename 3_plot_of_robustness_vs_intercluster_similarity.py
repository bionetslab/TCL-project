from melc import SelectHowManyClustersHungarianAlgorithm
import argparse

def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('adata_pickle_path', type=str, help='Please specify path to anndata pickle data!')
    parser.add_argument('--method', type=str, default='top-down', help='Please specify cluster agreement method between clusterings. The two options are ["top-down", "bottom-up"]. The hierarchical clustering method is agglomerative, meaning bottom-up.')
    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    SelectHowManyClustersHungarianAlgorithm.plot_of_robustness_vs_intercluster_similarity(args.adata_pickle_path, args.method)
    