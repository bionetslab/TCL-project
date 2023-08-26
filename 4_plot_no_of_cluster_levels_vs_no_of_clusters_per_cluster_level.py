from melc import plot_no_of_cluster_levels_vs_no_of_clusters_per_cluster_level
import argparse


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('adata_pickle_path', type=str, help='Please specify path to anndata pickle data!')
    parser.add_argument('clusterings_patientLevel_dict_path', type=str, help='Please specify path to clustering combinations calculated in a previous step. Note: Be careful to select the appropriate clustering path corresponding to the method you are interested in (for example, if you have calculated clusterings for both the "top-down" and "bottom-up approaches", please be sure to set the path accordingly).')
    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    plot_no_of_cluster_levels_vs_no_of_clusters_per_cluster_level.plot_no_of_cluster_levels_vs_no_of_clusters_per_cluster_level(args.adata_pickle_path, args.clusterings_patientLevel_dict_path)
    