from melc import DAsquidpy_Cluster_Cooccurrence_Across_Spatial_Dimensions_revised
import argparse


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('adata_pickle_path', type=str, help='Please specify path to anndata pickle data!')
    parser.add_argument('dependent_variable_name', type=str, help='Please specify the name of the dependent variable. Note: The aforementioned variable should be present under observation metadata (obsm) of the anndata onject containing gene/ protein expression.')
    parser.add_argument('--interval', default=50, type=int, help='Please specify the no. of intervlas/ bins that you want the cluster co-occurence across spatial dimensions to be computed for. Note: default value is 50 if not specified otherwise.')

    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    DAsquidpy_Cluster_Cooccurrence_Across_Spatial_Dimensions_revised.dasquidpy_cluster_cooccurrence_across_spatial_dimensions(args.adata_pickle_path, args.dependent_variable_name, args.interval)