from melc import DifferentialAnalysis_MultipleProteinProfiles
import argparse


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('adata_pickle_path', type=str, help='Please specify path to anndata pickle data!')
    parser.add_argument('dependent_variable_name', type=str, help='Please specify the name of the dependent variable. Note: The aforementioned variable should be present under observation metadata (obsm) of the anndata onject containing gene/ protein expression.')
    parser.add_argument('--N', type=int, default=2, help='Please specify number of proteins to be analyzed within protein profile.')
    parser.add_argument('--threshold', type=float, default=0.5, help='Please specify minimum value over which protein expressions are considered for further analyses.')
    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    DifferentialAnalysis_MultipleProteinProfiles.da_multiple_protein_profiles(args.adata_pickle_path, args.dependent_variable_name, args.N, args.threshold)
    