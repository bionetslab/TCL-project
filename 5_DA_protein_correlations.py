from melc import DifferentialAnalysis_ProteinCorrelation
import argparse


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('adata_pickle_path', type=str, help='Please specify path to anndata pickle data!')
    parser.add_argument('dependent_variable_name', type=str, help='Please specify the name of the dependent variable. Note: The aforementioned variable should be present under observation metadata (obsm) of the anndata onject containing gene/ protein expression.')
    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    DifferentialAnalysis_ProteinCorrelation.da_protein_correlation(args.adata_pickle_path, args.dependent_variable_name)
    